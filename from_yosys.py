# Warning: This file was written by ChatGPT at 4am. Proceed with caution.

from circt_compile import Circuit, DNAGenerator, IDGenerator

import json
import re
from collections.abc import Iterable

def _normalize_cell_type(t: str) -> str:
    """Normalize yosys cell type strings to canonical names: AND, NOT, DFF_P."""
    if t is None:
        return ""
    tt = re.sub(r'[\$_]', '', str(t)).upper()
    if 'AND' in tt:
        return 'AND'
    if 'NOT' in tt:
        return 'NOT'
    if 'DFF' in tt:
        return 'DFF_P'
    return tt

def parse_yosys_json_to_circuit(yosys_source, port_map):
    """
    Parse Yosys JSON into a Circuit instance using only AND / NOT / DFF_P cells.

    Behavior:
      - Does NOT allocate implicit inputs. Unknown numeric nets are treated as constant 0.
      - Yosys string bits are treated as constants ("1" -> 1, others -> 0).
      - If the user provides a mapping for an input port bit in port_map, we create c._in(index)
        even if the yosys bit is a string constant (we honor explicit mappings).
      - Output ports are connected with c.out(driver_cell, out_index).
    Args:
      yosys_source: path to .json file OR a Python dict already loaded.
      port_map: mapping port name -> int or list[int].
    Returns:
      (c, num_registers)
    """
    # load json if a path string provided
    if isinstance(yosys_source, str):
        with open(yosys_source, 'r') as fh:
            data = json.load(fh)
    elif isinstance(yosys_source, dict):
        data = yosys_source
    else:
        raise ValueError("yosys_source must be a path or dict")

    modules = data.get("modules", {})
    if not modules:
        raise ValueError("No 'modules' found in yosys JSON")
    top_name = data.get("top")
    if top_name and top_name in modules:
        top = modules[top_name]
    else:
        top_name = next(iter(modules))
        top = modules[top_name]

    c = Circuit()

    # map: net id (int) or string constant token -> driving cell index
    net_map = {}

    # helper constants cached under these keys
    def const_one():
        if 'CONST1' not in net_map:
            net_map['CONST1'] = c.one()
        return net_map['CONST1']
    def const_zero():
        if 'CONST0' not in net_map:
            net_map['CONST0'] = c._not(const_one())
        return net_map['CONST0']

    # resolve a yosys net *value* (int net id or string constant)
    def resolve_net_source_value(netval):
        if isinstance(netval, str):
            return const_one() if netval == "1" else const_zero()
        if isinstance(netval, int):
            # numeric net id: return driver if already present; else unresolved
            return net_map.get(netval, None)
        return None

    # helper that returns a driver cell index for a netval, using constants if necessary.
    # returns a tuple (driver_cell_index, resolved_bool)
    def get_driver_for(netval):
        if isinstance(netval, str):
            return (const_one() if netval == "1" else const_zero(), True)
        if isinstance(netval, int):
            if netval in net_map:
                return (net_map[netval], True)
            else:
                return (None, False)
        return (const_zero(), True)

    # ---- process ports to honour explicit input mappings ----
    module_ports = top.get("ports", {})
    # iterate ports; create c._in for mapped input bits, otherwise leave nets unresolved (they become 0)
    for pname, pinfo in module_ports.items():
        direction = pinfo.get("direction") if isinstance(pinfo, dict) else None
        bits = pinfo.get("bits") if isinstance(pinfo, dict) else None
        if bits is None:
            continue

        mapped = port_map.get(pname)
        if mapped is not None:
            if isinstance(mapped, Iterable) and not isinstance(mapped, (str, bytes)):
                mapped_list = list(mapped)
            else:
                mapped_list = [mapped]
        else:
            mapped_list = None

        if direction == "input":
            for i, bit in enumerate(bits):
                if mapped_list and i < len(mapped_list):
                    idx = mapped_list[i]
                    # if yosys gives an int net id, map that net id to an input cell
                    if isinstance(bit, int):
                        net_map[bit] = c._in(idx)
                    else:
                        # yosys bit is a string constant but user provided mapping -> honor mapping
                        net_map[bit] = c._in(idx)
                else:
                    # no mapping: do not create _in; leave unresolved (will be treated as const0 unless it's "1")
                    if isinstance(bit, str):
                        # string constants are known
                        net_map[bit] = (const_one() if bit == "1" else const_zero())
                    # integer bits left absent from net_map (unresolved)

        else:
            # for outputs and others, do nothing now; outputs handled later
            # but if yosys indicates constant string bits, record them
            for bit in bits:
                if isinstance(bit, str) and bit not in net_map:
                    net_map[bit] = (const_one() if bit == "1" else const_zero())

    # ---- collect cell descriptions for dependency processing ----
    cells = top.get("cells", {})
    # We'll keep a dict cellname -> parsed info
    parsed_cells = {}
    for cellname, cell in cells.items():
        ctype = _normalize_cell_type(cell.get("type") or cell.get("cell_type"))
        if ctype not in ("AND", "NOT", "DFF_P"):
            raise ValueError(f"Unsupported cell type '{cell.get('type')}' in cell {cellname}")

        connections = cell.get("connections") or cell.get("ports") or {}
        conn = {}
        if isinstance(connections, dict):
            for p, v in connections.items():
                if isinstance(v, list):
                    conn[p] = v
                else:
                    conn[p] = [v]
        parsed_cells[cellname] = {
            "type": ctype,
            "conn": conn,
            "processed": False
        }

    reg_count = 0
    # For DFFs we will record pending stores for which D driver isn't resolved yet:
    pending_stores = []   # list of tuples (reg_idx, D_netval) where D_netval is int or str

    progress = True
    while progress:
        progress = False
        for name, info in parsed_cells.items():
            if info["processed"]:
                continue
            ctype = info["type"]
            conn = info["conn"]

            if ctype == "DFF_P":
                # we can always allocate register index and set Q read cell immediately.
                D_list = conn.get("D") or conn.get("d") or []
                Q_list = conn.get("Q") or conn.get("q") or []
                if not D_list or not Q_list:
                    raise ValueError(f"DFF_P cell {name} missing D/Q connections: {conn}")
                D_net = D_list[0]
                Q_net = Q_list[0]
                reg_idx = reg_count
                reg_count += 1
                # Q driven by read(reg_idx)
                net_map[Q_net] = c.read(reg_idx)
                # record pending store - we will resolve D later when possible
                pending_stores.append((reg_idx, D_net))
                info["processed"] = True
                progress = True
                # continue - D resolution deferred

            elif ctype == "NOT":
                A_list = conn.get("A") or conn.get("a") or conn.get("IN") or conn.get("in") or []
                Y_list = conn.get("Y") or conn.get("y") or conn.get("OUT") or conn.get("out") or []
                if not A_list or not Y_list:
                    raise ValueError(f"NOT cell {name} missing A/Y connections: {conn}")
                A_net = A_list[0]
                # check if A is resolvable
                a_driver, a_ok = get_driver_for(A_net)
                if not a_ok:
                    # can't instantiate yet
                    continue
                # instantiate
                not_sig = c._not(a_driver)
                net_map[Y_list[0]] = not_sig
                info["processed"] = True
                progress = True

            elif ctype == "AND":
                A_list = conn.get("A") or conn.get("a") or []
                B_list = conn.get("B") or conn.get("b") or []
                Y_list = conn.get("Y") or conn.get("y") or []
                if not A_list or not B_list or not Y_list:
                    raise ValueError(f"AND cell {name} missing A/B/Y connections: {conn}")
                A_net = A_list[0]
                B_net = B_list[0]
                a_driver, a_ok = get_driver_for(A_net)
                b_driver, b_ok = get_driver_for(B_net)
                if not (a_ok and b_ok):
                    # inputs not yet resolvable; defer
                    continue
                and_sig = c._and(a_driver, b_driver)
                net_map[Y_list[0]] = and_sig
                info["processed"] = True
                progress = True

    # End of iterative instantiation pass.
    # Any remaining unprocessed combinational cells are unresolved (due to missing inputs).
    # Per your rule: do NOT allocate implicit inputs; set their outputs to 0.
    for name, info in parsed_cells.items():
        if info["processed"]:
            continue
        ctype = info["type"]
        conn = info["conn"]
        # set output to zero for NOT/AND that couldn't be processed
        if ctype in ("NOT", "AND"):
            Y_list = conn.get("Y") or conn.get("y") or conn.get("OUT") or conn.get("out") or []
            if Y_list:
                net_map[Y_list[0]] = const_zero()
        elif ctype == "DFF_P":
            # This should not happen: DFFs were always marked processed because we always allocated Q.
            # But just in case, ensure Q set to read and register allocated.
            D_list = conn.get("D") or conn.get("d") or []
            Q_list = conn.get("Q") or conn.get("q") or []
            if Q_list:
                reg_idx = reg_count
                reg_count += 1
                net_map[Q_list[0]] = c.read(reg_idx)
                pending_stores.append((reg_idx, D_list[0] if D_list else None))

    # Now resolve pending stores for DFFs (use const zero for unresolved D nets)
    for reg_idx, D_net in pending_stores:
        driver, ok = get_driver_for(D_net)
        if not ok or driver is None:
            driver = const_zero()
        c.store(driver, reg_idx)

    # Finally hook up outputs using port_map (use c.out as you requested)
    for pname, pinfo in module_ports.items():
        direction = pinfo.get("direction") if isinstance(pinfo, dict) else None
        if direction != "output":
            continue
        bits = pinfo.get("bits", [])
        if pname not in port_map:
            raise ValueError(f"Output port '{pname}' not present in provided port_map")
        mapped = port_map[pname]
        if isinstance(mapped, Iterable) and not isinstance(mapped, (str, bytes)):
            mapped_list = list(mapped)
        else:
            mapped_list = [mapped]
        if len(mapped_list) < len(bits):
            raise ValueError(f"Output port '{pname}' mapping length mismatch: expected {len(bits)} indices")
        for i, bit in enumerate(bits):
            out_idx = mapped_list[i]
            # driver: if present in net_map use it; otherwise treat as constant per yosys rule
            if bit in net_map and net_map[bit] is not None:
                driver = net_map[bit]
            else:
                # if bit is a string const use it; otherwise numeric unmapped -> zero
                driver = const_one() if isinstance(bit, str) and bit == "1" else const_zero()
            c.out(driver, out_idx)

    return c, reg_count
