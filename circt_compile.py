from cellscript import *

REGS_OFFSET = 1000

class Circuit:
    def __init__(self):
        self.cells = []

    def one(self):
        self.cells.append(("1", set(), set()))
        return len(self.cells) - 1

    def _not(self, i):
        self.cells[i][2].add(len(self.cells))
        self.cells.append(("not", {i}, set()))
        return len(self.cells) - 1

    def nor(self, i, j):
        self.cells[i][2].add(len(self.cells))
        self.cells[j][2].add(len(self.cells))
        self.cells.append(("nor", {i, j}, set()))
        return len(self.cells) - 1

    def _in(self, idx):
        # assert False
        self.cells.append(("in", set(), set(), idx))
        return len(self.cells) - 1

    def out(self, i, idx):
        return self.out_not(self._not(i), idx)

    def out_not(self, i, idx):
        self.cells[i][2].add(len(self.cells))
        self.cells.append(("out_not", {i}, set(), idx))
        return len(self.cells) - 1

    def store_not(self, i, idx):
        return self.out_not(i, idx + REGS_OFFSET)

    def store(self, i, idx):
        return self.out(i, idx + REGS_OFFSET)

    def read(self, idx):
        return self._in(idx + REGS_OFFSET)

    def _and(self, i, j):
        return self.nor(self._not(i), self._not(j))

    def _or(self, i, j):
        return self._not(self.nor(i, j))

    def to_tree(self, i):
        return (self.cells[i][0], [self.to_tree(i) for i in self.cells[i][1]], None if len(self.cells[i]) < 4 else self.cells[i][3])

    def to_outs_forest(self):
        forest = []
        for c, cell in enumerate(self.cells):
            if cell[0] == "out_not":
                forest.append(self.to_tree(c))
        return forest

    def depth(self, i):
        return 1 + max(self.depth(j) for j in self.cells[i][1])

    def simulate(self, i, ins):
        if self.cells[i][0] == "in":
            return ins[self.cells[i][3]]
        if self.cells[i][0] == "1":
            return 1
        params = [self.simulate(i, ins) for i in self.cells[i][1]]
        if self.cells[i][0] == "not":
            return 1 - params[0]
        if self.cells[i][0] == "nor":
            return (1 - params[0]) | (1 - params[1])
        if self.cells[i][0] == "out_not":
            return (1 - params[0])

    def simulate_outs(self, ins):
        forest = {}
        for c, cell in enumerate(self.cells):
            if cell[0] == "out_not":
                forest[cell[3]] = self.simulate(c, ins)
        return forest


class IDGenerator:
    def __init__(self, alpha = "GATC"):
        self.next_id = 2
        self.alpha = alpha

    def get_id(self):
        idxs = self.alpha[1:] # doesn't include Gs so IDs can't overlap
        while True:
            n = self.next_id
            self.next_id += 1
            dna = idxs[n % 3]
            while n != 0:
                n //= 3
                dna += idxs[n % 3]
            if GENE_PROMOTER not in dna:
                break
        return self.alpha + dna + self.alpha # prevent accidentally colliding with any gene promoters


class DNAGenerator:
    def __init__(self, ins, stateins, outs, stateouts):
        self.dna = ""
        self.inputs = ins
        self.stateins = stateins
        self.outputs = outs
        self.stateouts = stateouts
        self.idgen = IDGenerator()

    def input(self, i):
        if i >= REGS_OFFSET:
            return self.stateins[i - REGS_OFFSET]
        print(i, self.inputs)
        return self.inputs[i]

    def output(self, i):
        if i >= REGS_OFFSET:
            return self.stateouts[i - REGS_OFFSET]
        return self.outputs[i]

    def from_tree(self, tree):
        if tree[0] == "in":
            post_decode = "".join(AMINO_ACID_TO_DNA[s] for s in self.input(tree[2]))
            return post_decode
        if tree[0] == "out_not":
            id = self.from_tree(tree[1][0])
            self.dna += f" GG{id}{GENE_PROMOTER}{mk_protein(self.output(tree[2]))}"
            return None
        rid = self.idgen.get_id()
        if rid in self.dna:
            # print("rid", rid, "(", self.idgen.next_id - 1, ") already in dna:")
            print(self.dna)
        assert rid not in self.dna
        # print("rid", self.idgen.next_id - 1)
        repressor = mk_represser_gene(rid)

        if tree[0] == "1":
            self.dna += f" GG{GENE_PROMOTER}{repressor}"
        elif tree[0] == "not":
            id = self.from_tree(tree[1][0])
            self.dna += f" GG{id}{GENE_PROMOTER}{repressor}"
        elif tree[0] == "nor":
            id1 = self.from_tree(tree[1][0])
            id2 = self.from_tree(tree[1][1])
            self.dna += f" GG{id1}{GENE_PROMOTER}{id2}{repressor}"
        return rid


