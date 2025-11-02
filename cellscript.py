import random

GENE_PROMOTER = "TATAAA"
START_CODON = "AUG" # TAC in DNA
END_CODONS = {"UAA", "UAG", "UGA"} # ATT, ATC, ACT in DNA
DNA_TO_MRNA = {"T": "A", "A": "U", "C": "G", "G": "C"}
MRNA_TO_DNA = {"A": "T", "U": "A", "G": "C", "C": "G"}
INTRON_START = "GU"
INTRON_END = "ACCCCCCCAG"
CODON_TABLE = {
    # Phenylalanine
    "UUU": "F", "UUC": "F",
    # Leucine
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    # Isoleucine
    "AUU": "I", "AUC": "I", "AUA": "I",
    # Methionine (Start)
    "AUG": "M",
    # Valine
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    # Serine
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    # Proline
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "UAU": "Y", "UAC": "Y",
    # Histidine
    "CAU": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAU": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic acid
    "GAU": "D", "GAC": "D",
    # Glutamic acid
    "GAA": "E", "GAG": "E",
    # Cysteine
    "UGU": "C", "UGC": "C",
    # Tryptophan
    "UGG": "W",
    # Arginine
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
CODON_TABLE_REV = {CODON_TABLE[v]: v for v in CODON_TABLE}
AMINO_ACID_TO_DNA = {
    **{aa:'A' for aa in "AVLIPFMW"},   # hydrophobic -> A
    **{aa:'T' for aa in "STYNQC"},      # polar -> T
    **{aa:'G' for aa in "KRH"},        # basic -> G
    **{aa:'C' for aa in "DE"},         # acidic -> C
}
DNA_TO_AMINO_ACID = {"A": "A", "T": "S", "G": "K", "C": "D"}
REPRESSER_OPCODE = "FLI"
REPRESSER_HEAD = 2
REPRESSER_TAIL = 2

def is_valid_protein(p):
    return all(s in CODON_TABLE_REV for s in p)

def transcribe_dna(dna):
    return "".join(DNA_TO_MRNA[s] if s in DNA_TO_MRNA else "" for s in dna)

def rev_transcribe(rna):
    return "".join(MRNA_TO_DNA[s] for s in rna)

def splice_rna(rna):
    exons = ""
    i = 0
    while i < len(rna):
        if rna[i:].startswith(INTRON_START):
            end = rna.find(INTRON_END, i)
            if end != -1:
                i = end + len(INTRON_END)
                continue
        exons += rna[i]
        i += 1
    return exons

def translate(mrna):
    ppc = ""
    for i in range(0, len(mrna), 3):
        if i + 3 > len(mrna):
            continue
        codon = mrna[i:i + 3]
        if codon in CODON_TABLE:
            ppc += CODON_TABLE[codon]
    return ppc

class RNAPolymerase:
    def __init__(self):
        self.ptr = 0

    def transcribe_one_gene(self, ms):
        promoter = ms.pattern_match_once(GENE_PROMOTER, self.ptr)
        if promoter == -1:
            self.ptr = 0
            return None, True
        self.ptr = promoter + 1

        mrna = transcribe_dna(ms.dna[promoter + len(GENE_PROMOTER):])
        mrna = splice_rna(mrna)
        start = mrna.find(START_CODON)
        if start == -1:
            return None, False
        end = start + 3
        while not any(mrna[end:].startswith(codon) for codon in END_CODONS):
            end += 3

        return mrna[start + 3:end], False

class MoleculeSet:
    def __init__(self, dna):
        self.dna = dna
        self.hidden_dna = [None for _ in dna]
        self.proteins = []

    def clear(self):
        self.proteins = []
        self.hidden_dna = [None for _ in self.dna]

    def pattern_match(self, pattern, start=0):
        for i in range(start, len(self.dna) - len(pattern), 1):
            if self.dna[i:].startswith(pattern) and all(self.hidden_dna[i + j] == None for j in range(len(pattern))):
                yield i

    def pattern_match_once(self, pattern, start=0):
        for i in range(start, len(self.dna) - len(pattern), 1):
            if self.dna[i:].startswith(pattern) and all(self.hidden_dna[i + j] == None for j in range(len(pattern))):
                return i
        return -1

    def count(self):
        counts = {}
        for prot in self.proteins:
            if prot in counts:
                counts[prot] += 1
            else:
                counts[prot] = 1
        return counts

    def add(self, protein, p):
        if random.random() > p:
            return False

        if protein.startswith(REPRESSER_OPCODE):
            pattern = "".join(AMINO_ACID_TO_DNA[s] for s in protein[len(REPRESSER_OPCODE):])
            sites = list(self.pattern_match(pattern))
            if len(sites) == 0:
                print(":()")
                self.proteins.append(protein)
                return False
            # site = random.choice(sites)
            site = sites[0]
            if len(sites) > 1:
                print("warning, multiple sites for", protein, "with pattern", pattern)
                for site in sites:
                    print(">", self.dna[site:site+100])

            for i in range(len(REPRESSER_OPCODE) + len(pattern) + REPRESSER_HEAD + REPRESSER_TAIL):
                j = site + i - REPRESSER_HEAD - len(REPRESSER_OPCODE)
                if j >= 0 and j < len(self.dna):
                    self.hidden_dna[j] = protein
            return True
        else:
            self.proteins.append(protein)
        return False

protein_to_dna = lambda p: rev_transcribe("".join(CODON_TABLE_REV[a] for a in p))
mk_protein = lambda p: f"TAC{protein_to_dna(p)}ATT"
mk_represser_gene = lambda idx: mk_protein(REPRESSER_OPCODE + "".join(DNA_TO_AMINO_ACID[s] for s in idx))

def main():
    from circt_compile import Circuit, DNAGenerator, IDGenerator
    from from_yosys import parse_yosys_json_to_circuit
    # c = Circuit()
    # c.out(c._not(c.read(0)), 1)
    # c.out(c._or(c._in(0), c._in(1)), 0)
    # c.out(c.one(), 1)
    # c.out(c.one(), 2)
    # c.out(c.one(), 3)
    # c.out(c.one(), 4)
    idgen = IDGenerator(alpha="YHQNK")
    pre_tag = "Q"*8
    IN_PROT = [f"{pre_tag}INPT{idgen.get_id()}" for _ in range(6)]
    # in[0] = input?
    # in[1..5] = character bits
    OUT_PROTS = [f"{pre_tag}RESLT{idgen.get_id()}" for _ in range(7)]
    # out[0] = exit
    # out[1] = display?
    # out[2..6] = character bits
    c, nregs = parse_yosys_json_to_circuit("./example.json", {
        "clk_i": 0,
        "exit_o": 0,             # output index for exit
        "display_o": 1,          # output index for display (1-bit)
        "char_o": [2,3,4,5]      # output indices for char_o 4-bit bus
    })
    # for x in c.cells:
    #     if x[0] == "in": print(x)
    # print(c.simulate_outs([0]))
    # exit(1)
    # nregs = 1

    # print(c.cells)
    STATE_PROT_IN = [f"{pre_tag}STATEIN{idgen.get_id()}" for _ in range(nregs)]
    STATE_PROT_OUT = [f"{pre_tag}STATERES{idgen.get_id()}" for _ in range(nregs)]

    gen = DNAGenerator(IN_PROT, STATE_PROT_IN, OUT_PROTS, STATE_PROT_OUT)
    for top in c.to_outs_forest():
        gen.from_tree(top)
    # print(gen.dna)

    dna = gen.dna.replace(" ", "")

    print("Simulating DNA:", dna)
    rp = RNAPolymerase()
    ms = MoleculeSet(dna)
    first = True
    cycles = 0
    while True:
        # print(ms.hidden_dna)
        # print(ms.proteins)
    # for i in range(10):
        # print(ms.dna)
        # print("".join("x" if ms.hidden_dna[x] != None else " " for x in range(len(dna))))

        gene, cycled = rp.transcribe_one_gene(ms)
        if gene != None:
            protein = translate(gene)
            print("+ protein:", protein, "from", gene)
            ms.add(protein, 1)
        elif first:
            print("no promoter sequence found, giving up")
            exit(1)
        first = False

        if cycled:
            cycles += 1

        if cycles > 0:
            cycles = 0
            counts = ms.count()
            print(counts,OUT_PROTS)
            results = [int(p in counts) for p in OUT_PROTS]
            state = [int(p in counts) for p in STATE_PROT_OUT]
            print("res", results)
            print("state", state)
            if results[1] != 0:
                bits = 0
                for i in range(5):
                    if results[2 + i] != 0:
                        bits = bits | (1 << i)
                print(chr(bits + ord('a')), end='')
            if results[0] != 0:
                exit()

            ms.clear()
            print(state)
            for p, s in zip(STATE_PROT_IN, state):
                if s != 0:
                    while ms.add(p, 1):
                        pass


if __name__ == "__main__":
    main()

