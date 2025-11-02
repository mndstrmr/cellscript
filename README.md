# Cellscript

A hackathon esolang project for Oxford Compsoc, 01/11/2025, written in one night.

- Code is written as DNA, which is a sequence of genes.
- Genes are headered by the gene promoter (TATAAA).
- The gene is scanned by RNA Polymerase, which transcribes it into RNA.
- In the process introns are spliced out, leaving behind exons (comments ignored by the compiler!)
- Somewhere near the beggining is a start codon (AUG). Codons are three tuples which map to amino acids via tRNA.
- The sequence of codons is translated continuously until a stop codon is found (UAA, UAG or UGA).
- Amino acids get joined into a polypeptide chain, which becomes a protein (though protein folding is not implemented, for obvious reasons)

Proteins are heavily involved with all functions in the cell, most interestingly the transcription process itself. By binding to operation sites in the DNA, specific proteins
can block genes from being expressed. This is called gene repression.

In our model, a protein who's amino acid sequence begins with FLI (Phenylalanine, Leucine, Isoleucine) can bind to the DNA. Binding may occur at any location on the DNA where the amino acid sequence
following the FLI matches the DNA. We say that hydrophobic AAs map to A, polar AAs to T, basic AAs to G and acidic to C.

When a protein binds to DNA, the region it binds to, and two base pairs in front and behind, become blocked from transcription. This means that a gene promoter can be skipped, and the
entire gene will fail to be expressed. We can identify the gene we wish to disable by attaching an additional unique sequence either immediately before the gene promoter, immediately after
(before the start codon), or both, like so:

```
[id 1] TATAAA [id 2] --- (any junk) --- AUG (protein encoded as DNA) --- UAA
```
Therefore, if a protein binds to id 1 or id 2, the gene will be expressed. The same is true if a protein binds to both. This makes a NOR gate for gene expression: `not (a or b)`.

Nor gates form singleton universal gate sets, and so we can build arbitrary circuits with them.

See hello.txt for a DNA sequence for a program that prints hello. It is compiled from hello.sv (SystemVerilog is a hardware definition language, used to
define circuits such as these). We use yosys to produce an aigmap (only AND and NOT gates), which we then export as JSON and load in Python. We provide state in the form of additional
inputs and outputs.

Owing to the limited time constraints, and the frankly apalling time at which I am writing this, the code is not nice, and probably barely works. The hello program takes too long to have
observed success, but it seems to be getting there. Smaller programs have been seen to work.

In practice, none of this is that biologically realistic, but it is interesting, and perhaps serves as some kind of proof that cells are pretty damn powerful.
