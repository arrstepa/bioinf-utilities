# bioinf-utilities

## What is it?
bioinf-utilities is a toolkit for basic processing of DNA and RNA sequences and fastq files.

## Requirements

- Python 3.8+
- Biopython==1.78

Install dependencies:

```bash
pip install biopython==1.78
```

## NB: `SeqUtils.GC()` was removed in BioPython 1.80 and replaced with `SeqUtils.gc_fraction()`, but here GC is used (with BioPython >1.8 it will cause incorrect GC count

## Main Features
- DNA/RNA/AA sequences processing
	- Amino acids: one can obtain amino acid content for sequence
	- Nucleic acids: capability to transcribe DNA to mRNA, get complements, reverse sequences, and obtain reverse complements.
- fastq files filtration via BioPython based on GC content, phred quality and length
- fasta files procession
	- package allows user to convert multiline fasta files into oneline format

## Class hierarchy

```
BiologicalSequence (ABC)
├── NucleicAcidSequence
│   ├── DNASequence
│   └── RNASequence
└── AminoAcidSequence
```

## Where to Get It
The source code is currently hosted on GitHub at: https://github.com/arrstepa/bioinf-utilities.git



## Authors:
Arina Stepanova
