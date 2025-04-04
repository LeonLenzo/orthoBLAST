Script: orthoblast.py

Description:
    - Extracts gene sequences from GFF and genome FASTA for two species.
    - Performs a one-way BLASTn from Species 1 genes → Species 2 genes.
    - Ranks all hits by percent identity and saves the top 100 matches.

Inputs:
    species1.fasta   → Genome FASTA for Species 1
    species1.gff     → GFF3 annotation for Species 1 (must include "gene" features)
    species2.fasta   → Genome FASTA for Species 2
    species2.gff     → GFF3 annotation for Species 2

Requirements:
    - Biopython
    - BCBio.GFF
    - pandas
    - BLAST+ tools installed (makeblastdb, blastn)

Output:
    - top_100_conserved_genes.csv
T
