#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import glob
from Bio import SeqIO
from BCBio import GFF

species1_fasta = "GCA_016801405.1_ASM1680140v1_genomic.fna"
species1_gff = "GCA_016801405.1_ASM1680140v1_genomic.gff"
species2_fasta = "GCA_003171515.3_CUR_PTRM4_2.2_genomic.fna"
species2_gff = "GCA_003171515.3_CUR_PTRM4_2.2_genomic.gff"

species1_genes = "species1_genes.fasta"
species2_genes = "species2_genes.fasta"
blast_output = "blast_1to2.tsv"
output_csv = "top_100_conserved_genes.csv"


def extract_gene_sequences(genome_fasta, gff_file, output_fasta):
    sequences = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    with open(gff_file) as gff_handle, open(output_fasta, "w") as out_fasta:
        for rec in GFF.parse(gff_handle, base_dict=sequences):
            for feature in rec.features:
                if feature.type == "gene":
                    qualifiers = feature.qualifiers
                    gene_id = (
                        qualifiers.get("ID", [None])[0]
                        or qualifiers.get("Name", [None])[0]
                        or qualifiers.get("gene", [None])[0]
                        or "unknown"
                    )
                    gene_seq = feature.extract(rec.seq)
                    SeqIO.write(SeqIO.SeqRecord(gene_seq, id=gene_id, description=""), out_fasta, "fasta")

def make_blast_db(fasta_file, db_name):
    subprocess.run(
        ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", db_name],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

def run_blast(query, db, output_file):
    subprocess.run([
        "blastn", "-query", query, "-db", db,
        "-out", output_file, "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
        "-max_target_seqs", "1", "-evalue", "1e-5"
    ], check=True, stderr=subprocess.DEVNULL)

def load_blast_table(blast_file):
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
    df = pd.read_csv(blast_file, sep="\t", names=cols, header=None)
    return df

def cleanup_temp_files():
    files_to_remove = glob.glob("species1_db*") + \
                      glob.glob("species2_db*") + [
                          species1_genes,
                          species2_genes,
                          blast_output
                      ]
    for f in files_to_remove:
        try:
            os.remove(f)
            #print(f"🧹 Removed: {f}")
        except Exception as e:
            print(f"⚠️ Could not remove {f}: {e}")

def main():
    # Step 1: Extract genes from each genome
    extract_gene_sequences(species1_fasta, species1_gff, species1_genes)
    extract_gene_sequences(species2_fasta, species2_gff, species2_genes)

    # Step 2: Build BLAST DB
    make_blast_db(species2_genes, "species2_db")

    # Step 3: Run BLAST (Species1 → Species2)
    run_blast(species1_genes, "species2_db", blast_output)

    # Step 4: Load BLAST hits and rank by percent identity
    df = load_blast_table(blast_output)
    top100 = df.sort_values(by="pident", ascending=False).head(200)

    # Step 5: Save output
    final_table = top100[[
        "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"
    ]].rename(columns={
        "qseqid": "Species1_Gene",
        "sseqid": "Species2_Gene",
        "pident": "Percent_Identity",
        "length": "Alignment_Length",
        "evalue": "E-value",
        "bitscore": "Bitscore",
    })
    final_table.to_csv(output_csv, index=False)
    print(f"✅ Top 200 matches by identity saved to: {output_csv}")

    # Step 6: Cleanup
    cleanup_temp_files()

if __name__ == "__main__":
    main()
