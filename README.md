# Manuscript code
Pipeline for Identification of metagenomic viral contigs

Identification of viral contigs in metagenome data

(Ubuntu)1.Use CLC Genomics Workbench to assemble metagenomic data and do read mapping (version 11.0.1, QIAGEN Bioinformatics, Denmark)

(Ubuntu)2.Use Prodigal to predict ORF(V2.6.3)

(Ubuntu)3.Use diamond to align Refseq protein database (DIAMOND v0.9.25.126)

(Ubuntu)4.Use Megan6 to assign taxonomy (Community version)

(Ubuntu)5.Use hmmsearch to align Viral Protein Family(VPF) model database (HMMER 3.1b2)

(Ubuntu)6.Use Python script(confirm-virus-like-contigs.py) to filter contigs with 5 more hits.(python3)

(Ubuntu)7.Use Python script(extract-contig-length.py) to filter the 5 kb contig length.(python3)

(Ubuntu)8.Use Python script(extract-orf-sequences.py) to extract the orf sequence.(python3)

(Ubuntu)9.Use orf from step 8 to align with Pfam database (HMMER 3.1b2)

(Ubuntu)10.Use orf from step 8 to align with KEGG database (DIAMOND v0.9.25.126)

(Ubuntu)11.Use Python script(further-confirm-whether-virus-contigs-with-naturepaper-methods.py) to check if it matches any of the three filters.(python3)

(Ubuntu)12.Use Python script(combine_different_steps.py) to combine virus contigs from two different methods.(python3)
