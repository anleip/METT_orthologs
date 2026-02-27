# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
## Goal is to view synteny between matching proteins
## First get ordering of proteins from gene locus, file downloaded from NCBI -> view annotated genomes
bu_annotation = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/GCF_044361425.1_annotation.tsv",delimiter='\t')
pv_annotation = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/GCF_000012825.1_annotation.tsv",delimiter='\t')
bu_annotation.columns = [col.strip().replace(" ", "_") for col in bu_annotation.columns]
pv_annotation.columns = [col.strip().replace(" ", "_") for col in pv_annotation.columns]
bu_annotation.columns = [f"bu_{col}" for col in bu_annotation.columns]
pv_annotation.columns = [f"pv_{col}" for col in pv_annotation.columns]
bu_annotation['bu_order'] = np.arange(1,len(bu_annotation)+1)
pv_annotation['pv_order'] = np.arange(1,len(pv_annotation)+1)
# %%
print(bu_annotation)
# %%
print(pv_annotation)
# %%
# Import and filter BLASTp output
blastp = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/pv_v_bu_blast_protein.crunch",delimiter="\t",names=["pv_Protein_accession","bu_Protein_accession","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
blastp_loose = blastp[(blastp['evalue']<1e-5) & (blastp['pident']>=30)].copy()
blastp_strict = blastp[(blastp['evalue']<1e-10) & (blastp['pident']>=30) & (blastp['length']>=70)].copy()
blastp_strict_35 = blastp[(blastp['evalue']<1e-10) & (blastp['pident']>=35) & (blastp['length']>=70)].copy()
blastp_strict_50 = blastp[(blastp['evalue']<1e-10) & (blastp['pident']>=50) & (blastp['length']>=70)].copy()
blastp_strict_70 = blastp[(blastp['evalue']<1e-10) & (blastp['pident']>=70) & (blastp['length']>=70)].copy()
blastp_strict_90 = blastp[(blastp['evalue']<1e-10) & (blastp['pident']>=90) & (blastp['length']>=70)].copy()

blastp_loose["logE"] = -np.log10(blastp_loose["evalue"].replace(0, 1e-300))
blastp_strict["logE"] = -np.log10(blastp_strict["evalue"].replace(0, 1e-300))
blastp_strict_35["logE"] = -np.log10(blastp_strict_35["evalue"].replace(0, 1e-300))
blastp_strict_50["logE"] = -np.log10(blastp_strict_50["evalue"].replace(0, 1e-300))
blastp_strict_70["logE"] = -np.log10(blastp_strict_70["evalue"].replace(0, 1e-300))
blastp_strict_90["logE"] = -np.log10(blastp_strict_90["evalue"].replace(0, 1e-300))

# %%
print(len(blastp_strict_35))
# Merge order into BLASTp output table
# %%
blastp_loose_ordered = blastp_loose.merge(bu_annotation,on='bu_Protein_accession').merge(pv_annotation,on='pv_Protein_accession')
blastp_strict_ordered = blastp_strict.merge(bu_annotation,on='bu_Protein_accession').merge(pv_annotation,on='pv_Protein_accession')
blastp_strict_35_ordered = blastp_strict_35.merge(bu_annotation,on='bu_Protein_accession').merge(pv_annotation,on='pv_Protein_accession')
blastp_strict_50_ordered = blastp_strict_50.merge(bu_annotation,on='bu_Protein_accession').merge(pv_annotation,on='pv_Protein_accession')
blastp_strict_70_ordered = blastp_strict_70.merge(bu_annotation,on='bu_Protein_accession').merge(pv_annotation,on='pv_Protein_accession')
blastp_strict_90_ordered = blastp_strict_90.merge(bu_annotation,on='bu_Protein_accession').merge(pv_annotation,on='pv_Protein_accession')

# %%
print(blastp_loose_ordered)
print(blastp_strict_ordered)
# %%
print(blastp_loose_ordered['pv_order'],blastp_loose_ordered['bu_order'])

# %%
plt.style.use('dark_background')
plt.figure(figsize=(7,7))
sc = plt.scatter(
    blastp_loose_ordered['pv_order'],
    blastp_loose_ordered['bu_order'],
    c=blastp_loose_ordered["logE"],     # colour by E-value
    s=5,
    cmap="viridis"
)
plt.xlabel("P. vulgatus protein order")
plt.ylabel("B. uniformis protein order")
plt.title('BLASTp PV v.s. BU (Loose)')
cbar = plt.colorbar(sc)
cbar.set_label("-log10(E-value)")
plt.show()

# %%
# plot function
def dot_plot(x,y,col,title):
    plt.style.use('dark_background')
    plt.figure(figsize=(7,7))
    sc = plt.scatter(
        x,
        y,
        c=col,     # colour by E-value
        s=5,
        cmap="viridis"
    )
    plt.xlabel("P. vulgatus protein order")
    plt.ylabel("B. uniformis protein order")
    plt.title(title)
    cbar = plt.colorbar(sc)
    cbar.set_label("-log10(E-value)")
    plt.show()
# %%
dot_plot(blastp_strict_35_ordered['pv_order'],
    blastp_strict_35_ordered['bu_order'],
    blastp_strict_35_ordered["logE"], "p identity >35%")
# %%
dot_plot(blastp_strict_50_ordered['pv_order'],
    blastp_strict_50_ordered['bu_order'],
    blastp_strict_50_ordered["logE"], "p identity >50%")
# %%
dot_plot(blastp_strict_70_ordered['pv_order'],
    blastp_strict_70_ordered['bu_order'],
    blastp_strict_70_ordered["logE"], "p identity >70%")
# %%
dot_plot(blastp_strict_90_ordered['pv_order'],
    blastp_strict_90_ordered['bu_order'],
    blastp_strict_90_ordered["logE"], "p identity >90%")
# %%
#
#
# Using essential list from the collaborators
lab_df = pd.read_csv('/nfs/research/jlees/anlei/analysis/METT_ortholog/essential_genes_uniformis_vulgatus_table_only.csv')
print(lab_df)
# %%
# Can't find common gene ID or locus tag to associate the essential gene list to gene order!!