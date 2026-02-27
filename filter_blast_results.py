# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
blast_default = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/pv_v_bu_blast_default.crunch",delimiter="\t",names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
# %%
print(blast_default)
# %%
blastp = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/pv_v_bu_blast_protein.crunch",delimiter="\t",names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
# %%
print(blastp)
# %%
print(len(blastp))
# %%
blastp_loose = blastp[(blastp['evalue']<1e-5) & (blastp['pident']>=30)]
blastp_strict = blastp[(blastp['evalue']<1e-10) & (blastp['pident']>=30) & (blastp['length']>=70)]

# %%
print(blastp_strict)
# %%
# the reverse blastp result
blastp_r = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/bu_v_pv_blast_protein.crunch",delimiter="\t",names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
# %%
blastp_r_strict = blastp_r[(blastp_r['evalue']<1e-10) & (blastp_r['pident']>=30) & (blastp_r['length']>=70)]
blastp_r_strict.columns = [f"r_{col}" for col in blastp_r_strict.columns]
print(blastp_r_strict)
# %%
print(blastp_r)
# %%
names = blastp_strict.merge(blastp_r_strict,left_on=['qseqid','sseqid'],right_on=['r_sseqid','r_qseqid'],validate="1:1")
# %%
print(names)
# %%
reciprocal = names[['qseqid','sseqid']].copy()
# %%
print(reciprocal)
# %%
print(len(blastp_strict))
# %%
print(len(blastp_r_strict))
print(blastp_r_strict['r_qseqid'][0])
# %%
print(names.iloc[:,[0,1,12,13]])

# %%
plt.hist(blastp_strict['length'],bins=100)
plt.title('matching length distribution')
plt.show()
# %%
plt.hist(blastp_strict['pident'],bins=100)
plt.show()
# %%
# Looking for intersecting matches
tokens = []; tokens_r = []
for i, row in blastp_strict.iterrows():
    tokens.append([row['qseqid'],row['sseqid']])
for j, row in blastp_r_strict.iterrows():
    tokens_r.append([row['r_sseqid'],row['r_qseqid']])
# %%
intersect_q = []; intersect_s = []
for i, t in enumerate(tokens):
    if t in tokens_r:
        intersect_q.append(blastp_strict['qseqid'].iloc[i])
        intersect_s.append(blastp_strict['sseqid'].iloc[i])

# %%
print(len(intersect_q))
# %%
# make interseting df!
#interstct_df = pd.DataFrame(data=[intersect_q,intersect_s],columns=['pv_Protein_accession','bu_Protein_accession'])
# %%
