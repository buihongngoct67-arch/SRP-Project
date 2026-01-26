import pandas as pd
from pathlib import Path
path = r"C:\Users\DELL\Documents\GitHub\SRP-Project\data\Mus musculus practice.tsv"
Path(path).exists()
mus_musculus_df = pd.read_csv(path, sep="\t")

print(mus_musculus_df.columns)
filtered_df = mus_musculus_df[(mus_musculus_df['FDR'] < 0.01) & (mus_musculus_df["Expression"] == "present")]

print(filtered_df.head())
non_info_uberon = ['UBERON:0001062', 'UBERON:0000465', 'UBERON:0000061', 'UBERON:0000413', 'UBERON:0000039', 'UBERON:0000586']
uberon_df = filtered_df[(filtered_df['Anatomical entity ID'].str.startswith('UBERON')) & (~filtered_df['Anatomical entity ID'].isin(non_info_uberon))].copy()

gene_counts = uberon_df.groupby('Gene name')['Anatomical entity name'].nunique()

specific_genes = gene_counts[gene_counts <= 3].index

result_table = uberon_df[uberon_df['Gene name'].isin(specific_genes)]
final_table = result_table[['Gene name', 'Anatomical entity name', 'Sex', 'Developmental stage name', 'Anatomical entity ID']]
print(final_table)
final_table.to_excel('G:\\Tsukuba\\Lab Animal Science\\ProjectPython\\filtered_genes_table.xlsx', index=False)


child_ID = result_table['Anatomical entity ID'].unique()
print(child_ID)

from collections import defaultdict
import networkx as nx
__all__ = ["dfs_edges", "dfs_tree", "dfs_predecessors", "dfs_successors", "descendants", "ancestors"]
def dfs_predecessors(G, source=None, depth_limit=None, )
