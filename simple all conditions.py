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

test_column_table = result_table[['Gene ID', 'Gene name', 'Anatomical entity ID',
       'Anatomical entity name', 'Developmental stage ID',
       'Developmental stage name', 'Sex', 'Strain', 'Expression',
       'Call quality', 'FDR', 'Expression score', 'Expression rank']]
print(test_column_table.head())
test_column_table.to_excel('G:\\Tsukuba\\Lab Animal Science\\ProjectPython\\test_column_table.xlsx', index=False)


child_ID = result_table['Anatomical entity ID'].unique()
print(child_ID)

import networkx as nx

obo_file_path = r"C:\Users\DELL\Documents\GitHub\SRP-Project\uberon.obo"

def parse_obo_to_networkx(obo_file_path):
    G = nx.DiGraph()
    current_id = None
    with open(obo_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if line == '[Term]':
                current_id = None
            elif line.startswith('id:'):
                current_id = line.split('id:')[1].strip()
                G.add_node(current_id)
            elif line.startswith('is_a:') and current_id in child_ID:
                parent_id = line.split('is_a:')[1].split('!')[0].strip()
                G.add_edge(parent_id, current_id)
            elif line.startswith("relationship: part_of") and current_id in child_ID:
                parent_id = line.split("relationship: part_of")[1].split('!')[0].strip()
                G.add_edge(parent_id, current_id)
    return G

def graph_to_parent_child_table(G):
    data = []
    for parent, child in G.edges():
        data.append({'Child': child, 'Parent': parent})
    return pd.DataFrame(data)

G = parse_obo_to_networkx(obo_file_path)
df = graph_to_parent_child_table(G)
print(df.head())
print(df)

        


