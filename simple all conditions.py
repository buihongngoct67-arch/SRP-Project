
import pandas as pd
path = r"G:\Tsukuba\Lab Animal Science\ProjectPython\Mus musculus practice.tsv"
mus_musculus_df = pd.read_csv(path, sep="\t")



print(mus_musculus_df.columns)
filtered_df = mus_musculus_df[(mus_musculus_df['FDR'] < 0.01) & (mus_musculus_df["Expression"] == "present")]

print(filtered_df.head())
non_info_uberon = ['UBERON:0001062', 'UBERON:0000465', 'UBERON:0000061', 'UBERON:0000413', 'UBERON:0000039', 'UBERON:0000586']
uberon_df = filtered_df[(filtered_df['Anatomical entity ID'].str.startswith('UBERON')) & (~filtered_df['Anatomical entity ID'].isin(non_info_uberon))].copy()

specific_genes_df = uberon_df[uberon_df['Anatomical entity name'].apply(len) < 4]['Gene name']
result_table = uberon_df[uberon_df['Gene name'].isin(specific_genes_df)]
final_table = result_table[['Gene name', 'Anatomical entity name', 'Sex', 'Developmental stage name']]
print(final_table)
final_table.to_csv("result_lt4.csv")