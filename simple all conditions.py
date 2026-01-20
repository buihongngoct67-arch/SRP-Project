
import pandas as pd
path = r"G:\Tsukuba\Lab Animal Science\ProjectPython\Mus musculus practice.tsv"
mus_musculus_df = pd.read_csv(path, sep="\t")
print(mus_musculus_df.head())
print(mus_musculus_df.columns)
filtered_df = mus_musculus_df[(mus_musculus_df['FDR'] < 0.01) & (mus_musculus_df["Expression"] == "present")]
non_info_uberon = ['UBERON:0001062', 'UBERON:0000465', 'UBERON:0000061', 'UBERON:0000413', 'UBERON:0000039', 'UBERON:0000586']
uberon_df = filtered_df[(filtered_df['Anatomical entity ID'].str.startswith('UBERON')) & (~filtered_df['Anatomical entity ID'].isin(non_info_uberon))].copy()
grouped_df = uberon_df.groupby('Gene name').agg({'Anatomical entity name': lambda x:', '.join(x.unique()), 'Sex':lambda x:', '.join(x.unique()), 'Developmental stage name': lambda x: ', '.join(x.unique()), 'Anatomical entity ID': lambda x:', '.join(x.unique())}).reset_index()
specific_genes = grouped_df[grouped_df['Anatomical entity name'].apply(len) > 3]['Gene name'].tolist()
result_table = grouped_df[grouped_df['Gene name'].isin(specific_genes)]
final_table = result_table[result_table['Gene name', 'Anatomical entity name', 'Sex', 'Developmental stage name']]
print(final_table)
final_table.to_csv("result.csv")