import pandas as pd

df = pd.read_csv(r'C:\Users\nazim\Desktop\My Stuff\patchr_db\watchdb.abatch.txt', sep='\t', dtype=str)
df['assay_method'] = df['assay_method'].astype(str)
longest = df['assay_method'].apply(len).sort_values(ascending=False).head(10)
print(longest)
