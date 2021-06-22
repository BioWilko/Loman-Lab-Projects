import pandas as pd
import sys
import random

location = sys.argv[1]
start = sys.argv[2] #ISO format (YYYY-MM-DD)
end = sys.argv[3] #ISO format (YYYY-MM-DD)
n_samples = int(sys.argv[4])

fasta_dir = "/cephfs/covid/bham/artifacts/published/latest/fasta/"
out_dir = sys.argv[5]

out_fasta = open("%s%s_%s-%s.fasta" %(out_dir, location.replace(" ", "_"), start, end), "w")

df = pd.read_csv("/cephfs/covid/bham/artifacts/published/latest/majora.metadata.tsv", sep="\t", index_col="central_sample_id", dtype="object")

df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")

mask = (df['collection_date'] > start) & (df['collection_date'] <= end)

df = df.loc[mask]

df = df.loc[df["adm2"] == location]

cog_ids = df.index.tolist()

subset = random.sample(cog_ids, n_samples)
subset_df = df[df.index.isin(subset)].copy()

for id in subset:
	try:
		run_group = subset_df.loc[id]["run_group"]
		fasta_path = "%s%s.%s.climb.fasta" %(fasta_dir, id, run_group)
		fasta = open(fasta_path)
		out_fasta.write(fasta.read())
		fasta.close()
	except:
		subset_df.drop(labels=id, inplace=True, errors="ignore")

subset_df.to_csv("%s%s_%s-%s.tsv" %(out_dir, location.replace(" ", "_"), start, end), sep="\t", index_label="central_sample_id")