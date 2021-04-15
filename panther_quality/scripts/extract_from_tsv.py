import pandas as pd
import sys
import os

nextflow_bin = os.environ['PATH'].split(os.pathsep)[0]
sys.path.insert(0, nextflow_bin)

metadata_df = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)

samp_id = sys.argv[2]

samp_df = metadata_df.loc[metadata_df['central_sample_id'] == str(samp_id)]

samples_filtered = samp_df[["central_sample_id", "barcode", "artic_primers", "artic_protocol", "joined_run_name", "library_name"]].copy()

for index, row in samples_filtered.iterrows():
	if len(row["barcode"])  == 1:
		samples_filtered.loc[index, "barcode"] = "0" + row["barcode"]

samples_filtered.to_csv(path_or_buf="./" + samp_id + ".csv", index=False)
