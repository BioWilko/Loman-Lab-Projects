import pandas as pd
import math
import sys

summary_path = sys.argv[1]
out_path = sys.argv[2]

df_sum = pd.read_csv(summary_path, delimiter="\t")


data = {"qc_pass":{"0":0}, "qc_fail":{"0":0}}

for index, row in df_sum.iterrows():
    read_time = row["start_time"] + row["duration"]
    read_mins = math.ceil(read_time / 60)
    if row["passes_filtering"]:
        if int(read_mins) in data["qc_pass"]:
            data["qc_pass"][int(read_mins)] += 1
        else:
            data["qc_pass"][int(read_mins)] = 0
    else:
        if int(read_mins) in data["qc_fail"]:
            data["qc_fail"][int(read_mins)] += 1
        else:
            data["qc_fail"][int(read_mins)] = 0

df_out = pd.DataFrame(data).to_csv(out_path)