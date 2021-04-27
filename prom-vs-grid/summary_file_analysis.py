import pandas as pd
import math

df_chunk = pd.read_csv("/data/BHAM-COVID19-20415-163-164/BHAM-COVID19-20415-163-164/20210302_1506_2-E1-H1_PAG50332_b95b5f4e/sequencing_summary_PAG50332_e57674a9.txt", delimiter="\t", chunksize=1000000)


data = {"qc_pass":{"0":0}, "qc_fail":{"0":0}}

for chunk in df_chunk:
    for index, row in chunk.iterrows():
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

df_out = pd.DataFrame(data).to_csv("~/projects/prom-vs-grid/seq_summaries/prom_163-4_rot.csv")
