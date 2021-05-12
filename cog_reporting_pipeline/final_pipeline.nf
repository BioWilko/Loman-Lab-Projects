data_dir = "/data/nick_home/"
metadata = "/data/nick_home/updated_metadata_heartlands.tsv"
aln2type_headers = "/data/homes/samw/projects/cog_reporting/aln2type_headers.csv"

out_dir = '/data/homes/samw/projects/cog_reporting/reports/'

originating_lab = "QE"
start_date = "2020-01-01"
end_date = "2021-05-12"

process generate_ids {
  output:
    file("ids.csv") into id_file_ch

  """
  #!/usr/bin/env python
  import pandas as pd
  import numpy as np

  start_date = "${start_date}"
  end_date = "${end_date}"

  df = pd.read_csv("${metadata}", sep='\t', keep_default_na=False)
  
  df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")

  mask = (df['collection_date'] > start_date) & (df['collection_date'] <= end_date)
  
  df = df.loc[mask]

  samples = df.query("originating_lab == '${originating_lab}'")

  uniques = samples["central_sample_id"].unique()
  pd.DataFrame(uniques).to_csv("ids.csv", index=False, header=False)
  """
}

process filter {
  input:
    file ids from id_file_ch
  output:
    file("filtered.csv") into filtered_ch

  """
  #!/usr/bin/env python
  import pandas as pd
  import os.path
  from os import path

  ids = open("${ids}", "r")
  lines = ids.readlines()

  exists = []

  for line in lines:
    id = line.strip()
    if os.path.isfile("/data/nick_home/" + id + ".nanopolish-indel.consensus.fasta"):
      exists.append(id)

  df = pd.DataFrame(exists).to_csv("filtered.csv", index=False, header=False)
  """
}

filtered_ch
  .splitText()
  .map{it -> it.trim()}
  .into {pango_id_ch;aln_id_ch;report_id_ch}

process aln2type {
  errorStrategy 'ignore'
  input:
    val samp_id from aln_id_ch
  output:
    file("${samp_id}_voc.csv") into voc_multi_ch
    file("${samp_id}_snp.csv") into snp_multi_ch
  
  """
  awk -F'.' '/^>/{print \$1; next}{print}' < /data/nick_home/${samp_id}.nanopolish-indel.muscle.out.fasta > ${samp_id}.fasta
  aln2type ./ ./ ${samp_id}_voc.csv MN908947 ${samp_id}.fasta /data/homes/samw/aln2type/variant_definitions/variant_yaml/*.yml
  aln2type ./ ./ ${samp_id}_snp.csv MN908947 ${samp_id}.fasta /data/homes/samw/aln2type/variant_definitions/SNP_reporting_defs/*.yml
  """
}

process combine_voc_out {
  input:
    file voc_summaries from voc_multi_ch.toList()
  output:
    file("voc_summary.csv") into voc_summary_ch
  script:
    voc_csv_str = voc_summaries.join(',')

  """
  cat ${aln2type_headers} > voc_summary.csv
  tail -n +2 -q {${voc_csv_str}} >> voc_summary.csv
  """
}

process combine_snp_out {
  input:
    file snp_summaries from snp_multi_ch.toList()
  output:
    file("snp_summary.csv") into snp_summary_ch
  script:
    snp_csv_str = snp_summaries.join(',')

  """
  cat ${aln2type_headers} > snp_summary.csv
  tail -n +2 -q {${snp_csv_str}} >> snp_summary.csv
  """
}

process cat_consensus {
  input:
    val samp_ids from pango_id_ch.toList()
  output:
    file("combined_consensus.fasta") into consensus_out_ch
  script:
    id_str = samp_ids.join(',')

  """
  cat /data/nick_home/{${id_str}}.nanopolish-indel.consensus.fasta > temp_consensus.fasta
  awk -F'.' '/^>/ {print \$1; next}{print}' < temp_consensus.fasta > combined_consensus.fasta
  """
}

process pangolin {
  conda '/data/homes/samw/miniconda3/envs/pangolin/'
  input:
    file combined from consensus_out_ch
  output:
    file("lineage_report.csv") into lineage_summary

  """
  pangolin ${combined}
  """
}

process parse_basic_info{
  input:
    file lineage_summary from lineage_summary
    val ids from report_id_ch.toList()
  output:
    file("basic_info.csv") into basic_info_ch
  script:
    id_str = ids.join('","')

  """
  #!/usr/bin/env python
  import pandas as pd

  data_path = "/data/nick_home/"

  pango_df = pd.read_csv("${lineage_summary}", index_col="taxon")
  metadata_df = pd.read_csv((data_path + "updated_metadata_heartlands.tsv"), sep="\t")

  samp_ids = ["${id_str}"]
  report_dict = {}

  for id in samp_ids:
    cov_df = pd.read_csv((data_path + id + ".cov.txt"), index_col="sample", sep="\t")
    if len(metadata_df.loc[metadata_df["central_sample_id"] == id]) == 1:
      samp_df = metadata_df.loc[metadata_df["central_sample_id"] == id]
      report_dict[id] = dict.fromkeys(["COG ID", "Sender ID","Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
      report_dict[id]["COG ID"] = id
      report_dict[id]["Sender ID"] = samp_df["sender_sample_id"].iloc[0].replace(",",".")
      report_dict[id]["Sample Collection Date"] = samp_df["collection_date"].iloc[0]      
      report_dict[id]["Genome Recovered?"] = "Y" if cov_df.at[id, "perc"] >= 50 else "N"
      report_dict[id]["Coverage (%)"] = cov_df.at[id, "perc"]
      report_dict[id]["Lineage"] = pango_df.at[id, "lineage"]
      report_dict[id]["Sender Ct Value"] = samp_df["ct_1_ct_value"].iloc[0]
      report_dict[id]["Incident Code"] = samp_df["collecting_org"].iloc[0]
      report_dict[id]["Uploaded to COG-UK?"] = "Y" if cov_df.at[id, "perc"] >= 50 else "N"
      report_dict[id]["Uploaded to GISAID?"] = "Y" if cov_df.at[id, "perc"] >= 90 else "N"
    else:
      try:
        repeat_info = metadata_df.loc[metadata_df["exclude"].isnull()]
        samp_df = repeat_info.loc[metadata_df["central_sample_id"] == id]
        report_dict[id] = dict.fromkeys(["COG ID", "Sender ID", "Sample Collection Date", "Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
        report_dict[id]["COG ID"] = id
        report_dict[id]["Sender ID"] = samp_df["sender_sample_id"].iloc[0].replace(",",".")
        report_dict[id]["Sample Collection Date"] = samp_df["collection_date"].iloc[0]
        report_dict[id]["Coverage (%)"] = cov_df.at[id, "perc"]
        report_dict[id]["Genome Recovered?"] = "Y" if cov_df.at[id, "perc"] >= 50 else "N"
        report_dict[id]["Lineage"] = pango_df.at[id, "lineage"]      
        report_dict[id]["Sender Ct Value"] = samp_df["ct_1_ct_value"].iloc[0]
        report_dict[id]["Incident Code"] = samp_df["collecting_org"].iloc[0]
        report_dict[id]["Uploaded to COG-UK?"] = "Y" if cov_df.at[id, "perc"] >= 50 else "N"
        report_dict[id]["Uploaded to GISAID?"] = "Y" if cov_df.at[id, "perc"] >= 90 else "N"
      except:
        report_dict[id] = dict.fromkeys(["COG ID", "Sender ID","Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
        report_dict[id]["COG ID"] = id
        report_dict[id]["Genome Recovered?"] = "RPT NOT FOUND"
        report_dict[id]["Coverage (%)"] = cov_df.at[id, "perc"]

  basic_df = pd.DataFrame.from_dict(report_dict, orient="index").to_csv("basic_info.csv", index=False)
  """
}

process parse_VOCs {
  input:
    file voc_summary from voc_summary_ch
  output:
    file("parsed_vocs.csv") into parsed_voc_ch

  """
  #!/usr/bin/env python
  import pandas as pd

  df = pd.read_csv("${voc_summary}", usecols=["sample_id", "phe-label", "unique-id", "status"])
  df = df[df.status != "low-qc"]

  df.rename({"sample_id":"COG ID", "phe-label":"PHE Label", "unique-id":"Unique ID", "status":"Status"}, inplace=True, axis=1)

  df.to_csv("parsed_vocs.csv", index=False)
  """
}

process parse_SNPs {
  input:
    file snp_summary from snp_summary_ch
  output:
    file("parsed_snps.csv") into parsed_snp_ch

  """
  #!/usr/bin/env python
  import pandas as pd

  df = pd.read_csv("${snp_summary}")

  group_df = df.groupby(["sample_id"])

  ids = df["sample_id"].unique()
  mutations = ["E484K", "N501Y", "K417T", "K417N"]

  dict = {}

  status_dict = {"confirmed":"Y", "low-qc":"X"}

  for samp_id in ids:
    samp_df = group_df.get_group(samp_id)
    dict[samp_id] = { i : "N" for i in mutations }
    for index, row in samp_df.iterrows():
      dict[samp_id][row["phe-label"]] = status_dict[row["status"]]

  pd.DataFrame.from_dict(dict, orient='index').to_csv("parsed_snps.csv", index_label="COG ID")
  """
}

process generate_report {
  conda '/data/homes/samw/miniconda3/envs/xlsxwriter'
  publishDir out_dir
  input:
    file basic_info from basic_info_ch
    file parsed_snps from parsed_snp_ch
    file parsed_vocs from parsed_voc_ch
  output:
    file("*report.xlsx") into result_ch

  """
  #!/usr/bin/env python
  
  import pandas as pd
  import os

  basic_df = pd.read_csv("${basic_info}", index_col="COG ID")
  snps_df = pd.read_csv("${parsed_snps}", index_col="COG ID")
  vocs_df = pd.read_csv("${parsed_vocs}", index_col="COG ID")

  voc_vui_df = vocs_df.join(snps_df, how="outer")

  report_df = basic_df.join(voc_vui_df, how="outer")
  dfs = {"Results":report_df}

  writer = pd.ExcelWriter("${originating_lab}_${start_date}_${end_date}_report.xlsx", engine='xlsxwriter')
  for sheetname, df in dfs.items():  # loop through `dict` of dataframes
      df.to_excel(writer, sheet_name=sheetname)  # send df to writer
      worksheet = writer.sheets[sheetname]  # pull worksheet object
      for idx, col in enumerate(df):  # loop through all columns
          series = df[col]
          max_len = max((
              series.astype(str).map(len).max(),  # len of largest item
              len(str(series.name))  # len of column name/header
              )) + 3  # adding a little extra space
          worksheet.set_column(idx, idx, max_len)  # set column width
  writer.save()
  """
}