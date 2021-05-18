import java.text.SimpleDateFormat

def printHelp() {
  log.info"""
  Mandatory arguments:
    --location          Set location code ["Heartlands", "QE", "UHCW", "Shrewsbury"]

  Optional arguments: (Ensure trailing "/" is included for all paths set this way)
    --start_date        Set start date for report in format YYYY-MM-DD (default 01-01-2020)
    --end_date          Set start end date for report in format YYYY-MM-DD (defaults to current date)
    --report_dir        Set report output directory (default "/data/homes/samw/projects/cog_reporting/reports/")
    --data_dir          Set directory containing ARTIC pipeline out (default "/data/nick_home/")
    --metadata_path     Set path to metadata file (default "/data/nick_home/updated_metadata_heartlands.tsv")
"""
}

def date = new Date()
def sdf = new SimpleDateFormat("yyyy-MM-dd")
end_date = sdf.format(date)
start_date = "2020-01-01"

out_dir = "/data/homes/samw/projects/cog_reporting/reports/"
data_dir = "/data/nick_home/"
metadata_path = "/data/nick_home/updated_metadata_heartlands.tsv"


if (params.help){
  printHelp()
  exit 0
}

if (params.start_date){
  start_date = params.start_date
}

if (params.end_date){
  end_date = params.end_date
}

if (params.report_dir){
  out_dir = params.report_dir
}

if (params.data_dir){
  data_dir = params.data_dir
}

if (params.metadata_path){
  metadata_path = params.metadata_path
}

aln2type_headers = "/data/homes/samw/projects/cog_reporting/aln2type_headers.csv"

originating_lab = params.location

process generate_ids {
  output:
    file("ids.csv") into id_file_ch

  """
  #!/usr/bin/env python
  import pandas as pd
  import numpy as np

  start_date = "${start_date}"
  end_date = "${end_date}"

  df = pd.read_csv("${metadata_path}", sep='\t', keep_default_na=False)

  for index, row in df.iterrows():
    if pd.isnull(row["collection_date"]) == True and pd.isnull(row["received_date"]) == False:
      df.at[index, "collection_date"] = df.at[index, "received_date"]

  df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")

  mask = (df['collection_date'] > start_date) & (df['collection_date'] <= end_date)
  
  df = df.loc[mask]

  samples = df.query("originating_lab == '${originating_lab}'")

  uniques = samples["central_sample_id"].unique()
  pd.DataFrame(uniques).to_csv("ids.csv", index=False, header=False)
  """
}

process flag_undated_samples {
  publishDir out_dir, overwrite:true
  output:
    file("${originating_lab}_${end_date}_undated.csv") optional true

  """
  #!/usr/bin/env python
  import pandas as pd

  df = pd.read_csv("${metadata_path}", sep="\t", keep_default_na=False)

  location_df = df.query("originating_lab == '${originating_lab}'")

  mask = (location_df['collection_date'].isnull() == True) & (location_df['received_date'].isnull() == True)

  undated_df = location_df.loc[mask]

  if undated_df["central_sample_id"].empty == False:
    undated_df["central_sample_id"].to_csv("${originating_lab}_${end_date}_undated.csv", index=False, header=False)
  """
  }

process filter {
  input:
    file ids from id_file_ch
  output:
    file("filtered.csv") into filtered_ch
    file("fails.csv") into fails_ch

  """
  #!/usr/bin/env python
  import pandas as pd
  import os.path
  from os import path

  ids = open("${ids}", "r")
  lines = ids.readlines()

  exists = []
  not_exist = []

  for line in lines:
    id = line.strip()
    if os.path.isfile("${data_dir}" + id + ".nanopolish-indel.consensus.fasta"):
      exists.append(id)
    else:
      not_exist.append(id)

  pd.DataFrame(exists).to_csv("filtered.csv", index=False, header=False)
  pd.DataFrame(not_exist).to_csv("fails.csv", index=False, header=False)
  """
}

filtered_ch
  .splitText()
  .map{it -> it.trim()}
  .into {pango_id_ch;voc_id_ch;snp_id_ch;report_id_ch}

fails_ch
  .splitText()
  .map{it -> it.trim()}
  .set {fails_ch}

process aln2type_voc {
  input:
    val samp_id from voc_id_ch
  output:
    file("${samp_id}_voc.csv") optional true into voc_multi_ch
  
  """
  awk -F'.' '/^>/{print \$1; next}{print}' < ${data_dir}${samp_id}.nanopolish-indel.muscle.out.fasta > ${samp_id}.fasta
  aln2type ./ ./ ${samp_id}_voc.csv MN908947 ${samp_id}.fasta --no_call_deletion /data/homes/samw/aln2type/variant_definitions/variant_yaml/*.yml
  """
}

process aln2type_snp {
  input:
    val samp_id from snp_id_ch
  output:
    file("${samp_id}_snp.csv") optional true into snp_multi_ch
  
  """
  awk -F'.' '/^>/{print \$1; next}{print}' < ${data_dir}${samp_id}.nanopolish-indel.muscle.out.fasta > ${samp_id}.fasta
  aln2type ./ ./ ${samp_id}_snp.csv MN908947 ${samp_id}.fasta --no_call_deletion --output_unclassified /data/homes/samw/aln2type/variant_definitions/SNP_reporting_defs/*.yml
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
  cat ${data_dir}{${id_str}}.nanopolish-indel.consensus.fasta > temp_consensus.fasta
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
    val fails from fails_ch.toList()
  output:
    file("basic_info.tsv") into basic_info_ch
  script:
    id_str = ids.join('","')
    fails_str = fails.join('","')

  """
  #!/usr/bin/env python
  import pandas as pd

  data_path = "${data_dir}"

  pango_df = pd.read_csv("${lineage_summary}", index_col="taxon")
  metadata_df = pd.read_csv(("${metadata_path}"), sep="\t")

  samp_ids = ["${id_str}","${fails_str}"]
  report_dict = {}

  for id in samp_ids:
    try:
      cov_df = pd.read_csv((data_path + id + ".cov.txt"), index_col="sample", sep="\t")
      if len(metadata_df.loc[metadata_df["central_sample_id"] == id]) == 1:
        samp_df = metadata_df.loc[metadata_df["central_sample_id"] == id]
        report_dict[id] = dict.fromkeys(["COG ID", "Sender ID", "Sender Suborg", "Sample Collection Date", "Sample Received Date", "Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
        report_dict[id]["COG ID"] = id
        report_dict[id]["Sender ID"] = samp_df["sender_sample_id"].iloc[0]
        report_dict[id]["Sender Suborg"] = samp_df["collecting_suborg"].iloc[0]
        report_dict[id]["Sample Collection Date"] = samp_df["collection_date"].iloc[0]
        report_dict[id]["Sample Received Date"] = samp_df["received_date"].iloc[0]
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
          report_dict[id] = dict.fromkeys(["COG ID", "Sender ID", "Collector Suborg", "Sample Collection Date", "Sample Received Date", "Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
          report_dict[id]["COG ID"] = id
          report_dict[id]["Sender ID"] = samp_df["sender_sample_id"].iloc[0]
          report_dict[id]["Sender Suborg"] = samp_df["collecting_suborg"].iloc[0]          
          report_dict[id]["Sample Collection Date"] = samp_df["collection_date"].iloc[0]
          report_dict[id]["Sample Received Date"] = samp_df["received_date"].iloc[0]
          report_dict[id]["Coverage (%)"] = cov_df.at[id, "perc"]
          report_dict[id]["Genome Recovered?"] = "Y" if cov_df.at[id, "perc"] >= 50 else "N"
          report_dict[id]["Lineage"] = pango_df.at[id, "lineage"]      
          report_dict[id]["Sender Ct Value"] = samp_df["ct_1_ct_value"].iloc[0]
          report_dict[id]["Incident Code"] = samp_df["collecting_org"].iloc[0]
          report_dict[id]["Uploaded to COG-UK?"] = "Y" if cov_df.at[id, "perc"] >= 50 else "N"
          report_dict[id]["Uploaded to GISAID?"] = "Y" if cov_df.at[id, "perc"] >= 90 else "N"
        except:
          report_dict[id] = dict.fromkeys(["COG ID", "Sender ID", "Collector Suborg","Sample Collection Date", "Sample Received Date", "Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
          report_dict[id]["COG ID"] = id
          report_dict[id]["Genome Recovered?"] = "RPT NOT FOUND"
          report_dict[id]["Coverage (%)"] = cov_df.at[id, "perc"]
    except:
      report_dict[id] = dict.fromkeys(["COG ID", "Sender ID", "Sender Suborg", "Sample Collection Date", "Sample Received Date", "Sender Ct Value", "Genome Recovered?", "Coverage (%)", "Lineage", "Incident Code", "Uploaded to COG-UK?", "Uploaded to GISAID?"])
      samp_df = metadata_df.loc[metadata_df["central_sample_id"] == id]
      report_dict[id]["COG ID"] = id
      report_dict[id]["Genome Recovered?"] = "N"
      report_dict[id]["Coverage (%)"] = "0"     
      report_dict[id]["Sample Collection Date"] = samp_df["collection_date"].iloc[0]
      report_dict[id]["Sample Received Date"] = samp_df["received_date"].iloc[0]
      report_dict[id]["Sender ID"] = samp_df["sender_sample_id"].iloc[0]
      report_dict[id]["Sender Suborg"] = samp_df["collecting_suborg"].iloc[0]

  basic_df = pd.DataFrame.from_dict(report_dict, orient="index").to_csv("basic_info.tsv", index=False, sep="\t")
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

  df = df[df["status"] != "low-qc"]
  df = df[df["status"] != "alt-low-qc"]
  df = df[df["phe-label"] != "E484K"]

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
  mutations = ["E484K", "N501Y", "K417T", "K417N", "L452R"]

  dictionary = {}

  status_dict = {"confirmed":"Y", "low-qc":"X"}

  unclassified = df.loc[df["phe-label"] == "unclassified"]

  for index, row in unclassified.iterrows():
    dictionary[row["sample_id"]] = dict.fromkeys(mutations, "N")

  for samp_id in ids:
    samp_df = group_df.get_group(samp_id)
    if samp_df["phe-label"].any() != "unclassified":
      dictionary[samp_id] = { i : "N" for i in mutations }
      for index, row in samp_df.iterrows():
        dictionary[samp_id][row["phe-label"]] = status_dict[row["status"]]

  pd.DataFrame.from_dict(dictionary, orient='index').to_csv("parsed_snps.csv", index_label="COG ID")
  """
}

process generate_report {
  conda '/data/homes/samw/miniconda3/envs/xlsxwriter'
  publishDir out_dir, overwrite: true
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
  import numpy as np

  basic_df = pd.read_csv("${basic_info}", index_col="COG ID", sep="\t")
  snps_df = pd.read_csv("${parsed_snps}", index_col="COG ID")
  vocs_df = pd.read_csv("${parsed_vocs}", index_col="COG ID")

  voc_vui_df = vocs_df.join(snps_df, how="outer")

  report_df = basic_df.join(voc_vui_df, how="outer")
  report_df.sort_values("Sample Collection Date", inplace=True)

  mutations = ["E484K", "N501Y", "K417T", "K417N", "L452R"]

  for index, row in report_df.iterrows():
    if row["E484K"] == "N" and row["N501Y"] == "N" and row["K417T"] == "N" and row["K417N"] == "N" and row["L452R"] == "N" and row["Coverage (%)"] <= 26:
      for i in mutations:
        report_df.at[index, i] = np.nan

  table_shape = "A1:U" + str(len(report_df))

  dfs = {"Results":report_df}

  writer = pd.ExcelWriter("${originating_lab}_${end_date}_report.xlsx", engine='xlsxwriter')
  for sheetname, df in dfs.items():  # loop through `dict` of dataframes
      df.to_excel(writer, sheet_name=sheetname)  # send df to writer
      workbook = writer.book
      worksheet = writer.sheets[sheetname]  # pull worksheet object
      negative = workbook.add_format({'bg_color':   '#FFC7CE', 'font_color': '#9C0006'})
      positive = workbook.add_format({'bg_color':   '#C6EFCE', 'font_color': '#006100'})
      error = workbook.add_format({'bg_color':      '#FF0000', 'font_color': '#FFFFFF'})
      worksheet.conditional_format("Q1:Q1048576", {"type":"cell", "criteria":"equal to", "value":'"N"', "format":negative})
      worksheet.conditional_format("R1:R1048576", {"type":"cell", "criteria":"equal to", "value":'"N"', "format":negative})
      worksheet.conditional_format("S1:S1048576", {"type":"cell", "criteria":"equal to", "value":'"N"', "format":negative})
      worksheet.conditional_format("T1:T1048576", {"type":"cell", "criteria":"equal to", "value":'"N"', "format":negative})
      worksheet.conditional_format("U1:U1048576", {"type":"cell", "criteria":"equal to", "value":'"N"', "format":negative})
      worksheet.conditional_format("Q1:Q1048576", {"type":"cell", "criteria":"equal to", "value":'"Y"', "format":positive})
      worksheet.conditional_format("R1:R1048576", {"type":"cell", "criteria":"equal to", "value":'"Y"', "format":positive})
      worksheet.conditional_format("S1:S1048576", {"type":"cell", "criteria":"equal to", "value":'"Y"', "format":positive})
      worksheet.conditional_format("T1:T1048576", {"type":"cell", "criteria":"equal to", "value":'"Y"', "format":positive})
      worksheet.conditional_format("U1:U1048576", {"type":"cell", "criteria":"equal to", "value":'"Y"', "format":positive})
      worksheet.conditional_format("Q1:Q1048576", {"type":"cell", "criteria":"equal to", "value":'"X"', "format":error})
      worksheet.conditional_format("R1:R1048576", {"type":"cell", "criteria":"equal to", "value":'"X"', "format":error})
      worksheet.conditional_format("S1:S1048576", {"type":"cell", "criteria":"equal to", "value":'"X"', "format":error})
      worksheet.conditional_format("T1:T1048576", {"type":"cell", "criteria":"equal to", "value":'"X"', "format":error})
      worksheet.conditional_format("U1:U1048576", {"type":"cell", "criteria":"equal to", "value":'"X"', "format":error})
      for idx, col in enumerate(df):  # loop through all columns
          series = df[col]
          max_len = max((
              series.astype(str).map(len).max(),  # len of largest item
              len(str(series.name))  # len of column name/header
              )) + 2  # adding a little extra space
          worksheet.set_column(idx, idx, max_len)  # set column width
      worksheet.autofilter(table_shape)
      worksheet.freeze_panes(1, 0)
  writer.save()
  """
}