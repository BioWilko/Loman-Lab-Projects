Channel
  .value("BHAM-COVID19-20415-186")
  .set {library_ch}

Channel
  .from(1..37)
  .set {barcode_ch}

metadata = '/data/nick_home/metadata.tsv'

reference = '/data/homes/samw/projects/variant_and_lineage/MN908947.3.fasta'

process extract_metadata {
  input:
    val lib from library_ch
    val barcode_list from barcode_ch
  output:
    file "metadata.csv" into metadata_csv

  """
  #!/usr/bin/env python

  import pandas as pd
  import os

  barcodes = [${barcode_list}]

  metadata_df = pd.read_csv("${metadata}", sep="\t", dtype=str, keep_default_na=False)

  library_df = metadata_df.loc[metadata_df['library_name'] == str("${lib}")].copy()

  list_of_lists = []

  for barcode in barcodes:
      barcode_df = library_df.loc[library_df['barcode'] == str(barcode)]
      barcode_list = barcode_df[["central_sample_id", "barcode", "artic_primers", "artic_protocol", "joined_run_name", "library_name"]].values.flatten().tolist()
      if len(barcode_list[1]) == 1:
          barcode_list[1] = str(0) + str(barcode_list[1])  #Adds leading zero (shouldn't be necesary but we are where we are)
      list_of_lists.append(barcode_list)

  pd.DataFrame(list_of_lists, columns=["central_sample_id", "barcode", "artic_primers", "artic_protocol", "joined_run_name", "library_name"]).to_csv(path_or_buf="./metadata.csv", index=False)
  """
}

metadata_csv
  .splitCsv(header: true)
  .map{row -> tuple(row.central_sample_id, row.barcode, row."artic-primers", row."artic-protocol", row.run_name, row.library_name)}
  .set {gplex_data_ch}

process guppyplex {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  errorStrategy 'ignore'

  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name from gplex_data_ch
  output:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, file("${library_name}_barcode${barcode}_.fastq") into gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory /cephfs/grid/bham/${library_name}/${library_name}/*/fastq_pass/barcode${barcode}/ --prefix ${library_name}_barcode${barcode}
  """
}

process minion {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, fastq_file from gplex_out_ch
  output:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, file("*.consensus.fasta") into alignment_in_ch
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, file("*.consensus.fasta") into pangolin_in_ch

  """
  artic minion --normalise 200 --threads 16 --scheme-directory $HOME/projects/panther_ext_quality/dataset/primer-schemes --read-file ${fastq_file} --fast5-directory /cephfs/grid/bham/${library_name}/${library_name}/*/fast5_pass/barcode${barcode}// --sequencing-summary /cephfs/grid/bham/${library_name}/${library_name}/*/sequencing_summary* nCoV-2019/V3 ${library_name}_barcode${barcode}_out
  """
}

process alignment {
  conda '/data/homes/samw/miniconda3/envs/mafft_env/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, fasta_file from alignment_in_ch
  output:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, file("${library_name}_barcode${barcode}_msa") into msa_out_ch

  """
  cat ${reference} ${fasta_file} > to_align.fasta
  mafft --auto to_align.fasta > ${library_name}_barcode${barcode}_msa
  """
}

process variant_call {
    publishDir '/data/homes/samw/projects/variant_and_lineage/outputs', mode: "copy", overwrite: true
    errorStrategy 'ignore'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, msa_file from msa_out_ch
  output:
    file("${library_name}_barcode${barcode}_type_summary.csv") into variant_out_ch

  """
  aln2type ./ ./ ${library_name}_barcode${barcode}_type_summary.csv MN908947.3 ${msa_file} /data/homes/samw/aln2type/variant_definitions/variant_yaml/*yml 
  """
}

process pango_lineage {
  publishDir '/data/homes/samw/projects/variant_and_lineage/outputs', mode: "copy", overwrite: true
  conda '/data/homes/samw/miniconda3/envs/pangolin/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, fasta_file from pangolin_in_ch
  output:
    file("${library_name}_barcode${barcode}_lineage.csv") into lineage_out_ch

  """
  pangolin ${fasta_file}
  mv lineage_report.csv ${library_name}_barcode${barcode}_lineage.csv
  """
}