#!/usr/bin/env nextflow

tsv_dir_ch = Channel.value("/data/nick_home/updated_metadata_heartlands.tsv")

project_dir = "~/projects/panther_ext_quality/scripts/"

Channel
  .fromPath("/data/homes/samw/projects/panther_ext_quality/sample_ids")
  .splitText()
  .set {raw_samp_id}

sample_ch = raw_samp_id.flatMap {it.trim()}

process extract_from_meta {
	input:
		path tsv_dir from tsv_dir_ch
    val samp_id from sample_ch

  output:
    file "${samp_id}.csv" into samp_data_ch

	"""
  python ${project_dir}extract_from_tsv.py ${tsv_dir} ${samp_id}
  """
}

samp_data_ch
  .splitCsv(header: true)
  .map{row -> tuple(row.central_sample_id, row.barcode, row."artic-primers", row."artic-protocol", row.joined_run_name, row.library_name)}
  .set {guppy_data_ch}

process guppyplex {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name from guppy_data_ch
  output:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, file("${central_sample_id}_${library_name}_.fastq") into guppyplex_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory /cephfs/grid/bham/${library_name}/${library_name}/${run_name}/fastq_pass/barcode${barcode}/ --prefix ${central_sample_id}_${library_name}
  """
}

process minion {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, reads_file from guppyplex_ch
  output:
    set central_sample_id, library_name, file("multiqc_data/multiqc_data.json") into multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 16 --scheme-directory $HOME/projects/panther_ext_quality/dataset/primer-schemes --read-file ${reads_file} --fast5-directory /cephfs/grid/bham/${library_name}/${library_name}/${run_name}/fast5_pass/barcode${barcode}// --sequencing-summary /cephfs/grid/bham/${library_name}/${library_name}/${run_name}/sequencing_summary* nCoV-2019/V3 ${central_sample_id}_${library_name}_out
  multiqc .
  """
}

process json_to_csv {
  input:
    set central_sample_id, library_name, file(json) from multiqc_out_ch
  output:
    file "${central_sample_id}_${library_name}.csv" into coverage_ch
  """
  #!/usr/bin/env python

  import pandas as pd
  import os
  import json

  data = json.load(open("./${json}"))
  coverage = data["report_plot_data"]["custom_data_linegraph"]["datasets"][0][0]["data"]
  pd.DataFrame({"${central_sample_id}_${library_name}":coverage}).to_csv(path_or_buf="${central_sample_id}_${library_name}.csv", index=False)
  """
}

process combine {
  publishDir '/nextflow_out/', mode: 'copy', overwrite: true
  input:
    file csv_list from coverage_ch.toList()
  output:
    file "combined.tsv" into out_ch

  """
  echo "" > combined.tsv
  seq 1 98 >> combined.tsv

  for filename in ${csv_list}; do paste combined.tsv \$filename > tmpout && mv tmpout combined.tsv; done
  """
}