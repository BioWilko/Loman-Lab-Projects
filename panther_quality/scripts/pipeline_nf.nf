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
  .into {guppy_data_ch}

process guppyplex {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name from guppy_data_ch
  output:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, file("${central_sample_id}_${library_name}*") into guppyplex_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory /cephfs/grid/bham/${library_name}/${library_name}/${run_name}/fastq_pass/barcode${barcode}/ --prefix ${central_sample_id}_${library_name}
  """
}

process minion {
  publishDir '$HOME/projects/panther_ext_quality/multiqc_out/'
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set central_sample_id, barcode, "artic-primers", "artic-protocol", run_name, library_name, reads_file from guppyplex_ch
  output:
    file "${central_sample_id}_${library_name}_out*" into minion_out_ch

  """
  artic minion --strict --normalise 200 --threads 16 --scheme-directory $HOME/projects/panther_ext_quality/dataset/primer-schemes --read-file ${reads_file} --fast5-directory /cephfs/grid/bham/${library_name}/${library_name}/${run_name}/fast5_pass/barcode${barcode}// --sequencing-summary /cephfs/grid/bham/${library_name}/${library_name}/${run_name}/sequencing_summary* nCoV-2019/V3 ${central_sample_id}_${library_name}_out
  multiqc .
  """
}