

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