Channel
  .value("BHAM-COVID19-20415_100_103_prom")
  .set { prom_lib_ch }

Channel
  .value("BHAM-COVID19-20415-100")
  .set { grid_1_lib_ch }

Channel
  .value("BHAM-COVID19-20415-103")
  .set { grid_2_lib_ch }

Channel
  .value("~/projects/prom-vs-grid/100-103-barcoder_out")
  .set { prom_dir_ch }

Channel
  .value("~/projects/prom-vs-grid/100_fastq_combined")
  .set { grid_1_dir_ch }

Channel
  .value("~/projects/prom-vs-grid/103_fastq_combined")
  .set { grid_2_dir_ch }

Channel
  .of('01', '02', '03', '04', '05', '06', '07', '08', '09', 10..96)
  .set { prom_barcode_ch }

Channel
  .of('01', '02', '03', '04', '05', '06', '07', '08', '09', 10..48)
  .set { grid_1_barcode_ch }

Channel
  .of(49..96)
  .set { grid_2_barcode_ch }

reference = '/data/homes/samw/projects/variant_and_lineage/MN908947.3.fasta'

process guppyplex_prom {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    val barcode from prom_barcode_ch
    val library from prom_lib_ch
    val dir from prom_dir_ch
  output:
    set dir, barcode, library, file("prom_${barcode}_.fastq") into prom_gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${dir}/barcode${barcode}/ --prefix prom_${barcode}
  """
}

process minion_prom {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    set dir, barcode, library, fastq_file from prom_gplex_out_ch
  output:
    set barcode, library, file("multiqc_data/multiqc_data.json") into prom_multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 4 --read-file ${fastq_file} --scheme-directory $HOME --fast5-directory /cephfs/grid/bham/${library}/*/*/fast5_pass/ --sequencing-summary /cephfs/grid/bham/${library}/*/*/sequencing_summary* nCoV-2019/V3 prom_${barcode}_out
  multiqc .
  """
}

process guppyplex_grid_1 {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    val barcode from grid_1_barcode_ch
    val library from grid_1_lib_ch
    val dir from grid_1_dir_ch
  output:
    set barcode, library, file("grid_${barcode}_.fastq") into grid_1_gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${dir}/barcode${barcode}/ --prefix grid_${barcode}
  """
}

process minion_grid_1 {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    set barcode, library, fastq_file from grid_1_gplex_out_ch
  output:
    set barcode, library, file("multiqc_data/multiqc_data.json") into grid_1_multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 4 --read-file ${fastq_file} --scheme-directory $HOME --fast5-directory ~/projects/prom-vs-grid/100-103_fast5/barcode${barcode} --sequencing-summary ~/projects/prom-vs-grid/100_seq_summary.txt nCoV-2019/V3 grid_${barcode}_out
  multiqc .
  """
}

process guppyplex_grid_2 {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    val barcode from grid_2_barcode_ch
    val library from grid_2_lib_ch
    val dir from grid_2_dir_ch
  output:
    set barcode, library, file("grid_${barcode}_.fastq") into grid_2_gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${dir}/barcode${barcode}/ --prefix grid_${barcode}
  """
}

process minion_grid_2 {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'
  
  input:
    set barcode, library, fastq_file from grid_2_gplex_out_ch
  output:
    set barcode, library, file("multiqc_data/multiqc_data.json") into grid_2_multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 4 --read-file ${fastq_file} --scheme-directory $HOME/ --fast5-directory ~/projects/prom-vs-grid/100-103_fast5/barcode${barcode} --sequencing-summary ~/projects/prom-vs-grid/103_seq_summary.txt nCoV-2019/V3 grid_${barcode}
  multiqc .
  """
}

process prom_json_to_csv {
  input:
    set barcode, library, json from prom_multiqc_out_ch
  output:
    file "prom_${barcode}.csv" into prom_coverage_ch

  """
  #!/usr/bin/env python

  import pandas as pd
  import os
  import json

  data = json.load(open("${json}"))
  coverage = data["report_plot_data"]["custom_data_linegraph"]["datasets"][0][0]["data"]
  pd.DataFrame({"prom_${barcode}":coverage}).to_csv(path_or_buf="prom_${barcode}.csv", index=False)
  """
}

process grid_1_json_to_csv {
  input:
    set barcode, library, json from grid_1_multiqc_out_ch
  output:
    file "grid_${barcode}.csv" into grid_1_coverage_ch

  """
  #!/usr/bin/env python

  import pandas as pd
  import os
  import json

  data = json.load(open("${json}"))
  coverage = data["report_plot_data"]["custom_data_linegraph"]["datasets"][0][0]["data"]
  pd.DataFrame({"grid_${barcode}":coverage}).to_csv(path_or_buf="grid_${barcode}.csv", index=False)
  """
}

process grid_2_json_to_csv {
  input:
    set barcode, library, json from grid_2_multiqc_out_ch
  output:
    file "grid_${barcode}.csv" into grid_2_coverage_ch

  """
  #!/usr/bin/env python

  import pandas as pd
  import os
  import json

  data = json.load(open("${json}"))
  coverage = data["report_plot_data"]["custom_data_linegraph"]["datasets"][0][0]["data"]
  pd.DataFrame({"grid_${barcode}":coverage}).to_csv(path_or_buf="grid_${barcode}.csv", index=False)
  """
}

process combine {
  publishDir 'nextflow_out/', mode: 'copy', overwrite: true
  input:
    file grid_1_csv_list from grid_1_coverage_ch.toList()
    file grid_2_csv_list from grid_2_coverage_ch.toList()
    file prom_csv_list from prom_coverage_ch.toList()
  output:
    file "grid-prom_combined.tsv" into out_ch

  """
  echo "" > grid-prom_combined.tsv
  seq 1 98 >> grid-prom_combined.tsv

  for filename in "${grid_1_csv_list} ${grid_2_csv_list} ${prom_csv_list}"; do paste grid-prom_combined.tsv \$filename > tmpout && mv tmpout grid-prom_combined.tsv; done
  """
}