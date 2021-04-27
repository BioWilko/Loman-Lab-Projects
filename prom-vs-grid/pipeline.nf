Channel
  .value("BHAM-COVID19-20415-163-164")
  .set { prom_lib_ch }

Channel
  .value("BHAM-COVID19-20415-163")
  .set { grid_163_lib_ch }

Channel
  .value("BHAM-COVID19-20415-164")
  .set { grid_164_lib_ch }

Channel
  .of('01', '02', '03', '04', '05', '06', '07', '08', '09', 10..15, 17..45, 49..91)
  .set { prom_barcode_ch }

Channel
  .of('01', '02', '03', '04', '05', '06', '07', '08', '09', 10..15, 17..45)
  .set { grid_163_barcode_ch }

Channel
  .of(49..91)
  .set { grid_164_barcode_ch }

reference = '/data/homes/samw/projects/variant_and_lineage/MN908947.3.fasta'

process guppyplex_prom {
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    val barcode from prom_barcode_ch
    val library from prom_lib_ch
  output:
    set barcode, library, file("prom_${barcode}_.fastq") into prom_gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ~/projects/prom-vs-grid/163-4_basecaller_out/barcode${barcode}/ --prefix prom_${barcode}
  """
}

process minion_prom {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set barcode, library, fastq_file from prom_gplex_out_ch
  output:
    set barcode, library, file("multiqc_data/multiqc_data.json") into prom_multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 4 --scheme-directory $HOME/projects/panther_ext_quality/dataset/primer-schemes --read-file ${fastq_file} --fast5-directory /data/BHAM-COVID19-20415-163-164/*/*/fast5_pass/ --sequencing-summary /data/BHAM-COVID19-20415-163-164/*/*/sequencing_summary* nCoV-2019/V3 prom_${barcode}_out
  multiqc .
  """
}

process guppyplex_163 {
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    val barcode from grid_163_barcode_ch
    val library from grid_163_lib_ch
  output:
    set barcode, library, file("163_${barcode}_.fastq") into grid_163_gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory /cephfs/grid/bham/${library}/*/*/fastq_pass/barcode${barcode}/ --prefix 163_${barcode}
  """
}

process minion_163 {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set barcode, library, fastq_file from grid_163_gplex_out_ch
  output:
    set barcode, library, file("multiqc_data/multiqc_data.json") into grid_163_multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 4 --scheme-directory $HOME/projects/panther_ext_quality/dataset/primer-schemes --read-file ${fastq_file} --fast5-directory /cephfs/grid/bham/${library}/*/*/fast5_pass/barcode${barcode}/ --sequencing-summary /cephfs/grid/bham/${library}/*/*/sequencing_summary* nCoV-2019/V3 163_${barcode}_out
  multiqc .
  """
}

process guppyplex_164 {
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    val barcode from grid_164_barcode_ch
    val library from grid_164_lib_ch
  output:
    set barcode, library, file("164_${barcode}_.fastq") into grid_164_gplex_out_ch

  """
  artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory /cephfs/grid/bham/${library}/*/*/fastq_pass/barcode${barcode}/ --prefix 164_${barcode}
  """
}

process minion_164 {
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set barcode, library, fastq_file from grid_164_gplex_out_ch
  output:
    set barcode, library, file("multiqc_data/multiqc_data.json") into grid_164_multiqc_out_ch

  """
  artic minion --strict --normalise 200 --threads 4 --scheme-directory $HOME/projects/panther_ext_quality/dataset/primer-schemes --read-file ${fastq_file} --fast5-directory /cephfs/grid/bham/${library}/*/*/fast5_pass/barcode${barcode}/ --sequencing-summary /cephfs/grid/bham/${library}/*/*/sequencing_summary* nCoV-2019/V3 164_${barcode}_out
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

process grid_163_json_to_csv {
  input:
    set barcode, library, json from grid_163_multiqc_out_ch
  output:
    file "grid_${barcode}.csv" into grid_163_coverage_ch

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

process grid_164_json_to_csv {
  input:
    set barcode, library, json from grid_164_multiqc_out_ch
  output:
    file "grid_${barcode}.csv" into grid_164_coverage_ch

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
    file grid_163_csv_list from grid_163_coverage_ch.toList()
    file grid_164_csv_list from grid_164_coverage_ch.toList()
    file prom_csv_list from prom_coverage_ch.toList()
  output:
    file "163-4_combined.tsv" into out_ch

  """
  echo "" > 163-4_combined.tsv
  seq 1 98 >> 163-4_combined.tsv

  for filename in "${grid_163_csv_list} ${grid_164_csv_list} ${prom_csv_list}"; do paste 163-4_combined.tsv \$filename > tmpout && mv tmpout 163-4_combined.tsv; done
  """
}
