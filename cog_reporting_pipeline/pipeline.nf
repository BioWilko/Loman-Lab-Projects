Channel
  .from("\"prom\"", "\"grid\"", "\"grid\"")
  .toList()
  .set {prefix_ch}

Channel
  .from("\"BHAM-COVID19-20415_100_103_prom\"", "\"BHAM-COVID19-20415-100\"", "\"BHAM-COVID19-20415-103\"")
  .toList()
  .set {library_ch}

Channel
  .from((1..98).toList(), (1..48).toList(), (49..96).toList())
  .toList()
  .set {barcode_ch}

Channel
  .from("\"/data/homes/samw/projects/prom-vs-grid/100-103-barcoder_out/\"", "\"/data/homes/samw/projects/prom-vs-grid/100_fastq_combined/\"", 
    "\"/data/homes/samw/projects/prom-vs-grid/103_fastq_combined/\"")
  .toList()
  .set {fastq_dir_ch}

Channel
  .from("\"/cephfs/grid/bham/BHAM-COVID19-20415_100_103_prom/BHAM-COVID19-20415_100_103_prom/20201118_1419_2-E9-H9_PAE17046_c4370812/fast5_pass/\"", 
    "\"/data/homes/samw/projects/prom-vs-grid/100-103_fast5/\"", "\"/data/homes/samw/projects/prom-vs-grid/100-103_fast5/\"")
  .toList()
  .set{fast5_dir_ch}

Channel
  .from("\"/cephfs/grid/bham/BHAM-COVID19-20415_100_103_prom/BHAM-COVID19-20415_100_103_prom/20201118_1419_2-E9-H9_PAE17046_c4370812/sequencing_summary_PAE17046_23ee0e08.txt\"",
    "\"/data/homes/samw/projects/prom-vs-grid/100_seq_summary.txt\"", "\"/data/homes/samw/projects/prom-vs-grid/103_seq_summary.txt\"")
  .toList()
  .set {summary_ch}

scheme_dir = "$HOME"

reference_fasta = "/data/homes/samw/projects/variant_and_lineage/MN908947.3.fasta"

aln2type_yml_dir = "/data/homes/samw/aln2type/variant_definitions/variant_yaml/*yml"

out_dir = "~/projects/test_out/"

process generate_tuples {
  input:
    set prefixes from prefix_ch.toList()
    set libraries from library_ch.toList()
    set barcodes from barcode_ch.toList()
    set fastq_dirs from fastq_dir_ch.toList()
    set fast5_dirs from fast5_dir_ch.toList()
    set seq_summaries from summary_ch.toList()
  output:
    file "tuples.csv" into tuple_csv_ch

  """
  #!/usr/bin/env python

  import csv
  import os
  from glob import glob

  prefixes = ${prefixes}
  libraries = ${libraries}
  barcodes = ${barcodes}
  fastq_dirs = ${fastq_dirs}
  fast5_dirs = ${fast5_dirs}
  seq_summaries = ${seq_summaries}

  csv_out = [["prefix", "library", "barcode", "fastq_dir", "fast5_dir", "seq_summary"]]

  for index, library in enumerate(libraries):
    for barcode in barcodes[index]:
      if int(barcode) < 10:
        leading_barcode = str(0) + str(barcode)
      else:
        leading_barcode = barcode
      sample_tuple = []
      sample_tuple.extend([prefixes[index], library, leading_barcode])
      if len(glob(fastq_dirs[index] + "*/")) != 0:
        sample_tuple.append(fastq_dirs[index] + "barcode" + str(leading_barcode) + "/")
      else:
        sample_tuple.append(fastq_dirs[index])
      if len(glob(fast5_dirs[index] + "*/")) != 0:
        sample_tuple.append(fast5_dirs[index] + "barcode" + str(leading_barcode) + "/")
      else:
        sample_tuple.append(fast5_dirs[index])
      sample_tuple.append(seq_summaries[index])
      csv_out.append(sample_tuple)

  with open("tuples.csv", "w", newline="") as f:
      writer = csv.writer(f)
      writer.writerows(csv_out)
  """
}

tuple_csv_ch
  .splitCsv(header: true)
  .map(row -> tuple(row.prefix, row.library, row.barcode, row.fastq_dir, row.fast5_dir, row.seq_summary))
  .set {gplex_in_ch}

process guppyplex {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'
  input:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary from gplex_in_ch
  output:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("${prefix}_${barcode}_.fastq") into gplex_out_ch

  """
    artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ${fastq_dir} --prefix ${prefix}_${barcode}
  """
}

process minion {
  errorStrategy 'ignore'
  conda '/data/homes/samw/miniconda3/envs/artic/'

  input:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, fastq_file from gplex_out_ch
  output:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("multiqc_data/multiqc_data.json") into multiqc_out_ch
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("*.consensus.fasta") into pangolin_in_ch
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("*.consensus.fasta") into mafft_in_ch


  """
  artic minion --strict --normalise 200 --threads 4 --read-file ${fastq_file} --scheme-directory ${scheme_dir} --fast5-directory ${fast5_dir} --sequencing-summary ${seq_summary} nCoV-2019/V3 ${prefix}_${barcode}_out
  multiqc .
  """
}

process alignment {
  conda '/data/homes/samw/miniconda3/envs/mafft_env/'
  input:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, consensus from mafft_in_ch
  output:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("${prefix}_${barcode}_msa") into aln2type_in_ch

  """
  cat ${reference_fasta} ${consensus} > to_align.fasta
  mafft --auto to_align.fasta > ${prefix}_${barcode}_msa
  """
}

process variant_call {
  publishDir "${out_dir}", mode: "copy", overwrite: true 
  errorStrategy 'ignore'
  input:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, msa_file from aln2type_in_ch
  output:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("${prefix}_${barcode}_type.csv") into variant_out_ch

  """
  aln2type ./ ./ ${prefix}_${barcode}_type.csv MN908947.3 ${msa_file} ${aln2type_yml_dir} 
  """
}

process pango_lineage {
  publishDir '${out_dir}', mode: "copy", overwrite: true
  conda '/data/homes/samw/miniconda3/envs/pangolin/'
  input:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, consensus from pangolin_in_ch
  output:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("${prefix}_${barcode}_lineage.csv") into lineage_out_ch

  """
  pangolin ${consensus}
  mv lineage_report.csv ${prefix}_${barcode}_lineage.csv
  """
}

process json_to_csv {
  input:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, json from multiqc_out_ch
  output:
    set prefix, library, barcode, fastq_dir, fast5_dir, seq_summary, file("${prefix}_${barcode}.csv") into coverage_ch

  """
  #!/usr/bin/env python

  import pandas as pd
  import os
  import json

  data = json.load(open("${json}"))
  coverage = data["report_plot_data"]["custom_data_linegraph"]["datasets"][0][0]["data"]
  pd.DataFrame({"${prefix}_${barcode}":coverage}).to_csv(path_or_buf="${prefix}_${barcode}.csv", index=False)
  """
}

process coverage_combine {
  publishDir '${out_dir}', mode: 'copy', overwrite: true
  input:
    file csv_list from coverage_ch.toList()
  output:
    file "coverage_combined.tsv" into out_ch

  """
  echo "Amplicons" > coverage_combined.tsv
  seq 1 98 >> coverage_combined.tsv

  for filename in "${csv_list}"; do paste coverage_combined.tsv \$filename > tmpout && mv tmpout coverage_combined.tsv; done
  """
}