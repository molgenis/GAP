def getBarcodes() {
barcodes = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map{ tuple(it.SentrixBarcode_A, it.Project, it.manifest) }
    .unique()
}

process echoFinalReport {
  label 'small'
  executor "local"
  echo true

  input:
  tuple val(barcode), val(Project), val(manifestFile)

  output:
  path "${barcode}.txt", emit: reports

  script:
  """
  echo "${barcode}" > "${barcode}.txt"
  """
}

process gtcToFinalReport {
  label 'small'
  module = ["$params.beadArrayVersion", "$params.gapVersion"]
  input:
  tuple val(barcode), val(Project), val(manifestFile)

  output:
  path "${barcode}.txt", emit: finalreports

  script:
  """

  mkdir -p $params.outdir

  $baseDir/bin/gtc_final_report.py \
  --manifest "$params.manifestDir/${manifestFile}" \
  --samplesheet "$params.samplesheet" \
  --gtc_directory "$params.gtcDir/$barcode" \
  --output_file "${barcode}.txt"
  """
}

process  mergeFinalReports() {
  label 'large'
  echo false
  publishDir "$params.outdir/finalreport", mode: 'copy'

  input:
    path ( report_list ).collect()

  output:
    path "finalreport.txt" , emit: report

  script:
  """
  mkdir -p "$params.outdir/tmp"
  mv $report_list "$params.outdir/tmp"
  mergeFinalReports.sh "$params.outdir/tmp/" 'finalreport.txt'
  rm -r "$params.outdir/tmp"
  """
}

process finalReportToOptical {
  label 'large'
  echo false

  input:
  path "finalreport.txt"

  output:
  path ('chr_*')

  script:
  """

  GS_to_Opticall_includeControls.sh -c "$params.recalculate" -i 'finalreport.txt' -c 'no' -o ./
"""
}

process OptiCall {
  label 'large'
  errorStrategy 'ignore'
  module = ["$params.opticallVersion"]
  publishDir "$params.outdir/opticall", mode: 'copy'
  echo false

  input:
  path ( optical_chr )

  output:
  path ('chr_*.calls'), emit: calls
  path ('chr_*.probs'), emit: probs

  script:
  """
  echo "${optical_chr}"

  opticall \
  -in "${optical_chr}" \
  -out ./"${optical_chr}"
  """
}

process OptiCallToGenSample {
  label 'large'
  echo false
  publishDir "$params.outdir/oxford_gen_sample", mode: 'copy'

  input:
  path ( calls )
  path ( probs )

  output:
  path ('chr_*.gen'), emit: gen
  path ('chr_*.sample'), emit: sample

  script:
  """
  echo "${calls}"
  echo "${probs}"

  optiCall_to_OxfordFile.sh "${calls}" "${probs}" ./

"""
}
