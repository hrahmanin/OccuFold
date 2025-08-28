nextflow.enable.dsl=2

// Stop if no region
if( !params.region ) error "Provide --region"

// Defaults
params.outdir = params.outdir ?: 'results'
params.peaks  = params.peaks  ?: null
params.motifs = params.motifs ?: null
params.use_predicted_occupancy = params.use_predicted_occupancy ?: true
// params.model_weights and params.ctcfpfm come from nextflow.config

workflow {
  PREPROCESS( Channel.value(params.region) )
  PREDICT( PREPROCESS.out )
}

process PREPROCESS {
  def safe = ((params.region as String).replaceAll('[_,]','')).replaceAll(/[^A-Za-z0-9._-]/,'_')

  tag safe
  publishDir "${params.outdir}/step1", mode: 'copy'

  input:
  val region

  output:
  path "${safe}.csv"

  script:
  def peaksOpt  = params.peaks  ? "--peaks ${params.peaks}"   : ""
  def motifsOpt = params.motifs ? "--motifs ${params.motifs}" : ""

  """
  mkdir -p .mplconfig
  export MPLCONFIGDIR="\$PWD/.mplconfig"

  if [ -n "${params.peaks}" ] && [ ! -r "${params.peaks}" ]; then
    echo "ERROR: --peaks not readable: ${params.peaks}" >&2
    exit 2
  fi

  REGION_CANON="\$(echo "${region}" | tr -d '_,')"

  python ${projectDir}/scripts/step1_region_snippet.py \
    --region "\${REGION_CANON}" \
    --outdir step1 \
    ${peaksOpt} ${motifsOpt}

  mv step1/*.csv "${safe}.csv"
  """
}

process PREDICT {
  publishDir "${params.outdir}/step2", mode: 'copy'

  input:
  path csv

  output:
  path "${file(csv).baseName}.occupancy.csv"

  script:
  def stem    = file(csv).baseName
  def useFlag = params.use_predicted_occupancy ? "--use-predicted-occupancy" : ""

  // ---- Resolve paths to ABSOLUTE under ${projectDir} if they’re relative
  def isAbs = { String p -> p != null && p.startsWith('/') }
  def weightsPath = params.model_weights ? (isAbs(params.model_weights) ? params.model_weights : "${projectDir}/${params.model_weights}") : null
  def pfmPath     = params.ctcfpfm      ? (isAbs(params.ctcfpfm)      ? params.ctcfpfm      : "${projectDir}/${params.ctcfpfm}")      : null

  def weightsOpt = weightsPath ? "--weights ${weightsPath}" : ""
  def pfmOpt     = pfmPath     ? "--ctcfpfm ${pfmPath}"     : ""
  def peaksOpt   = (!params.use_predicted_occupancy && params.peaks) ? "--peaks ${params.peaks}" : ""

  """
  mkdir -p .mplconfig
  export MPLCONFIGDIR="\$PWD/.mplconfig"

  # Fail fast if aux files aren’t visible inside the task
  [ -z "${weightsPath}" ] || [ -r "${weightsPath}" ] || { echo "Missing weights: ${weightsPath}" >&2; exit 2; }
  [ -z "${pfmPath}" ]     || [ -r "${pfmPath}" ]     || { echo "Missing PFM: ${pfmPath}" >&2; exit 2; }

  python ${projectDir}/scripts/step2_predict_occupancy.py \
    --input ${csv} \
    --outdir step2 \
    ${useFlag} ${weightsOpt} ${pfmOpt} ${peaksOpt}

  mv step2/*.occupancy.csv "${stem}.occupancy.csv"
  """
}
