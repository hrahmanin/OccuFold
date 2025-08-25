nextflow.enable.dsl=2

// Hard-stop if no region
if( !params.region ) error "Provide --region"

// Defaults (kept here so the file is self-contained)
params.outdir = params.outdir ?: 'results'
params.peaks  = params.peaks  ?: null
params.motifs = params.motifs ?: null

workflow {
  PREPROCESS( Channel.value(params.region) )
}

process PREPROCESS {
  // Safe name for outputs/tags
  def safe = (params.region as String).replaceAll(/[^A-Za-z0-9._-]/,'_')

  tag safe
  publishDir "${params.outdir}/step1", mode: 'copy'

  input:
  val region

  // Collect the CSV we write
  output:
  path "${safe}.csv"

  script:
  // Optional CLI args
  def peaksOpt  = params.peaks  ? "--peaks ${params.peaks}"   : ""
  def motifsOpt = params.motifs ? "--motifs ${params.motifs}" : ""

  """
  # Writable cache for Matplotlib (avoids home/.config)
  mkdir -p .mplconfig
  export MPLCONFIGDIR="\$PWD/.mplconfig"

  # If peaks specified, ensure readable (fail fast with clear message)
  if [ -n "${params.peaks}" ] && [ ! -r "${params.peaks}" ]; then
    echo "ERROR: --peaks not readable: ${params.peaks}" >&2
    exit 2
  fi

  # Convert region like 'chr1:1_000_000-1_050_000' -> 'chr1:1000000-1050000'
  REGION_CANON="\$(echo "${region}" | tr -d '_,')"

  # Run the host script (repo is bind-mounted and PYTHONPATH is set by config)
  python ${projectDir}/scripts/step1_region_snippet.py \
    --region "\${REGION_CANON}" \
    --outdir step1 \
    ${peaksOpt} ${motifsOpt}

  # Standardize output name
  mv step1/*.csv "${safe}.csv"
  """
}
