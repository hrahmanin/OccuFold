nextflow.enable.dsl=2
if( !params.region ) throw new IllegalArgumentException("Provide --region")
params.outdir = params.outdir ?: "results"
params.peaks  = params.peaks  ?: null
params.motifs = params.motifs ?: null

workflow {
  PREPROCESS( Channel.value(params.region) )
}

process PREPROCESS {
  def safe = params.region.replaceAll('[^A-Za-z0-9._-]','_')
  tag "${safe}"
  publishDir "${params.outdir}/step1", mode: 'copy'
  container "$projectDir/occufold_tools.sif"
  env { PYTHONPATH = "${projectDir}:${System.getenv('PYTHONPATH') ?: ''}" }
  input:  val region
  output: path "${safe}.csv"
  script:
  def peaksOpt  = params.peaks  ? "--peaks ${params.peaks}"   : ""
  def motifsOpt = params.motifs ? "--motifs ${params.motifs}" : ""
  """
  python ${projectDir}/scripts/step1_region_snippet.py \
    --region "${region}" \
    --outdir step1 \
    ${peaksOpt} ${motifsOpt}
  mv step1/*.csv "${safe}.csv"
  """
}
