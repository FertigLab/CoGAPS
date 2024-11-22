process COGAPS {
  tag "$meta.id"
  label 'process_medium'
  label 'process_long'
  container 'ghcr.io/fertiglab/cogaps:master'

  input:
    tuple val(meta), path(dgCMatrix)
  output:
    tuple val(meta), path("${prefix}/cogapsResult.rds"), emit: cogapsResult
    path  "versions.yml",                                emit: versions

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  touch "${prefix}/cogapsResult.rds"
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  Rscript -e 'library("CoGAPS");
      sparse <- readRDS("$dgCMatrix");
      data <- as.matrix(sparse);
      #avoid errors with distributed params
      dist_param <- NULL;
      if(!("$params.distributed"=="null")){
        dist_param <- "$params.distributed"};
      params <- CogapsParams(seed=42,
                             nIterations = $params.niterations,
                             nPatterns = $params.npatterns,
                             sparseOptimization = as.logical($params.sparse),
                             distributed=dist_param);
      if (!(is.null(dist_param))){
        params <- setDistributedParams(params, nSets = $params.nsets);
      };
      cogapsResult <- CoGAPS(data = data, params = params, nThreads = $params.nthreads);
      saveRDS(cogapsResult, file = "${prefix}/cogapsResult.rds")'

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}
