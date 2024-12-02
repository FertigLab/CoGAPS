process COGAPS {
  tag "$prefix"
  label 'process_medium'
  label 'process_long'
  container 'ghcr.io/fertiglab/cogaps:master'

  input:
    tuple val(meta), path(dgCMatrix), val(cparams)

  output:
    tuple val(meta), path("${prefix}/cogapsResult.rds"), emit: cogapsResult
    path  "versions.yml",                                emit: versions

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}/${cparams.niterations}-${cparams.npatterns}-${cparams.sparse}-${cparams.distributed}"
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
  prefix = task.ext.prefix ?: "${meta.id}/${cparams.niterations}-${cparams.npatterns}-${cparams.sparse}-${cparams.distributed}"
  """
  mkdir -p "${prefix}"
  Rscript -e 'library("CoGAPS");
      sparse <- readRDS("$dgCMatrix");
      data <- as.matrix(sparse);
      #avoid errors with distributed params
      dist_param <- NULL;
      if(!("$cparams.distributed"=="null")){
        dist_param <- "$cparams.distributed"};
      params <- CogapsParams(seed=42,
                             nIterations = $cparams.niterations,
                             nPatterns = $cparams.npatterns,
                             sparseOptimization = as.logical($cparams.sparse),
                             distributed=dist_param);
      if (!(is.null(dist_param))){
        params <- setDistributedParams(params, nSets = $cparams.nsets);
      };
      cogapsResult <- CoGAPS(data = data, params = params, nThreads = $cparams.nthreads);
      saveRDS(cogapsResult, file = "${prefix}/cogapsResult.rds")'

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}

process PREPROCESS {
  tag "$meta.id"
  label 'process_low'
  container 'docker.io/satijalab/seurat:5.0.0'

  input:
      tuple val(meta), path(data) 
  output:
      tuple val(meta), path("${prefix}/dgCMatrix.rds"), emit: dgCMatrix
      path "versions.yml"                             , emit: versions

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"

  """
  mkdir "${prefix}"
  touch "${prefix}/dgCMatrix.rds"
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"

  Rscript -e 'res <- Seurat::Read10X("$data/raw_feature_bc_matrix/");
              res <- Seurat::NormalizeData(res);
              saveRDS(res, file="${prefix}/dgCMatrix.rds")';

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}