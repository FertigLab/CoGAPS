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
      cogapsResult <- CoGAPS(data = data, params = params, nThreads = $cparams.nthreads,
                             outputFrequency = floor($cparams.niterations/10));
      saveRDS(cogapsResult, file = "${prefix}/cogapsResult.rds")'

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}

process COGAPS_TENX2DGC {
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

  Rscript -e 'res <- Seurat::Read10X("$data/filtered_feature_bc_matrix/");
              res <- Seurat::NormalizeData(res);
              saveRDS(res, file="${prefix}/dgCMatrix.rds")';

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}

process COGAPS_ADATA2DGC {
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
        hdf5r: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  Rscript -e 'message("Reading", "$data");
              f <- hdf5r::h5file(filename = "$data", mode="r");
              enctype <- hdf5r::h5attributes(f[["X/"]])[["encoding-type"]];

              i <- hdf5r::readDataSet(f[["X/indices"]]);
              p <- hdf5r::readDataSet(f[["X/indptr"]]);
              x <- hdf5r::readDataSet(f[["X/data"]]);
              var <- hdf5r::readDataSet(f[["var/_index"]]);
              obs <- hdf5r::readDataSet(f[["obs/_index"]]);

              message("Got", enctype, " ", length(var), "x", length(obs));

              if(enctype=="csr_matrix"){
                dimnames <- list(var, obs)
                transpose <- FALSE
              } else if (enctype=="csc_matrix"){
                dimnames <- list(obs, var)
                transpose <- TRUE
              } else {
                stop("Unknown encoding type")
              };
              message("Creating dgCMatrix");
              res <- Matrix::sparseMatrix(i=i, p=p, x=x, dims=lengths(dimnames),
                                          dimnames=dimnames, index1=FALSE, repr="C");

              if(transpose){
                res <- Matrix::t(res)
              }; 

              message("Normalizing data");
              res <- Seurat::NormalizeData(res);
              message("Saving dgCMatrix");
              saveRDS(res, file="${prefix}/dgCMatrix.rds")';

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hdf5r: \$(Rscript -e 'print(packageVersion("hdf5r"))' | awk '{print \$2}')
        Matrix: \$(Rscript -e 'print(packageVersion("Matrix"))' | awk '{print \$2}')
        Seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}