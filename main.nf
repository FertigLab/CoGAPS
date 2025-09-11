process COGAPS {
  tag "$prefix"
  label 'process_medium'
  label 'process_high_memory'
  label 'process_long'
  container 'ghcr.io/fertiglab/cogaps:nextflow'

  input:
    tuple val(meta), path(dgCMatrix), val(cparams)

  output:
    tuple val(meta), path("${prefix}/cogapsResult.rds"), emit: cogapsResult
    path  "versions.yml",                                emit: versions

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
        nsets <- $cparams.nsets;
        allow_cpus <- as.numeric($task.cpus);
        if( allow_cpus < 2){
          stop("Error: distributed mode requires at least 2 cpus")
        }
        if (nsets > allow_cpus){
          message("Warning: nsets is greater than available cpus. Setting nsets to ", allow_cpus);
        } 
        params <- setDistributedParams(params, nSets = min(nsets,allow_cpus));
      };
      cogapsResult <- CoGAPS(data = data, params = params, nThreads = $cparams.nthreads,
                             outputFrequency = 100);
      saveRDS(cogapsResult, file = "${prefix}/cogapsResult.rds")'

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

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


  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"

  Rscript -e 'res <- Seurat::Read10X("$data/filtered_feature_bc_matrix/");
              saveRDS(res, file="${prefix}/dgCMatrix.rds")';

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

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
}

process COGAPS_ADATA2DGC {
  tag "$meta.id"
  label 'process_medium'
  container 'docker.io/satijalab/seurat:5.0.0'

  input:
      tuple val(meta), path(data) 
  output:
      tuple val(meta), path("${prefix}/dgCMatrix.rds"), emit: dgCMatrix
      path "versions.yml"                             , emit: versions

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
}

process COGAPS_PREPROCESS {
  tag "$prefix"
  label 'process_medium'
  container 'ghcr.io/fertiglab/cogaps:nextflow'

  input:
    tuple val(meta), path(dgCMatrix)

  output:
    tuple val(meta), path("${prefix}/dgCMatrix.rds"),    emit: dgCMatrix
    path  "versions.yml",                                emit: versions

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p "${prefix}"
  Rscript -e 'library("Matrix");
      library(sparseMatrixStats)
      sparse <- readRDS("$dgCMatrix");

      #drop rows with > 95% zero counts
      message("filtering rows with >95% zeros");
      nz <- rowSums(sparse != 0);
      sparse <- sparse[nz > 0.05 * ncol(sparse),];
      message("filtered to ", nrow(sparse), " columns of ", length(nz));

      #drop columns with > 95% zero counts
      message("filtering columns with >95% zeros");
      nz <- colSums(sparse != 0);
      sparse <- sparse[,nz > 0.05 * nrow(sparse)];
      message("filtered to ", ncol(sparse), " rows of ", length(nz));

      #select top N genes
      message("finding top ", ${params.n_top_genes}, " genes");
      vars <- rowVars(sparse);
      ngenes <- min(length(vars),${params.n_top_genes});
      top_genes <- order(vars, decreasing=TRUE)[1:ngenes];
      sparse <- sparse[top_genes,];
      message("selected top ", length(top_genes), " genes of ", length(vars));

      saveRDS(sparse, file = "${prefix}/dgCMatrix.rds")'

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  touch "${prefix}/dgCMatrix.rds"
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}


//example workflow
workflow {
  //example channel with data folders, for example
  ch_adata = Channel.fromPath("${params.input}/**.h5ad")
    .map { tuple([id:it.getName().replace('.', '-')], it)}

  ch_rds = Channel.fromPath("${params.input}/**.rds")
    .map { tuple([id:it.getName().replace('.', '-')], it)}

  //make a channel with desired pattern number
  def patterns = params.npatterns.split(',').collect { it.toInteger() }
  ch_patterns = Channel.from(patterns)

  //example channel with cparams
  ch_fixed_params = Channel.of([niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1])

  ch_cparams = ch_patterns
    .combine(ch_fixed_params)
    .map { tuple([id:it[0].toString(), npatterns:it[0], niterations:it[1].niterations, sparse:it[1].sparse, distributed:it[1].distributed, nsets:it[1].nsets, nthreads:it[1].nthreads]) }

  // convert adata to dgCMatrix
  COGAPS_ADATA2DGC(ch_adata)

  // preprocess dgCMatrix
  ch_preprocess = COGAPS_ADATA2DGC.out.dgCMatrix
    .map { tuple(it[0], it[1]) }

  ch_preprocess = ch_preprocess.mix(ch_rds)
  
  COGAPS_PREPROCESS(ch_preprocess)

  // ch_cogaps_input of converted adatas and rdses
  ch_input = COGAPS_PREPROCESS.out.dgCMatrix
    .map { tuple(it[0], it[1]) }

  // combine the two channels as input to CoGAPS
  ch_input = ch_input.combine(ch_cparams)

  COGAPS(ch_input)
}

//example:
//nextflow run main.nf --input tests/nextflow --outdir out -c nextflow.config -profile docker --max_memory 10GB --max_cpus 8

