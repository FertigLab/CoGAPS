process COGAPS {
  tag "$prefix"
  label 'process_high'
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
      #select top 5K genes
      message("finding top ", ${params.n_top_genes}, " genes");
      vars <- apply(sparse, 1, var);
      top_genes <- order(vars, decreasing=TRUE)[1:${params.n_top_genes}];
      sparse <- sparse[top_genes,];
      message("selected top ", length(top_genes), " genes of ", length(vars));
   
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
  label 'process_medium'
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


//example workflow
workflow {
  //example channel with data folders, for example
  ch_adata = Channel.fromPath("${params.input}/**.h5ad")
    .map { tuple([id:it.getName().replace('.', '-')], it)}

  ch_rds = Channel.fromPath("${params.input}/**.rds")
    .map { tuple([id:it.getName().replace('.', '-')], it)}

  //example channel with cparams
  ch_cparams = Channel.of([npatterns: 5, niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1],
                          [npatterns: 10, niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1],
                          [npatterns: 15, niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1],
                          [npatterns: 20, niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1],
                         )

  // convert adata to dgCMatrix
  COGAPS_ADATA2DGC(ch_adata)

  // ch_cogaps_input of converted adatas and rdses
  ch_input = COGAPS_ADATA2DGC.out.dgCMatrix
  ch_input = ch_input.mix(ch_rds)

  // combine the two channels as input to CoGAPS
  ch_input = ch_input.combine(ch_cparams)

  COGAPS(ch_input)
}

//example:
//nextflow run main.nf --input tests/nextflow --outdir out -c nextflow.config -profile docker 
