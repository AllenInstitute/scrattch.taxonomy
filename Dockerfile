FROM rocker/tidyverse:4.2

MAINTAINER Nelson Johansen nelson.johansen@alleninstitute.org

RUN export GITHUB_PAT=1000

## Would have liked to do these next 2 RUN under FROM:python3.8 but artifact passing wasn't immediatly clear. build-essential libxml2-dev python-is-python3
RUN apt-get update
RUN apt-get install -y wget python3-dev python3-venv python3-pip
RUN pip3 install anndata==0.8.0 numpy

RUN R -e 'install.packages("reticulate")'
RUN R -e 'install.packages("anndata", update=TRUE)'
RUN R -e 'install.packages("ontologyIndex")'  
RUN R -e 'install.packages("stringdist")' 
RUN R -e 'install.packages("R.utils")'

RUN R -e 'install.packages("BiocManager", update=FALSE)' 
RUN R -e 'BiocManager::install(c( "AnnotationDbi", "data.table", "GO.db", \
                                  "impute", "limma", "preprocessCore", "xml2", "rols"), dependencies=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c( "munsell", "rhdf5", "dplyr", \
                                  "optparse", "foreach", "doParallel", "futile.logger", \
                                  "ggplot2", "WGCNA"), dependencies=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c( "randomForest", "LaplacesDemon", "reshape2", \
                                  "feather", "future", "tibble", "dendextend", \
                                  "Matrix", "MatrixExtra"), dependencies=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c( "mgcv", "edgeR", "caret", \
                                  "ggbeeswarm", "pvclust", \
                                  "cowplot" ), dependencies=NA, update=TRUE)'
RUN R -e 'BiocManager::install(c("bigstatsr", "umap"), dependencies=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c("beachmat", "BiocNeighbors"), dependencies=NA, update=TRUE)' 
# RUN R -e 'BiocManager::install("scOntoMatch", update=FALSE)'  # Uncomment if new code fails. Installation of "ontologyIndex" SHOULD suffice.

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.5.tar.gz", repos=NULL, type="source")'
RUN R -e 'install.packages("https://cloud.r-project.org/src/contrib/profmem_0.6.0.tar.gz", repos=NULL, type="source")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.8-0.tar.gz", repos=NULL, type="source")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_1.1.0.tar.gz", repos=NULL, type="source")'

RUN R -e 'install.packages("reticulate")'
RUN R -e 'install.packages("arrow")'
RUN R -e 'install.packages("anndata", update=TRUE)'
RUN R -e 'install.packages("jsonlite", update=TRUE)'


## Seurat setup
RUN apt-get update && \
    apt-get install -y libgeos-dev libglpk-dev
RUN R -e 'install.packages("sp")'
RUN R -e 'install.packages("spam")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz", repos=NULL, type="source",)'
RUN R -e 'remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))' # Version 5 screws up Seurat mapping
RUN R -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))' # Version 5 screws up Seurat mapping

## Remote installs
RUN R -e 'install.packages("remotes", update=TRUE)'
RUN R -e 'remotes::install_github("krlmlr/bindrcpp")'
RUN R -e 'remotes::install_github("igraph/rigraph")'
RUN R -e 'remotes::install_github("i-cyto/Rphenograph")' 
RUN R -e 'remotes::install_github("PavlidisLab/patchSeqQC")'
RUN R -e 'remotes::install_github("cysouw/qlcMatrix")'

## Allen Institute R installs
RUN R -e 'remotes::install_github("AllenInstitute/CCN")'
RUN R -e 'remotes::install_github("AllenInstitute/patchseqtools")' 
RUN R -e 'remotes::install_github("AllenInstitute/tasic2016data")'
RUN R -e 'remotes::install_github("AllenInstitute/hodge2019data")'
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.io")' 
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.vis")' 
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.hicat")' 
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.bigcat")'
RUN R -e 'remotes::install_github("AllenInstitute/mfishtools")'

## set python virtual environment
ENV VIRTUAL_ENV=/pyenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# cell_type_mapper install from GitHub
RUN git clone -b update/uns/to/precomp/stats/params --single-branch https://github.com/AllenInstitute/cell_type_mapper.git
RUN pip install -r ./cell_type_mapper/requirements.txt
RUN pip install -e ./cell_type_mapper
RUN pip install anndata==0.8.0 numpy==1.26.4 

## scrattch-taxonomy install from local source
COPY scrattch.taxonomy_1.2.tar.gz ./scrattch.taxonomy_1.2.tar.gz
RUN R -e 'install.packages("scrattch.taxonomy_1.2.tar.gz", repos=NULL, type="source")'

## scrattch-mapping install from local source
COPY scrattch.mapping_1.2.tar.gz ./scrattch.mapping_1.2.tar.gz
RUN R -e 'install.packages("scrattch.mapping_1.2.tar.gz", repos=NULL, type="source")'

# ## scrattch-patchseq install from local source
COPY scrattch.patchseq_1.2.tar.gz ./scrattch.patchseq_1.2.tar.gz
RUN R -e 'install.packages("scrattch.patchseq_1.2.tar.gz", repos=NULL, type="source")'

## 
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_1.1.0.tar.gz", repos=NULL, type="source")'

## Clean up
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/downloaded_packages

## Strip binary installed libraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
RUN strip /usr/local/lib/R/site-library/*/libs/*.so

##