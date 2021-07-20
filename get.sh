# Script to retrieve the source data from zenodo
# $ pip install zenodo_get
# $ bash get.sh <your-folder>

DOI="10.5281/zenodo.4813082"

# define an env variable SOURCE_DATA in export.sh

PATH_SOURCE_DATA=$PWD/source-data

mkdir -p ${PATH_SOURCE_DATA}

zenodo_get ${DOI} -o ${PATH_SOURCE_DATA}
