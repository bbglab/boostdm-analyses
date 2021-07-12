# Script to retrieve the source data from zenodo
# $ pip install zenodo_get
# $ source get.sh <your-folder>

DOI="10.5281/zenodo.4813082"

# define an env variable SOURCE_DATA in export.sh

SOURCE_DATA=$1

bash -c "mkdir -p ${SOURCE_DATA} && export SOURCE_DATA=${SOURCE_DATA}"

zenodo_get ${DOI} -o ${SOURCE_DATA}
