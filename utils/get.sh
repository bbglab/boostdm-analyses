# source get.sh <target data folder>

DOI="10.5281/zenodo.4813082"

# define the env variable SOURCE_DATA in export.sh

SOURCE_DATA=$1

bash -c "mkdir -p ${SOURCE_DATA} && export SOURCE_DATA=${SOURCE_DATA}"

# source export.sh

echo "source data: ${SOURCE_DATA}"

# pip install zenodo_get

zenodo_get ${DOI} -o ${SOURCE_DATA}
