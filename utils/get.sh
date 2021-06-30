# this script requires
# pip install zenodo_get

DOI="10.5281/zenodo.4813082"

# define the env variable SOURCE_DATA in export.sh

source export.sh

echo $SOURCE_DATA

zenodo_get ${DOI} -o ${SOURCE_DATA}
