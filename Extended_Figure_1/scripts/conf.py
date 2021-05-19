import os
import matplotlib as mpl
import matplotlib.pyplot as plt

os.environ['INTOGEN_DATASETS'] = "/workspace/projects/boostdm/nature-release/boostdm-output/intogen"
intogen_data = os.environ['INTOGEN_DATASETS']

os.environ['DRIVERS_PATH'] = "/workspace/projects/boostdm/nature-release/boostdm-output/intogen/drivers.tsv"
drivers_path = os.environ['DRIVERS_PATH']

os.environ['COHORTS_PATH'] = "/workspace/projects/boostdm/nature-release/boostdm-output/intogen/cohorts.tsv"
cohorts_path = os.environ['COHORTS_PATH']

# PFAM info
PFAM_files = '/workspace/projects/boostdm/nature-release/boostdm-output/intogen/pfam_biomart.tsv.gz'
PFAM_info = '/workspace/projects/boostdm/nature-release/boostdm-output/intogen/pfam_names.info.csv'

# CDS coordinates
path_coord = '/workspace/projects/boostdm/nature-release/boostdm-output/intogen/cds_biomart.tsv'

# oncotree
oncotree_path = os.path.join(intogen_data, "oncotree", "tree_cancer_types.json")

# output_boostdm
output_boostdm = "/workspace/projects/boostdm/nature-release/boostdm-output/"

# colors
dict_colors_role = {"Act": "#a6611a",
                    "LoF": "#018571",
                    "ambiguous": "#f5f5f5",#808080
                    "Amb": "#f5f5f5"}
all_observed_mutations = "/workspace/projects/boostdm/nature-release/boostdm-output/discovery/mutations.tsv"

def config_params(font_size=10):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


config_params()
