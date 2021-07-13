import os
import matplotlib as mpl
import matplotlib.pyplot as plt

base_path = os.environ['PATH_SOURCE_DATA']

# output_boostdm
output_boostdm = os.path.join(base_path,"boostdm-output")

os.environ['INTOGEN_DATASETS'] = os.path.join(output_boostdm,"intogen")
intogen_data = os.environ['INTOGEN_DATASETS']

os.environ['DRIVERS_PATH'] = os.path.join(intogen_data,"drivers.tsv")
drivers_path = os.environ['DRIVERS_PATH']

os.environ['COHORTS_PATH'] = os.path.join(intogen_data,"cohorts.tsv")
cohorts_path = os.environ['COHORTS_PATH']

# PFAM info
PFAM_files = os.path.join(intogen_data,"pfam_biomart.tsv.gz")
PFAM_info = os.path.join(intogen_data,"pfam_names.info.csv")

# CDS coordinates
path_coord = os.path.join(intogen_data,"cds_biomart.tsv") 

# oncotree
oncotree_path = os.path.join(intogen_data, "oncotree", "tree_cancer_types.json")

# colors
dict_colors_role = {"Act": "#a6611a",
                    "LoF": "#018571",
                    "ambiguous": "#f5f5f5",#808080
                    "Amb": "#f5f5f5"}


# colors
dict_colors_role = {"Act": "#a6611a",
                    "LoF": "#018571",
                    "ambiguous": "#f5f5f5",#808080
                    "Amb": "#f5f5f5"}


degrons_path = os.path.join(intogen_data,"degrons","degron_instances.tsv")
uniprot_conv = os.path.join(intogen_data,"degrons","uniprot_transcript.tsv")


all_observed_mutations = os.path.join(output_boostdm,"discovery","mutations.tsv")

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
