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


degrons_path = "/workspace/projects/boostdm/nature-release/boostdm-output/intogen/degrons/degron_instances.tsv"
uniprot_conv = "/workspace/projects/boostdm/nature-release/boostdm-output/intogen/degrons/uniprot_transcript.tsv"


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

# naming conventions

features = ["shap_CLUSTL_SCORE","shap_CLUSTL_cat_1","shap_CLUSTL_cat_2","shap_HotMaps_cat_1",
            "shap_HotMaps_cat_2","shap_smRegions_cat_1","shap_smRegions_cat_2", "shap_csqn_type_synonymous","shap_csqn_type_missense","shap_csqn_type_nonsense",
            "shap_csqn_type_splicing","shap_PhyloP",    
            "shap_Acetylation","shap_Phosphorylation",
            "shap_Regulatory_Site","shap_Ubiquitination","shap_Methylation", 
           "shap_nmd"]

name_features = {"shap_nmd":"NMD","shap_CLUSTL_SCORE":"l. cluster score","shap_CLUSTL_cat_1":"l. cluster ttype","shap_CLUSTL_cat_2":"l. cluster","shap_HotMaps_cat_1":"3D cluster ttype",
                "shap_HotMaps_cat_2":"3D cluster","shap_smRegions_cat_1":"sig. domain ttype","shap_smRegions_cat_2":"sig. domain other ttype","shap_role_Act":"Act","shap_role_LoF":"LoF","shap_Acetylation":"Acetylation",
                "shap_O-GlcNAc":"O-GlcNAc","shap_Methylation":"Methylation","shap_Phosphorylation":"Phospho site","shap_Regulatory_Site":"Reg. site","shap_Sumoylation":"SUMO","shap_Ubiquitination":"ub site",
                "shap_PhyloP":"PhyloP","shap_csqn_type_synonymous":"synonymous","shap_csqn_type_missense":"missense","shap_csqn_type_stop_gained":"nonsense","shap_csqn_type_stop_lost":"stop lost",
                "shap_csqn_type_splicing":"splice variant","shap_csqn_type_nonsense":"nonsense mutation"}