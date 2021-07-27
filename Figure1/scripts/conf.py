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
PFAM_files = os.path.join(intogen_data, "pfam_biomart.tsv.gz")
PFAM_info = os.path.join(intogen_data, "pfam_names.info.csv")

# CDS coordinates
path_coord = os.path.join(intogen_data, "cds_biomart.tsv")
# oncotree
oncotree_path = os.path.join(intogen_data, "oncotree", "tree_cancer_types.json")

# colors
dict_colors_role = {"Act": "#a6611a",
                    "LoF": "#018571",
                    "ambiguous": "#f5f5f5",
                    "Amb": "#f5f5f5"}

degrons_path = os.path.join(intogen_data,"degrons", "degron_instances.tsv")
uniprot_conv = os.path.join(intogen_data,"degrons", "uniprot_transcript.tsv")
all_observed_mutations = os.path.join(output_boostdm, "discovery", "mutations.tsv") # mutations.tsv

# naming conventions
features = ["shap_CLUSTL_SCORE","shap_CLUSTL_cat_1","shap_CLUSTL_cat_2","shap_HotMaps_cat_1",
            "shap_HotMaps_cat_2","shap_smRegions_cat_1","shap_smRegions_cat_2", "shap_csqn_type_synonymous",
            "shap_csqn_type_missense","shap_csqn_type_nonsense",
            "shap_csqn_type_splicing","shap_PhyloP",    
            "shap_Acetylation","shap_Phosphorylation",
            "shap_Regulatory_Site","shap_Ubiquitination","shap_Methylation", 
            "shap_nmd"]
name_features = {"shap_nmd":"NMD","shap_CLUSTL_SCORE":"l. cluster score","shap_CLUSTL_cat_1":"l. cluster ttype",
                 "shap_CLUSTL_cat_2":"l. cluster","shap_HotMaps_cat_1":"3D cluster ttype",
                 "shap_HotMaps_cat_2":"3D cluster","shap_smRegions_cat_1":"sig. domain ttype",
                 "shap_smRegions_cat_2":"sig. domain other ttype","shap_role_Act":"Act","shap_role_LoF":"LoF",
                 "shap_Acetylation":"Acetylation","shap_O-GlcNAc":"O-GlcNAc","shap_Methylation":"Methylation",
                 "shap_Phosphorylation":"Phospho site","shap_Regulatory_Site":"Reg. site","shap_Sumoylation":"SUMO",
                 "shap_Ubiquitination":"ub site","shap_PhyloP":"PhyloP","shap_csqn_type_synonymous":"synonymous",
                 "shap_csqn_type_missense":"missense","shap_csqn_type_stop_gained":"nonsense",
                 "shap_csqn_type_stop_lost":"stop lost","shap_csqn_type_splicing":"splice variant",
                 "shap_csqn_type_nonsense":"nonsense mutation"}

# Selected tumor types for extended figure 8 and 9
selecte_ttypes=['LGG', 'RCCC', 'PAAD', 'BLCA', 'OS', 'SCLC', 'HNSC', 'UCEC', 'HGG', 'OV', 'COREAD', 'BRCA', 'UCS',
                'PRAD', 'HC', 'GBM', 'CESC', 'LUAD', 'LUSC', 'ST', 'ESCA', 'UVM', 'DLBCL', 'CH', 'AML', 'MBL', 'CLL',
                'PIA', 'MM', 'NB', 'ALL', 'CM', 'SBCC', 'NHLY', 'THYM', 'THCA', 'RPC', 'S', 'ACC', 'RCH', 'PCPG',
                'PAIS', 'NSPH', 'MESO', 'LUNG_CANCER', 'LEIS', 'EWS', 'MGCT']
# path samples
samples_intogen = os.path.join(output_boostdm,"discovery","samples_anonymized.json")
dictionary_dndscv = os.path.join(base_path,"dndscv","dictionary_excess_dndscv.json")
# path clinvar and oncokb
clinvar_muts = os.path.join(base_path,"external","all_clinvar_hg38_with_origin.toboostdm.tsv.gz")
oncokb_output = os.path.join(base_path,"external","output_oncokb.txt")


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
