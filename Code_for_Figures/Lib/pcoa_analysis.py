
import re
import pandas as pd

def load_species_tree(infile_tree = "/dartfs/rc/lab/R/RossB/SalernoP/metagenomics_data_backup/PCoA/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk"):
    # A species tree provided by the metaphlan repository
    # https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils
    
    with open(infile_tree, 'r') as f:
        tree_str = next(f).rstrip("\n")

    tree_str = re.sub("GCA_[0-9]+\|","",tree_str)
    return tree_str

def remove_missing_species(df_metaphlan, tree_str):
    # Remove species that are not in the tree and will screw up the beta diversity calculation
    
    missing = [x for x in df_metaphlan.index if x not in tree_str]
    print("total missing {}".format(len(missing)))
    return df_metaphlan.loc[[x not in missing for x in df_metaphlan.index],:]

def combine_dfs(df_list):
    # Merges a list of metaphlan dataframes
    # Fills in the missing species with zeros
    df_full = df_list[0]
    for i in range(1, len(df_list)):
        df_full = df_full.merge(df_list[i],left_index=True,right_index=True,how="outer").fillna(0)
    return df_full
