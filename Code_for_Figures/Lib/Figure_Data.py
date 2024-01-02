import sys
sys.path.append("/home/averster/Documents/RossLab/Lib/")
from bioinformatics_parsing import parse_metagenome

import pandas as pd
import numpy as np
from pathlib import Path

from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()


days_to_months = 365 / 12

def load_hayden(type="Healthy", infile="/home/averster/Documents/RossLab/Data/Hayden_metaphlan_combined_s.xlsx", rel_ab=True):
    df_hayden = parse_metagenome.load_combined_metaphlan_file(Path(infile))
    if rel_ab:
        df_hayden = df_hayden.apply(lambda x: x / np.sum(x),0)
    df_meta_hayden = pd.read_excel("/home/averster/Documents/RossLab/Data/Hayden_metadata.xlsx", engine="openpyxl")
    if type == "Healthy":
        df_meta_hayden = df_meta_hayden.query('~Indv.str.contains("B")', engine="python")
    elif type == "CF":
        df_meta_hayden = df_meta_hayden.query('Indv.str.contains("B")', engine="python")
    else:
        raise Exception("Invalid type")

    # Processing Hayden
    df_model_hayden = df_hayden.T
    df_model_hayden["SRR"] = df_model_hayden.index
    df_model_hayden = df_model_hayden.merge(df_meta_hayden[["SRR","month","Indv"]].rename(columns={"month":"age","Indv":"Individual"}), on = "SRR", how = "inner")
    df_model_hayden.index = df_model_hayden["SRR"]
    df_model_hayden.index.name = None
    df_model_hayden = df_model_hayden.drop(columns = ["SRR"])
    return df_model_hayden


def load_kostic(infile="/home/averster/Documents/RossLab/Data/Kostic_metaphlan_combined_s_relative_abundance.xlsx", rel_ab=True):
    infile_meta = "/home/averster/Documents/RossLab/Data/Kostic_metadata.xlsx"
    df_meta_kostic = pd.read_excel(infile_meta, engine='openpyxl')

    df_kostic = parse_metagenome.load_combined_metaphlan_file(Path(infile))
    if rel_ab:
        df_kostic = df_kostic.apply(lambda x: x / np.sum(x),0)
    # Processing Kostic

    df_meta_kostic_sub = df_meta_kostic[["Gid_shotgun","Age_at_Collection","Subject_ID"]].rename(columns={"Age_at_Collection":"age", "Subject_ID":"Individual"})

    df_model_kostic = df_kostic.T
    df_model_kostic["Gid_shotgun"] = df_model_kostic.index
    df_model_kostic = df_model_kostic.merge(df_meta_kostic_sub, on = "Gid_shotgun")
    df_model_kostic.index = df_model_kostic["Gid_shotgun"]
    df_model_kostic = df_model_kostic.drop(columns=["Gid_shotgun"])

    df_model_kostic["age"] = df_model_kostic["age"] / days_to_months
    df_model_kostic["Individual"] = df_model_kostic["Individual"].astype("str")
    return df_model_kostic


def load_vatanen(infile="/home/averster/Documents/RossLab/Data/Vatanen_metaphlan_combined_s_relative_abundance.xlsx", rel_ab=True):
    infile_meta = "/home/averster/Documents/RossLab/Data/VatanenMetadata.csv"
    df_meta_vatanen = pd.read_csv(infile_meta)

    df_vatanen = parse_metagenome.load_combined_metaphlan_file(Path(infile))
    if rel_ab:
        df_vatanen = df_vatanen.apply(lambda x: x / np.sum(x),0)
    # Processing Vatanen
    df_model_vatanen = df_vatanen.T

    df_model_vatanen["gid_wgs"] = df_model_vatanen.index
    df_model_vatanen = df_model_vatanen.merge(df_meta_vatanen[["age_at_collection","gid_wgs","subjectID"]].dropna(), how = "left")
    df_model_vatanen.index = df_model_vatanen["gid_wgs"]
    df_model_vatanen = df_model_vatanen.drop(columns = ["gid_wgs"])

    df_model_vatanen["age"] = df_model_vatanen["age_at_collection"] / days_to_months
    df_model_vatanen = df_model_vatanen.drop(columns = ["age_at_collection"])
    df_model_vatanen = df_model_vatanen.rename(columns={"subjectID":"Individual"})
    df_model_vatanen["Individual"] = df_model_vatanen["Individual"].astype("str")
    return df_model_vatanen

def load_yassour(infile="/home/averster/Documents/RossLab/Data/Yassour_metaphlan_combined_s_relative_abundance.xlsx", rel_ab=True):
    df_yassour = parse_metagenome.load_combined_metaphlan_file(Path(infile))
    if rel_ab:
        df_yassour = df_yassour.apply(lambda x: x / np.sum(x),0)
    # Processing Yassour

    df_model_yassour = df_yassour.T
    df_model_yassour["age"] = [float(x.split("_")[1]) for x in df_model_yassour.index]
    df_model_yassour["Individual"] = [x.split("_")[0] for x in df_model_yassour.index]
    df_model_yassour = df_model_yassour.fillna(0)
    df_model_yassour["Individual"] = df_model_yassour["Individual"].astype("str")
    return df_model_yassour

def load_dartmouth(use_counts=False): 
    if use_counts:
        infile = "/home/averster/Documents/RossLab/Data/dart_CF_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx"
    else:
        infile = "/home/averster/Documents/RossLab/Data/dart_CF_metaphlan_combined_s_relative_abundance.xlsx"

    df_dart =  parse_metagenome.load_combined_metaphlan_file(Path(infile))
    if not use_counts:
        df_dart = df_dart.apply(lambda x: x / np.sum(x),0)

    infile_metadata = "/home/averster/Documents/RossLab/Data/CF_MetaData_v3.xlsx"

    df_meta_dart = pd.read_excel(infile_metadata, engine='openpyxl')

    df_model_dart = df_dart.T
    df_model_dart["Directory"] = df_model_dart.index
    df_model_dart = df_model_dart.merge(df_meta_dart[["Age_At_Collection_Date_Months","Directory", "Individual"]], how="left").rename(columns = {"Age_At_Collection_Date_Months":"age"})
    df_model_dart.index = df_model_dart["Directory"]
    df_model_dart = df_model_dart.drop(columns="Directory")
    df_model_dart = df_model_dart.query('age > 0') #Drops the two controls
    if use_counts:
        df_model_dart = df_model_dart.drop(columns=["k__Bacteria"])

    return df_model_dart


def load_diabimmune_data(use_rarified=False, use_counts=False):

    if use_counts:
        indir = Path("/home/averster/Documents/RossLab/Data/")
        df = load_vatanen(indir / "Vatanen_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx")
        ro.globalenv['df_model_vatanen'] = load_vatanen(indir / "Vatanen_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx", False)
        ro.globalenv['df_model_yassour'] = load_yassour(indir / "Yassour_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx", False)
        ro.globalenv['df_model_kostic'] = load_kostic(indir / "Kostic_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx", False)
    elif use_rarified:
        indir = Path("/home/averster/Documents/RossLab/Data/")
        ro.globalenv['df_model_vatanen'] = load_vatanen(indir / "VatanenRarified_metaphlan_combined_s.xlsx")
        ro.globalenv['df_model_yassour'] = load_yassour(indir / "YassourRarified_metaphlan_combined_s.xlsx")
        ro.globalenv['df_model_kostic'] = load_kostic(indir / "KosticRarified_metaphlan_combined_s.xlsx")
    else:
        ro.globalenv['df_model_vatanen'] = load_vatanen()
        ro.globalenv['df_model_yassour'] = load_yassour()
        ro.globalenv['df_model_kostic'] = load_kostic()

    ro.r('''

        library('tidyverse')
        source("/home/averster/Documents/RossLab/CF_Paper/CF_Metagenomics_Manuscript/av_code_for_figures/Lib/InfantAgeModels.R")
        df_list <- add_species_zeros(list(df_model_vatanen, df_model_yassour, df_model_kostic), FALSE)
        df_model_vatanen_full <- df_list[[1]]
        df_model_yassour_full <- df_list[[2]]
        df_model_kostic_full <- df_list[[3]]
        df_model_diabimmune <- rbind(df_model_vatanen_full, df_model_yassour_full)
        df_model_diabimmune <- rbind( df_model_diabimmune, df_model_kostic_full)
    ''')
    df_model_diabimmune = ro.globalenv['df_model_diabimmune']
    return df_model_diabimmune


def load_TEDDY(use_counts=False):
    if use_counts:
        infile = "/home/averster/Documents/RossLab/Data/TEDDY_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx"
    else:
        infile = "/home/averster/Documents/RossLab/Data/TEDDY_metaphlan_combined_s_relative_abundance.xlsx"
    df_teddy = parse_metagenome.load_combined_metaphlan_file(infile)
    
    indir = Path("/home/averster/Documents/RossLab/Data/")
    infile_meta = indir / "TEDDY_10_Percent_Rand_Selected_220113.txt"
    df_teddy_meta = pd.read_csv(infile_meta, sep=" ")

    if not use_counts:
        df_teddy = df_teddy.apply(lambda x: x / np.sum(x),0)

    df_model_teddy = df_teddy.T
    df_model_teddy["BioSample"] = df_model_teddy.index
    df_model_teddy = df_model_teddy.merge(df_teddy_meta[["BioSample","dbGaP_Subject_ID","collinterval"]].rename(columns={"collinterval":"age","dbGaP_Subject_ID":"Individual"}), on="BioSample", how="inner")
    df_model_teddy.index = df_model_teddy["BioSample"]
    df_model_teddy.index.name = None
    df_model_teddy = df_model_teddy.drop(columns = ["BioSample"])
    if use_counts:
        df_model_teddy = df_model_teddy.drop(columns=["k__Bacteria"])
    return df_model_teddy


if __name__ == "__main__":
    df_model_diabimmune = load_diabimmune_data()
    print(df_model_diabimmune)
