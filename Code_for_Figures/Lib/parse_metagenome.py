
import numpy as np
import pandas as pd
from pathlib import Path
import logging
import re
import tqdm

def load_single_metaphlan(infile, tax_level="s"):
    df_meta = pd.read_csv(infile, sep = "\t", header = 3)

    if df_meta.index[0] == 0:
        cols = list(df_meta.columns)
        cols[0] = cols[0].replace("#","")
        df_meta.columns = cols
    else:
        # Ocasionally the metaphlan file does not have "additional_species" and it loads the clade_name as the index
        # We need to correct that
        assert df_meta.index[0] == "k__Bacteria"
        df_meta.insert(0, "clade_name", df_meta.index)
        df_meta.columns = ["clade_name", "NCBI_tax_id", "relative_abundance", "additional_species"]

    locs_oi = [i for i, x in enumerate(df_meta["clade_name"]) if re.search("{}__[A-Za-z_]+$".format(tax_level), x)]
    expr_dict = {}
    for (i, dat) in df_meta.iloc[locs_oi,:].iterrows():
        expr_dict[dat["clade_name"]] = dat["relative_abundance"]
    return expr_dict


def load_multi_metaphlan(indir, tax_level="s", simply_name = False):
    indir = Path(indir)
    master_dict = {}
    for infile in tqdm.tqdm(indir.glob("*/*_metaphlan.txt"), total=len(list(indir.glob("*/*_metaphlan.txt")))):
        sample = infile.parts[-2]
        try:
            master_dict[sample] = load_single_metaphlan(infile, tax_level)
        except:
            logging.info("problem with {}".format(sample))
    df_metaphlan = pd.DataFrame(master_dict).fillna(0)
    if simply_name:
        df_metaphlan.index = [x.split("|")[-1] for x in df_metaphlan.index]
    return df_metaphlan


def load_library_size(infile_library_size):
    df_library_size = pd.read_csv(infile_library_size, sep = ",")
    df_library_size.index = df_library_size["name"]
    library_size = df_library_size["reads"].to_dict()
    return library_size

def abundance_length_library_normalized(infile_counts, library_size, size_hash):
        
    # create the abundance
    df_abundance = pd.read_csv(infile_counts)
    df_abundance.index = df_abundance["Tag"]
    df_abundance = df_abundance.drop(columns = ["Tag","total","skipped"])
    
    cols_drop = [] # columns that don't have a species identification
    
    library_size_vector = np.array([library_size[x] for x in df_abundance.index])
    for col in df_abundance:
        if col not in size_hash:
            cols_drop.append(col)
            continue
        df_abundance[col] = df_abundance[col] / size_hash[col] / library_size_vector
    df_abundance = df_abundance.drop(columns=cols_drop)
    return df_abundance


def hash_sequence_length(infile_seqs):
    size_hash = {}
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f, "fasta"):
            size_hash[seq.id] = len(seq)
    return size_hash


def gene_metagenome_abundance(infile_counts, infile_library_size, infile_seqs):
    size_hash = hash_sequence_length(infile_seqs)
    library_size = load_library_size(infile_library_size)
    df_abundance = abundance_length_library_normalized(infile_counts, library_size, size_hash)
    return df_abundance


def genus_metagenome_abundance(infile_counts, infile_library_size, infile_seqs = "/home/averster/Documents/RossLab/MetaphlanFiles/g__Bacteroides.fasta"):
    # Load up the seqs and record which apply to which species
    sp_marker_hash = defaultdict(list)
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f, "fasta"):
            m = re.search("s__([A-Za-z_]+)",seq.description)
            if not m:
                # This occurs with some weird markers for sub-species
                # Eg. 1263044__GCA_000433375.1_00023 for s__Bacteroides_faecis_CAG_32
                continue
            sp = m.group(1)
            sp_marker_hash[sp].append(seq.id)
    
    size_hash = hash_sequence_length(infile_seqs)
    library_size = load_library_size(infile_library_size)
    df_abundance = abundance_length_library_normalized(infile_counts, library_size, size_hash)
    
    # average the sp_markers
    # Drop absent markers
    df_abundance_sp = pd.DataFrame()
    for sp in sp_marker_hash:
        df_abundance_sp["{}_mean".format(sp)] = df_abundance[sp_marker_hash[sp]].apply(lambda x: np.mean(x), 1)
        df_abundance_sp["{}_sd".format(sp)] = df_abundance[sp_marker_hash[sp]].apply(lambda x: np.std(x), 1)

    df_abundance_sp = df_abundance_sp.fillna(0)
    
    return df_abundance_sp


def load_combined_metaphlan_file(infile_excel, simplify_taxa = True):
    df = pd.read_excel(infile_excel, engine='openpyxl')
    if simplify_taxa:
        df.index = [x.split("|")[-1] for x in df["clade"]]
    else:
        df.index = df["clade"]
    df = df.drop("clade", axis=1)
    df.columns = [re.sub("_metaphlan.+","",x) for x in df.columns]

    return df



def merge_two_metaphlan_tables(df, df2):
    return df.merge(df2, how="outer", left_index=True, right_index=True).fillna(0)
