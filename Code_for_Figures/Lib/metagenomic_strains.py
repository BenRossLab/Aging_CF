
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter
from bioinformatics_parsing import parse_sequences, parse_blast, parse_phylogeny
import re
from scipy.stats import chi2_contingency, fisher_exact
import pickle

from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def load_fpmarker_counts(infile):
    df_Fp = pd.read_csv(infile)
    df_Fp = df_Fp.rename(columns={"Tag":"BioSample"})
    df_Fp.index = df_Fp["BioSample"]
    df_Fp = df_Fp.drop(columns=["BioSample"])
    return df_Fp

def load_clusters(infile):
    df_clusters = pd.read_csv(infile, sep=",")
    df_clusters = df_clusters.rename(columns={"samples":"BioSample"})
    return df_clusters


def all_clades(df_Fp):
    return set([re.search("clade([A-Z]{1})",x).group(1) for x in df_Fp.columns if re.search("clade([A-Z]{1})",x) ])

def identify_clades(df_Fp, clades_oi, min_reads = 25, min_frac_markers = 0.1):
    results = {}
    for clade in clades_oi:
        df_Fp_clade = df_Fp[[x for x in df_Fp.columns if "clade{}".format(clade) in x]]
        detected = df_Fp_clade.apply(lambda x: (np.sum(x >= min_reads) / len(x)) >= min_frac_markers,1)
        results[clade] = detected

    df_clade_identified = pd.DataFrame(results)
    return df_clade_identified

def summarize_clades_per_cluster(df_results, clades_oi):
    results = {}
    for (clade, dat) in df_results.groupby("cluster_ordered"):
        results[clade] = dat[clades_oi].apply(lambda x: np.sum(x) / len(x),0)
    df_plot = pd.DataFrame(results)
    return df_plot


def load_cluster_assignments(infile_assignments, infile_assignments_cf, infile_map):
    df_cluster_assignments = pd.read_csv(infile_assignments)
    df_cluster_assignments_cf = pd.read_csv(infile_assignments_cf)
    df_cluster_assignments_map = pd.read_csv(infile_map)

    df_cluster_assignments = df_cluster_assignments.merge(df_cluster_assignments_map)
    df_cluster_assignments_cf = df_cluster_assignments_cf.merge(df_cluster_assignments_map)
    return df_cluster_assignments, df_cluster_assignments_cf


def combine_data(df_Fp, df_model = None, df_assignments = None, clades_oi = ["A","B","C","D","E","F","G","H","I","J","K","L"], min_reads=10):
    
    df_results = identify_clades(df_Fp, clades_oi, min_reads=min_reads)
    df_plot = df_results[clades_oi]

    df_plot["n_strains"] = df_plot.loc[:,clades_oi[0]:clades_oi[-1]].sum(1)

    if df_model is not None:
        df_plot = df_plot.merge(df_model[["age","s__Faecalibacterium_prausnitzii"]], left_index=True, right_index=True)
    
    if df_assignments is not None:
        df_plot = df_plot.merge(df_assignments, left_index=True, right_on="samples")
    
    return df_plot


def summarize_clades(df_plot, group_vars = ["assignment","type"], clades_oi = ["A","B","C","D","E","F","G","H","I","J","K","L"], min_cluster = 10):
    
    results = []
    for (group, dat) in df_plot.query('n_strains > 0').groupby(group_vars):
        if dat.shape[0] < min_cluster:
            continue
        if len(group_vars) == 1:
            group_header = [group]
        else:
            group_header = list(group)
        results.append(group_header + dat[clades_oi].apply(lambda x: np.sum(x) / len(x),0).values.tolist())

    df_plot_summary = pd.DataFrame(results, columns = group_vars + clades_oi)
    df_plot_summary = df_plot_summary.rename(columns={"assignment":"category"})
    return df_plot_summary


def calculate_pvals(df_plot, group_var="assignment", clades_oi = ["A","B","C","D","E","F","G","H","I","J","K","L"]):
    lines = []
    for clade_oi in clades_oi:
        r = []
        for (i, dat) in df_plot.query('(n_strains > 0)').groupby(group_var):
            r.append([dat[clade_oi].sum(), dat.shape[0]])
        table = pd.DataFrame(r)
        if table.shape[0] > 2:
            pval = chi2_contingency(table)[1]
        else:
            pval = fisher_exact(table)[1]
        lines.append([clade_oi, pval])
        
    df_pvals = pd.DataFrame(lines, columns = ["clade","pval"])
    return df_pvals


def load_marker_count(infile_fp_markers):
    fp_markers = parse_sequences.load_seq_list_from_fasta(infile_fp_markers)
    clade_markers = []
    for seq in fp_markers:
        clade_markers.append(seq.id.split("_")[0])

    clade_marker_counts = Counter(clade_markers)
    return clade_marker_counts


def blast_for_markers(clade_marker_counts, infile_fp_markers, indir):

    df_prokka = parse_blast.find_all_prokka_blast(indir)
    df_prokka["name"] = df_prokka["name"].apply(parse_phylogeny.get_gca_from_ncbi_filename)

    lines = []
    for (infile_db, df_out) in parse_blast.blast_list_of_databases(infile_fp_markers, df_prokka["blastdb"].values):
        blast_counts = Counter([x.split("_")[0] for x in df_out["Query"].unique()])
        for key in blast_counts:
            lines.append([Path(infile_db).name,key, blast_counts[key] / clade_marker_counts[key]])

    df_blast = pd.DataFrame(lines, columns=["infile","clade","fraction"]).pivot(index="infile",columns="clade", values="fraction").fillna(0)
    return df_blast


def define_dicts(df_blast):
    color_dict = {
        "cladeA": "#fa7980",
        "cladeB": "#2a8000",
        "cladeC": "#d1c100",
        "cladeD": "#d53709",
        "cladeE": "#689d90",
        "cladeF": "#e6f74d",
        "cladeG": "#fe15b5",
        "cladeH": "#00e757",
        "cladeI": "#eba428",
        "cladeJ": "#1bd7ee",
        "cladeK": "#3e4fc7",
        "cladeL": "#805a40"}

    clade_dict = {}

    for (genome, dat) in df_blast.iterrows():
        clades = dat[dat>0.95].index
        if len(clades) == 1:
            clade_dict[parse_phylogeny.get_gca_from_ncbi_filename(genome)] = clades[0]
    return clade_dict, color_dict


def write_tree(df_blast, color_dict, clade_dict, outfile_1, outfile_2, infile = "/home/averster/Documents/RossLab/Genomes/Fp_core50_tree_full/RAxML_bestTree.Fp_core50_tree_full"):
    clade_dict, color_dict = define_dicts(df_blast)
    df_colors = pd.DataFrame({"clade":color_dict.keys(), "colors":color_dict.values()})
    df_markers = pd.DataFrame({"clade": clade_dict.keys(), "FpGroup":clade_dict.values()})

    ro.globalenv['df_colors'] = df_colors
    ro.globalenv['df_markers'] = df_markers
    ro.globalenv['infile'] = infile
    ro.globalenv['outfile_1'] = outfile_1
    ro.globalenv['outfile_2'] = outfile_2

    ro.r('''

        library('tidyverse')
        library('ggtree')
        library('treeio')
        library('ggrepel')
        library('phangorn')
        library('extrafont')


        t <- read.tree(infile)

        nodes_oi <- t$tip.label[str_detect(t$tip.label,"SRR") | str_detect(t$tip.label,"021721")]
        node <- MRCA(t, nodes_oi)
        node <- Ancestors(t, node, type="parent")
        t <- root(t, node=node)

        p2 <- ggtree(t)

        for (i in 1:nrow(df_colors)) {
            clade_oi <- df_colors$clade[i]
            color_oi <- df_colors$colors[i]
            nodes_oi <- (df_markers %>% filter(FpGroup == clade_oi) %>% pull(clade))
            if (length(nodes_oi) < 2)
                next
            node <- MRCA(t, nodes_oi)
            p2 <- collapse(p2, node, 'max',fill=color_oi, color="black", alpha=1.0)
            
        }

        for (i in 1:nrow(df_colors)) {

            clade_oi <- df_colors$clade[i]
            color_oi <- df_colors$colors[i]
            nodes_oi <- (df_markers %>% filter(FpGroup == clade_oi) %>% pull(clade))
            if (length(nodes_oi) >= 2)
                next
            
            nodes_oi_n <- which(t$tip.label%in%nodes_oi)
            p2 <- p2 + geom_point2(aes_string(subset=paste0("node%in%c(",paste(nodes_oi_n,collapse = ","),")")), shape=21, size=3, fill=color_oi)
        }

        nodes_oi <- t$tip.label[str_detect(t$tip.label,"SRR")]
        node <- MRCA(t, nodes_oi)
        p2 <- collapse(p2, node, 'max',fill="#805a40", color="black", alpha=1.0)

        p2$data$label[str_detect(p2$data$label,"GCA")] <- NA
        p2$data$label[str_detect(p2$data$label,"SRR")] <- "Cat"
        p2$data$label[str_detect(p2$data$label,"021721")] <- "CF"

        pdf(outfile_1, height=4.8, width=4)
        print(p2)    #+ geom_text_repel(aes(label=label), max.overlaps=50)
        dev.off()


        t <- read.tree(infile)

        nodes_oi <- t$tip.label[str_detect(t$tip.label,"SRR") | str_detect(t$tip.label,"021721")]
        node <- MRCA(t, nodes_oi)
        node <- Ancestors(t, node, type="parent")
        t <- tree_subset(t, node, levels_back=0)

        p2 <- ggtree(t)

        i <- 12
        clade_oi <- df_colors$clade[i]
        color_oi <- df_colors$colors[i]
        stopifnot(clade_oi == "cladeL")
        nodes_oi <- (df_markers %>% filter(FpGroup == clade_oi) %>% pull(clade))

        nodes_oi_n <- which(t$tip.label%in%nodes_oi)
        p2 <- p2 + geom_point2(aes_string(subset=paste0("node%in%c(",paste(nodes_oi_n,collapse = ","),")")), shape=21, size=8, fill=color_oi)


        p2$data$label[str_detect(p2$data$label,"GCA")] <- NA
        p2$data$label[str_detect(p2$data$label,"SRR")] <- "Cat"
        p2$data$label[str_detect(p2$data$label,"021721")] <- "CF"

        print("GGtree will fail to actually use Arial")
        print("In Inkscape, go Extensions->Font->Replace Font. Replace Helvetica with Arial.")
        pdf(outfile_2,  height=4.8, width=4, family="Arial")
        print(p2 + geom_tiplab(size=10/.pt))
        dev.off()


    ''')


if __name__ == "__main__":

    infile_fp_markers = "/home/averster/Documents/RossLab/Sequences/faecalibacterium_12clades_markers_FINAL_DeFilippis2020.fna"
    indir = Path("/home/averster/Documents/RossLab/Genomes/Fp_Genbank")
    clade_marker_counts = load_marker_count(infile_fp_markers)

    infile_blast_save = "df_blast_save.p"
    if Path(infile_blast_save).exists():
        df_blast = pickle.load(open(infile_blast_save,"rb"))
    else:
        df_blast = blast_for_markers(clade_marker_counts, infile_fp_markers, indir)
        pickle.dump(df_blast, open(infile_blast_save,"wb"))

    clade_dict, color_dict = define_dicts(df_blast)

    assert df_blast.query('cladeL > 0.95').shape[0] > 0, "We should be finding at least one cladeL genome"

    write_tree(df_blast, color_dict, clade_dict, "Figure_Fp_Tree_full.pdf", "Figure_Fp_Tree_zoom.pdf")