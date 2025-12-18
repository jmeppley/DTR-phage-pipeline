from pathlib import Path
import os
from glob import glob
import pandas as pd
import numpy as np
from snakemake.logging import logger

configfile: 'config.yml'
SCRIPT_DIR = srcdir('scripts')

###################################
# INPUT VARS                      #
###################################
SAMPLE        = config['sample'] 
STYPE         = config['stype']
K             = config['KMER_CALC']['k']
DATABASE_NAME = ["nr_euk"] # can also point to dbs like "progenomes" and "mar"
TAX_RANK      = ["0","1","2","3","4","5","6"]


###################################
# OUTPUT VARS                     #
###################################
OUTPUT_ROOT  = Path(config['output_root'])
VERSION      = config['version']
RACON_ROUNDS = config['RACON']['repeat']




###################################
# DTR FINDING FILES               #
###################################
DTR_READS_DIR                = OUTPUT_ROOT / SAMPLE / STYPE / VERSION / 'dtr_reads'
DTR_READS_PREFIX             = '{}/output'.format(str(DTR_READS_DIR))
DTR_READS_TMP_DIR            = DTR_READS_DIR / 'aln_tmp'
DTR_READS_FASTA              = '{}.dtr.fasta'.format(str(DTR_READS_PREFIX))
DTR_READS_STATS              = '{}.dtr.stats.tsv'.format(str(DTR_READS_PREFIX))
DTR_READS_HIST               = '{}.dtr.hist.png'.format(str(DTR_READS_PREFIX))

FILTERED_FASTA = DTR_READS_FASTA

###################################
# K-mer UAMP files and plots      #
###################################
KMER_BIN_ROOT                     = OUTPUT_ROOT /  SAMPLE / STYPE / VERSION / 'kmer_binning'
KMER_BINS_MEMBERSHIP              = KMER_BIN_ROOT / 'bin_membership.tsv'
KMER_BINS_LIST                    = KMER_BIN_ROOT / 'bin_list.txt'
KMER_BIN_STATS                    = KMER_BIN_ROOT / 'bin_stats.csv'
KMER_FREQS_TMP                    = KMER_BIN_ROOT / 'kmer_comp.tmp'
KMER_FREQS                        = KMER_BIN_ROOT / 'kmer_comp.tsv'
KMER_FREQS_UMAP                   = KMER_BIN_ROOT / 'kmer_comp.umap.tsv'
KMER_FREQS_UMAP_TAX               = KMER_BIN_ROOT / 'kmer_comp.umap.{database}.{rank}.png'
KMER_FREQS_UMAP_QSCORE            = KMER_BIN_ROOT / 'kmer_comp.umap.qscore.png'
KMER_FREQS_UMAP_GC                = KMER_BIN_ROOT / 'kmer_comp.umap.gc.png'
KMER_FREQS_UMAP_READLENGTH        = KMER_BIN_ROOT / 'kmer_comp.umap.readlength.png'
KMER_FREQS_UMAP_BINS_PLOT         = KMER_BIN_ROOT / 'kmer_comp.umap.bins.png'
KMER_FREQS_UMAP_BINS_COORDS       = KMER_BIN_ROOT / 'kmer_comp.umap.bins.tsv'


###################################
# BIN ANALYSIS FILES              #
###################################
BINS_ROOT                   = KMER_BIN_ROOT / 'bins'
BIN_DIR                     = BINS_ROOT / '{bin_id}'
BIN_READLIST                = BIN_DIR / 'read_list.txt'
BIN_FASTA                   = BIN_DIR / '{bin_id}.reads.fa'
BIN_GENOME_SIZE             = BIN_DIR / 'genome_size.txt'
SKIP_BINS                   = config['SKIP_BINS'][SAMPLE][STYPE][VERSION]
BINNED_ANALYSIS_ROOT        = KMER_BIN_ROOT / "refine_bins"

#####################################
# Alignment clustering  FILES       #
#####################################
ALN_CLUST_DIR            = BINNED_ANALYSIS_ROOT / 'alignments'
ALN_CLUST_PAF            = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.ava.paf'
ALN_CLUST_OUTPUT_PREFIX  = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust'
ALN_CLUST_OUTPUT_HEATMAP = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust.heatmap.png'
ALN_CLUST_OUTPUT_INFO    = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust.info.csv'
ALN_CLUST_READS_COMBO    = ALN_CLUST_DIR / 'all_bins.clust.info.csv'


###################################
# Separate cluster-specific reads #
###################################
BIN_CLUSTER_ROOT                = BINNED_ANALYSIS_ROOT / 'align_cluster_reads'
BIN_CLUSTER_DIR                 = BIN_CLUSTER_ROOT / '{bin_clust_id}'
BIN_CLUSTER_READS_INFO          = BIN_CLUSTER_DIR / '{bin_clust_id}.readinfo.csv'
BIN_CLUSTER_READS_LIST          = BIN_CLUSTER_DIR / 'readlist.csv'
BIN_CLUSTER_READS_FASTA         = BIN_CLUSTER_DIR / '{bin_clust_id}.reads.fasta'
BIN_CLUSTER_REF_READ_LIST       = BIN_CLUSTER_DIR / '{bin_clust_id}.ref_readlist.csv'
BIN_CLUSTER_POL_READS_LIST      = BIN_CLUSTER_DIR / '{bin_clust_id}.pol_readlist.csv'
BIN_CLUSTER_REF_READ_FASTA      = BIN_CLUSTER_DIR / '{bin_clust_id}.ref_read.fa'
BIN_CLUSTER_POL_READS_FASTQ     = BIN_CLUSTER_DIR / '{bin_clust_id}.pol_reads.fq'
BIN_CLUSTER_POL_READS_FASTA     = BIN_CLUSTER_DIR / '{bin_clust_id}.pol_reads.fa'


####################################
# Polish cluster reads using Racon #
####################################
POLISH_DIR                       = BINNED_ANALYSIS_ROOT / 'align_cluster_polishing'
RACON_DIR                        = POLISH_DIR / 'racon'
BIN_CLUSTER_RACON_POLISHED_FASTA = RACON_DIR / '{bin_clust_id}' / '{{bin_clust_id}}.ref_read.racon_{repeats}x.fasta'.format(repeats=RACON_ROUNDS)


#####################################
# Polish cluster reads using Medaka #
#####################################
MEDAKA_DIR                                    = POLISH_DIR / 'medaka'
BIN_CLUSTER_POLISHED_REF_TMP                  = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.tmp.fasta'
BIN_CLUSTER_POLISHED_REF                      = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.fasta'
BIN_CLUSTER_POLISHED_POL_VS_REF_PAF           = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.paf'
BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS       = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.strands.summary.tsv'
BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.strands.reads.tsv'
BIN_CLUSTER_POLISHED_REF_PRODIGAL             = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.fasta'
BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT         = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.txt'
BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS       = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.stats.txt'

#########################################
# Check for fixed or cyclic permutation #
#########################################
DTR_ALIGN_PREFIX                 = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read')
DTR_ALIGN_COORD_PLOT             = '{}.dtr.aligns.png'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_TSV                    = '{}.dtr.aligns.tsv'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_CYC_PERM_TSV           = str(POLISH_DIR / 'polished.cyclic_permut.stats.tsv')

######################################
# Combine Medaka polished references #
######################################
ALL_POL_PREFIX                   = str(POLISH_DIR / 'polished')
ALL_POL_UNTRIMMED                = '{}.untrimmed.fasta'.format(ALL_POL_PREFIX)
ALL_POL_CDS_SUMMARY              = '{}.cds.summary.tsv'.format(ALL_POL_PREFIX)
ALL_POL_STRANDS                  = '{}.pol_strands.tsv'.format(ALL_POL_PREFIX)
ALL_POL_STRAND_ANNOTS            = '{}.pol_strands.reads.tsv'.format(ALL_POL_PREFIX)
ALL_POL                          = '{}.seqs.fasta'.format(ALL_POL_PREFIX)
ALL_POL_UNIQ                     = '{}.seqs.unique.fasta'.format(ALL_POL_PREFIX)
ALL_POL_DTR_STATS                = '{}.dtr.stats.tsv'.format(ALL_POL_PREFIX)
ALL_POL_DTR_GC_STATS             = '{}.dtr.gc.tsv'.format(ALL_POL_PREFIX)
ALL_POL_NUCMER_PREFIX            = '{}.nuc'.format(ALL_POL_PREFIX)
ALL_POL_NUCMER_DELTA             = '{}.nuc.delta'.format(ALL_POL_PREFIX)
ALL_POL_NUCMER_COORDS            = '{}.nuc.coords'.format(ALL_POL_PREFIX)
ALL_POL_STATS                    = '{}.stats.tsv'.format(ALL_POL_PREFIX)
ALL_POL_STATS_UNIQ               = '{}.stats.unique.tsv'.format(ALL_POL_PREFIX)
ALL_POL_CDS_PLOT_UNIQ_ALL        = '{}.unique.cds.all.png'.format(ALL_POL_PREFIX)
ALL_POL_CDS_PLOT_UNIQ_DTR_NPOL10 = '{}.unique.cds.dtr_npol10.png'.format(ALL_POL_PREFIX)

wildcard_constraints:
    rank = '\d+',


##### load rules #####

include: 'rules/align_clusters.smk'
include: 'rules/polish.smk'
include: 'rules/qc_genomes.smk'
include: 'rules/annotate.smk'
include: 'rules/dedup.smk'
include: 'rules/dtr_align.smk'

#############################################
# Required steps for assembly-free analysis #
#############################################

# 1. all_kmer_count_and_bin
# 2. all_kaiju              (optional for taxonomic annotation of UMAP plots)
# 3. all_populate_kmer_bins
# 4. all_alignment_clusters
# 5. all_polish_and_annotate
# 6. all_combine_dedup_summarize
# 7. all_linear_concatemer_reads

# minimal kmer bins

rule all:
    input:
        SUMMARY_PLOT,
        KMER_FREQS_UMAP_QSCORE,
        KMER_FREQS_UMAP_GC,
        KMER_FREQS_UMAP_READLENGTH,
        KMER_FREQS_UMAP_BINS_PLOT,
        KAIJU_RESULTS_KRONA_HTML,
        expand(str(KMER_FREQS_UMAP_TAX), database=DATABASE_NAME, rank=TAX_RANK),
        KMER_BIN_STATS,
        lambda w: expand_template_from_bins(w, ALN_CLUST_OUTPUT_HEATMAP),
        ALN_CLUST_READS_COMBO,
        lambda w: expand_template_from_bin_clusters(w, BIN_CLUSTER_REF_READ_FASTA),
        lambda w: expand_template_from_bin_clusters(w, BIN_CLUSTER_POL_READS_FASTA),
        lambda w: expand_template_from_bin_clusters(w, DTR_ALIGN_COORD_PLOT),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS),
        ALL_POL_CDS_PLOT_UNIQ_ALL,
        ALL_POL_CDS_PLOT_UNIQ_DTR_NPOL10,
        ALL_POL,
        ALL_POL_UNIQ,
        ALL_POL_STATS


rule call_umap_freq_map_bins:
    """ Determines the bins """
    input: KMER_FREQS_UMAP
    output: KMER_BINS_MEMBERSHIP
    conda: '../envs/umap.yml'
    params:
        min_cluster = config['UMAP']['bin_min_reads'],
    shell:
        'python {SCRIPT_DIR}/run_hdbscan.py -o {output} -c {params.min_cluster} {input}'

checkpoint generate_bins_and_stats:
    """ 
    collects all bins stats into KMER_BIN_STATS,
    and, for each bin_id, generates in BINS_ROOT:
     * BIN_READLIST
     * BIN_GENOME_SIZE
    """
    input: 
        bin_membership=KMER_BINS_MEMBERSHIP,
    output: 
        all_stats=KMER_BIN_STATS,
        bins_dir=directory(BINS_ROOT)
    run:
        df_reads = pd.read_csv(input.bin_membership, sep='\t')

        # aggregate read stats into bin stats
        df_bins = df_reads.groupby('bin_id').agg(bases=pd.NamedAgg('length', sum),
                                                 genomesize=pd.NamedAgg('length', np.mean),
                                                 coverage=pd.NamedAgg('length', len),
                                                 readlist=pd.NamedAgg('read', set),
                                                )
        df_bins['rl_gs_est'] = True

        # write bin stats (without readlist) to stats file
        df_bins[['bases','genomesize','coverage','rl_gs_est']] \
               .sort_index() \
               .to_csv(output.all_stats, sep='\t')

        # write bin specific files
        for bin_id, row in df_bins.iterrows():
            #
            # create bin specific subdir of BINS_ROOT
            bin_dir = str(BIN_DIR).format(bin_id=bin_id)
            os.makedirs(bin_dir, exist_ok=True)
            #
            # save genome size to file
            gsize_fn = str(BIN_GENOME_SIZE).format(bin_id=bin_id)
            with open(gsize_fn, 'w') as fh:
                fh.write(f'{row["genomesize"]}\n')
            #
            # save read list to file
            rlist_fn = str(BIN_READLIST).format(bin_id=bin_id)
            with open(rlist_fn, 'w') as fh:
                fh.write('{}\n'.format('\n'.join(row['readlist'])))

def expand_template_from_bins(wildcards, template):
    # get dir through checkpoints to throw Exception if checkpoint is pending
    checkpoint_dir = checkpoints.generate_bins_and_stats.get(**wildcards).output
    # get bins from files
    bins, = glob_wildcards(BIN_READLIST)
    # skips from config
    bins = [b for b in bins if b not in SKIP_BINS]
    # expand template
    return expand(str(template), bin_id=bins)

rule kmer_binned_fasta:
    input:
        read_list=BIN_READLIST,
        reads_fasta=FILTERED_FASTA,
    output: BIN_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.read_list} -o {output} {input.reads_fasta}'
