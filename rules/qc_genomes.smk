############################
# QC and summarize genomes #
############################

rule combine_all_polished_ref_reads:
    input: 
        dirname=MEDAKA_DIR,
    output: temp(ALL_POL_UNTRIMMED)
    shell:
        'cat {input.dirname}/*/*.medaka.fasta > {output}'

rule remove_adapters_from_ref_reads:
    input: ALL_POL_UNTRIMMED
    output: ALL_POL
    params:
        switches=config['PORECHOP']['switches']
    threads: config['max_threads']
    conda: '../envs/porechop.yml'
    shell:
        'porechop {params.switches} -t {threads} -i {input} -v 1 -o {output}'

rule find_all_DTR_genomes:
    input: ALL_POL
    output:
        stats = temp(ALL_POL_DTR_STATS),
    conda: '../envs/dtr.yml'
    params:
        prefix = '{}'.format(ALL_POL_PREFIX),
        ovlp   = config['DTR_FIND']['ovlp'],
        chunksize = config['DTR_FIND']['chunksize']
    threads: config['max_threads']
    shell:
        'python {SCRIPT_DIR}/find_dtr_all_seqs.py --no_fasta --no_hist -t {threads} '
        '-p {params.prefix} -o {params.ovlp} -c {params.chunksize} {input}'

rule aggregate_prodigal_statistics:
    output:
        table=temp(ALL_POL_CDS_SUMMARY),
    params:
        pol_dir=MEDAKA_DIR
    conda: '../envs/clustering.yml'
    shell:
        'python {SCRIPT_DIR}/combine_cds_summary.py -o {output.table} {params.pol_dir}'

rule build_bin_cluster_summary_table:
    input:
        dtr_table=ALL_POL_DTR_STATS,
        cds_table=ALL_POL_CDS_SUMMARY
    output: temp(ALL_POL_DTR_GC_STATS)
    params:
        bins_dir=BINS_DIR,
        clust_dir=BIN_CLUSTER_DIR,
        prefix='{}{}'.format(SAMPLE,STYPE)
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/calculate_genomes_gc_contents.py -p {params.prefix} '
        '-o {output} {input.dtr_table} {input.cds_table} {params.bins_dir} {params.clust_dir}'

rule combine_bin_cluster_strand_counts_into_table:
    input: 
        pol_dir=MEDAKA_DIR
    output:
        counts=temp(ALL_POL_STRANDS),
        annots=temp(ALL_POL_STRAND_ANNOTS)
    run:
        import os
        import pandas as pd
        from glob import glob
        count_fns = glob(os.path.join(input.pol_dir, '*', '*.ref_read.strands.summary.tsv'))
        count_dfs = [pd.read_csv(fn, sep='\t') for fn in count_fns]
        count_df  = pd.concat(count_dfs)
        count_df['bin']     = count_df['clust_id'].map(lambda x: int(x.split('_')[0]))
        count_df['cluster'] = count_df['clust_id'].map(lambda x: int(x.split('_')[1]))
        count_df            = count_df.sort_values(['bin','cluster']).drop(['bin','cluster'], axis=1).round({'frac_pos': 2, 'frac_neg': 2})
        count_df.to_csv(output.counts, sep='\t', index=False)
        annot_fns = glob(os.path.join(input.pol_dir, '*', '*.ref_read.strands.reads.tsv'))
        annot_dfs = [pd.read_csv(fn, sep='\t') for fn in annot_fns]
        annot_df  = pd.concat(annot_dfs)
        annot_df.to_csv(output.annots, sep='\t', index=False)

rule combine_dtr_aligns:
    input:
        medaka_dir = MEDAKA_DIR
    output: 
        cyc_perm_stats = temp(DTR_ALIGN_CYC_PERM_TSV)
    run:
        import os
        import pandas as pd
        from glob import glob
        fns = glob(os.path.join(input.medaka_dir, '*', '*.dtr.aligns.tsv'))
        df = pd.concat([pd.read_csv(fn, sep='\t') for fn in fns])
        df['dtr_len'] = df['tend'] - df['tstart']
        df['left_dist'] = df['tend']
        df['right_dist'] = df['tlen'] - df['tstart']
        for_df = []
        for group,df_ in df.groupby('bin_clust_id'):
            cyc_perm = False
            # If the distance from the genome end to the inner DTR boundary is more than 2x the median DTR length, keep.
            # If > 25% of the reads pass that filter, label the cluster as being cyclically permuted.
            df_filt = df_[df_['left_dist'] > (2 * df_['dtr_len'].median())]
            df_filt = df_filt[df_filt['right_dist'] > (2 * df_filt['dtr_len'].median())]
            if df_filt.shape[0]>(0.25*df_.shape[0]):
                cyc_perm = True
            for_df.append( (df_.loc[0,'bin_clust_id'], cyc_perm))
        df_out = pd.DataFrame.from_records(for_df, columns=['clust_id','cyc_perm'])
        df_out.to_csv(output.cyc_perm_stats, sep='\t', index=False)

rule combine_all_draft_stats:
    input:
        dtr_gc_stats=ALL_POL_DTR_GC_STATS,
        cds_stats=ALL_POL_CDS_SUMMARY,
        strand_stats=ALL_POL_STRANDS,
        cyc_perm_stats = DTR_ALIGN_CYC_PERM_TSV
    output: ALL_POL_STATS
    run:
        import pandas as pd
        df_gc     = pd.read_csv(input.dtr_gc_stats, sep='\t').set_index('clust_id')
        df_cds    = pd.read_csv(input.cds_stats, sep='\t').set_index('clust_id')
        df_strand = pd.read_csv(input.strand_stats, sep='\t').set_index('clust_id')
        df_cyc    = pd.read_csv(input.cyc_perm_stats, sep='\t').set_index('clust_id')
        df_cds = df_cds.drop(['read','ref_len'], axis=1)
        df_strand = df_strand.drop('read', axis=1)

        df = df_cds.join(df_strand, how='left')
        df = df.join(df_cyc, how='left')
        df = df.join(df_gc, how='left').reset_index().sort_values('clust_id')
        df.to_csv(output[0], sep='\t', index=False)