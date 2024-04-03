#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 22:56:57 2024

@author: aysun
"""
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from utils import cluster_utils

# Initial rough clustering based on MASH distance
mash_df = pd.read_csv('output/mashdist_allvall.tsv', sep='\t', index_col=0, header=0)
ph_df = pd.read_csv('output/phylosift_results.tsv', sep='\t', index_col=0, header=0) 
stat_df = pd.read_csv('output/checkm_asmstats.tsv', sep='\t', index_col=0, header=0)    

# Prepare scripts to run fastANI and refine clusters base on ANI
genome_list = stat_df.index.to_list()
D = mash_df.loc[genome_list, genome_list].to_numpy()
d = squareform(D)
Z = linkage(d, 'single')
c_df = cluster_utils.cluster_genomes(Z, D, genome_list, th_dif=0.05, use_dist=True)
fastani_file = 'run_fastani_withincluster.sh'
ff = open(fastani_file, 'w')
ff.write('#!/bin/bash\n')
for c, row in c_df.iterrows():
    if row.cluster_size == 1: 
        continue
    output_file = 'output/fastani_input_withincluster_{}.txt'.format(c)
    with open(output_file, 'w') as f:
        for g in row.genomes:
            f.write('data/genome/{}.genome.fa\n'.format(g))
    cmd = 'fastANI --ql {} --rl {} -t 16 -o fastani_output_withincluster_{}.out'.format(output_file, 
                                                                                        output_file, 
                                                                                        c)
    ff.write(cmd + '\n')

# 0. Make a high quality set of representative genomes
entero_df = pd.read_csv('data/enterococcus_metadata.tsv')
stat_df = pd.read_csv('output/enterococcus_asmstats.tsv', sep='\t', index_col=0, header=0)    
pass_phylo = ph_df[(ph_df.apply(lambda x: sum([1 for g in x[1:] if g>0]), axis=1) > 30) 
                   & (ph_df.flag < 2)].index.values
fail_phylo = ph_df.index.difference(pass_phylo).values
sort_columns = ['num_comp', 'contamination', 'phylosift_flag', 'completeness','n50']
for species in entero_df.species.unique():
    all_genomes = entero_df[entero_df.species == species].index.values
    if len(all_genomes) == 1:
        repgenome = all_genomes[0]
    else:
        genomes = set(all_genomes.copy()).difference(fail_phylo)
        if len(genomes) == 0: 
            genomes = all_genomes.copy()
        repgenome = stat_df.loc[genomes].sort_values(by=sort_columns, 
                                                     ascending=[True, True, False, False, False])
        repgenome = repgenome.iloc[0].name
    entero_df.loc[all_genomes, 'repgenome'] = repgenome
    entero_df.loc[all_genomes, 'repgenome_num_comp'] = stat_df.loc[repgenome, 'num_comp']
    entero_df.loc[all_genomes, 'repgenome_n50'] = stat_df.loc[repgenome, 'n50']
    entero_df.loc[all_genomes, 'repgenome_dataset'] = stat_df.loc[repgenome, 'dataset']
    
    repgenome_path = 'data/genome/{}.fna'.format(repgenome)              
    entero_df.loc[all_genomes, 'rep_genome_path'] = repgenome_path
    
entero_df.to_csv('output/enterococcus_metadata_repgenomes.tsv', sep='\t', index=True, header=True)    