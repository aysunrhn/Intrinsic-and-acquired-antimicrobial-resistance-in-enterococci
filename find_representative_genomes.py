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
mash_df = pd.read_csv('output/mashdist-allvall.tsv', sep='\t', index_col=0, header=0)
stat_df = pd.read_csv('output/checkm_asmstats.tsv', sep='\t', index_col=0, header=0)    

genome_list = stat_df.index.to_list()
D = mash_df.loc[genome_list, genome_list].to_numpy()
d = squareform(D)
Z = linkage(d, 'single')
c_df = cluster_utils.cluster_genomes(Z, D, genome_list, th_dif=0.05, use_dist=True)

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