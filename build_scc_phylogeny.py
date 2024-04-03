#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:57:32 2024

@author: aysun
"""
import os
import pandas as pd
from utils import orthofinder_utils

#%% Load and process OrthoFinder results
entero_df = pd.read_csv('output/enterococcus_metadata_repgenomes.tsv', sep='\t', 
                        index_col=0, header=0)
genome_list = entero_df.repgenome.unique()

orthofinder_dir = 'output/orthofinder/out'
ortho_count = pd.read_csv(os.path.join(orthofinder_dir, 'Orthogroups/Orthogroups.GeneCount.tsv'), 
                          sep='\t', index_col=0, header=0)
# Extract SCC genes and align them
scc_list = orthofinder_utils.extract_scc_orthogroups(orthofinder_dir)

align_script = os.path.join(orthofinder_dir, 'align_scc_repgenomes.sh')
orthofinder_utils.align_scc_orthogroups(orthofinder_dir, list(scc_list), align_script)

glist = ['_'.join(g.split('_')[:2]) if g.startswith('GC') else g for g in genome_list]
fix_prefix = {g1: g2 for g1, g2 in zip(glist, genome_list)} # a small fix for some genome names
orthofinder_utils.concat_scc_alignments(orthofinder_dir, fix_prefix, list(genome_list), 
                                        msa_file='SCCorthogroups_aligned.fasta')

