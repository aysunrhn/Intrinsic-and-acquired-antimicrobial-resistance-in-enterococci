#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 23:02:16 2024

@author: aysun
"""

def cluster_genomes(Z, D, genome_list, th_dif=0.05, use_dist=True):
    """
    Cluster genomes based on their distance or similarity, uses hierchical
    clustering from scipy

    Parameters
    ----------
    Z : np.array
        Z matrix produced by running 'linkage' function in scipy
    D : np.array
        Distance or similarity matrix
    genome_list : list
        List of genomes
    th_dif : float, optional
        Distance threshold to define clusters. The default is 0.05.
    use_dist : bool, optional
        Parameter to treat matrices as distance values. The default is True.

    Returns
    -------
    c_df : pd.DataFrame
        A dataframe mapping each cluster to its members, cluster size and min/max ANI

    """
    import pandas as pd
    import numpy as np
    from scipy.cluster.hierarchy import fcluster
    
    if isinstance(genome_list, list): 
        genome_list = np.array(genome_list)
    
    # Treat matrix entries as distances
    if use_dist:  
        clusters = fcluster(Z, t=th_dif, criterion='distance')
        c_dict = {c: {'genomes': [], 'cluster_size': 0, 'min_dist': 0, 'max_dist': 1} 
                  for c in range(1, max(clusters) + 1)}
        for m in set(clusters):
            idx = np.where(clusters == m)[0]
            genomes = genome_list[idx]
            min_dist = 1
            max_dist = 0
            for i in idx:
                for j in idx:
                   if i == j: 
                       continue
                   d = D[i, j]
                   if d < min_dist: 
                       min_dist = d
                   if d > max_dist: 
                       max_dist = d
            if min_dist == 1: 
                max_dist = 1
            if len(genomes) == 1: # only 1 genome in the cluster
                min_dist = 0
                max_dist = 0
            c_dict[m]['genomes'] = genomes
            c_dict[m]['cluster_size'] = len(genomes)
            c_dict[m]['min_dist'] = min_dist
            c_dict[m]['max_dist'] = max_dist
            
    # Treat matrix entries as similarities
    else:
        th_dif = 1 - th_dif           
        clusters = fcluster(Z, t=th_dif, criterion='distance')
        c_dict = {c: {'genomes': [], 'cluster_size': 0, 'min_sim': 0, 'max_sim': 1} 
                  for c in range(1, max(clusters) + 1)}
        for m in set(clusters):
            idx = np.where(clusters == m)[0]
            genomes = genome_list[idx]
            min_sim = 1
            max_sim = 0
            for i in idx:
                for j in idx:
                   if i == j: 
                       continue
                   d = 1 - D[i, j]
                   if d < min_sim: 
                       min_sim = d
                   if d > max_sim: 
                       max_sim= d
            if min_sim == 1: 
                max_sim= 1
            if len(genomes) == 1: 
                min_sim = 1
                max_sim = 1
            c_dict[m]['genomes'] = genomes
            c_dict[m]['cluster_size'] = len(genomes)
            c_dict[m]['min_sim'] = min_sim
            c_dict[m]['max_sim'] = max_sim
            
    c_df = pd.DataFrame(c_dict.values(), index=c_dict.keys())
    
    return c_df