#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:08:49 2024

@author: aysun
"""
def clean_alignment(alignment_file, genome_list, gene):
    """
    Helper function to clean protein/gene alignment files and change the sequence
    IDs to indicate genome names

    Parameters
    ----------
    alignment_file : str
        Path to the alignment file, should be in fasta format
    genome_list : list
        List of genomes that the genes come from
    gene : str
        Name of the gene we aligned

    Returns
    -------
    None.
    Writes a new, cleaned alignment file in the same directory

    """
    import os
    from Bio import AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    if not os.path.isfile(alignment_file):
        print("ABORT!! -- Alignment file not found!")
        return
    msa = AlignIO.read(alignment_file, 'fasta')
    newmsa = []
    seqlen = msa.get_alignment_length()
    for genome in genome_list:
        rec = [r for r in msa if r.name == genome + '_' + gene]
        if not rec:
            seq = Seq(''.join(['-']*seqlen))
            rec = SeqRecord(seq, id='{}_{}'.format(genome, gene), name='', description='')
        else:
            rec = rec[0]
        newmsa.append(rec)
    newmsa = MultipleSeqAlignment(newmsa)
    AlignIO.write(newmsa, gene+'.cleaned.aln', format='fasta')

def extract_scc_orthogroups(orthofinder_dir):
    """
    Extract SCC orthogroups after running OrthoFinder
    
    Parameters
    ----------
    orthofinder_dir : str
        Directory where all the OrthoFinder results are stored

    Returns
    -------
    core_list : list
        A list of SCC orthogroups

    """
    import os
    import pandas as pd
    
    sc_file = os.path.join(orthofinder_dir, 'Orthogroups', 'Orthogroups_SingleCopyOrthologues.txt')
    ortho_file = os.path.join(orthofinder_dir, 'Orthogroups', 'Orthogroups.tsv')
    with open(sc_file, 'r') as f:
        sc_list = [line.strip() for line in f]
    df = pd.read_csv(ortho_file, sep='\t', index_col=0, header=0)
    num_genomes = len(df.columns)
    genome_count = df.apply(lambda x: num_genomes - sum(x.isna()), axis=1)
    core_list = genome_count[genome_count == num_genomes].index.intersection(sc_list)
    
    return core_list
    
def align_scc_orthogroups(orthofinder_dir, genome_list=None, scc_list=None, 
                          align_script=None, align_prot=True, g2path=None):
    """
    Align SCC orthogroup sequences, assumes the muscle binary is in the path

    Parameters
    ----------
    orthofinder_dir : str
        Directory where all the OrthoFinder results are stored
    genome_list : list, optional
        List of genomes that the genes belong to. The default is None.
    scc_list : list, optional
        List of SCC genes. The default is None.
    align_script : str, optional
        Path to write out alignment script, if you don't want to run the alignment
        within python. The default is None.
    align_prot : bool, optional
        Parameter to align protein sequences (instead of nucleotides). The default is True.
    g2path : dict, optional
        Dictionary mapping genomes to the data path where assemblies are stored. 
        The default is None.

    Returns
    -------
    None.
    Either writes an alignment file you can run through cli, or just produces SCC
    alignments
    """
    import os
    import subprocess
    import pandas as pd
    from Bio import SeqIO
    
    if not scc_list: 
        scc_list = extract_scc_orthogroups(orthofinder_dir)
    if not genome_list: 
        ortho_file = os.path.join(orthofinder_dir, 'Orthogroups', 'Orthogroups.tsv')
        genome_list = pd.read_csv(ortho_file, sep='\t', index_col=0, header=0).columns
        
    # Just to make NCBI acc IDs look cleaner
    glist = ['_'.join(g.split('_')[:2]) if g.startswith('GC') else g for g in genome_list]
    if align_script:
        out_handle = open(align_script, 'w+')
        out_handle.write('#!/bin/bash\n')
    
    # Align SCC proteins
    if align_prot:
        for i, ortho in enumerate(scc_list):
            print("Processing orthogroup {} -- {}".format(i+1, ortho))
            ortho_file = os.path.join(orthofinder_dir, 'Single_Copy_Orthologue_Sequences', 
                                      ortho + '.fa')
            in_file = os.path.join(orthofinder_dir, 'Single_Copy_Orthologue_Sequences', 
                                   ortho + '.input.fa')
            rec_list = [rec for rec in SeqIO.parse(ortho_file,'fasta') 
                        if '_'.join(rec.id.split('_')[:-1]) in genome_list]
            SeqIO.write(rec_list, in_file,'fasta')
            alignment_file = os.path.join(orthofinder_dir, 'Single_Copy_Orthologue_Sequences', 
                                          ortho + '.aln')
            cmd = ('muscle -in {} '.format(in_file) + '-out {}'.format(alignment_file))
            if not align_script:
                print("Aligning orthogroup {}".format(ortho))
                output = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                        universal_newlines=True)
                if output.returncode == 0:
                    clean_alignment(alignment_file, genome_list, ortho)
                else:
                    print("Ooops something went wrong!!")
            else:
                out_handle.write(cmd + '\n')
    # Align SCC nucleotides
    else:
        for i, ortho in enumerate(scc_list):
            print("Processing orthogroup {} -- {}".format(i+1, ortho))
            ortho_file = os.path.join(orthofinder_dir, 'Orthogroup_Sequences', 
                                      ortho + '.fa')
            in_file = os.path.join(orthofinder_dir, 'Single_Copy_Orthologue_Sequences', 
                                   ortho + '.input.fa')
            rec_list = []
            for rec in SeqIO.parse(ortho_file,'fasta'):
                genome = '_'.join(rec.id.split('_')[:-1])
                cds = rec.id
                if genome in glist:
                    cdspath = g2path[genome]
                    addrec = [rec for rec in SeqIO.parse(cdspath, 'fasta') if rec.id == cds]
                    if addrec: 
                        addrec = addrec[0]
                        addrec.id = addrec.id.split()[0]
                        rec_list.append(addrec)
            SeqIO.write(rec_list, in_file,'fasta')
            alignment_file = os.path.join(orthofinder_dir, 'Single_Copy_Orthologue_Sequences', 
                                          ortho + '.aln')
            cmd = ('muscle -in {} '.format(in_file) + '-out {}'.format(alignment_file))
            if not align_script:
                print("Aligning orthogroup {}".format(ortho))
                output = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
                                        universal_newlines=True)
                if output.returncode == 0:
                    clean_alignment(alignment_file, genome_list, ortho)
                else:
                    print("Ooops something went wrong!!")
            else:
                out_handle.write(cmd + '\n')    
        
    if align_script: 
        out_handle.close()
    
def concat_scc_alignments(orthofinder_dir, fix_prefix, scc_list=None, genome_list=None,
                          align_prot=False, msa_file='SCCorthogroups_aligned.fasta',
                          partitions_file='SCCorthogroups-partition_file.txt'):
    """
    Concatenate SCC orthogroup alignments in the Orthofinder results directory to 
    obtain a single fasta alignment file of all genomes    

    Parameters
    ----------
    orthofinder_dir : str
        Directory where all the OrthoFinder results are stored
    fix_prefix : dict
        A dictionary, if you want to fix the prefixes of sequence IDs in the alignment
    scc_list : list, optional
        List of SCC genes. The default is None.
    genome_list : list, optional
        List of genomes that the genes belong to. The default is None.
    align_prot : bool, optional
        Parameter to align protein sequences (instead of nucleotides). The default is True.
    msa_file : str, optional
        Path to the final alignment file to write. The default is 'SCCorthogroups_aligned.fasta'.
    partitions_file : str, optional
        Path to the partitions file, will be input to IQTREE. The default is 
        'SCCorthogroups-partition_file.txt'.

    Returns
    -------
    None.
    Writes a MSA file
    """
    import os 
    import pandas as pd
    from Bio import SeqIO, AlignIO
    from Bio.SeqRecord import SeqRecord
    
    if not genome_list:
        ortho_file = os.path.join(orthofinder_dir, 'Orthogroups', 'Orthogroups.tsv')
        genome_list = pd.read_csv(ortho_file, sep='\t', index_col=0, header=0).columns
    if not scc_list: 
        scc_list = extract_scc_orthogroups(orthofinder_dir)
    genome2Seq = {g: '' for g in genome_list}
    start = 1
    partID = 1
    if partitions_file: 
        fpart = open(os.path.join(orthofinder_dir, partitions_file), 'w')
        alignment_type = 'DNA' if not align_prot else 'PROT'
    for i, ortho in enumerate(scc_list):
        print("Processing orthogroup {} -- {}".format(i+1, ortho))
        alignment_file = os.path.join(orthofinder_dir, 'Single_Copy_Orthologue_Sequences', 
                                      ortho + '.aln')
        msa = AlignIO.read(alignment_file, 'fasta')
        for rec in msa:
            genome = '_'.join(rec.id.split('_')[:-1])
            if genome.endswith('pr') and 'pr_' in rec.id: genome = genome.replace('pr', '')
            genome = fix_prefix.get(genome, genome)
            if not genome in genome_list:
                print("This genome {} isn't supposed to be here!! ABORT!".format(genome))
                continue
            genome2Seq[genome] = genome2Seq[genome] + rec.seq
        gene_len = len(str(rec.seq))
        stop = start + gene_len - 1
        if partitions_file: 
            fpart.write('{}, part{} = {}-{}\n'.format(alignment_type, partID, start, stop))
        start += gene_len
        partID += 1

    all_rec = [SeqRecord(seq, id=genome, name='', description='') for genome, seq in genome2Seq.items()]
    SeqIO.write(all_rec, os.path.join(orthofinder_dir, msa_file), format='fasta')
    if partitions_file: 
        fpart.close()