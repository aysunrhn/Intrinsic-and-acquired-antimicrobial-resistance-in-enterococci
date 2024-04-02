#!/bin/bash
# Bash script to download metadata about genomes assemblies of a particular species from the NCBI Assembly Database, removing anomalous entries
# Takes 2 positional arguments: species name and the output file path to write the metadata retrieved
# ./download_ncbi_assembly_metadata.sh enterococcus enterococcus_ncbi_assembly_metadata.tsv

esearch -query "$1 AND latest[filter] AND all[filter] NOT anomalous[filter]" -db assembly | efetch -format docsum | xtract -pattern DocumentSummary -element \
    GbUid,Organism,SpeciesName,Taxid,AssemblyName,AssemblyAccession,LastMajorReleaseAccession,AssemblyType,AssemblyStatus,BioSampleAccn,Coverage,ContigN50,ScaffoldN50,PartialGenomeRepresentation,AsmUpdateDate,SubmissionDate,SubmitterOrganization,FtpPath_GenBank,FtpPath_Stats_rpt -block Synonym -sep ";" -element Genbank,RefSeq,Similarity -block RS_BioProjects -subset Bioproj \
    -element BioprojectAccn -block GB_BioProjects -subset Bioproj -element BioprojectAccn -group Biosource -block InfraspeciesList -subset Infraspecie -sep "_" \
    -element Sub_type,Sub_value >> "$2"
