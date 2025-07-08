# Determining Phylogenetically Averaged Motif (PAM) score representations for the non-coding regions, promoters and 3'UTRs, of Caenorhabditis Elegans

## Table of Contents

1. Getting the promoter and 3'UTR sequences for Caenorhabtidis Elegans and orthologous species
2. Getting the Position Weight Matrices (PWMs) for Transcription Factors (TFs) and RNA-Binding Proteins (RBPs) of Caenorhabtidis Elegans
3. Building the PAM score representations of promoters and 3'UTRs

## 1. Getting the promoter and 3'UTR sequences for Caenorhabtidis Elegans and orthologous species

The WormBase gene ids for all protein-coding genes of Caenorhabtidis Elegans (C. Elegans) were extracted from 
the most recent gtf file, which can be accessed through this link: [C. Elegans GTF file](https://asia.ensembl.org/Caenorhabditis_elegans/Info/Index) 
which contained a total of 19985 gene ids for protein-coding genes. This list of gene ids can be accessed through the 
`celegans_gene_ids.xlsx` file.

To get the promoter and 3'UTR sequences, Ensembl's REST API was used and documented one-to-one orthologs for each C. Elegans gene 
were saved from the metazoan species available.

This is a step-by-step guide about how it was done, which can be done by running the `fetching_promoters_and_3UTRs.py` script 
(warning that it takes a long time to run due to the volume of fetch requests required). All the fetch requests start with the 
following Ensemble REST API specific identifier https://rest.ensembl.org:

### 1.1 Fetch the one-to-one orthologs of each of C. Elegans' protein-coding genes from the metazoan species in Ensembl
For each C. Elegans' gene id, the orthologs from the metazoan compara were fetched using the following url:
/homology/id/caenorhabditis_elegans/{gene_id}?type=orthologues;compara=metazoan;content-type=application/json. This would be the 
url for the gene id WBGene00000001 for example: https://rest.ensembl.org/homology/id/caenorhabditis_elegans/WBGene00000001?type=orthologues;compara=metazoa;content-type=application/json.
Only the one-to-one orthologs were considered and for each one-to-one ortholog, the orthologous species name and gene id were saved 
which are required to be able to get the sequence information later.

### 1.2 Fetch the gene information for each C. Elegans and orthologous metazoan species
For each C. Elegans gene as well as each metazoan species' one-to-one orthologous genes, information about each gene's 
transcripts are fetched through the url: /lookup/id/{gene_id}?expand=1;species={species_name}. This would be the url for 
the gene id WBGene00000031 of C. Elegans for example: https://rest.ensembl.org/lookup/id/WBGene00000031?expand=1;species=caenorhabditis_elegans. 
Since a protein-coding gene can have multiple transcripts, only the Ensemble canonical transcripts were considered and for 
each gene id's canonical transcript, the following information was saved to json files:
- start and end coordinates of the gene (transcript region) and canonical transcript (translation region). For WBGene00000031 for example,
  - For the transcript region, the start and end coordinates are 5192995 and 5194608
  - For the coding region, the start and end coordinates are 5193139 and 5194567
- strand (+1 or -1) (for example, it is -1 for WBGene00000031)
- Region name of the sequence the gene is found in (for example, it is V for WBGene00000031)

### 1.3 Fetch the promoter and 3'UTR sequences for each gene information saved
The length of the promoters and 3'UTRs can be specified in the script, but for the current work, promoter sequences were 
determined as 600bp upstream of the transcription start site while 3'UTR sequences as 200bp downstream of the translation end site. 

For each gene (C. Elegans and orthologs), the gene information saved previously in step 1.2 were used to extract these upstream 
and downstream regions. The start and end coordinated for the promoter and 3'UTR sequences were determined using the start and 
end coordinates of the transcript and translation regions as well as the strandness. The promoter and 3'UTR sequences were then 
fetched using the following url: /sequence/region/{species}/{chromosome}:{seq_start}..{seq_end}:{strand}. 
This would be the url for the promoter sequence of the gene id WBGene00000031 for C. Elegans for example: 
https://rest.ensembl.org/sequence/region/caenorhabditis_elegans/V:5194609..5195208:-1

The final output files were saved as json files into two folders, one for promoter sequences and one for 3'UTR sequences. 
Each json file was for a C. Elegans gene where the keys are the species name (C. Elegans and orthologs) and the 
values are strings for either the promoter or 3'UTR sequences associated with each species for the gene.

## 2. Getting the Position Weight Matrices (PWMs) for Transcription Factors (TFs) and RNA-Binding Proteins (RBPs) of Caenorhabtidis Elegans

Position Frequency Matrices (PFMs) for RNA Binding Proteins (RBPs) were downloaded from CisBP-RNA database (https://cisbp-rna.ccbr.utoronto.ca/bulk.php) 
on July 8, 2025 for a total of 60 RBP motifs to be scanned for binding to 3'UTRs. PFMs for Transcription factors (TFs) were 
downloaded from Jaspar at this [link](https://jaspar.elixir.no/search?q=&collection=CORE&tax_group=all&tax_id=6239&type=all&class=all&family=all&version=latest), 
where the advanced options for searching are specified, for a total of 94 TFs to be scanned for binding to promoter sequences. All PFMs were converted to 
PWMs following the method outlined in Alam et al. 2023 paper, and was done using the `collecting_motif_data.py` script.

The matrix information for each PWM was saved as a dictionary where the keys are the motif names and the values are the 
matrices as numpy arrays of shape: 4 nucleotides x length of PWM. These dictionaries were saved as pickle files, one for 
the PWMs of TFs `PWMs_of_TF_for_caenorhabditis_elegans.pkl` and one for the PWMs of RBPs `PWMs_of_RBP_for_caenorhabditis_elegans.pkl`.

# 3. Building the PAM score representations of promoters and 3'UTRs


