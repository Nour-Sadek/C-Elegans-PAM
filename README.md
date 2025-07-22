# Determining Phylogenetically Averaged Motif (PAM) score representations for the non-coding regions, promoters and 3'UTRs, of Caenorhabditis Elegans

The project aims to replicate what was previously done for Saccharomyces Cerevisiae and human promoters and 3'UTRs (Alam et al. 2023), 
but for the corresponding non-coding regions for Caenorhabditis Elegans genes, as the resources to test the hypotheses generated 
by the promoter and 3'UTR representations are available for that model organism.

Briefly, the authors mapped promoter and 3'UTR sequences to a feature space/representation that would be informative about 
the function of the gene associated with those non-coding regions. They used Position Weight Matrices (PWMs) of Transcription 
Factors (TFs) for promoters and RNA-Binding Proteins (RBPs) for 3'UTRs. PWMs can be considered as equivalent to a kernel of 
a 1d convolutional layer of shape 4 nucleotides x length of the motif where the values for each nucleotide at each position
is the probability of that nucleotide being in that position so that the DNA sequence matches the consensus sequence 
needed for confidence binding of the RBP/TF to that sequence. To determine the confidence of the specific TF/RBP binding, 
its corresponding PWM is scanned over the sequence of interest and the output of the scan at each position gives a confidence 
score where higher scores means higher probability of binding and lower scores means lower probability of binding. 

The aim is to scan the PWMs of consensus TF and RBP motifs over promoter and 3'UTR sequences respectively where a high score 
after convolutional scanning could be interpreted as high probability of that motif binding to that region, which could be used 
to infer function of the gene associated with that non-coding region. However, there are a lot of false positives in terms of 
sequence matching, and so phylogenetic averaging was introduced where function, more so than sequence, would be conserved 
across homologs, and so scanning the PWM over the non-coding region and its orthologs, then averaging the scores over the 
orthologues - thus calculating a Phylogenetically Averaged Motif, or PAM score - is thought to potentially allow the interpretation 
of the averaged high match score as real biological signal rather than false positive. 

There is a representation vector for every non-coding region and PWM motif combo where that representation is built by 
averaging the maximum PWM score over all orthologous sequences of the gene that corresponds to that non-coding region. 

The goal, after clustering and visualization of the representation heatmap, is to inspect clusters enriched for various GO 
categories, and test hypotheses regarding the potential function of unknown genes based on how its representation clustered 
with the representation of other genes.


## Table of Contents

1. Collection of the promoter and 3'UTR sequences for Caenorhabtidis Elegans and homologous species
2. Collection of the Position Weight Matrices (PWMs) for Transcription Factors (TFs) and RNA-Binding Proteins (RBPs) of Caenorhabtidis Elegans
3. Building of the PAM score representations of promoters and 3'UTRs of Caenorhabtidis Elegans genes
4. Showcase of the output heatmaps
5. Go Enrichment Analysis
6. References

## 1. Data collection of the promoter and 3'UTR sequences for Caenorhabtidis Elegans and homologous species

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
Only the one-to-one orthologs were considered and for each one-to-one ortholog, the homologous species name and gene id were saved 
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
end coordinates of the transcript and translation regions as well as the strand-ness. The promoter and 3'UTR sequences were then 
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

## 3. Building the PAM score representations of promoters and 3'UTRs of Caenorhabtidis Elegans genes

The same workflow that was done in Alam et al. 2023 paper was implemented here but PyTorch was used instead of Keras/tensorflow. 
The implementation and generation of the PAM score representations can be found in the `building_the_representation.py` script.
One minor modification was that the PWMs were not padded to the same length as that of the longest PWM (which allows them to be 
all loaded together as one convolutional layer), but rather each PWM was scanned individually over every sequence. This was done 
because padding was thought to prevent the scanning of the start and end of the non-coding regions, so this was ensures that the 
entire length of the sequences get scanned, but this trades off the speed of the convolutional scanning.

### 3.1 How the sequences were filtered

On top of including sequences that only had the four letter A, G, T, and C, different additional cutoffs to which genes are incorporated 
in the convolutional scanning were implemented and compared:
- The first scanning was done over genes that had a minimum of 10 homologous species from the metazoan compara. This lead 
to the final inclusion of 
  - 5596 genes and a total of 322077 sequences for promoters
  - 5705 genes and a total of 372990 sequences for 3'UTRs
- The second scanning was done over genes that had a minimum of 5 homologous nematoda species. The nematoda 
species are: `["ascaris_summ", "brugia_malayi", "caenorhabditis_brenneri", "caenorhabditis_briggsae", 
"caenorhabditis_elegans", "caenorhabditis_japonica", "caenorhabditis_remanei", "loa_loa", "necator_americanus", 
"onchocerca_volvulus", "pristionchus_pacificus", "strongyloides_ratti", "trichinella_spiralis", "trichuris_muris"]`. This lead 
to the final inclusion of
  - 8233 genes and a total of 67938 sequences for promoters
  - 8467 genes and a total of 70623 sequences for 3'UTRs

These two different cutoffs were to decide whether including more evolutionarily distant species would lead to better PAM 
score representations than including fewer, but evolutionarily closer species.

There is also a third scanning that was done over genes that had orthologous sequences for all the 5 caenorhabditis 
species `["caenorhabditis_brenneri", "caenorhabditis_briggsae", "caenorhabditis_elegans", "caenorhabditis_japonica", 
"caenorhabditis_remanei"]` only, and that is because these four species are the closest evolutionarily to caenorhabditis 
elegans, so it is interesting to see how much the representations will differ. This lead to the final inclusion of:
- 3515 genes and a total of 17575 sequences for promoters
- 3908 genes and a total of 19540 sequences for 3'UTRs

### 3.2 Workflow of the convolutional scanning and normalization

Briefly, for each gene, both the forward and reverse complement of PWMs were scanned over the homologous promoter sequences, 
with only the forward PWM being scanned over 3'UTR sequences considering that RNA is single-stranded, and for promoter sequences, 
the max score at every position between the forward and reverse scans was chosen. And since more confident PWM binding is 
associated with higher convolutional scores, a max pool over the positional scores was done to get the max PWM score. After that, 
the scores for all the orthologous sequences were averaged to get the Phylogenetically Averaged Motif (PAM) score for that 
gene's non-coding region for that PWM.

In order to allow for comparisons of PAM scores between PWMs of different lengths, normalization of those PAM scores was done. 
Briefly, each PWM was scrambled 100 times and similarly scanned over the sequences. For these scrambled PWMs, the speed of 
convolutional layers was exploited where the 100 scrambled PWMs were loaded as filters in a single convolutional layer. The 
z-score between the real PAM score of the original PWM and the distribution of the 100 scrambled PWMs was calculated, and 
those z-scores were used in building the final representation.

### 3.3 Visualization of the promoter and 3'UTR PAM score representations

The output of the script essentially generates a dataframe populated with PAM scores between the promoter/3'UTR sequences
(rows) and the motifs (columns). Cluster 3 software was used to cluster the PAM scores. The pythonic package Bio.Cluster was used 
where, similarly to Alam et al. 2023 paper, the columns were median centered, followed by performing hierarchical clustering 
on both genes and motifs using un-centered correlation, average linkage. This provides multiple output files, which are read 
by the software Java TreeView to allow for proper visualization of the PAM scores representation heatmap.

### 4. Showcase of the output heatmaps

These are the output heatmaps from Java TreeView where the heatmap of clustered promoter representations is on the left 
and of clustered 3'UTR representations is on the right.

For when a minimum of 10 orthologous sequences from the metazoan compara per gene was the cut-off:

<img width="1038" height="1060" alt="Image" src="https://github.com/user-attachments/assets/055c2340-141d-43cf-8f67-ab2c1b71e7ef" />

For when a minimum of 5 orthologous sequences from the nemotoda species per gene was the cut-off:

<img width="988" height="1069" alt="Image" src="https://github.com/user-attachments/assets/64bd67b8-fc51-4a06-854b-850e88d4c537" />

For when orthologous sequences from the 5 caenorhabditis species per gene was the cut-off:

<img width="998" height="1069" alt="Image" src="https://github.com/user-attachments/assets/18fd183b-53b8-4ff0-969c-f439bb2f0def" />

## 6. References

Alam A, et al. Functional similarity of non-coding regions is revealed in phylogenetic average motif score representations. biorxiv. 2023. doi: https://doi.org/10.1101/2023.04.09.536185

Alok J. Saldanha, Java Treeview—extensible visualization of microarray data, Bioinformatics, Volume 20. https://jtreeview.sourceforge.net/

M. J. L. de Hoon, S. Imoto, J. Nolan, and S. Miyano: Open Source Clustering Software. Bioinformatics, 20 (9): 1453--1454 (2004). http://bonsai.hgc.jp/~mdehoon/software/cluster/

Ray D, et al. A compendium of RNA-binding motifs for decoding gene regulation. Nature. 2013 Jul 11;499(7457):172-7. doi: 10.1038/nature12311. 
PubMed PMID: 23846655. https://cisbp-rna.ccbr.utoronto.ca/

Rauluseviciute I, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. 
Nucleic Acids Res. 2024 Jan 5;52(D1):D174-D182.; doi: 10.1093/nar/gkad1059. https://jaspar.elixir.no/

Sarah CD, et al. Ensembl 2025. Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D948–D957, https://doi.org/10.1093/nar/gkae1071
