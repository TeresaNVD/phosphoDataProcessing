# Final project: scripts description

Author: Teresa Nunez de Villavicencio Diaz

Contact: tnunezde@uwo.ca

Supervisor: David W. Litchfield (Biochemistry)

Professor: Art Poon

The project folder (folder: final_project or ./) contains four python scripts (folder: ./scripts) that are integrated in an analysis pipeline by a bash script (folder: ./).

Python scripts:
>1. exp_summary.py: summarizes information like the proteins, gene names, and sequences found in the experiment after filtering contaminants and reverse hits.
>2. phosphosite.py: generates the information on the phosphosites found and assign to the phosphosites known kinase-substrate relationships.
>3. motif.py: scans sequences for kinase motifs.
>4. function_seq.py: retrieves the sequence and function of the proteins identified in UniProt database.

Bash script:
>run\_python_scripts.sh: Integrates all of the python scripts and provide snippets of the results. Should be executed at the root of the project folder otherwise specified by the user.

To save the output of the bash script to a file run on the command line at the root of the project folder:

```bash
./run_python_scripts.sh > output.txt
```

Requirements:
>Biopython, Weblogo ([Weblogo installation page](http://weblogo.threeplusone.com/manual.html#CLI "Title"))

---

## <center>Experiment summary (1)</center>

## Script name (1)

`exp_summary.py`

## Objective (1)

The Python script: experiment summary (exp_summary), allows to summarize important information extracted from the Phospho(STY)Sites.txt table generated with the MaxQuant software as part of the computational proteomics pipeline for the analysis of phosphoproteomics data generated here at UWO in the group of David Litchfield. MaxQuant is a freely available program widely used in proteomics groups around the world.

[MaxQuant page](http://www.biochem.mpg.de/5111795/maxquant "Title")

[MaxQuant documentation page](http://www.coxdocs.org/doku.php?id=maxquant:start "Title")

## Input (1)

Typically the raw files generated by the mass spectrometers are analyzed using a preferred computational pipeline for generating information on peptide and protein identification and quantification. 
In our laboratory raw files are analyzed using MaxQuant, which returns tabular information on the identified/quantified peptides and proteins. Since our laboratory is interested in the dynamics of phosphosites in response to
treatments with kinase inhibitors we make use of proteomics to assess the changes in the phosphoproteome of cells upon treatment. As a result, an important part is to understand the MaxQuant output contained in the Phospho(STY)Sites.txt table.

**Input:** Phospho(STY)Sites.txt

**Content:** A tab separated file containing information on the phosphopeptides (including the phosphosite) identified and quantified.

**Path:**  /final_project/data/

**Parameters:**  HERE

## Output (1)

Depending on the number of raw files analyzed in a given run using MaxQuant the Phospho(STY)Sites.txt file will contain more or less columns with information specific for each raw file making the size of the file and its readability in Excel challenging most of the time. In this regard the current script allows for computational access to this file and for summarizing the information in it.

1. Lists all the raw files linked to a given MaxQuant run and thus to a given Phospho(STY)Sites.txt output. Important for understanding what information is being analyzed specially if little information was written down.
2. Filters and lists the entries (rows) of contaminants and reverse hits (FDR control). These are entries that we are not interested on in the biological sense but provide quality information on the experiment being analyzed.
3. Lists all the gene names associated to phosphopeptides identified in the experiment.
4. Lists all the proteins (format: UniProt accessions) associated to phosphopeptides identified in the experiment.
5. Lists all the amino acids modifications associated to phosphopeptides identified in the experiment (e.g., Met oxidation and Asn/Gln deamidation).
6. Lists all the sequences containing phosphopeptides identified in the experiment (uses MaxQuant sequence window by default: 31 amino acids centered at the 16 residue).
7. Formats all the sequences containing phosphopeptides identified in the experiment using the PhosphoSitePlus (15 amino acids centered at the 8 residue) and Perseus short sequence (13 amino acids centered at the 7 residue) formats.

Lists: means the script outputs this information to a new file (unique elements only).

Filters: creates a filtered version of the Phospho(STY)Sites.txt file.

## Use (1)

The extraction of such information is essential for the analysis of the data in the biological context and constitute the entry for many functional analysis tools such as: enrichment tools (list of proteins and genes) and motif discovery and sequence logo representation tools (list of sequences).

---

## <center>Phosphosite (2)</center>

## Script name (2)

`phosphosite.py`

## Objective (2)

The information regarding the phosphosites identified and quantified in a given MaxQuant run is found in different columns across the Phospho(STY)Sites.txt file. The purpose of this script is to bring together this information and output it in a format that can be use to search the data for kinase-substrate information using known post-translational modifications (PTM) databases and prediction tools.

In addition the kinase-substrate relationships listed in PhosphoSitePlus are mined to retrieved known associations found in the experiment.

[Example #1 of phosphosite format](https://research.bioinformatics.udel.edu/iptmnet/search_bulk "Title")

[Example #2 of phosphosite format](http://networkin.info/index_ht.shtml "Title")

## Input (2)

This scripts takes to files as input the Phospho(STY)Sites_filtered.txt generated by (1) and a file containing kinase-substrate relationship in classical formatting more specifically Kinase_Substrate_Dataset.txt downloaded from PhosphoSitePlus.

**Input:** Phospho(STY)Sites.txt or contaminant and reverse filtered Phospho(STY)Sites_filtered.txt. Kinase_Substrate_Dataset.txt downloaded from PhosphoSitePlus

**Content:** A tab separated file containing information on the phosphopeptides (including the phosphosite) identified and quantified.

**Path:**  /final_project/data/

**Parameters:**  HERE

## Output (2)

1. Converts the columns containing information on the protein indentifier, the position of the modified site within the protein, and the amino acid at that position to three possible site formats defined by ptm_mode variable:
* PhosphoSitePlus or PSP (PTM database): 'psp'; e.g., Q9Y4H2   T779-p
* iPTMnet (PTM database): 'iptmnet'; e.g., Q9Y4H2   T   779
* NetworKIN (kinase-substrate relationship prediction): HERE
* Classic (also used in PhosphoSitePlus and in general for reporting phosphosite information): e.g., Q9Y4H2   T779
2. Lists the converted phosphosites.
3. Maps the classic formatted phosphosites to the kinase-substrate relationships information downloaded from PhosphoSitePlus.
4. Lists the known kinase-substrate relationships found (kinase_susbstrate.tsv).

[PhosphoSitePlus or PSP page](https://www.phosphosite.org/homeAction.action "Title")

[iPTMnet page](https://research.bioinformatics.udel.edu/iptmnet/ "Title")

[NetworKIN page](http://networkin.info/index_ht.shtml "Title")

## Use (2)

The generation of phosphosite information is essential for finding the kinases and pathways being modulated in a given experiment. The list of phosphosites outputted by this script will serve as the query for the iPTMnet database which is a major repository of PTM information (not only phosphorylation but acetylation, methylation, etc) as well as the input for prediction tools that use kinases consensus motifs and tools like STRING to estimate functional associations. To give an estimate of such kinase-substrate relationships represented the script matches the phosphosites found in the experiment to the information in PhosphoSitePlus one of the sources integrated in iPTMnet.

---

## <center>Motif (3)</center>

## Script name (3)

`motif.py`

## Objective (3)

Although known kinase-substrate information is preferred to understand the underlying biology of a phosphoproteomics experiment the information in the PTM databases or even the literature is often incomplete and centered around the proteins studied the most. Usually then the researcher has to rely on the mapping of know kinase consensus motifs rather than to look for individual substrates. In this regard, huge effort has been put to systematically annotate known kinase motifs and one of the main sources is the HPRD database. Thus using the motifs provided in this database any list of sequences can be scanned for matches and generate hypothesis as to which kinases are modulated in the experimental setting used. Then, the scientists can perform a basic motif search using regex and scripting and also use kinase-substrate prediction tools (essentially the same but occasionally including other biological data like protein interaction to the prediction) with an existing tool and compare the results. 

The lists of sequences centered at the modified residue are frequently plotted in a sequence logo to illustrate the variation or regularity at the different positions, which works as a graphic way of pointing towards the kinase motifs represented.

[HPRD motif page](http://hprd.org/serine_motifs "Title")

[Weblogo API page](http://weblogo.threeplusone.com/manual.html#API "Title")

## Input (3)

This script takes two files one providing the motifs and the associated kinases and the other providing the sequences, which will be the list of unique sequences found in the experiment ()different sequence windows can be used.

**Input:** HPRD_motifs_Jul2018.tsv and sequence file generated with exp_summary.py script (1); different sequence windows can be used PSP or MAxQuant short sequence (see script 1)

**Content:** A tab separated file containing two columns: motif and kinase, in that order (generated from HPRD_motifs_Jul2018.tsv using the batch code below) and a one column file with the sequences centered at the modified site.

**Path:**  /final_project/data/ ; /final_project/results/

**Parameters:**  HERE

## Output (3)

Since the HPRD motifs are written in natural language first they are converted to regex expressions after that the motif scan is performed.

1. Converts the motifs downloaded from HPRD to regex.
2. Lists the converted motifs.
3. Reads the list of sequences.
4. Lists the matched sequences (kinase_motif_scan.tsv).
5. Generates a sequence logo to visualize the frequency of amino acids at each posiiton of the list of sequence passed as input (uses Weblogo).
6. Lists all the sequences that matched the CK2 kinase motif as our group focuses in the study of this particular kinase.

## Use (3)

The motif scan is complementary to the phosphosite.py analysis and allowing the researcher to expand the results obtained from the known kinase-substrate analysis performed with phosphosite.py and the additional analyses performed using the outputted lists of sites and sequences (see script 1). The sequence logo generator allows for a graphic representation of the motifs that may be enriched.

---

## <center>Function and sequence (4)</center>

## Script name (4)

`function_seq.py`

## Objective (4)

So far we have performed functional analysis at the phosphosite level but also we would want to know as much as possible about the proteins that contain the phosphosites. In this regard people frequently turn to UniProt as it is a protein-centric database. For programmatic database access UniProt uses RESTful URLs or can be accessed using packages in Biopython.

[UniProt page](https://www.uniprot.org/ "Title")

[UniProt programmatic access page](https://www.uniprot.org/help/programmatic_access "Title")

[Bio.SwissProt page](http://biopython.org/DIST/docs/api/Bio.SwissProt.Record-class.html "Title")

## Input (4)

**Input:**  A Uniprot accession (e.g., P68400) or identifier (e.g., CSK21_HUMAN). The input is isoform sensitive e.g., P68400-1 is different from P68400-2

**Content:** -

**Path:**  -

**Parameters:**  HERE

## Output (4)

1. Protein sequence in fasta format
2. Lists protein functional information.

## Use (4)

The proteins functional information can be used to understand the functions of the proteins being targeted for example one can extract the UniProt keywords and for instance generate a word-cloud to visualize the frequent termns found. 

---

## <center>Output files</center>

Brief description of the content of each input and output file.

Input files in ./data:

| File | Description |
| ------------- | ------------- |
| Phospho(STY)Sites.txt  | Phosphoproteomics experiment result  |
| Kinase\_Substrate_Dataset.txt  | Known kinase-susbtarte relationships from PhosphoSitePlus  |
| HPRD\_motifs_Jul2018.tsv  | HPRD motifs from HPRD database  |

Output files in ./data:

| File | Description |
| ------------- | ------------- |
| Phospho(STY)Sites_filtered.txt  | Input file filtered from contaminats and reverse hits  |

Output files in ./results:

| File | Description |
| ------------- | ------------- |
| header_info.tsv  | Column names  |
| raw_files.tsv  | Raw files name  |
| contaminants.tsv  | Protein name of contaminants  |
| reverse_hits.tsv  | Reverse hits (FDR control)  |
| gene_names.tsv  | Gene names of proteins identified  |
| proteins.tsv  | UniProt accession of proteins identified  |
| modifications.tsv  | Type of PTMs found in experiment  |
| sequences.tsv  | Sequence window of 31 amino acids centered at the modified residue |
| sequences.fasta | Fasta sequence of proteins identified |
| sequence6.tsv  | Sequence window of 13 amino acids centered at the modified residue |
| sequence7.tsv  | Sequence window of 15 amino acids centered at the modified residue |
| sequence6.svg  | Weblogo for sequence window of 13 amino acids centered at the modified residue |
| sequence7.svg  | Weblogo for sequence window of 15 amino acids centered at the modified residue |
| ck2\_hit_sequence6.tsv  | Same as above but only for CK2 motif matches |
| ck2\_hit_sequence7.tsv  | Same as above but only for CK2 motif matches |
| ck2\_hit_sequence6.svg  | Same as above but only for CK2 motif matches |
| ck2\_hit_sequence7.svg  | Same as above but only for CK2 motif matches |
| HPRD\_motifs\_Jul2018\_as\_regex.tsv  | HPRD motifs as regex |
| kinase\_motif\_scan\_sequence6.tsv  | HPRD motifs scan for sequence6 |
| kinase\_motif\_scan\_sequence7.tsv  | HPRD motifs scan for sequence7 |
| kinase\_substrate.tsv  | known kinase-substrate relationships |
| sites_classical.tsv  | formatted phosphosites |
| sites_iptmnet.tsv  | formatted phosphosites |
| sites_psp.tsv  | formatted phosphosites |
| sites_networkin.tsv  | formatted phosphosites |
| swp_info.tsv  | functional information of proteins identified |

