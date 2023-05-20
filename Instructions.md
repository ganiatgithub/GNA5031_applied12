# GNA5031 Applied Session 11

**Case study on the antibiotic resistome**

## Learning Objectives:

At the conclusion of this session, students should be able to:

•	Examine different genomic sequence files and understand their content
•	Query an antimicrobial resistance gene (ARG) database to identify likely ARGs within genomes and plasmids
•	Use BLAST to identify an organism by its 16S rRNA genome sequence
•	Explain the benefits and implications of searching chromosomes compared to plasmids for ARGs

Software required: `DIAMOND`, `seqkit`, text editing software (`BBedit` for Mac, `Notepad++` for PC), `Microsoft excel`
 
 
## Introduction

The antibiotic resistome refers to the collection of all antimicrobial resistance genes (ARGs) in both pathogenic and non-pathogenic bacteria. To understand the resistome of a particular environment or setting, metagenomic sequencing can be used to examine the genomes of all microbial community members and profile the ARGs they contain. However, depending on the depth of sequencing and the success of metagenome assembly efforts, with this approach it may be challenging to link ARGs with the organisms that have them, to verify whether the presence of such genes corresponds to an antibiotic resistance phenotype, and to understand transmission dynamics. 

In this case scenario, a collaborator has been studying antimicrobial resistant bacteria in agricultural runoff from a dairy farm. Recent infections in the cows have not resolved despite the farmer’s use of antibiotics, and there is concern that the farm’s runoff may be a hotspot for antimicrobial resistance that may be entering nearby waterways.
Using a range of antibiotic selection agar plates to culture bacteria from the runoff, your collaborator has isolated a number of bacterial strains that grow in the presence of antibiotics. Three of these isolates exhibited high levels of resistance and your collaborator wishes to investigate them further. They have done some initial analyses and provided you with three files for each isolate:

1.	They initially conducted 16S rRNA sequencing of each isolate they found, to taxonomically identify them. They have provided you with the 16S rRNA sequence for each of the three concerning isolates in nucleotide format. `organismx_16S.fna`
2.	They have sequenced the genomes of these isolates and have provided you with the set of predicted protein sequences from each genome. `organismx_genomic.faa`
3.	To examine plasmid-borne resistance, they have isolated and purified plasmids from these isolates, and provided their sequences in nucleotide format. `organismx_plasmid.fna`

In this workshop you will analyse this data to identify the ARGs in these isolates that are likely to confer their resistance phenotypes, and consider the implications of this resistance, so that you may provide direction to your collaborator.

## Part 1: Identify the isolates and search their genomes for ARGs
Let’s start by checking the files you have been given. Use head <file>, less <file> and seqkit stats <file>, replacing <file> with the filename, to understand their contents.
It would be helpful to know which organisms we’re dealing with, so let’s check those 16S rRNA sequences to identify what these organisms are.
Head to NCBI BLAST and select BLASTn. Copy the sequence into the query box, and change the search to BLAST against rRNA/ITS databases – specifically 16S ribosomal RNA sequences. You can leave the rest of the options at the default settings, and run the search.
Examine the results. 

### Login to virtual machine
[Instructions in detail](https://docs.google.com/document/d/1WBYDpS5utSvHylmRrgBNzXmZwQ74kTkv/edit)

```
ssh -x gnii0001@gna5031s1-gnii0001-01.rep.monash.edu # currently na
# next type in your passwords
```

### Obtain data and software

```
conda activate gna5031

conda install -c bioconda diamond

diamond help # test if diamond has been installed

conda install -c bioconda seqkit
seqkit -v # test if seqkit has been instaled

git clone https://github.com/ganiatgithub/GNA5031_applied12.git # obtain all information needed for this session.
```

### Exercise 1
1.	How many protein sequences do you have for each genome, and how large are the plasmids?
<model answer>

2.	What is the taxonomic identity of each of your isolates based on their 16S rRNA gene sequence? How similar is it to known reference strains?
<model answer>

Now we know what species these isolates are, we want to search their genes to determine if any of them are ARGs. For this we use DIAMOND, a sequence alignment tool which works similarly to BLAST but runs faster when searching a lot of sequences.

A copy of the Comprehensive Antibiotic Resistance Database (CARD) has been provided to you. To use it with DIAMOND, we need to make a DIAMOND database – a version of the database that can be recognised and searched by DIAMOND.

```
diamond makedb --in protein_fasta_protein_homolog_model.fasta -d CARD.dmnd
```

Then, we can run DIAMOND using the new CARD file as a database. We are using the blastp function from DIAMOND, because we are comparing a protein query (the proteins from the isolate) with a protein database (CARD). So we don’t have to run it three times, we’ll run it in a loop on all three isolate files:
```
for isolate in *.faa
do
name=”$(basename -- $isolate | sed ‘s/.faa//’)”
diamond blastp --db CARD.dmnd --query $isolate --out “$isolate”_CARD_results.txt –outfmt 6 --max-target-seqs 1 --max-hsps 1 --id 50
done
```
**What’s happening?**

This command begins by running a loop on all three isolate files: for each one (which ends in .faa), do the following:

What it does for each file is define the “name” variable, which is simply the sample name. The basename command strips any directory names that might come before the filename, then sed removes the .faa extension. This leaves us with the isolate name in the $name variable, which we use to name the output files.
Then, the diamond blastp command takes our new CARD.dmnd database and our isolate proteins, and searches those proteins in the database. The --out flag allows us to name the output file (which uses the isolate name, coming from the $name variable we defined above). We add a few additional options to the end to control our output:
•	--outfmt 6: DIAMOND has several options for its output, this is a tabular format 
•	--max-target-seqs 1: DIAMOND will output only one best hit when a protein matches something in the database
•	--max-hsps-1: DIAMOND will output only one best ‘high scoring pair’ per alignment. This prevents matches appearing twice if the query protein happens to align equally well in more than one place on the same database protein.
•	--id 50: This controls the minimum percentage identity between the query and database proteins. For proteins, 50% is considered a generous match to allow ARGs less similar to the reference sequences to be picked up.
The loop is closed with done, and you should have one DIAMOND results file for each isolate, named appropriately.
Have a look at each of these files with cat <file>

### Exercise 2
1.	How many ARGs have been identified in each genome, and how similar are they to known reference ARG sequences?
<model answer>

2.	Summarise the ARGs for each isolate and the class of antibiotics they confer resistance to. Which antibiotics would you recommend your collaborator use to continue the culture of these isolates?
<model answer>

3.	Based on your results, it is clear that the genotype may not match the phenotype. Why might this be the case? How would you suggest your collaborator determine whether the isolates are indeed resistant to these antibiotics? 
<model answer>
 
## Part 2: Identifying plasmid-borne ARGs and making recommendations
During your analysis you notice something strange – even though the [species name] isolate came from antibiotic selection, it doesn’t seem to contain any ARGs. What’s going on? Let’s do some forensic bioinformatics.

First, let’s rule out a simple mistake – perhaps your collaborator got the files mixed up and they’ve given you the protein sequences from the wrong organism. Let’s check some of the proteins from the isolate to confirm that it is what we think it is.

Take the first few proteins in the file:
`head -n 4 isolate.faa`

Copy these and check them in NCBI BLAST. This time, use Protein BLAST and leave all of the settings at their default.
The closest match should correspond with the 16S rRNA sequence file that identifies this isolate as [species name].
Why can’t we find a gene?
It’s likely that the resistant phenotype in this isolate comes from an ARG on the plasmid. Let’s check all the plasmids for ARGs.
Using the CARD database we made, let’s query it with the plasmid sequences the collaborator provided. These are nucleotide sequences (you can see this when you check the file with head), so this time we use the blastx command, which takes a DNA sequence, translates to a protein sequence in all 6 open reading frames, and aligns them to the protein database.

```
for plasmid in *plasmid.fna
do
name=”$(basename -- $plasmid | sed ‘s/.fna//’)”
diamond blastx --db CARD.dmnd --query $plasmid --out “$isolate”_plasmid_results.txt –outfmt 6 --max-target-seqs 1 --max-hsps 1 --id 50
done
```

This command works exactly the same as the blastp command before.
Use less to view the results for each plasmid.
### Exercise 3

1.	What ARGs are present on the plasmids, and what’s the taxonomy associated with the reference sequence that they most closely match?
2.	Why are there inconsistencies between the taxonomy of the isolate and the ARGs found on the plasmids?
3.	If you had metagenomic data in addition to these isolates, how would you make use of both datasets? What further analyses could be done to assist your collaborator in understanding where the antimicrobial resistance is coming from, and how they may recommend the farmer to approach treatment for the animals?
 
 
# Misc
[where to look for plasmids and genomes](https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/Pseudomonas%20aeruginosa)
[Plasmid sequence of Pseudomonas aeruginosa HS18-89](https://www.ncbi.nlm.nih.gov/nuccore/CP084322.1?report=fasta)
[Genbank site for Pseudomonas aeruginosa HS18-89](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/751/785/GCA_028751785.1_ASM2875178v1/)