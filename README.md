# MetaPro Metatranscriptomics Practical Lab

**This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

**This tutorial was produced by Billy Taj (billy.taj@sickkids.ca), Mobolaji Adeolu (adeolum@mcmaster.ca), John Parkinson (john.parkinson@utoronto.ca) & Xuejian Xiong (xuejian@sickkids.ca)**

## Overview

This tutorial will take you through the MetaPro pipeline for processing metatranscriptomic data. The pipeline, developed by the Parkinson lab, consists of various steps which are as follows:

1.  Remove low-quality reads that MetaPro has low confidence in processing.
2.  Remove adapter sequences, which are added during library preparation and sequencing steps, and trim low quality bases and sequencing reads.
3.  Remove duplicate reads to reduce processing time for following steps.
4.  Remove vector contamination (reads derived from cloning vectors, spike-ins, and primers).
5.  Remove host reads (if exploring a microbiome in which the host is an issue).
6.  Remove abundant rRNA sequences which typically dominate metatranscriptomic datasets despite the use of rRNA removal kits.
7.  Add duplicated reads, removed in step 2, back to the data set to improve quality of assemblies.
8.  Classify reads to known taxonomic groups and visualize the taxonomic composition of your dataset.
9.  Assemble the reads into contigs to improve annotation quality.
10.  Annotate reads to known genes.
11.  Map identified genes to the swiss-prot database to identify enzyme function
12.  Generate normalized expression values associated with each gene.
13. Visualize the results using KEGG metabolic pathways as scaffolds in Cytoscape.

The MetaPro metatranscriptomic pipeline includes existing bioinformatic tools and a series of Python scripts that handle the orchestrated invocation of the bioinformatics tools, read-conflict resolution, file format conversion, and output parsing. We will go through these steps to illustrate the complexity of the process and the underlying tools and scripts.  The operation of the pipeline is designed to run all steps automatically in succession.  However, for the purposes of the tutorial, a single-step tutorial-mode has been included in MetaPro's design.  One can now call the MetaPro pipeline to perform each step in isolation, for learning purposes.


New, faster, and/or more accurate tools are being developed all the time, and it is worth bearing in mind that any pipelines need to be flexible to incorporate these tools as they get adopted as standards by the community. For example, over the past two years, our lab has transitioned from Trimmomatic to AdapterRemoval, and from BLAST to DIAMOND.
Note:  This workshop was designed for use with DIAMOND v0.826.  Newer versions of DIAMOND will be incompatible with the pre-compiled database files we have made as part of this exercise.  
To illustrate the process we are going to use sequence reads generated from the contents of the colon of a mouse. These are 150 bp single-end reads. Paired-end reads can also be used, and are often preferred because they can improve annotation quality when there is enough overlap in the read pairs to improve the effective average read length. Working with paired-end data involves an additional data processing step (merging of overlapping reads) produces more files during data processing (files for merged/singleton reads, forward reads, and reverse reads), but the structure of a pipeline for paired-end data is similar to the pipeline described here and can be readily adapted.

Rather than use the entire set of 25 million read, which might take several days to process on a desktop, the tutorial will take you through processing a subset of 100,000 single-ended reads.

**Note**
The purpose of this tutorial is to demonstrate MetaPro's various steps.  The pipeline is fully capable of running all of the steps without the intervention of the user beyond an initial call to the program.  If you wish to simply use MetaPro, do not use the --tutorial option.
This tutorial also assumes that the pipeline files are contained in a directory

## Preliminaries

### Install Docker
MetaPro operates in a containerized environment, controlled by Docker.  The pipeline and its dependent programs exist within the Docker container image.  Please install Docker to continue.

```
https://www.docker.com/products/docker-desktop
```
If you are running the tutorial on a computing cluster environment, your admin may already have Singularity installed.  For the purposes of this tutorial, Docker is equivalent to Singularity.
Next, pull the MetaPro docker image
```
docker pull parkinsonlab/metapro:develop
```
OR
```
singularity pull docker://parkinsonlab/metapro:develop
```

Docker and Singularity maintain different access modes to use their containers.  
1) Scripted-mode: where the user calls software within the container and run. 
2) Interactive-mode: where the user can enter into the container use it like an operating systems

For the purposes of the tutorial, we will be using Docker in interactive mode.
To use Docker in interactive mode:

```
docker run -it -v <a folder in your directory>:<an equivalent folder to mount to in the container instance> <the docker image>
```

an example would be:
```
docker run -it -v C:\Users\Billy\Documents\MetaPro_tutorial:/MetaPro_docker_tutorial parkinsonlab/metapro:develop
```

an equivalent singularity command would be:
```
singularity shell -B <directories to bind-mount> <path to the singularity image.sif>
```

The MetaPro Docker container runs Ubuntu Linux 18.04

### Precomputed Files
MetaPro's tools may take a long time to run if the user does not have the necessary computing resources.  Therefore, we will provide pre-computed files so that the user is not forced to run the software. 

```
wget https://github.com/ParkinsonLab/MetaPro_tutorial/releases/download/1.0/tutorial_files.tar.gz
tar -xzvf tutorial_files.tar.gz
```
 
### Input files

Our data set consists of 150 bp single-end Illumina reads generated from mouse colon contents. To inspect its contents:

```
tar -xvf tutorial_files.tar.gz 
less mouse1.fastq
```

**Notes**:

-   Type `q` to exit `less`.

### Databases, and licenses:
This tutorial relies on a few external databases and libraries to perform the filtering tasks associated with MetaPro.  
We have assembled the smaller databases in our precomputed files package  

[The UniVec Core database](https://ftp.ncbi.nih.gov/pub/UniVec/)  
[A mouse host sequence database](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/)  In this tutorial, we will use one from Ensembl  

There are optional databases that are mentioned in this tutorial.  Due to the size of these references, they are not required to be present for this tutorial, but if one were to use MetaPro, it is highly suggested that they are obtained:

-   [The ChocoPhlan Pangenome Database](http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/)
-   [The NCBI Non-redundant (NR) Protein Database](ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz)
-   [The GB version of the nucleotide accession2taxid table](https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)
-   [The Centrifuge NT database](https://ccb.jhu.edu/software/centrifuge/manual.shtml#nt-database)
    -   To complete the centrifuge database install, the required utilities are placed in /pipeline_tools/centrifuge
-   [The Kaiju Database](https://github.com/bioinformatics-centre/kaiju)
    -   To complete the kaiju database install, the required utilities are placed in /pipeline_tools/kaiju
    -   MetaPro relies on the full database. (makeDB.sh -r)
-   [The NCBI Taxdump database](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
-   [The Swiss-Prot database (fasta)](https://www.uniprot.org/downloads)
-   [The PRIAM Database](http://priam.prabi.fr/REL_JAN18/Distribution.zip)
    
These optional databases require indexing prior to use.  
-   MetaPro requires 2 versions of ChocoPhlan:
    -   One with all of sequences in a single file (for BWA)
    -   One with all of the sequences separated. (for BLAT)


-   the commands used to build the indexed databases are as follows (You don't need to do these!)
    To unpack chocphlan and combine all of the sequences:
    -   `tar -xvf chocophlan.tar.gz && cd chocophlan`
    -   `for i in $(ls | grep ".gz"); do gunzip $i; done'
    -   `for i in $(ls | grep ".m8"); do cat $i >> chocophlan_full.fasta`
    -   `mv chocophlan_full.fasta ..`
     
    To index:    
        -   To prepare ChocoPhlAn for BWA:
            -   `bwa index -a bwtsw <path to chocophlan_full.fasta>`
            -   `samtools faidx <path to chocophlan_full.fasta>`
        -   To prepare NR for DIAMOND:        
            -   `diamond makedb -p 8 --in <path to nr> -d <path to nr>`
        -   To prepare Kaiju:
            -   `/pipeline_tools/kaiju/makeDB.sh -r <a suitable destination for the Kaiju DB>`
        -   To 
            

Additionally, a license from MetaGeneMark is required to run to the contig assembly step
-   [MetaGeneMark](http://exon.gatech.edu/Genemark/license_download.cgi)

### Configuration file:
MetaPro controls many of its features with a Configuration file, a copy has been provided for you in the precomputed files, but it needs to be altered to include the path to the databases.


```
[Databases]
database_path: /project/j/jparkin/Lab_Databases #(a shortcut to help lay out paths)
UniVec_Core: %(database_path)s/univec_core/UniVec_Core.fasta #(the UniVec core database)
Adapter: %(database_path)s/Trimmomatic_adapters/TruSeq3-PE-2.fa #(The adapters database)
Host: %(database_path)s/Mouse_cds/Mouse_cds.fasta #(The host database)
Rfam: %(database_path)s/Rfam/Rfam.cm #(The Infernal rRNA database)
DNA_DB: %(database_path)s/ChocoPhlAn/ChocoPhlAn.fasta #(The BWA database.  It assumes the index is in the same directory)
DNA_DB_Split: %(database_path)s/ChocoPhlAn/ChocoPhlAn_split/ #(The split database, for BLAT)
Prot_DB: %(database_path)s/nr_01_22_2020/nr #(it assumes there's a file named "nr.dmnd")
Prot_DB_reads: %(database_path)s/nr_01_22_2020/nr #(The NR database file which contains the sequences)
accession2taxid: %(database_path)s/accession2taxid/accession2taxid #(the GB version of the nucleotide accession2taxid table)
WEVOTEDB: %(database_path)s/WEVOTE_db/ #(refer to the unpacked taxdump.tar.gz)
nodes: %(database_path)s/WEVOTE_db/nodes.dmp #(a file within the NCBI taxdump package)
names: %(database_path)s/WEVOTE_db/names.dmp #(a file within the NCBI taxdump package)
Kaiju_db: %(database_path)s/kaiju_db/kaiju_db_nr.fmi #(This one is made from indexing the NCBI NR database)
Centrifuge_db: %(database_path)s/centrifuge_db/nt #(The Centrifuge NT database)
SWISS_PROT: %(database_path)s/swiss_prot_db/swiss_prot_db #(diamond-indexed swissprot sequences)
SWISS_PROT_map: /pipeline/custom_databases/SwissProt_EC_Mapping.tsv (This is in the container: /pipeline/custom_databases)
PriamDB: %(database_path)s/PRIAM_db/ (The PRIAM Database)
DetectDB: %(database_path)s/DETECTv2 (
EC_pathway: /pipeline/custom_databases/EC_pathway.txt #(This is in the container:  /pipeline/custom_databases) 
path_to_superpath: /pipeline/custom_databases/pathway_to_superpathway.csv #(This is in the container:  /pipeline/custom_databases)
MetaGeneMark_model: /pipeline_tools/mgm/MetaGeneMark_v1.mod #(This is in the container already.  do not alter)
```

### Checking read quality with FastQC

```
/pipeline_tools/FastQC/fastqc mouse1.fastq
unzip mouse1_fastqc.zip
```

The FastQC report is generated in a HTML file, `mouse1_fastqc.html`. You'll also find a zip file which includes data files used to generate the report.

To read the HTML report file, open it in a web browser (such as FireFox).  Then you can go through the report and find the following information:

-   Basic Statistics: Basic information of the mouse RNA-seq data, e.g. the total number of reads, read length, GC content.
-   Per base sequence quality: An overview of the range of quality values across all bases at each position.
-   Per Base Sequence Content: A plot showing nucleotide bias across sequence length.
-   Adapter Content: Provides information on the level of adapter contamination in your sequence sample.

## Processing the Reads

### Step 1: Remove adapter sequences, trim low quality sequences, and remove duplicate reads. 
For the first step, MetaPro removes adaptor sequences, trims low-quality reads, and removes duplicate reads in one pass.
```
The format is:
read1='<path to input sequence>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality
```
The commands would look like:
```
read1=/home/billy/mpro_tutorial/mouse1.fastq
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality
```


**Notes**:
All MetaPro steps share the same file directory scheme:
- data: where the interim files are placed for each run.  This includes intermediate steps.
- final results: where the end-phase deliverables are placed, assuming the pipeline will continue running.
- All of MetaPro's commands are generated in separate shellscripts in each folder.


In this Quality-filtering stage, MetaPro will do several actions:
-filter reads below a quality score of 75
-filter reads below a minimum length of 30 bp
-remove adapters
-remove duplicate reads within the dataset


<!-- ***Question 1: How many low quality sequences have been removed?*** -->

Checking read quality with FastQC:

```
cd <output_directory>/quality_filter/data/4_quality_filter/
/pipeline_tools/FastQC/fastqc singletons_hq.fastq
unzip singletons_hq_fastqc.zip
```

Compare with the previous report to see changes in the following sections:

-   Basic Statistics
-   Per base sequence quality




**Read quality filtering**

AdapterRemoval, which was used to remove the adapters and trim low quality bases in the reads, uses a sliding window method to remove contigous regions of low quality bases in reads. However, it is worthwhile to impose an overall read quality threshold to ensure that all reads being used in our analyses are of sufficiently error-free. For this we use the tool VSEARCH which can be found at this [website](https://github.com/torognes/vsearch) (when processing paired-end data, this step should come **after** the read merging step):

```
Example command only, MetaPro already automatically perfoms this task
vsearch --fastq_filter mouse1_trim.fastq --fastq_maxee 2.0 --fastqout mouse1_qual.fastq
```



**Notes**:

-   The command line parameters are:
    -   `--fastq_filter ` Instructs VSEARCH to use the quality filtering algorithm to remove low quality reads
    -   `--fastq_maxee 2.0` The expected error threshold. Set at 1. Any reads with quality scores that suggest that the average expected number of errors in the read are greater than 1 will be filtered.
    -   `--fastqout` Indicates the output file contain the quality filtered reads

Checking read quality with FastQC:

```
fastqc mouse1_qual.fastq
firefox mouse1_qual_fastqc.html
```

Compare with the previous reports to see changes in the following sections:

-   Basic Statistics
-   Per base sequence quality
-   Per sequence quality

<!-- ***Question 2: How has the per read sequence quality curve changed?*** -->

### Step 2. Remove duplicate reads

To significantly reduce the amount of computating time required for identification and filtering of rRNA reads, we perform a dereplication step to remove duplicated reads using the software tool CD-HIT which can be obtained from this [website](https://github.com/weizhongli/cdhit).

```
Example command only.  MetaPro already calls this command as part of the Quality-filtering step
/usr/local/prg/cd-hit-v4.6.7-2017-0501/cd-hit-auxtools/cd-hit-dup -i mouse1_qual.fastq -o mouse1_unique.fastq
```

**Notes**:

-   The command line parameters are:
    -   `-i`: The input fasta or fastq file.
    -   `-o`: The output file containing dereplicated sequences, where a unique representative sequence is used to represent each set of sequences with multiple replicates.
-   A second output file `mouse1_unique.fastq.clstr` is created which shows exactly which replicated sequences are represented by each unique sequence in the dereplicated file and a third, empty, output file, `mouse1_unique.fastq2.clstr` is also created which is only used for paired-end reads.

<!-- ***Question 3: Can you find how many unique reads there are?*** -->

While the number of replicated reads in this small dataset is relatively low, with larger datasets, this step can reduce file size by as much as 50-80%

### Step 3. Remove vector contamination

To identify and filter reads from sources of vector, adapter, linker, and primer contamination we use the Burrows Wheeler aligner (BWA) and the BLAST-like alignment tool (BLAT) to search against a database of cow sequences. As a reference database for identifying contaminating vector and adapter sequences we rely on the UniVec\_Core dataset which is a fasta file of known vectors and common sequencing adapters, linkers, and PCR Primers derived from the NCBI Univec Database. 

Now we call MetaPro to perform the vector filtering
```
read1='<path to your quality filter final results mouse.fastq>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector
```
The commands would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/quality_filter/final_results/singletons.fastq
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector
```

MetaPro will automatically run the following:

```
bwa index -a bwtsw UniVec_Core
samtools faidx UniVec_Core
makeblastdb -in UniVec_Core -dbtype nucl

bwa mem -t 4 UniVec_Core mouse1_unique.fastq > mouse1_univec_bwa.sam
samtools view -bS mouse1_univec_bwa.sam > mouse1_univec_bwa.bam
samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam
samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam
vsearch --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta
blat -noHead -minIdentity=90 -minScore=65  UniVec_Core mouse1_univec_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_univec.blatout
python3 /pipeline/Scripts/read_BLAT_Filter_v3.py single high mouse1_univec_bwa.fastq mouse1_univec.blatout mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq
```

**Notes**:

-   The commands to the following tasks:
    -   `bwa index, samtools faidx, and makeblastdb`: Index the UniVec core database for BWA and BLAT 
    -   `bwa mem`: Generates alignments of reads to the vector contaminant database
    -   `samtools view`: Converts the .sam output of bwa into .bam for the following steps
    -   `samtools fastq`: Generates fastq outputs for all reads that mapped to the vector contaminant database (`-F 4`) and all reads that did not map to the vector contaminant database (`-f 4`)

<!-- ***Question 4: Can you find how many reads BWA mapped to the vector database?*** -->

Now we want to perform additional alignments for the reads with BLAT to filter out any remaining reads that align to our vector contamination database. However, BLAT only accepts fasta files so we have to convert our reads from fastq to fasta. This can be done using VSEARCH.

```
Example call.  MetaPro automatically runs this
vsearch --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta

```

**Notes**:

-   The VSEARCH command used, `--fastq_filter`, is the same as the command used to filter low quality reads in Step 1. However, here we give no filter criteria so all input reads are passed to the output fasta file.

Now we can use BLAT to perform additional alignments for the reads against our vector contamination database.

```
Example call.  MetaPro automatically runs this

```

**Notes**:

-   The command line parameters are:
    -   `-noHead`: Suppresses .psl header (so it's just a tab-separated file).
    -   `-minIdentity`: Sets minimum sequence identity is 90%.
    -   `-minScore`: Sets minimum score is 65. This is the matches minus the mismatches minus some sort of gap penalty.
    -   `-fine`: For high-quality mRNAs.
    -   `-q`: Query type is RNA sequence.
    -   `-t`: Database type is DNA sequence.

Lastly, a python script is used to filter the reads that BLAT does not confidently align to any sequences from our vector contamination database.

```
python3 /pipeline/Scripts/read_BLAT_Filter_v3.py single high mouse1_univec_bwa.fastq mouse1_univec.blatout mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq

```

**Notes**:

The argument structure for this script is:
`read_BLAT_Filter_v3.py <operating mode: either "single" or "paired"> <filter stringency.  to handle paired-read conflicts.  "low" or "high"> <Input_Reads.fq> <BLAT_Output_File> <Unmapped_Reads_Output> <Mapped_Reads_Output>`

Here, BLAT does not identify any additional sequences which align to the vector contaminant database. However, we have found that BLAT is often able find alignments not identified by BWA, particularly when searching against a database consisting of whole genomes.

some alignments to vector contaminants missed by BWA in large multi-million read datasets.
In handling, paired-ended data, some alignments may yield paired-data to be broken.  To handle this, the filter stringency option decides how to resolve this discrepancy.  
Low filter stringency will only remove reads where both pairs aligned to a vector.
High filter stringency will remove reads where either pair aligned to a vector.


### Step 4. Remove host reads

To identify and filter host reads (here, reads of mouse origin) we repeat the steps above using a database of mouse DNA sequences. For our purposes we use a [mouse genome database](ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz) downloaded from Ensembl.

The MetaPro call is:
```
read1='<path to your vector filter final results mouse.fastq>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host
```

The commands would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/vector_read_filter/final_results/singletons.fastq
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host
```

This call will perform the following steps:
- Prepare the host database for alignment (BWA + BLAT)
- perform alignment using BWA
- convert the unaligned reads from BWA to a format for BLAT
- perform alignment of the unaligned reads using BLAT
- Run a script to remove the host reads from the input sample 


```
bwa index -a bwtsw mouse_cds.fa
samtools faidx mouse_cds.fa
makeblastdb -in mouse_cds.fa -dbtype nucl
bwa mem -t 4 mouse_cds.fa mouse1_univec_blat.fastq > mouse1_mouse_bwa.sam
samtools view -bS mouse1_mouse_bwa.sam > mouse1_mouse_bwa.bam
samtools fastq -n -F 4 -0 mouse1_mouse_bwa_contaminats.fastq mouse1_mouse_bwa.bam
samtools fastq -n -f 4 -0 mouse1_mouse_bwa.fastq mouse1_mouse_bwa.bam
vsearch --fastq_filter mouse1_mouse_bwa.fastq --fastaout mouse1_mouse_bwa.fasta
blat -noHead -minIdentity=90 -minScore=65  mouse_cds.fa mouse1_mouse_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_mouse.blatout
./1_BLAT_Filter.py mouse1_mouse_bwa.fastq mouse1_mouse.blatout mouse1_mouse_blat.fastq mouse1_mouse_blat_contaminats.fastq
```

<!-- ***Question 5: How many reads did BWA and BLAT align to the mouse host sequence database?*** -->

***Optional:*** In your own future analyses you can choose to complete steps 3 and 4 simultaneously by combining the vector contamination database and the host sequence database using `cat UniVec_Core mouse_cds.fa > contaminants.fa`. However, doing these steps together makes it difficult to tell how much of your reads came specifically from your host organism.

### Step 5. Remove abundant rRNA sequences

rRNA genes tend to be highly expressed in all samples and must therefore be screened out to avoid lengthy downstream processing times for the assembly and annotation steps. MetaPro uses [Barrnap] (https://github.com/tseemann/barrnap) and [Infernal] (http://infernal.janelia.org/).
You could use sequence similarity tools such as BWA or BLAST for this step, but we find Infernal, albeit slower, is more sensitive as it relies on a database of covariance models (CMs) describing rRNA sequence profiles based on the Rfam database. Due to the reliance on CMs, Infernal, can take as much as 4 hours for ~100,000 reads on a single core.  In an effort to shrink the computing time, we leverage a computing cluster's multiple cores.
Here, MetaPro demonstrates the case for automation.
MetaPro subdivides the input data, coordinates the concurrent processes, and collects the results into one single file after all of the scanning has been complete.


MetaPro will perform the following:
- subdivide the input data in user-defined chunksizes.
- Each chunk is then run independently:
- - run each chunk through Barrnap
- - using the results of Barrnap, filter the data chunk into mRNA, and leftover data for further scanning.
- - run each leftover chunk through Infernal.
- - filter the Barrnap leftover chunk using the Infernal results, to get mRNA, and "other"
- collect all of the data pieces (Barrnap mRNA, Infernal mRNA) into mRNA, and "other"

By running things this way, the rRNA step takes 4 minutes(as recorded with a 40-core computing node with 200 GB RAM, and an rRNA chunksize of 1000 reads), but it requires significant computing power, memory, and storage space, not available on a typical desktop PC.
If you were to run this on your own, you will need the RFam database.
The call to MetaPro would be:
```
-NOTE: example form.

read1='<path to your host filter final results mouse.fastq>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA
```

The command would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/host_read_filter/final_results/singletons.fastq
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA
```

instead, we have provided you with the results.
``` 
tutorial_files/mouse1_run/rRNA_filter/final_results.
```

Here, we only remove a few thousand reads than map to rRNA, but in some datasets rRNA may represent up to 80% of the sequenced reads.

<!-- ***Question 6: How many rRNA sequences were identified? How many reads are now remaining?*** -->


### Step 6. Rereplication / duplicate repopulation

After removing contaminants, host sequences, and rRNA, we need to replace the previously removed replicate reads back in our data set.

```
read1='<path to your rRNA filter final results mouse.fastq>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop
```

The command would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/rRNA_filter/final_results/mRNA/singletons.fastq
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop
```

Now that we have filtered vectors, adapters, linkers, primers, host sequences, and rRNA, check read quality with FastQC:

```
fastqc mouse1_mRNA.fastq
firefox mouse1_mRNA_fastqc.html
```

<!-- ***Question 8: How many total contaminant, host, and rRNA reads were filtered out?*** -->


### Step 7. Contig assembly

We have now gathered the putative mRNA scripts, and here we assemble the mRNA into contigs. 
Previous studies have shown that assembling reads into larger contigs significantly increases our ability to annotate them to known genes through sequence similarity searches. Here we will apply the SPAdes genome assemblers' transcript assembly algorithm to our set of putative mRNA reads.
The typical workload of MetaPro ranges anywhere from 40-million to 100-million reads, and more.  To make this process as efficient as possible, MetaPro assembles the data into contigs in an effort to shrink the number of reads to annotate.  

```
Example only: do not run.
read1='<path to your rereplicated mRNA final results mouse.fastq>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs
```
The command would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/duplicate_repopulation/final_results/singletons.fastq
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs
```

The expected output will be in:
```
tutorial_files/mouse1_run/assemble_contigs/final_results
```


**Notes**:
In this step, MetaPro does the following:
-   Assemble the reads into contigs using SPAdes
-   Use MetaGeneMark to predict the genes in these contigs
-   Use BWA to align the mRNA reads against these split-contigs to find out which reads were consumed by the process.
-   Produce a relational map of split-contig and their constituent reads.

[SPAdes](https://cab.spbu.ru/software/spades/) assembles long contigs, but MetaPro requires that each contig only represent 1 gene.  Thus the need to disassemble them, using [MetaGeneMark](http://exon.gatech.edu/Genemark/meta_gmhmmp.cgi)  MetaGeneMark requires the user to register and obtain a free license.  Thus, we have provided the results in the precomputed files package.

<!--
***Question 9: How many assemblies did SPAdes produce?  
Hint: try using the command`tail mouse1_contigs.fasta`***
-->
<!--
***Question 10: How many reads were not used in contig assembly? How many reads were used in contig assembly? How many contigs did we generate?***
-->



### Step 8. Annotate reads to known genes/proteins

Here we will attempt to infer the specific genes our putative mRNA reads originated from. In our pipeline we rely on a tiered set of sequence similarity searches of decreasing accuracy - BWA, BLAT, and DIAMOND. While BWA provides high stringency, sequence diversity that occurs at the nucleotide level results in few matches observed for these processes. Nonetheless it is quick. To avoid the problems of diversity that occur at the level of nucleotide, particularly in the absence of reference microbial genomes, we use a cascaded method involving two other tools: BLAT, and DIAMOND. BLAT provides a more sensitive alignment, along with quality scores to rank the matches.  DIAMOND is used to provide more sensitive peptide-based searches, which are less prone to sequence changes between strains.

Since BWA and BLAT utilize nucleotide searches, we rely on the [ChocoPhlan pangenome database] (http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/) that we obtained from The Huttenhower lab, which contains over 10000 organisms in separate .ffn files.  We create a merged copy of these sequences, and index it for BWA to use.  We leave it in its separated state for BLAT to use.  

For DIAMOND searches we use the [Non-Redundant (NR) protein database] (ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz) from the NCBI.

This is a computationally intensive step.  We employ our subdivision strategy here, similar to our design for rRNA removal, as seen below
-   The mRNA read data (contigs, remaining singletons, and remaining paired reads if applicable) are split into chunksizes (GA_chunksize in the configuration)
-   Each chunk is sent through BWA to be aligned against the ChocoPhlan database  
-   A map of genes to constituent reads is formed for each chunk.
-   Reads not annotated by BWA are isolated, and are sent through to BLAT.
-   Another map of genes to constituent BLAT-annotated reads are formed for each chunk.
-   Reads not annotated by BLAT are isolated, and are sent through to DIAMOND.
-   A 3rd batch of maps of annotations to constituent reads are formed for the DIAMOND reads for each chunk.

In all, this leaves us with (split into chunks):
-   a batch of BWA-annotated gene-to-read maps
-   a batch of BLAT-annotated gene-to-read maps
-   a batch of DIAMOND-annotated protein-to-read maps
-   a batch of reads unannotated by BWA, BLAT, and DIAMOND.

These batches of files are then sent to a custom script that will perform all of the final merging:
-   collect and merge every map to form one map.  
-   collect and merge the unannotated DIAMOND reads
-   collect all of the genes that were found in the reads, and convert them into proteins.  Then merge them with the proteins found in DIAMOND.  This step is for downstream analysis.

The command is:
```
Example only: do not run. Heavily computationally intensive
read1='<path to your unassembled singletons.fastq>'
contig='<path to your contigs.fasta>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA
```

The command would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/home/billy/mpro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA
```

The outputs are located in:
```
tutorial_files/GA_FINAL_MERGE/final_results
```


**Notes**:
-   The gene annotation step will create 4 different subdirectories: GA_BWA, GA_BLAT, GA_DIAMOND, and GA_FINAL_MERGE
-   MetaPro will take the top-quality hits of each tool to count towards annotation:
    -   In BWA: the CIGAR string is decoded.  Any hit with less than 90% match is rejected 
    -   In BLAT and DIAMOND: the only reads that pass annotation are reads with all 3 conditions satisfied:
        -   A sequence identity score of greater than 85
        -   An alignment length of greater than 65%
        -   A bitscore of over 60
-   MetaPro makes a number of extra considerations to account for multi-mapped reads:
    -   Every read scanned by BWA, BLAT, and DIAMOND is affixed with a quality score suffix taken from the aligner's report
    -   In BWA: the alignment score is used.
    -   In BLAT and DIAMOND: the bitscore is used.
    -   Should there be a case where a read is annotated to multiple genes or proteins, each tool will use their quality scores to rank the hit
        -   Alignment information is iterated through, from the top of the file down to the bottom.
        -   The hit with the higher score is used for all cases.
        -   If the scores are tied, then the incumbent hit is used.
    -   This read-ranking is done on each chunk when the gene-to-read maps are formed.
    -   There are further considerations made for paired-end annotation:
        -   MetaPro views paired-end data as 2 copies of the same read
        -   If there are disagreements between the forward and reverse read's annotation, the quality scores are used to resolve the conflict.
        -   If both reads agree on the same gene or protein, the read is counted only once. 
        -   This paired-read conflict resolution is performed in the GA_FINAL_MERGE step.
    
-   Unless you are running this tutorial on a computing cluster, most systems do not have enough memory to handle indexing or searching large databases like `ChocoPhlan` (19GB) and `nr` (>60GB). The descriptions in this section are purely for your information. Please use our precomputed gene, protein, and read mapping files from the tar file `tar -xzf tutorial_files.tar.gz`

### Step 9. Taxonomic Classification

Now that we have putative mRNA transcripts, we can begin to infer the origins of our mRNA reads. Firstly, we will attempt to use a reference based short read classifier to infer the taxonomic orgin of our reads. Here we will use [Kaiju] (https://github.com/bioinformatics-centre/kaiju), [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml), and our Gene Annotation results to generate taxonomic classifications for our reads based on a reference database. 
Kaiju can classify prokaryotic reads at speeds of millions of reads per minute using the proGenomes database on a system with less than 16GB of RAM (~13GB). Using the entire NCBI nr database as a reference takes ~43GB. Similarly fast classification tools require >100GB of RAM to classify reads against large databases. 

However, Kaiju still takes too much memory for the systems in the workshop so we have precompiled the classifications, `mouse1_classification.tsv`, in the tar file `tutorial_files.tar.gz`.
Centrifuge is a lightweight rapid microbial classification engine.  It uses methods similar to BWA and the Ferrgina-Manzini (FM) index  to make quick work of assigning taxomony.

The ChocoPhlan Pangenome Database contains taxonomic information that MetaPro extracts.  Kaiju, Centrifuge, and the extracted taxa are combined using [WEVOTE](https://github.com/aametwally/WEVOTE).  WEVOTE is the Weighted Voting Taxonomic Identification system.  It performs consensus merging of various taxa results and reconciles the taxa identification from various sources.  
MetaPro uses this to settle on one confident taxon amongst Kaiju, Centrifuge, and the ChocoPhlan database choices.

The command to run this step is:
```
Example only: do not run.  Also MetaPro assumes the Gene Annotation step has completed
read1='<path to your unassembled singletons.fastq>'
contig='<path to your contigs.fasta>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA
```
The command would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/home/billy/mpro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA
```



Instead, we have provided the results here:
```
tutorial_files/mouse1_run/taxonomic_annotation/final_results
```

We can use [Krona] (https://github.com/marbl/Krona/wiki) to generate a hierarchical multi-layered pie chart summary of the taxonomic composition of our dataset.

To use Krona, the export of MetaPro's taxonomic annotations need to be appended

```
python3 /pipeline/Scripts/alter_taxa_for_krona.py tutorial_files/taxonomic_classifications.tsv tutorial_files/mouse1_classification.tsv
/pipeline_tools/kaiju/kaiju2krona -t nodes.dmp -n names.dmp -i mouse1_classification.tsv -o mouse1_classification_Krona.txt
/pipeline_tools/KronaTools/scripts/ImportText.pl -o mouse1_classification.html mouse1_classification_Krona.txt
```

We can then view this pie chart representation of our dataset using a web browser

<!-- 
***Question 10: What is the most abundant family in our dataset? What is the most abundant phylum?  
Hint: Try decreasing the `Max depth` value on the top left of the screen and/or double clicking on spcific taxa.***
-->



### Step 10. Enzyme Function Annotation

To help interpret our metatranscriptomic datasets from a functional perspective, we rely on mapping our data to functional networks such as metabolic pathways and maps of protein complexes. Here we will use the KEGG carbohydrate metabolism pathway.

MetaPro uses 3 tools to produce its enzyme annotations: [PRIAM](http://priam.prabi.fr/), DIAMOND with the SwissProt Database, and [Detect](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq266)


We need to match our annotated genes the enzymes in the KEGG pathway. To do this, we will use DIAMOND to identify homologs of our genes/proteins from the SWISS-PROT database that have assigned enzyme functions. Diamond is a relatively coarse and straight forward way to annotate enzyme function by homology. 
We also use more robust methods for enzymatic function annotation, such as our own probability density based enzyme function annotation tool, DETECT, as well as PRIAM 
MetaPro combines the predictions of all 3 tools to give 2 answers, a lower-confidence set of enzyme predictions, and a higher-confidence set of predictions.
PRIAM is incredibly resource-intensive, and slow to run.  For the sake of brevity, the results have been provided for this tutorial.


The command is:
```
Example only: do not run.  Also, this command assumes that MetaPro has performed the Gene annotation phase
read1='<path to your unassembled singletons.fastq>'
contig='<path to your contigs.fasta>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC
```

The command would look like:
```
read1=/home/billy/mpro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/home/billy/mpro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
config=/home/billy/mpro_tutorial/config_mouse.ini
output=/home/billy/mpro_tutorial/mouse1_run
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC
```

The results are provided in:
```
tutorial_files/enzyme_annotation/final_results
```

**Notes**:
-   MetaPro's high-confidence and low-confidence are determined by the following:
    -   DIAMOND: Low-confidence hits are ones with an e-value of 1e-5 or smaller.  High-confidence hits are ones with an e-value of 1e-10 or smaller
    -   PRIAM: Low-confidence hits are ones with e-values lower than 1e-5.  High-confidence hits are ones where their probability value is 0.5 or higher
    -   DETECT: There is no separation.
    
-   MetaPro reconciles the 3 annotations in the following manner:
    -   The Enzymes predicted by DETECT are taken, followed by the annotations that agree between PRIAM and DIAMOND
    -   In the event of multiple enzymes being annotated to the same protein:
        -   Every enzyme annotation comes with a probability score.
        -   MetaPro includes a enzyme co-occurence database (compiled from the ENZYME database) that contains pairs of enzymes known to exist together
        -   Using this co-occurence database, MetaPro filters invalid predictions.
            -   In cases where more-than-2 enzymes are annotated to a protein, the top 2 enzymes are taken, based on the probability score.
            -   If a pair of enzymes do not exist in the database, the enzyme with the higher probability score is declared the proper annotation.
            -   Otherwise, the annotation is declared as a pair of enzymes.

<!--
***Question 14: How many high-confidence unique enzyme functions were identified in our dataset?***
-->
### Step 11. Generate annotation outputs

We have removed low quality bases/reads, vectors, adapters, linkers, primers, host sequences, and rRNA sequences and annotated reads to the best of our ability - now lets summarize our findings. We do this by looking at the relative expression of each of our genes in our microbiome.

MetaPro generates many output files:
-   A gene expression table of RPKMs in terms of the 20-most prevalent taxa in the sample
-   A Cytoscape-compatible network file
-   An enzyme superpathway heatmap to visualize the distribution of enzymes found.
-   A read count summary table, to track where the reads went.
-   A gene/protein-to-read map of all genes and proteins identified by MetaPro, followed by its constituent reads
-   A histogram of read quality
-   A summary of all taxa identified, followed by the number of reads associated with those taxa.

The command to generate the outputs is:
```
Example only: do not run.  Also, this command assumes that MetaPro has performed the gene, taxa, and enzyme annotations.
read1='<path to your unassembled singletons.fastq>'
contig='<path to your contigs.fasta>'
config='<path to config file>'
output='<path to output folder>'
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial output
```



**Notes**:

<!--
***Question 15: have a look at the `mouse1_RPKM.txt` file. What are the most highly expressed genes? Which phylum appears most active?***
-->
### Step 12. Visualize the results using a KEGG Pathway as a scaffold in Cytoscape.

To visualize our processed microbiome dataset in the context of the carbohydrate metabolism pathways, we use the network visualization tool - Cytoscape together with the enhancedGraphics and KEGGscape plugins. Some useful commands for loading in networks, node attributes and changing visual properties are provided below (there are many cytoscape tutorials available online).
Note: Use cytoscape 3.7.2, instead of the latest version
```
https://github.com/cytoscape/cytoscape/releases/3.7.2/
```


**Download the metabolic pathway**

First, download the carbohydrate metabolism pathways from KEGG using the following commands:

```
wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml
```

You can find other [pathways on KEGG] (http://www.genome.jp/kegg-bin/get_htext?htext=br08901.keg) which can also be imported into Cytoscape by selecting the `Download KGML` option on the top of the page for each pathway.

**Install the Cytoscape plugins**

-   Select `Apps` -> `App Manager`
-   Search for `enhancedGraphics`
-   Select `enhancedGraphics` in the middle column then click `Install` in the bottom right
-   Search for `KEGGScape`
-   Select `KEGGScape` in the middle column then click `Install` in the bottom right

**Import an XML from KEGG into Cytoscape**

-   Select `File` -> `Import` -> `Network` -> `File...`
-   Select the XML file, `ec00010.xml` or `ec00500.xml` and click `Open`
-   Check `Import pathway details from KEGG Database` box then select `OK`

**Loading a node attribute text file (.txt) - this will map attributes to nodes in your network which you can subsequently visualize**

-   Select `File` -> `Import` -> `Table` -> `File...`
-   Select the `mouse1_cytoscape.txt` file and click `Open`
-   Change the `Key Column for network` from `shared name` to `KEGG_NODE_LABEL`
-   Click OK

**Visualizing your node attributes**

-   In the left `Control Panel` select the `Style` tab
-   Check the `Lock node width and height` box
-   Click the left-most box by the `Size` panel and change the default node size to 20.0
-   Click the blank box immediately to the right of the box you clicked to change the default size, change the `Column` field to `RPKM` and the `Mapping Type` field to `Continuous Mapping`
-   Click the left-most box by the `Image/Chart 1` panel, switch to the `Charts` tab, Click the doughnut ring icon, and press the `>>` "add all" button between the two column fields before clicking apply (make sure to remove overall RPKM from the fields that are added to the doughnut ring)
-   If you do not see the `Image/Chart 1` panel, select `Properties` -> `Paint` -> `Custom Paint 1` -> `Image/Chart 1` from the to left corner of the control panel
-   To improve the visualization you can modify colour properties under `Image/Chart 1` -> `Charts` -> `Options`, or modify other properties such as Label Font Size, Label Position, Fill Color, Node location, and edge properties

**Notes**:

-   A cytoscape file with node attributes precalculated is provided for your convenience, `tar -xzf tutorial_files.tar.gz Example.cys`, feel free to open it and play with different visualizations and different layouts - compare the circular layouts with the spring embedded layouts for example. If you want to go back to the original layout then you will have to reload the file.
-   Cytoscape can be temperamental. If you don't see pie charts for the nodes, they appear as blank circles, you can show these manually. Under the 'properties' panel on the left, there is an entry labeled 'Custom Graphics 1'. Double click the empty box on the left (this is for default behavior) - this will pop up a new window with a choice of 'Images' 'Charts' and 'Gradients' - select 'Charts', choose the chart type you want (pie chart or donut for example) and select the different bacterial taxa by moving them from "Available Columns" to "Selected Columns". Finally click on 'Apply' in bottom right of window (may not be visible until you move the window).

**Visualization Questions:**

- Which genes are most highly expressed in these two systems?
- Which taxa are responsible for most gene expression?
- Can you identify sub-systems (groups of interacting genes) that display anomalous taxonomic profiles?
- Think about how you might interpret these findings; for example are certain taxa responsible for a specific set of genes that operate together to fulfill a key function?
- Can you use the gene annotations to identify the functions of these genes through online searches?
- Think about the implications of sequence homology searches, what may be some caveats associated with interpreting these datasets?
