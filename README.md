# 1 About
## 1.1 Motivation
The prediction and prioritization of large-scale CNV harmfulness is a
predominant problem in medical genomics.  Our method aims to automate
the logic used by a molecular geneticist to simplify harmful variant
identification.  To our knowledge, our method is the first to combine
functional context and phenotype to discover clinically harmful genes
(and CNVs) in patients with a variety of disorders.

## 1.2 Results
We have compared several feature and gene weighing systems along with
different machine learning approaches at the gene and CNV levels of
classification.  When combining the best performing methodologies and
parameters, we achieved an accuracy of 94%, with 91% precision and 98%
recall on our dataset of Agilent CGH 180k Microarray CNVs.

## 1.3 Availability
Our software will be made freely available at [https://github.com/compbio-UofT/cnv-prioritization](https://github.com/compbio-UofT/cnv-prioritization).

# 2 Overview
1. Data
    1. Genomic positions of exons
    2. CNVs
    3. Gene-gene interaction network
    4. Gene Ontology
    5. Human Phenotype Ontology
    6. Database of Genomic Variants
2. Code
    1. External
        1. Go-perl
            1. Data::Stag
            2. Map2slim
        2. Bedtools
        3. Weka
    2. Python Preqrequisites
        1. Networkx
        2. Numpy
        3. Pandas
        4. Matplotlib
    3. Preprocessing
        1. parse\_mania\_by\_weight.py
        2. enst\_to\_go.py
        3. parsecnvs.py
        4. enst\_to\_hpo.py
        5. dgv\_overlap.py
    4. Main scripts
        1. merge\_cnv\_gene\_go.py
        2. preprocess\_weka\_similarity\_121006c.py
        3. randomize\_for\_weka.py
        4. eval.py
    5. Called by eval.py
        1. weka\_remove.py
        2. calc\_results\_on\_cnvs.py
        3. noharm\_merge.py
        4. dgv\_arff.py
        5. parse\_j48graph.py

# 3 Usage
## 3.1 Prerequisites
### 3.1.1 Weka

Weka is required to run the machine learning algorithms used by cnv-prioritization.

    wget http://prdownloads.sourceforge.net/weka/weka-3-6-9.zip
    unzip weka-3-6-9.zip

### 3.1.2 Bedtools

Bedtools is required to annotate bedfiles (genes, go terms, etc).

    wget http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz
    tar zxf BEDTools.v2.17.0.tar.gz
    cd bedtools-2.1    7.0
    make
    
### 3.1.3 Go-perl

Go-perl is used to create go-slim mappings with map2slim.

#### 3.1.3.1 Preqrequisites for Go-perl

    wget http://search.cpan.org/CPAN/authors/id/C/CM/CMUNGALL/Data-Stag-0.11.tar.gz
    tar zxf Data-Stag-0.11.tar.gz
    cd Data-Stag-0.11 
    perl Makefile.PL PREFIX=~
    make install
    PERL5LIB=$PERL5LIB:~/share/perl5
    
#### 3.1.3.2 Go-perl

    wget http://search.cpan.org/CPAN/authors/id/C/CM/CMUNGALL/go-perl-0.14.tar.gz
    tar zxf go-perl-0.14.tar.gz
    cd go-perl-0.14
    perl Makefile.PL PREFIX=~
    make install
    
## 3.2 Preprocessing
### 3.2.1 Exon Annotations
Format: chrom, start, end, gene

Quick Download

    wget https://www.dropbox.com/s/d3n4uh21jwrxpfl/refseqgeneshg18_ensg_121218_exons.bed

Create from original sources

	Source:
    UCSC table browser file of ensembl genes (ENSG)

	Url:
    http://genome.ucsc.edu/cgi-bin/hgTables

	Settings:
    Clade - mammal
    Genome - human
    Assembly - Mar. 2006 (NCBI36/hg18)
    Group - Genes and Gene Prediction Tracs
    Track - Ensembl Genes
    Table - knownGene
    Region - Genome
    Format - all fields from selected data
    Output file - temp.ucsc.genes

	Transform:

    less temp.ucsc.genes | tail -n +2 | awk '{print $3"\t"$5"\t"$6"\t"$13}' | > refseqgeneshg18_ensg_120213.bed
    tail -n +2 refseqgeneshg18.ensg.121008.full | ruby -lane '$F[9].split(",").zip($F[10].split(",")).each{|x,y| print [$F[2],x,y,$F[12]].join("\t")}' > refseqgeneshg18_ensg_121218_exons.bed

### 3.2.2 CNV data
#### 3.2.2.1 Dataset 1

Quick download

    wget https://www.dropbox.com/s/a4uc78q8mekcj9r/db180k120227.csv
    
Format

	sample_id, array_design, genome_build, gender, phenotype, chr, cytoband, start, stop, size, max_start, max_stop, max_size, num_probes, type, classification, final_classification, inheritance, inheritance_coverage, p_value

#### 3.2.2.2 Dataset 2

Quick download

    wget https://www.dropbox.com/s/qptchg7qvhirnsz/Jan16.csv

Format

	Date Scanned, Lab no., Sex, Clinical Indication, Protocol, Genome Build, chromosome, start, stop, band, size, max start, max stop, Max size, No. probes, Gain/Loss, Classification,    Inheritance, Log2 ratio, ISCN

### 3.2.3 Gene-gene network (GeneMania)
Gene mania is used for gene-gene interactions.

Quick download

    wget https://www.dropbox.com/s/7m1y2urc8z1n344/mania_gt_0.0006.txt

Create from original sources

    wget http://genemania.org/data/current/precombined/Homo_sapiens.COMBINED.tgz
    tar zxf Homo_sapiens.COMBINED.tgz
    cd Homo_sapiens.COMBINED
    parse_mania_by_weight.py COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt mania_gt_0.0006.txt 0.0006 1 17474
    
    
Format
    
    gene, <gene:distance>*30

For each gene, list its closest 30 neighbours by random walk, the distance is the dijstra distance
This step takes a long time, on average, 20 min per 100 genes.  To speed this up, paralleize the generation using the <start> and <length> parameters.

    e.g. parse_mania_by_weight.py COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt mania_gt_0.0006.txt 0.0006 1 100
    e.g. parse_mania_by_weight.py COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt mania_gt_0.0006.txt 0.0006 101 100
    e.g. parse_mania_by_weight.py COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt mania_gt_0.0006.txt 0.0006 ... ... etc

### 3.2.4 Go slim

#### 3.2.4.1 Create go slim mapping

Quick download

    wget https://www.dropbox.com/s/yc7rj7r65eczn2z/goslim_mappings.txt

Create from original sources

    wget http://www.geneontology.org/GO_slims/goslim_generic.obo
    wget http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo
    map2slim goslim_generic.obo gene_ontology_ext.obo -outmap goslim_mappings.txt
    
#### 3.2.4.2 Gene id mappings ENSG/entrez to GO terms

Quick download

    wget https://www.dropbox.com/s/qhvqj11xy55g802/ensg_or_entrez_to_go.txt
    
Create from original sources

    wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz
    gunzip gene_association.goa_human.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
    gunzip idmapping.dat.gz
    less idmapping.dat | grep -P "\sENSG\d" > idmapping_filtered_ensg.dat
    enst_to_go.py gene_association.goa_human idmapping_filtered_ensg.dat goslim_mappings.txt > ensg_or_entrez_to_go.txt

#### 3.2.4.3 TODO: Remove this

    touch go_parents.txt
    touch patient_filter.txt

### 3.2.5 DGV data
DGV is used to annoate control regions.

Create from original sources

    wget http://dgvbeta.tcag.ca/dgv/docs/NCBI36_hg18_2012-11-23.txt
    less NCBI36_hg18_2012-11-23.txt | ruby -lane 'print (["chr"+$F[1],$F[2],$F[3]] + []).join("\t")' | tail -n +2 > NCBI36_hg18_2012-11-23.bed
    

### 3.2.6 Mapping from between GO/Genes/HPO
Quick download

    wget https://www.dropbox.com/s/rx2lqkobld7hn9q/go_and_ensg_to_hpo.txt
    
Create from original sources
    
    wget http://compbio.charite.de/hudson/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/genes_to_phenotype.txt
    less idmapping.dat | grep GeneID > idmapping_filtered_entrez.dat
    enst_to_hpo.py genes_to_phenotype.txt idmapping_filtered_ensg.dat idmapping_filtered_entrez.dat ensg_or_entrez_to_go.txt > go_and_ensg_to_hpo.txt
    
## 3.3 Annotate CNVs

### 3.3.1 Format input data (CNVs)

    for f in db180k120227.csv Jan16.csv; do
        parsecnvs.py                                                       `# parse cnvs` \
            ${constants}/${f}                                              `# input file` \
            ${f}                                                           `# source name` \
            "De novo|Likely Clinically Significant|Clinically Significant" `# case list ` \
            "CNV (seen in normal individuals)"                             `# ctrl list` \
            '*'                                                            `# pheno list` \
            >> dbAll.gff                                                   `# output file`
    done

### 3.3.2 Intersect CNVs with gene annotations
    bedtools-2.17.0/bin/intersectBed                          `# intersect bed` \
        -wao                                                  `# all lines in case/ctrl, 1 line per each overlap` \
        -a ${constants}/dbAll.gff                             `# case/ctrl file` \
        -b ${constants}/refseqgeneshg18_ensg_121218_exons.bed `# refseq genes hg18` \
        > case_control_refseq_exons.gff                       `# output file, with genes`

### 3.3.3 Intersect CNVs with DGV
    bedtools-2.17.0/bin/intersectBed               `# intersect bed` \
        -wao                                       `# all lines in case/ctrl, duplicate lines for each overlap` \
        -a ${constants}/dbAll.gff                  `# case/ctrl file` \
        -b ${constants}/NCBI36_hg18_2012-11-23.bed `# refseq genes hg18` \
        > cnvs_w_dgv.bed                           `# output file, with genes`
    gzip cnvs_w_dgv.bed

### 3.3.4 Calculate amout of DGV overlap for each CNV
    dgv_overlap.py cnvs_w_dgv.bed.gz | gzip > cnvs_w_dgv_overlap.txt.gz
    
## 3.4 Run

### 3.4.1 Setup variables
    constants=<directory where prerequisite files are located>
    scripts=<directory where code is loacated>
    top_x=10                                                    `# number of neighbours to include`
    dijkstra_distance_cutoff=200                                `# do not use neighbours above this distance`
    f_list[4]='x**(1/float(2))'                                 `# neighbour weighting function list`
    function_number=4                                           `# index of chosen function`
    function=${f_list[$function_number]}                        `# chosen function`
    out_file=weka_out_wo_sim_${function_number}_${top_x}        `# output file`

### 3.4.2 Format for weka
Input

	gff file of CNVs
	
Output

	Arff file
	Instances are genes
	Features are go-terms for the genes and their closest X neighbours
	
Command

    merge_cnv_gene_go.py                                   `# merge with refseq with go terms for weka` \
        ${constants}/mania_gt_0.0006.txt                           `# gene network` \
        ${constants}/ensg_or_entrez_to_go.txt                      `# refseq <-> GO` \
        ${constants}/test.gff                 `# case/ctrl, with genes` \
        weka_header.txt                                            `# weka header` \
        "${constants}/goslim_mappings.txt"                         `# go slim file` \
        -1                                                         `# evidence code filter level; 0=don't filter` \
        'y'                                                        `# filter by go slim` \
        "${top_x}"                                                 `# take to X nodes` \
        "${function}"                                              `# neighbour weight function` \
                                                                   `# [x, 2**x, x**2, sqrt(x)]]` \
        "${constants}/go_and_ensg_to_hpo.txt"                      `# go to hpo file` \
        "${dijkstra_distance_cutoff}"                              `# dijkstra distance cutoff` \
        ${out_file}.gz                                             `# weka body`
    gunzip -c ${out_file}.gz > ${out_file}.out

### 3.4.3 Calculate similarity scores
TODO: Test

    pushd /tmp/
    wget http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
    popd
    java -jar ontology-tools-0.5-2012-12-18.jar -l -Q ${out_file} -R go_and_ensg_to_hpo.txt -o ${out_file}.out -e "EXP,IDA,IPI,IMP,IGI,IEP,ISS,ISO,ISA,ISM,IBA,IBD,IKR,IRD,TAS,IC,IEA,IGC,RCA" >     ontology-tools.log
    gzip ${out_file}.out

### 3.4.4 Annotate CNVs with similarity scores
~2m

    weighted_gene_duplication=sim
    file_postfix=${function_number}_${top_x}_${weighted_gene_duplication}
    weighted_gene_duplication_full=weighted_gene_duplication_${weighted_gene_duplication}
    sim_out_file=weka_w_sim_${file_postfix}.pickle
    balance=bt_patient
    similarity_rank_cutoff=30

    preprocess_weka_similarity_121006c.py     `# ` \
        ${constants}/${out_file}.gz           `# ` \
        ${constants}/${out_file}.out.gz       `# ` \
        by_sim                                `# ` \
        remove                                `# ` \
        "${function}"                         `# ` \
        "${weighted_gene_duplication_full}"   `# ` \
        pluszero                              `# ` \
        ${sim_out_file}                       `# ` \
        1                                     `# ` \
        ${similarity_rank_cutoff}             `# ` \
        ${balance}

### 3.4.5 Run Weka
Create randomized set

Evalulate

    iteration=1
    folds=4
    remaining=new
    sim_dir=sim_job_3
    sim_file=../${sim_out_file}
    classifier_type_gene=rf
    classifier_type_cnv=nb
    gtrain=0,1
    gtest=2,2
    ctest=3,3
    SGE_TASK_ID=1
    array_id=0

    eval.py "SGE_TASK_ID=${SGE_TASK_ID};array_id=${array_id};skip=${skip};iteration=${iteration};function_number=${function_number};weighted_gene_duplication=${weighted_gene_duplication};top_x=${top_x};balance=${balance};folds=${folds};sim_file=${sim_file};remaining=${remaining};iteration=${iteration};classifier_type_gene=${classifier_type_gene};classifier_type_cnv=${classifier_type_cnv};gtrain=${gtrain};gtest=${gtest};ctest=${ctest}"
