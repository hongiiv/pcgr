## Cancer Predisposition Sequencing Reporter (CPSR)

### Overview

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. The workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven), utilizing the Docker environment for encapsulation of code and software dependencies. While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

*CPSR* accepts a query file with raw germline variant calls encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. analyzing SNVs/InDels). The software performs extensive variant annotation and produces an interactive HTML report, in which the user can investigate three main sets of variants identified in the query set:

1. Germline variants in a selected set of [configurable cancer predisposition genes](predisposition.md), that are **previously reported** as pathogenic or likely pathogenic in ClinVar (with no conflicting interpretations)

2. **Unclassified variants** constitute the set of germline variants within the configurable cancer predisposition gene list that are either:
	* Registered as *variant of uncertain significance (VUS)* in [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), or
	* *Is a novel protein-coding variant* (i.e. not reported in [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), and not found in [gnomAD](http://gnomad.broadinstitute.org/) or [1000 Genomes Project](http://www.internationalgenome.org/) user-defined population datasets), or
	* *Is a rare protein-coding variant* (e.g. minor allele frequency (MAF) < 0.001 in user-defined [gnomAD](http://gnomad.broadinstitute.org/) or [1000 Genomes Project](http://www.internationalgenome.org/) population datasets)
		* *The upper MAF threshold (e.g. 0.001) for listing of unclassified variants can be configured by the user*


3. Variants overlapping with previously identified hits in genome-wide association studies (GWAS) of cancer phenotypes (i.e. low to moderate risk conferring alleles), using [NHGRI-EBI Catalog of published genome-wide association studies](https://www.ebi.ac.uk/gwas/) as the underlying source.

The (**classified** and **unclassified**) variant sets can be interactively explored and ranked further through different types of filters (associated phenotypes, genes, variant consequences, population MAF etc.). Importantly, the unclassified variants are assigned and ranked according to a *pathogenicity score*, which is based on the aggregation of scores according to previously established [ACMG criteria](https://www.ncbi.nlm.nih.gov/pubmed/25741868) and also cancer-specific criteria, as outlined and specified in several previous studies ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052); [Maxwell et al., *Am J Hum Genet.*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395); [Amendola et al., *Am J Hum Genet.*,  2016](https://www.ncbi.nlm.nih.gov/pubmed/27181684)). See also [*Related work*](https://github.com/sigven/cpsr#related-work) below).

##### Cancer predisposition genes

We have compiled a comprehensive list of genes that are implicated in cancer predisposition and cancer syndromes. Three different sources were combined:
* A list of 152 genes that were curated and established within TCGA’s pan-cancer study ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052))
* A list of 107 protein-coding genes that has been manually curated in COSMIC’s [Cancer Gene Census v86](https://cancer.sanger.ac.uk/census),
* A list of 148 protein-coding genes established by experts within the Norwegian Cancer Genomics Consortium (http://cancergenomics.no)

The combination of the three sources resulted in a non-redundant set of 209 protein-coding genes of relevance for predisposition to tumor development. We want to make it explicit that this list of 209 genes is by no means regarded as an international consensus, but should rather be subject to continuous update by the international community that carry expertise on genetic risk factors for cancer.


### Example report

* [Cancer predisposition sequencing report](http://folk.uio.no/sigven/example.cpsr.grch37.html)

### Annotation resources included in _cpsr - 0.3.0_

* [VEP v94](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor (GENCODE version 28/19 (grch38/grch37) as the gene reference dataset), includes [gnomAD r2](http://gnomad.broadinstitute.org/), [dbSNP build 151/150](http://www.ncbi.nlm.nih.gov/SNP/), [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
* [dBNSFP v3.5](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (August 2017)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (November 2018)
* [DisGeNET](http://www.disgenet.org) - Database of gene-disease associations (v5.0, May 2017)
* [UniProt/SwissProt KnowledgeBase 2018_010](http://www.uniprot.org) - Resource on protein sequence and functional information (November 2018)
* [Pfam v32](http://pfam.xfam.org) - Database of protein families and domains (September 2018)
* [CancerMine v6](https://zenodo.org/record/1336650#.W9QMdRMzaL4) - Literature-derived database of tumor suppressor genes/proto-oncogenes (November 2018)
* [NHGRI-EBI GWAS catalog](https://www.ebi.ac.uk/gwas//) - GWAS catalog for cancer phenotypes (October 29th 2018)

### Documentation

*IN PROGRESS*

### News
* *November 19th 2018*: **0.3.0 pre-release**
  * Bug fixing and bundle update
* *November 12th 2018*: **0.2.1 pre-release**
  * Improved ACMG classification transparency
* *November 6th 2018*: **0.2.0 pre-release**
  * Adjustments of ACMG classification criteria
	* Mechanisms of disease for cancer susceptibility genes (GoF vs. LoF) retrieved from [Maxwell et al., *Am J Hum Genet*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395)
	* Exceptions for HFE/SERPINA1 wrt. high population MAF (*BA1*)
	* Threshold for genes with "primarily truncations" set to 90% pathogenic truncations (*BP1*)
	* Consider only pathogenic variants (not likely pathogenic) when checking for novel peptide changes at pathogenic loci (*PS1*/*PM5*)
* *October 27th 2018*: **0.1.1 pre-release**
	* Added documentation of ACMG evidence items in report output
	* GWAS hits are optionable to include
* *October 5th 2018*: **0.1.0 pre-release**
	* Initial release of CPSR - reporting of germline variants for cancer predisposition

### Getting started

#### STEP 0: Install PCGR (dev version)

Make sure you have a working installation of PCGR (**dev version**) and the accompanying dev data bundle(s) (walk through [steps 0-2](https://github.com/sigven/pcgr#getting-started)).

#### STEP 1: Download the latest release

Download the [pre-release](https://github.com/sigven/cpsr/releases/tag/v0.3.0) of *cpsr* (run script and configuration file)

#### STEP 2: Configuration

A few elements of the workflow can be figured using the *cpsr* configuration file, encoded in [TOML](https://github.com/toml-lang/toml) (an easy to read file format).

The initial step of the workflow performs [VCF validation](https://github.com/EBIvariation/vcf-validator) on the input VCF file. This procedure is very strict, and often causes the workflow to return an error due to various violations of the VCF specification. If the user trusts that the most critical parts of the input VCF is properly encoded,  a setting in the configuration file (`vcf_validation = false`) can be used to turn off VCF validation.

An exhaustive, predefined list of 209 cancer predisposition/syndrome genes can also be configured.

#### STEP 3: Run example

Run the workflow with **cpsr.py**, which takes the following arguments and options:

	usage: cpsr.py [-h] [--input_vcf INPUT_VCF] [--force_overwrite] [--version]
			[--basic] [--docker-uid DOCKER_USER_ID] [--no-docker]
			pcgr_base_dir output_dir {grch37,grch38} configuration_file
			sample_id

	Cancer Predisposition Sequencing Reporter (CPSR) - report of cancer-predisposing
	germline variants

	positional arguments:
	pcgr_base_dir         Directory that contains the PCGR data bundle
				    directory, e.g. ~/pcgr-dev
	output_dir            Output directory
	{grch37,grch38}       Genome assembly build: grch37 or grch38
	configuration_file    Configuration file (TOML format)
	sample_id             Sample identifier - prefix for output files

	optional arguments:
	-h, --help            show this help message and exit
	--input_vcf INPUT_VCF
				    VCF input file with somatic query variants
				    (SNVs/InDels). (default: None)
	--force_overwrite     By default, the script will fail with an error if any
				    output file already exists. You can force the
				    overwrite of existing result files by using this flag
				    (default: False)
	--version             show program's version number and exit
	--basic               Run functional variant annotation on VCF through
				    VEP/vcfanno, omit report generation (STEP 4) (default:
				    False)
	--docker-uid DOCKER_USER_ID
				    Docker user ID. Default is the host system user ID. If
				    you are experiencing permission errors, try setting
				    this up to root (`--docker-uid root`) (default: None)
	--no-docker           Run the CPSR workflow in a non-Docker mode (see
				    install_no_docker/ folder for instructions (default:
				    False)




The *cpsr* software bundle contains an example VCF file. It also contains a configuration file (*cpsr.toml*).

Analysis of the example VCF can be performed by the following command:

`python ~/cpsr-0.3.0/cpsr.py --input_vcf ~/cpsr-0.3.0/example.vcf.gz`
` ~/pcgr-dev ~/cpsr-0.3.0 grch37 ~/cpsr-0.3.0/cpsr.toml example`

Note that the example command also refers to the PCGR directory (*pcgr-dev*), which contains the data bundle that are necessary for both *PCGR* and *CPSR*.

This command will run the Docker-based *cpsr* workflow and produce the following output files in the _cpsr_ folder:

  1. __example.cpsr.grch37.pass.vcf.gz (.tbi)__ - Bgzipped VCF file with functional/clinical annotations
  2. __example.cpsr.grch37.pass.tsv.gz__ - Compressed TSV file (generated with [vcf2tsv](https://github.com/sigven/vcf2tsv)) with functional/clinical annotations
  3. __example.cpsr.grch37.html__ - Interactive HTML report with clinically relevant variants in cancer predisposition genes organized into tiers
  4. __example.cpsr.grch37.json.gz__ - Compressed JSON dump of HTML report content
  5. __example.cpsr.snvs_indels.tiers.grch37.tsv__ - TSV file with most important annotations of tier-structured SNVs/InDels

### Related work

* [CharGer - Characterization of Germline variants](https://github.com/ding-lab/CharGer)

### Contact

sigven@ifi.uio.no
