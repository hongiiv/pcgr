pcgrr/R/acmg.R
pcg_report <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = NULL, pcgr_data = pcgr_data)
sample_calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
pcg_report_tumor_only <- pcgrr::generate_report_data_tumor_only(sample_calls, pcgr_data, pcgr_version,
                                                                            sample_name, pcgr_config, genome_seq, genome_assembly = assembly)

blogdown::serve_site()


2019-04-30 01:19:26 [INFO] Extending annotation descriptions related to UniprotKB/SwissProt protein features
2019-04-30 01:19:28 [INFO] Extending annotation descriptions related to Database of Curated Mutations (DoCM)
2019-04-30 01:19:28 [INFO] Extending annotation descriptions related to KEGG pathways
2019-04-30 01:19:29 [INFO] Extending annotation descriptions related to ClinVar

docker run -dit --name uta_20170629 -p 50827:5432 biocommons/uta:uta_20170629

docker run --rm -it  -e "DATA_DATE=2019-03-27" -e "UTA_DB_URL=postgresql://anonymous@0.0.0.0:50827/uta/uta_20170629" -e "HGVS_SEQREPO_DIR=/files/resources/seq_repo/latest" --network host -v /Users/hongchangbum/git/brca-exchange/pipeline/pipeline_running/monthly_releases/data_release_2019-03-27/resources:/files/resources -v /Users/hongchangbum/git/brca-exchange/pipeline/pipeline_running/monthly_releases/data_release_2019-03-27/brca_out:/files/data -v /Users/hongchangbum/git/brca-exchange/pipeline/pipeline_running/luigi_pipeline_credentials.cfg:/opt/luigi_pipeline_credentials.cfg -v /Users/hongchangbum/git/brca-exchange/pipeline/pipeline_running/previous_releases/latest_release.tar.gz:/files/previous_release.tar.gz -v /Users/hongchangbum/git/brca-exchange/pipeline/pipeline_running/monthly_releases/data_release_2019-03-27/release_notes_data_release_2019-03-27.txt:/files/release_notes.txt -v /Users/hongchangbum/git/brca-exchange/pipeline/pipeline_running/monthly_releases/data_release_2019-03-27/code:/opt/brca-exchange -v /var/run/docker.sock:/var/run/docker.sock brcachallenge/brca-exchange-pipeline:data_release_2019-03-27 /bin/bash


import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.variantmapper
hgvsparser = hgvs.parser.Parser()


hdp = hgvs.dataproviders.uta.connect()



import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)

hgvs_g = 'NC_000007.13:g.36561662C>T'
hgvs_c = 'NM_001637.3:c.1582G>A'
var_c = am.g_to_c(var_g, 'NM_001637.3')

----
docker run --rm -u root -t -i -v /Users/hongchangbum/git/pcgr/data:/data -v /Users/hongchangbum/git/pcgr/examples:/examples -v /Users/hongchangbum/git/pcgr/test_out:/test_out ngenebio/pcgr:dev /bin/bash


gene_summary<-as.data.frame(read.csv(file='~/Desktop/gene_summary.csv',header=TRUE, sep=','))
pcgr_data[['gene_summary']] <- gene_summary
save(pcgr_data,file='pcgr_data.rda')

install.packages('/pcgrr', repos = NULL, type="source")


pcgr_data$gene_summary <- pcgr_data$gene_summary %>% mutate(ENTREZ_ID = as.character(ENTREZ_ID))

View(pcgr_data$pcgr_all_annotation_columns)
View(pcgr_data$tier1_tags_display)

pcgr_data$pcgr_all_annotation_columns<-append(pcgr_data$pcgr_all_annotation_columns,'GENE_BACKGROUND')
pcgr_data$tier1_tags_display<-append(pcgr_data$tier1_tags_display,'GENE_BACKGROUND')

pcgr_data$pcgr_all_annotation_columns<-append(pcgr_data$pcgr_all_annotation_columns,'BIOMARKER_DESCRIPTION')
pcgr_data$tier1_tags_display<-append(pcgr_data$tier1_tags_display,'BIOMARKER_DESCRIPTION')



crosstalk::bscols(
  DT::datatable(variants_tier1_prognostic_shared, escape=F,extensions=c("Buttons","Responsive"), width = "100%",options=list(buttons = c('csv','excel'),dom = 'Bfrtip')) %>%
  DT::formatStyle('EVIDENCE_LEVEL', backgroundColor = DT::styleEqual(c('A: Validated','A: FDA/NCCN/ELN guidelines','B: Clinical evidence','B1: Clinical evidence: late trials','B2: Clinical evidence: early trials','C: Case study','D: Preclinical evidence','E: Indirect evidence'), c("#009E73","#009E73","#56B4E9", "#56B4E9","#56B4E9","#0072B2","#E69F00", "#F0E442")))
)


<a href='http://pfam.xfam.org/family/PF07714' target="_blank">Protein tyrosine kinase</a>
<a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04010' target="_blank">MAPK signaling pathway</a>, <a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04014' target="_blank">Ras signaling pathway</a>, <a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04015' target="_blank">Rap1 signaling pathway</a>, <a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04072' target="_blank">Phospholipase D signaling pathway</a>, <a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04151' target="_blank">PI3K-Akt signaling pathway</a>, <a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04640' target="_blank">Hematopoietic cell lineage</a>, <a href='http://www.genome.jp/dbget-bin/www_bget?path:hsa04916' target="_blank">Melanogenesis</a>



Error in match.arg(theme, themes()) :
  'arg' should be one of “default”, “cerulean”, “journal”, “flatly”, “readable”, “spacelab”, “united”, “cosmo”, “lumen”, “paper”, “sandstone”, “simplex”, “yeti”
Calls: <Anonymous> ... discover_rmd_resources -> render -> <Anonymous> -> base -> match.arg
Execution halted

table.dataTable.dtr-inline.collapsed>tbody>tr>td:first-child:before, table.dataTable.dtr-inline.collapsed>tbody>tr>th:first-child:before {
    /* top: 9px; */

table.dataTable.dtr-inline.collapsed>tbody>tr>th:first-child:before {
   left: 4px;
    height: 14px;
    width: 14px;
    display: block;
    position: absolute;
    color: white;
    border: 2px solid white;
    border-radius: 14px;
    box-shadow: 0 0 3px #444;
    box-sizing: content-box;
    text-align: center;
    font-family: 'Courier New', Courier, monospace;
    line-height: 14px;
    content: '+';
    background-color: #f04124;
}


----
New Bioinformatics Portal Provides One-Stop Shop for Multi-Omics Data Analysis
Apr 22, 2019 | Uduak Grace Thomas
Premium
NEW YORK (GenomeWeb) – Seeking to help biologists easily analyze multi-omics data from genomic, proteomic, transcriptomic, and other kinds of studies, scientists from Sanford Burnham Prebys (SBP) Medical Discovery Institute, the Genomics Institute of the Novartis Research Foundation (GNF), and the University of California San Diego have developed Metascape, an open-access, web-based portal that automatically pulls information from various open-source repositories.

biomarker_hits_snv_indels_specific <- pcgrr::get_clinical_associations_snv_indel(pcg_report_snv_indel[['variant_set']][['all']],pcgr_data, pcgr_config,tumor_type_specificity = 'specific_tumortype',biomarker_mapping_stringency = biomarker_mapping_stringency)

pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']] <- biomarker_hits_snv_indels_specific$clinical_evidence_item
pcg_report_snv_indel[['clinical_evidence_item']][['any_tumortype']] <- biomarker_hits_snv_indels_any$clinical_evidence_item

pcgrr/R/biomarkers.R
get_clinical_associations_snv_indel <- function(sample_calls, pcgr_data, pcgr_config, tumor_type_specificity = 'any_tumortype', biomarker_mapping_stringency = 1){


pcgrr/R/report.R
init_pcg_report <- function(config = NULL, sample_name = 'SampleX', pcgr_version = '0.6.0', genome_assembly = 'grch37', class = NULL, pcgr_data = NULL, type = 'somatic'){



annotating data/example_maf.txt...
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=CUL1&alteration=Y466S&tumorType=cancer&consequence=missense_variant&proteinStart=466&proteinEnd=466
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=AKT3&alteration=E182*&tumorType=LUAD&consequence=stop_gained&proteinStart=182&proteinEnd=182
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=PIK3CA&alteration=E542K&tumorType=GBM&consequence=missense_variant&proteinStart=542&proteinEnd=542
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=FGFR3&alteration=V271M&tumorType=LUAD&consequence=missense_variant&proteinStart=271&proteinEnd=271
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=EGFR&alteration=H304Y&tumorType=GBM&consequence=missense_variant&proteinStart=304&proteinEnd=304
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=PTEN&alteration=C136R&tumorType=GBM&consequence=missense_variant&proteinStart=136&proteinEnd=136
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=FGFR2&alteration=Q212K&tumorType=GBM&consequence=missense_variant&proteinStart=212&proteinEnd=212
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=ATM&alteration=L2890R&tumorType=LUAD&consequence=missense_variant&proteinStart=2890&proteinEnd=2890

http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=BRAF&alteration=V600E&tumorType=&consequence=missense_variant&proteinStart=600&proteinEnd=600

http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=BRAF&alteration=V600E&tumorType=MEL&consequence=missense_variant&proteinStart=600&proteinEnd=600

select a.id, a.entrez_gene_id, b.hugo_symbol, b.curatedRefSeq, b.curatedIsoform, a.description, a.evidence_type, b.name, b.oncogene from evidence a left outer join gene b on a.entrez_gene_id = b.entrez_gene_id where a.evidence_type='GENE_SUMMARY';


Separate by comma. Evidence type includes GENE_SUMMARY, GENE_BACKGROUND, MUTATION_SUMMARY, ONCOGENIC, MUTATION_EFFECT, VUS, PROGNOSTIC_IMPLICATION, DIAGNOSTIC_IMPLICATION, TUMOR_TYPE_SUMMARY, STANDARD_THERAPEUTIC_IMPLICATIONS_FOR_DRUG_SENSITIVITY, STANDARD_THERAPEUTIC_IMPLICATIONS_FOR_DRUG_RESISTANCE, INVESTIGATIONAL_THERAPEUTIC_IMPLICATIONS_DRUG_SENSITIVITY, INVESTIGATIONAL_THERAPEUTIC_IMPLICATIONS_DRUG_RESISTANCE

BRAF, an intracellular kinase, is frequently mutated in melanoma, thyroid and lung cancers among others.

BRAF is a serine/threonine kinase that plays a key role in the regulation of the mitogen-activated protein kinase (MAPK) cascade (PMID: 15520807), which under physiologic conditions regulates the expression of genes involved in cellular functions, including proliferation (PMID: 24202393). Genetic alterations in BRAF are found in a large percentage of melanomas, thyroid cancers and histiocytic neoplasms as well as a small fraction of lung and colorectal cancers. The most common BRAF point mutation is V600E, which deregulates the protein's kinase activity leading to constitutive BRAF activation, as BRAF V600E can signal as a monomer independently of RAS or upstream activation (PMID: 20179705). Other BRAF mutations have been found that affect the protein's propensity to dimerize (PMID: 16858395, 26343582, 12068308). The product of these alterations is a BRAF kinase that can activate MAPK signaling in an unregulated manner and, in some instances, is directly responsible for cancer growth (PMID: 15520807). Inhibitors of mutant BRAF, including vemurafenib and dabrafenib, are FDA-approved for the treatment of late-stage or unresectable melanoma.

http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=KRAS&alteration=G12C&tumorType=LUAD&consequence=missense_variant&proteinStart=12&proteinEnd=12
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=RB1&alteration=Q702*&tumorType=GBM&consequence=stop_gained&proteinStart=702&proteinEnd=702
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=TP53&alteration=R248Q&tumorType=GBM&consequence=missense_variant&proteinStart=248&proteinEnd=248
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=NF1&alteration=X1445_splice&tumorType=GBM&consequence=splice_region_variant
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=STK11&alteration=H168R&tumorType=LUAD&consequence=missense_variant&proteinStart=168&proteinEnd=168
done!
annotating data/example_fusions.txt...
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=NOT-A-GENE&alteration=Deletion&tumorType=
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=MLL2-intragenic&alteration=&tumorType=GBM&consequence=fusion&alterationType=structural_variant
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=ALK-EML4&alteration=&tumorType=LUAD&consequence=fusion&alterationType=structural_variant
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=EGFR-intragenic&alteration=&tumorType=GBM&consequence=fusion&alterationType=structural_variant
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=TMPRSS2-ERG&alteration=&tumorType=GBM&consequence=fusion&alterationType=structural_variant
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=TMPRSS2-ERG&alteration=&tumorType=LUAD&consequence=fusion&alterationType=structural_variant
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=ERBB2-intragenic&alteration=&tumorType=GBM&consequence=fusion&alterationType=structural_variant
done!
annotating data/example_cna.txt...
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=MET&alteration=Amplification&tumorType=LUAD
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=MET&alteration=Amplification&tumorType=GBM
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=ERBB2&alteration=Amplification&tumorType=LUAD
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=CDK4&alteration=Deletion&tumorType=LUAD
http://oncokb.org/legacy-api/indicator.json?source=cbioportal&hugoSymbol=CDK4&alteration=Amplification&tumorType=GBM
done!
annotating data/example_clinical.txt...
done!
annotating data/example_clinical.oncokb.txt...
/Users/hongchangbum/miniconda2/lib/python2.7/site-packages/matplotlib/cbook/deprecation.py:107: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
  warnings.warn(message, mplDeprecation, stacklevel=1)
done!



Query query = new Query(id, queryType, entrezGeneId, hugoSymbol, alteration, alterationType, svType, tumorType, consequence, proteinStart, proteinEnd, hgvs);
        Set<LevelOfEvidence> levelOfEvidences = levels == null ? LevelUtils.getPublicAndOtherIndicationLevels() : LevelUtils.parseStringLevelOfEvidences(levels);
        IndicatorQueryResp resp = IndicatorUtils.processQuery(query, geneStatus, levelOfEvidences, source, highestLevelOnly, null);


====
root@319fa5c6ddc5:/# ./pcgr.py --input_vcf data/examples/validation_test.vcf.gz test_out/ grch37 examples/pcgr_conf.BRCA.toml CancerExample
usage: pcgr.py [options] <PCGR_DIR> <OUTPUT_DIR> <GENOME_ASSEMBLY> <CONFIG_FILE> <SAMPLE_ID>
pcgr.py: error: argument genome_assembly: invalid choice: 'examples/pcgr_conf.BRCA.toml' (choose from 'grch37', 'grch38')
root@319fa5c6ddc5:/# ./pcgr.py --input_vcf data/examples/validation_test.vcf.gz . test_out grch37 examples/pcgr_conf.BRCA.toml Blood_Cancer_Example --no-docker --force_overwrite


View(sample_calls %>% select(GENE_HOGO_SYMBOL,HGVSp ,CONSEQUENCE ,CHROM, REF, ALT, CBMDB_ID, CIVIC_ID,CIVIC_ID_2))



2019-05-02 05:56:31 - pcgr-validate-config - WARNING - Prediction of MSI status is not perfomed in tumor-only mode (vcf_tumor_only = true)
2019-05-02 05:56:31 - pcgr-validate-config - WARNING - Estimation of mutational burden is not performed in tumor-only mode (vcf_tumor_only = true)
2019-05-02 05:56:31 - pcgr-validate-config - WARNING - Estimation of mutational signatures is not perfomed in tumor-only mode (vcf_tumor_only = true)
2019-05-02 05:56:31 - pcgr-validate-input - INFO - STEP 0: Validate input data
2019-05-02 05:56:31 - pcgr-validate-input - INFO - pcgr_validate_input.py / /data/examples/validation_test.vcf.gz None /examples/pcgr_conf.BRCA.toml grch37 --output_dir /test_out
2019-05-02 05:56:32 - pcgr-validate-input - INFO - Skipping validation of VCF file - as defined in configuration file (vcf_validation = false)
2019-05-02 05:56:32 - pcgr-validate-input - INFO - Checking if existing INFO tags of query VCF file coincide with PCGR INFO tags
2019-05-02 05:56:32 - pcgr-validate-input - INFO - No query VCF INFO tags coincide with PCGR INFO tags
2019-05-02 05:56:32 - pcgr-validate-input - INFO - Finished

2019-05-02 05:56:32 - pcgr-vep - INFO - STEP 1: Basic variant annotation with Variant Effect Predictor (94, GENCODE release 19, grch37)
2019-05-02 05:56:32 - pcgr-vep - INFO - vep --input_file /test_out/validation_test.pcgr_ready.vcf.gz --output_file /test_out/validation_test.pcgr_ready.vep.vcf --vcf --check_ref --flag_pick_allele --pick_order canonical,appris,biotype,ccds,rank,tsl,length --force_overwrite --species homo_sapiens --assembly GRCh37 --offline --fork 4 --hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad --variant_class --regulatory --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --total_length --allele_number --no_stats --no_escape --xref_refseq --dir /data/grch37/.vep --cache_version 94 --fasta /data/grch37/.vep/homo_sapiens/94_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
2019-05-02 05:56:35 - pcgr-vep - INFO - Finished

2019-05-02 05:56:35 - pcgr-vcfanno - INFO - STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar, dbNSFP, UniProtKB, cancerhotspots.org, CiVIC, CBMDB, DoCM, TCGA, ICGC-PCAWG, IntoGen_drivers, BRCA_Common_SNP)
2019-05-02 05:56:35 - pcgr-vcfanno - INFO - pcgr_vcfanno.py --num_processes 4 --brca_comm_snp --dbnsfp --docm --clinvar --icgc --civic --cbmdb --intogen_driver_mut --tcga --uniprot --cancer_hotspots --pcgr_onco_xref /test_out/validation_test.pcgr_ready.vep.vcf.gz /test_out/validation_test.pcgr_ready.vep.vcfanno.vcf /data/grch37
cat: /data/grch37/brca_comm_snp/brca_comm_snp.vcfanno.vcf_info_tags.txt: No such file or directory
vcfanno -p=4 /test_out/validation_test.pcgr_ready.vep.vcfanno.vcf.tmp.conf.toml /test_out/validation_test.pcgr_ready.vep.vcf.gz > /test_out/validation_test.pcgr_ready.vep.vcfanno.vcf.tmp.unsorted.1 2> /test_out/validation_test.pcgr_ready.vep.vcfanno.log
2019-05-02 05:56:41 - pcgr-vcfanno - INFO - Finished

2019-05-02 05:56:41 - pcgr-summarise - INFO - STEP 3: Cancer gene annotations with pcgr-summarise
2019-05-02 05:56:41 - pcgr-summarise - INFO - pcgr_summarise.py /test_out/validation_test.pcgr_ready.vep.vcfanno.vcf.gz /data/grch37
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Completed summary of functional annotations for 1 variants on chromosome 4
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Completed summary of functional annotations for 3 variants on chromosome 7
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Completed summary of functional annotations for 1 variants on chromosome 9
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Completed summary of functional annotations for 6 variants on chromosome 13
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Completed summary of functional annotations for 5 variants on chromosome 17
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Completed summary of functional annotations for 1 variants on chromosome X
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - bgzip -f /test_out/validation_test.pcgr_ready.vep.vcfanno.annotated.vcf
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - tabix -f -p vcf /test_out/validation_test.pcgr_ready.vep.vcfanno.annotated.vcf.gz
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Number of non-PASS/REJECTED variant calls: 0
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - Number of PASSed variant calls: 17
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - write_pass_vcf
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - bgzip -f /test_out/validation_test.pcgr_ready.vep.vcfanno.annotated.pass.vcf
2019-05-02 05:56:42 - pcgr-gene-annotate - INFO - tabix -f -p vcf /test_out/validation_test.pcgr_ready.vep.vcfanno.annotated.pass.vcf.gz
2019-05-02 05:56:42 - pcgr-summarise - INFO - Converting VCF to TSV with https://github.com/sigven/vcf2tsv
2019-05-02 05:56:42 - pcgr-summarise - INFO - vcf2tsv.py /test_out/Blood_Cancer_Example.pcgr_acmg.grch37.pass.vcf.gz --compress /test_out/Blood_Cancer_Example.pcgr_acmg.grch37.pass.tsv
2019-05-02 05:56:42 - pcgr-summarise - INFO - Finished

2019-05-02 05:56:42 - pcgr-writer - INFO - STEP 4: Generation of output files - variant interpretation report for precision oncology
/pcgr.R /test_out /test_out/Blood_Cancer_Example.pcgr_acmg.grch37.pass.tsv.gz None Blood_Cancer_Example /examples/pcgr_conf.BRCA.toml dev grch37 /
[1] "/test_out/Blood_Cancer_Example.pcgr_acmg.grch37.pass.tsv.gz"

rlogging::message(paste0("Number of SNVs: ", n_snvs))


2019-05-02 05:57:09 [INFO] Excluding 0 variants from non-nuclear chromosomes/scaffolds
2019-05-02 05:57:09 [INFO] Number of PASS variants: 17
2019-05-02 05:57:09 [INFO] Number of SNVs: 9
2019-05-02 05:57:09 [INFO] Number of deletions: 1
2019-05-02 05:57:09 [INFO] Number of insertions: 4
2019-05-02 05:57:09 [INFO] Number of block substitutions: 2
2019-05-02 05:57:09 [INFO] Extending annotation descriptions related to UniprotKB/SwissProt protein features
2019-05-02 05:57:11 [INFO] Extending annotation descriptions related to Gene summary & background
2019-05-02 05:57:11 [INFO] Extending annotation descriptions related to Database of Curated Mutations (DoCM)
2019-05-02 05:57:12 [INFO] Extending annotation descriptions related to KEGG pathways
2019-05-02 05:57:12 [INFO] Extending annotation descriptions related to ClinVar
2019-05-02 05:57:14 [INFO] Total sample calls (Blood_Cancer_Example): 17
2019-05-02 05:57:14 [INFO] Excluding coinciding germline variants in 1000 Genomes Project populations
2019-05-02 05:57:14 [INFO] Total sample calls remaining: 17
2019-05-02 05:57:14 [INFO] Excluding coinciding germline variants in any population in the genome aggregation database (gnomAD)
2019-05-02 05:57:14 [INFO] Total sample calls remaining: 17
2019-05-02 05:57:14 [INFO] ------
2019-05-02 05:57:14 [INFO] Generating data for tiered cancer genome report - germline-filtered callset tier model pcgr_acmg'
2019-05-02 05:57:14 [INFO] Number of protein-coding variants: 14
2019-05-02 05:57:14 [INFO] Looking up SNV/InDel biomarkers for precision oncology - Blood_Cancer_NOS
2019-05-02 05:57:15 [INFO] 7 clinical evidence item(s) found .. (2 unique variant(s)), mapping = exact
2019-05-02 05:57:15 [INFO] Underlying variant(s):
2019-05-02 05:57:15 [INFO] KIT missense_variant missense_variant:ENST00000288135.5:c.2447A>T:exon17:p.D816V 4:g.55599321A>T
2019-05-02 05:57:15 [INFO] JAK2 missense_variant missense_variant:ENST00000381652.3:c.1849G>T:exon14:p.V617F 9:g.5073770G>T
2019-05-02 05:57:15 [INFO] 0 clinical evidence item(s) found .. mapping = codon
2019-05-02 05:57:15 [INFO] 2 clinical evidence item(s) found .. (2 unique variant(s)), mapping = exon
2019-05-02 05:57:15 [INFO] Underlying variant(s):
2019-05-02 05:57:15 [INFO] FLT3 intron_variant NA 13:g.28608140A>ACATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGG
2019-05-02 05:57:15 [INFO] FLT3 inframe_insertion inframe_insertion:ENST00000241453.7:c.1751_1798dup:exon14:p.Y599_D600insASDNEYFYVDFREYEY 13:g.28608257T>TCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGG
2019-05-02 05:57:16 [INFO] Looking up SNV/InDel biomarkers for precision oncology - any tumortype
2019-05-02 05:57:16 [INFO] 39 clinical evidence item(s) found .. (4 unique variant(s)), mapping = exact
2019-05-02 05:57:16 [INFO] Underlying variant(s):
2019-05-02 05:57:16 [INFO] EGFR missense_variant missense_variant:ENST00000275493.2:c.2369C>T:exon20:p.T790M 7:g.55249071C>T
2019-05-02 05:57:16 [INFO] KIT missense_variant missense_variant:ENST00000288135.5:c.2447A>T:exon17:p.D816V 4:g.55599321A>T
2019-05-02 05:57:16 [INFO] JAK2 missense_variant missense_variant:ENST00000381652.3:c.1849G>T:exon14:p.V617F 9:g.5073770G>T
2019-05-02 05:57:16 [INFO] BRCA1 stop_gained stop_gained:ENST00000471181.2:c.4327C>T:exon12:p.R1443X 17:g.41234451G>A
2019-05-02 05:57:16 [INFO] 0 clinical evidence item(s) found .. mapping = codon
2019-05-02 05:57:16 [INFO] 2 clinical evidence item(s) found .. (2 unique variant(s)), mapping = exon
2019-05-02 05:57:16 [INFO] Underlying variant(s):
2019-05-02 05:57:16 [INFO] FLT3 intron_variant NA 13:g.28608140A>ACATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGG
2019-05-02 05:57:16 [INFO] FLT3 inframe_insertion inframe_insertion:ENST00000241453.7:c.1751_1798dup:exon14:p.Y599_D600insASDNEYFYVDFREYEY 13:g.28608257T>TCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGG
2019-05-02 05:57:18 [INFO] Generating tiered set of result variants for output in tab-separated values (TSV) file
2019-05-02 05:57:18 [INFO] Number of noncoding/silent variants: 2
2019-05-02 05:57:18 [INFO] ------
2019-05-02 05:57:18 [INFO] ------
2019-05-02 05:57:18 [INFO] Assigning elements to PCGR value boxes
2019-05-02 05:57:18 [INFO] ------
2019-05-02 05:57:18 [INFO] Writing JSON file with report contents
2019-05-02 05:57:19 [INFO] ------
2019-05-02 05:57:19 [INFO] Rendering HTML report with rmarkdown
2019-05-02 05:57:31 - pcgr-writer - INFO - Finished



sample_calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)



The IDH2 (isocitrate dehydrogenase 2) protein is an enzyme that catalyzes the oxidative decarboxylation of isocitrate to α-ketoglutarate (α-KG) in the tricarboxylic acid (TCA) cycle. IDH2 utilizes NADP(+) as an electron acceptor and it is expressed in the mitochondria, where it plays a role in cell metabolism and energy production via TCA cycle. Cancer-associated mutations in the catalytic site of IDH2 confer a gain-of-function of neomorphic enzymatic activity allowing the mutant enzyme to convert α-KG to the “oncometabolite” D-2-hydroxyglutarate (2-HG) (PMID: 20171147). 2-HG promotes tumor development by inhibiting a variety of enzymes that require α-KG as a substrate, including enzymes involved in DNA demethylation, histone demethylation, adaptation to hypoxia and collagen maturation (PMID: 23630074). IDH2 mutations have been identified in hematologic malignancies, particularly in acute myeloid leukemia (AML), as well as in solid tumors such as gliomas and cholangiocarcinomas.

<a href='http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=12600' target="_blank">COSM12600</a>, <a href='http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=29117' target="_blank">COSM29117</a>


  if (!("COSMIC" %in% colnames(vcf_data_df))){
    cosmic_annotation_links <- pcgrr::generate_annotation_link(vcf_data_df,
                                                               vardb = "COSMIC",
                                                      group_by_var = "VAR_ID",
                                                      url_prefix = "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=",
                                                      link_key_var = "COSMIC_MUTATION_ID",
                                                      link_display_var = "COSMIC_MUTATION_ID")

2019-05-03 00:27:50 - pcgr-summarise - INFO - pcgr_summarise.py /test_out/validation_test.pcgr_ready.vep.vcfanno.vcf.gz /data/grch37
2019-05-03 00:27:50 - pcgr-gene-annotate - INFO - {'ENST00000288135': {'CANCER_SUSCEPTIBILITY_CUI': 'C0238198', 'ENTREZ_ID': '3815', 'UNIPROT_ACC': 'P10721', 'ONCOGENE': '1', 'CHEMBL_COMPOUND_ID': 'CHEMBL941&CHEMBL535&CHEMBL1421&CHEMBL1336&CHEMBL19461

if csq_fields[vep_csq_fields2index['PICK']] == "1"



2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'missense_variant', 'MODERATE', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000309486', 'protein_coding', '9/22', '', 'ENST00000309486.4:c.605T>A', 'ENSP00000310938.4:p.Leu202His', '1633/7114', '605/4704', '202/1567', 'L/H', 'cTc/cAc', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000310938', '', 'Q9UE29&Q9NQR3&Q7KYU6&Q3YB53&Q3YB50&Q3YB49&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V4Z8&G4V4Z7&E7EP70&C9IZW4&C4PFY7', 'UPI000014170B', 'NM_007297.3', 'hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0&Pfam_domain:PF12820&PIRSF_domain:PIRSF001734', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'missense_variant', 'MODERATE', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000346315', 'protein_coding', '10/19', '', 'ENST00000346315.3:c.1493T>A', 'ENSP00000246907.4:p.Leu498His', '1687/6451', '1493/4875', '498/1624', 'L/H', 'cTc/cAc', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000246907', '', 'Q9NQR3&Q92897&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3B891&K4K7V3&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&E7EWN5&E7EMP0&C9IZW4&C4PFY7', 'UPI000014170D', '', 'PIRSF_domain:PIRSF001734&Pfam_domain:PF12820&hmmpanther:PTHR13763:SF0&hmmpanther:PTHR13763', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000351666', 'protein_coding', '', '6/18', 'ENST00000351666.3:c.548-3006T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000338007', '', 'G3XAC3&Q9UE29&Q9NQR3&Q92897&Q7KYU6&K7EPC7&K4JXS7&G1UI37&C4PFY7', 'UPI000014172D', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000352993', 'protein_coding', '', '9/21', 'ENST00000352993.3:c.670+1808T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000312236', 'P38398', 'Q9UE29&Q9NQR3&Q92897&Q7KYU6&K7EPC7&K4JXS7&G1UI37&C4PFY7', 'UPI000013ECD3', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'missense_variant', 'MODERATE', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000354071', 'protein_coding', '10/18', '', 'ENST00000354071.3:c.1493T>A', 'ENSP00000326002.6:p.Leu498His', '1725/6411', '1493/4797', '498/1598', 'L/H', 'cTc/cAc', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000326002', '', 'Q9NQR3&Q92897&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3B891&K4K7V3&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&F8W8H7&E7EWN5&C9IZW4', 'UPI000013CC06', '', 'PIRSF_domain:PIRSF001734&Pfam_domain:PF12820&hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'missense_variant', 'MODERATE', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000357654', 'protein_coding', '10/23', '', 'ENST00000357654.3:c.1493T>A', 'ENSP00000350283.3:p.Leu498His', '1612/7094', '1493/5592', '498/1863', 'L/H', 'cTc/cAc', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', 'CCDS11453.1', 'ENSP00000350283', 'P38398', 'Q9UE29&Q9NQR3&Q92897&Q7KYU6&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3LRH8&Q3B891&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&E9PFZ0&E7EWN5&E7EP70&C9IZW4&C4PFY7', 'UPI0000126AC8', 'NM_007294.3', 'hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0&Pfam_domain:PF12820&PIRSF_domain:PIRSF001734', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000468300', 'protein_coding', '', '10/21', 'ENST00000468300.1:c.787+706T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', 'CCDS11455.2', 'ENSP00000417148', 'P38398', 'Q9UE29&Q9NQR3&Q92897&K4JXS7&G1UI37&C4PFY7', 'UPI000037834E', 'NM_007299.3', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -

========flag_pick_allele
vep --input_file /test_out/validation_test.pcgr_ready.vcf.gz --output_file /test_out/validation_test.pcgr_ready.vep.vcf --vcf --check_ref --flag_pick_allele --pick_order canonical,appris,biotype,ccds,rank,tsl,length --force_overwrite --species homo_sapiens --assembly GRCh37 --offline --fork 4 --hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad --variant_class --regulatory --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --total_length --allele_number --no_stats --no_escape --xref_refseq --dir /data/grch37/.vep --cache_version 94 --fasta /data/grch37/.vep/homo_sapiens/94_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz



2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'missense_variant', 'MODERATE', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000471181', 'protein_coding', '10/24', '', 'ENST00000471181.2:c.1493T>A', 'ENSP00000418960.2:p.Leu498His', '1725/5936', '1493/5655', '498/1884', 'L/H', 'cTc/cAc', '', '1', '', '-1', '', '1', 'SNV', 'HGNC', '1100', 'YES', '', 'CCDS11456.2', 'ENSP00000418960', '', 'Q9UE29&Q9NQR3&Q92897&Q7KYU6&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3B891&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&E9PFC7&E7EWN5&C9IZW4&C4PFY7', 'UPI0000E0360B', 'NM_007300.3', 'hmmpanther:PTHR13763:SF0&hmmpanther:PTHR13763&Pfam_domain:PF12820&PIRSF_domain:PIRSF001734', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - 1

2019-05-03 01:37:42 - pcgr-gene-annotate - INFO - {



'vep_csq_fields2index': {'SAS_AF_1KG': 41, 'CANONICAL': 26, 'IMPACT': 2, 'AFR_AF_1KG': 37, 'TREMBL': 31, 'NFE_AF_GNOMAD': 48, 'Amino_acids': 15, 'AMR_AF_1KG': 38, 'SAS_AF_GNOMAD': 50, 'MOTIF_NAME': 54, 'ENSP': 29, 'CCDS': 28, 'AFR_AF_GNOMAD': 43, 'Feature_type': 5, 'STRAND': 20, 'Gene': 4, 'PICK': 22, 'INTRON': 9, 'APPRIS': 27, 'DOMAINS': 34, 'EXON': 8, 'FIN_AF_GNOMAD': 47, 'Protein_position': 14, 'HIGH_INF_POS': 56, 'ALLELE_NUM': 18, 'ASJ_AF_GNOMAD': 45, 'Codons': 16, 'EAS_AF_1KG': 39, 'GLOBAL_AF_GNOMAD': 42, 'VARIANT_CLASS': 23, 'OTH_AF_GNOMAD': 49, 'SWISSPROT': 30, 'FLAGS': 21, 'EUR_AF_1KG': 40, 'BIOTYPE': 7, 'DISTANCE': 19, 'cDNA_position': 12, 'Existing_variation': 17, 'AMR_AF_GNOMAD': 44, 'GLOBAL_AF_1KG': 36, 'PHENO': 53, 'UNIPARC': 32, 'MOTIF_POS': 55, 'Feature': 6, 'SYMBOL_SOURCE': 24, 'HGVS_OFFSET': 35, 'HGVSp': 11, 'CDS_position': 13, 'EAS_AF_GNOMAD': 46, 'Consequence': 1, 'HGVSc': 10, 'SYMBOL': 3}, 

'dbnsfp_prediction_algorithms': ['SIFT', 'LRT', 'Mutationtaster', 'MutationAssessor', 'FATHMM', 'FATHMM_MKL_coding', 'PROVEAN', 'M-CAP', 'MutPred', 'metaSVM', 'metaLR', 'GERP_rs', 'splice_site_ada', 'splice_site_rf'], 


'vep_csq_index2fields': {1: 'Consequence', 2: 'IMPACT', 3: 'SYMBOL', 4: 'Gene', 5: 'Feature_type', 6: 'Feature', 7: 'BIOTYPE', 8: 'EXON', 9: 'INTRON', 10: 'HGVSc', 11: 'HGVSp', 12: 'cDNA_position', 13: 'CDS_position', 14: 'Protein_position', 15: 'Amino_acids', 16: 'Codons', 17: 'Existing_variation', 18: 'ALLELE_NUM', 19: 'DISTANCE', 20: 'STRAND', 21: 'FLAGS', 22: 'PICK', 23: 'VARIANT_CLASS', 24: 'SYMBOL_SOURCE', 26: 'CANONICAL', 27: 'APPRIS', 28: 'CCDS', 29: 'ENSP', 30: 'SWISSPROT', 31: 'TREMBL', 32: 'UNIPARC', 34: 'DOMAINS', 35: 'HGVS_OFFSET', 36: 'GLOBAL_AF_1KG', 37: 'AFR_AF_1KG', 38: 'AMR_AF_1KG', 39: 'EAS_AF_1KG', 40: 'EUR_AF_1KG', 41: 'SAS_AF_1KG', 42: 'GLOBAL_AF_GNOMAD', 43: 'AFR_AF_GNOMAD', 44: 'AMR_AF_GNOMAD', 45: 'ASJ_AF_GNOMAD', 46: 'EAS_AF_GNOMAD', 47: 'FIN_AF_GNOMAD', 48: 'NFE_AF_GNOMAD', 49: 'OTH_AF_GNOMAD', 50: 'SAS_AF_GNOMAD', 53: 'PHENO', 54: 'MOTIF_NAME', 55: 'MOTIF_POS', 56: 'HIGH_INF_POS'}}


vep_csq_fields2index
2019-05-03 01:37:42 - pcgr-gene-annotate - INFO - {'SAS_AF_1KG': 41, 'CANONICAL': 26, 'IMPACT': 2, 'AFR_AF_1KG': 37, 'TREMBL': 31, 'NFE_AF_GNOMAD': 48, 'Amino_acids': 15, 'AMR_AF_1KG': 38, 'SAS_AF_GNOMAD': 50, 'MOTIF_NAME': 54, 'ENSP': 29, 'CCDS': 28, 'AFR_AF_GNOMAD': 43, 'Feature_type': 5, 'STRAND': 20, 'Gene': 4, 'PICK': 22, 'INTRON': 9, 'APPRIS': 27, 'DOMAINS': 34, 'EXON': 8, 'FIN_AF_GNOMAD': 47, 'Protein_position': 14, 'HIGH_INF_POS': 56, 'ALLELE_NUM': 18, 'ASJ_AF_GNOMAD': 45, 'Codons': 16, 'EAS_AF_1KG': 39, 'GLOBAL_AF_GNOMAD': 42, 'VARIANT_CLASS': 23, 'OTH_AF_GNOMAD': 49, 'SWISSPROT': 30, 'FLAGS': 21, 'EUR_AF_1KG': 40, 'BIOTYPE': 7, 'DISTANCE': 19, 'cDNA_position': 12, 'Existing_variation': 17, 'AMR_AF_GNOMAD': 44, 'GLOBAL_AF_1KG': 36, 'PHENO': 53, 'UNIPARC': 32, 'MOTIF_POS': 55, 'Feature': 6, 'SYMBOL_SOURCE': 24, 'HGVS_OFFSET': 35, 'HGVSp': 11, 'CDS_position': 13, 'EAS_AF_GNOMAD': 46, 'Consequence': 1, 'HGVSc': 10, 'SYMBOL': 3}

=======
Question: How to select predefine transcript when annotate VCF
 --filter "Feature is YOUR_TRANSCRIPT_ID"


2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000491747', 'protein_coding', '', '10/22', 'ENST00000491747.2:c.787+706T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', 'CCDS11454.2', 'ENSP00000420705', 'P38398', 'Q9UE29&Q9NQR3&Q92897&Q7KYU6&K7EPC7&K4JXS7&G1UI37&C4PFY7&B4DES0', 'UPI000013C85E', 'NM_007298.3', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'missense_variant', 'MODERATE', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000493795', 'protein_coding', '9/22', '', 'ENST00000493795.1:c.1352T>A', 'ENSP00000418775.1:p.Leu451His', '1584/5732', '1352/5451', '451/1816', 'L/H', 'cTc/cAc', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', 'CCDS11459.2', 'ENSP00000418775', '', 'Q9UE29&Q9NQR3&Q7KYU6&Q4EW25&Q3YB53&Q3YB50&Q3YB49&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&E9PFZ0&E7EP70&C9IZW4&C4PFY7', 'UPI000013C860', '', 'Pfam_domain:PF12820&PIRSF_domain:PIRSF001734&hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000586385', 'protein_coding', '', '1/7', 'ENST00000586385.1:c.5-30087T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000465818', '', 'Q9NQR3&C6YB45&C4PFY7', 'UPI0001B065E7', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000591534', 'protein_coding', '', '1/10', 'ENST00000591534.1:c.-43-19517T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000467329', '', 'Q9UE29&Q9NQR3&Q7KYU6&K7EPC7&K4JXS7&C4PFY7', 'UPI0002840F63', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO - ['T', 'intron_variant', 'MODIFIER', 'BRCA1', 'ENSG00000012048', 'Transcript', 'ENST00000591849', 'protein_coding', '', '1/4', 'ENST00000591849.1:c.-99+31233T>A', '', '', '', '', '', '', '', '1', '', '-1', '', '', 'SNV', 'HGNC', '1100', '', '', '', 'ENSP00000465347', '', 'K7EJW3', 'UPI0002840F64', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
2019-05-03 00:55:21 - pcgr-gene-annotate - INFO -



17	41246055	.	A	T	10033	.	AB=0.522523;ABP=5.94472;ANN=T|missense_variant|MODERATE|BRCA1|BRCA1|transcript|NM_007300.3|protein_coding|10/24|c.1493T>A|p.Leu498His|1725/7270|1493/5655|498/1884||,T|missense_variant|MODERATE|BRCA1|BRCA1|transcript|NM_007297.3|protein_coding|9/22|c.1352T>A|p.Leu451His|1633/7115|1352/5451|451/1816||,T|missense_variant|MODERATE|BRCA1|BRCA1|transcript|NM_007294.3|protein_coding|10/23|c.1493T>A|p.Leu498His|1725/7207|1493/5592|498/1863||,T|intron_variant|MODIFIER|BRCA1|BRCA1|transcript|NM_007298.3|protein_coding|9/21|c.787+706T>A||||||,T|intron_variant|MODIFIER|BRCA1|BRCA1|transcript|NM_007299.3|protein_coding|10/21|c.787+706T>A||||||,T|non_coding_transcript_exon_variant|MODIFIER|BRCA1|BRCA1|transcript|NR_027676.1|pseudogene|10/23|n.1629T>A||||||;AO=348;CIGAR=1X;DP=666;DPB=666;DPRA=0;EPP=7.22845;EPPR=5.22274;GTI=0;HGVS_p=L498H,L451H,L498H,.,.,.;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=2028.64;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=13533;QR=12164;RO=318;RPL=219;RPP=53.5532;RPPR=48.9253;RPR=129;RUN=1;SAF=170;SAP=3.40965;SAR=178;SRF=154;SRP=3.69315;SRR=164;TYPE=snp;freebayes_AC=1;freebayes_AF=0.5;freebayes_AN=2;technology.illumina=1;dbNSFP_PROVEAN_score=-5.28,-5.28,-5.32,-5.08,-5.21,-5.21;dbNSFP_SiPhy_29way_logOdds=1.0301;dbNSFP_GERP___RS=-0.912;dbNSFP_Polyphen2_HVAR_score=0.521,0.521,0.986,0.986,0.976,0.967;dbNSFP_MutationAssessor_score=3.155;dbNSFP_DANN_score=0.97988827224212338;dbNSFP_MetaSVM_pred=T;dbNSFP_integrated_confidence_value=0;dbNSFP_VEST3_score=0.291,0.23,0.679,0.614,0.294,0.294,0.291,0.699,0.66,0.279,0.291,0.7,0.279;dbNSFP_Interpro_domain=BRCA1,_serine-rich_domain;dbNSFP_FATHMM_pred=D,D,D,.,.,.;dbNSFP_integrated_fitCons_score=0.6512;dbNSFP_CADD_raw=2.410025;dbNSFP_SIFT_score=0.003,0.003,0.003,0.002,0.002,0.002;dbNSFP_MetaLR_score=0.8727;dbNSFP_phyloP100way_vertebrate=0.586000;dbNSFP_Polyphen2_HDIV_score=0.832,0.832,1,0.999,0.999,0.998;dbNSFP_LRT_pred=N;dbNSFP_PROVEAN_pred=D,.,D,.,.,.;dbNSFP_fathmm_MKL_coding_pred=N;dbNSFP_phastCons100way_vertebrate=0.050000;dbNSFP_fathmm_MKL_coding_score=0.04369;dbNSFP_aapos=498,451,498,498,498,472;dbNSFP_CADD_phred=18.89;dbNSFP_MetaSVM_score=-0.0807;dbNSFP_Polyphen2_HDIV_pred=P,P,D,D,D,D;dbNSFP_FATHMM_score=-3.03,-3.03,-3.03,-3.03,-3.03,-3.03;dbNSFP_MutationTaster_score=1,1,1,1,1,1,1,1,1,1,1,1,1,1;dbNSFP_MutationTaster_pred=N,N,N,N,N,N,N,N,N,N,N,N,N,N;dbNSFP_MutationAssessor_pred=M;dbNSFP_Polyphen2_HVAR_pred=P,P,D,D,D,D;dbNSFP_MetaLR_pred=D;dbNSFP_SIFT_pred=D,.,D,.,.,.;dbNSFP_LRT_score=0.411127;dbNSFP_phastCons20way_mammalian=0.834000;dbNSFP_phyloP20way_mammalian=1.199000;ROI=BRCA1_15;CSQ=T|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000309486|protein_coding|9/22||ENST00000309486.4:c.605T>A|ENSP00000310938.4:p.Leu202His|1633/7114|605/4704|202/1567|L/H|cTc/cAc||1||-1|||SNV|HGNC|1100||||ENSP00000310938||Q9UE29&Q9NQR3&Q7KYU6&Q3YB53&Q3YB50&Q3YB49&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V4Z8&G4V4Z7&E7EP70&C9IZW4&C4PFY7|UPI000014170B|NM_007297.3|hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0&Pfam_domain:PF12820&PIRSF_domain:PIRSF001734|||||||||||||||||||||||,T|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000346315|protein_coding|10/19||ENST00000346315.3:c.1493T>A|ENSP00000246907.4:p.Leu498His|1687/6451|1493/4875|498/1624|L/H|cTc/cAc||1||-1|||SNV|HGNC|1100||||ENSP00000246907||Q9NQR3&Q92897&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3B891&K4K7V3&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&E7EWN5&E7EMP0&C9IZW4&C4PFY7|UPI000014170D||PIRSF_domain:PIRSF001734&Pfam_domain:PF12820&hmmpanther:PTHR13763:SF0&hmmpanther:PTHR13763|||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000351666|protein_coding||6/18|ENST00000351666.3:c.548-3006T>A||||||||1||-1|||SNV|HGNC|1100||||ENSP00000338007||G3XAC3&Q9UE29&Q9NQR3&Q92897&Q7KYU6&K7EPC7&K4JXS7&G1UI37&C4PFY7|UPI000014172D|||||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000352993|protein_coding||9/21|ENST00000352993.3:c.670+1808T>A||||||||1||-1|||SNV|HGNC|1100||||ENSP00000312236|P38398|Q9UE29&Q9NQR3&Q92897&Q7KYU6&K7EPC7&K4JXS7&G1UI37&C4PFY7|UPI000013ECD3|||||||||||||||||||||||||,T|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000354071|protein_coding|10/18||ENST00000354071.3:c.1493T>A|ENSP00000326002.6:p.Leu498His|1725/6411|1493/4797|498/1598|L/H|cTc/cAc||1||-1|||SNV|HGNC|1100||||ENSP00000326002||Q9NQR3&Q92897&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3B891&K4K7V3&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&F8W8H7&E7EWN5&C9IZW4|UPI000013CC06||PIRSF_domain:PIRSF001734&Pfam_domain:PF12820&hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0|||||||||||||||||||||||,T|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000357654|protein_coding|10/23||ENST00000357654.3:c.1493T>A|ENSP00000350283.3:p.Leu498His|1612/7094|1493/5592|498/1863|L/H|cTc/cAc||1||-1|||SNV|HGNC|1100|||CCDS11453.1|ENSP00000350283|P38398|Q9UE29&Q9NQR3&Q92897&Q7KYU6&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3LRH8&Q3B891&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&E9PFZ0&E7EWN5&E7EP70&C9IZW4&C4PFY7|UPI0000126AC8|NM_007294.3|hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0&Pfam_domain:PF12820&PIRSF_domain:PIRSF001734|||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000468300|protein_coding||10/21|ENST00000468300.1:c.787+706T>A||||||||1||-1|||SNV|HGNC|1100|||CCDS11455.2|ENSP00000417148|P38398|Q9UE29&Q9NQR3&Q92897&K4JXS7&G1UI37&C4PFY7|UPI000037834E|NM_007299.3||||||||||||||||||||||||,T|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000471181|protein_coding|10/24||ENST00000471181.2:c.1493T>A|ENSP00000418960.2:p.Leu498His|1725/5936|1493/5655|498/1884|L/H|cTc/cAc||1||-1||1|SNV|HGNC|1100|YES||CCDS11456.2|ENSP00000418960||Q9UE29&Q9NQR3&Q92897&Q7KYU6&Q4EW25&Q3YB53&Q3YB50&Q3YB49&Q3B891&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&G1UI37&E9PFC7&E7EWN5&C9IZW4&C4PFY7|UPI0000E0360B|NM_007300.3|hmmpanther:PTHR13763:SF0&hmmpanther:PTHR13763&Pfam_domain:PF12820&PIRSF_domain:PIRSF001734|||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000491747|protein_coding||10/22|ENST00000491747.2:c.787+706T>A||||||||1||-1|||SNV|HGNC|1100|||CCDS11454.2|ENSP00000420705|P38398|Q9UE29&Q9NQR3&Q92897&Q7KYU6&K7EPC7&K4JXS7&G1UI37&C4PFY7&B4DES0|UPI000013C85E|NM_007298.3||||||||||||||||||||||||,T|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000493795|protein_coding|9/22||ENST00000493795.1:c.1352T>A|ENSP00000418775.1:p.Leu451His|1584/5732|1352/5451|451/1816|L/H|cTc/cAc||1||-1|||SNV|HGNC|1100|||CCDS11459.2|ENSP00000418775||Q9UE29&Q9NQR3&Q7KYU6&Q4EW25&Q3YB53&Q3YB50&Q3YB49&K7EPC7&K4K7V3&K4JXS7&K4JUB1&G4V503&G4V502&G4V500&G4V4Z8&G4V4Z7&E9PFZ0&E7EP70&C9IZW4&C4PFY7|UPI000013C860||Pfam_domain:PF12820&PIRSF_domain:PIRSF001734&hmmpanther:PTHR13763&hmmpanther:PTHR13763:SF0|||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000586385|protein_coding||1/7|ENST00000586385.1:c.5-30087T>A||||||||1||-1|||SNV|HGNC|1100||||ENSP00000465818||Q9NQR3&C6YB45&C4PFY7|UPI0001B065E7|||||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000591534|protein_coding||1/10|ENST00000591534.1:c.-43-19517T>A||||||||1||-1|||SNV|HGNC|1100||||ENSP00000467329||Q9UE29&Q9NQR3&Q7KYU6&K7EPC7&K4JXS7&C4PFY7|UPI0002840F63|||||||||||||||||||||||||,T|intron_variant|MODIFIER|BRCA1|ENSG00000012048|Transcript|ENST00000591849|protein_coding||1/4|ENST00000591849.1:c.-99+31233T>A||||||||1||-1|||SNV|HGNC|1100||||ENSP00000465347||K7EJW3|UPI0002840F64|||||||||||||||||||||||||


docker pull sigven/pcgr:dev

docker run --rm -u root -t -i -v /Users/hongchangbum/git/pcgr/data:/data -v /Users/hongchangbum/git/pcgr/examples:/examples -v /Users/hongchangbum/git/pcgr/test_out:/test_out sigven/pcgr:dev /bin/bash

install.packages('/pcgrr', repos = NULL, type="source")

./pcgr.py --input_vcf data/examples/validation_test.vcf.gz . test_out grch37 examples/pcgr_conf.BRCA.toml Blood_Cancer_Example --no-docker --force_overwrite

newRow <- data.frame(genesymbol='FLT3', variant_name='FLT3-ITD',biomarker_description='FLT3-ITD (internal tandem duplications) frequently occur in patients with hematologic malignancies such as chronic myelogenous leukemia, acute myeloid leukemia (AML) and myelodysplastic syndrome, but particularly in cytogenetically normal AML (CN-AML).', refbuild='GRCh37',evidence_type='Prognostic',evidence_level='A:Validated',evidence_description='Meta-analysis of studies involving cytogentically normal younger (<60) patients showed reduced overall and relapse-free survival for patients with FLT3-ITD.', evidence_id='EID100002',evidence_direction='Supports',clinical_significance='Poor Outcome', variant_origin='Somatic Mutation', status='accepted', pubmed_html_link='<a href="http://www.ncbi.nlm.nih.gov/pubmed/24801015">', cancer_type='Acute Myeloid Leukemia',disease_ontology_id='DOID:9119', therapeutic_context=NA,rating='3', alteration_type='MUT', mapping_category='exon', mapping_rank='1',eitem_exon='14',eitem_codon=NA,eitem_consequence='inframe_insertion')

pcgr_data$civic_biomarkers <- rbind(pcgr_data$civic_biomarkers, newRow)

save(pcgr_data, file='pcgr_data.rda')
save(pcgr_data, file='~/git/pcgr/data/grch37/rda/pcgr_data.rda')

newRow <- data.frame(genesymbol='FLT3', variant_name='FLT3-ITD',biomarker_description='FLT3-ITD (internal tandem duplications) frequently occur in patients with hematologic malignancies such as chronic myelogenous leukemia, acute myeloid leukemia (AML) and myelodysplastic syndrome, but particularly in cytogenetically normal AML (CN-AML).', refbuild='GRCh37',evidence_type='Prognostic',evidence_level='A:Validated',evidence_description='Meta-analysis of studies involving cytogentically normal younger (<60) patients showed reduced overall and relapse-free survival for patients with FLT3-ITD.', evidence_id='EID100001',evidence_direction='Supports',clinical_significance='Poor Outcome', variant_origin='Somatic Mutation', status='accepted', pubmed_html_link='<a href="http://www.ncbi.nlm.nih.gov/pubmed/24801015">', cancer_type='Acute Myeloid Leukemia',disease_ontology_id='DOID:9119', therapeutic_context=NA,rating='3', alteration_type='MUT', mapping_category='exon', mapping_rank='1',eitem_exon='14',eitem_codon=NA,eitem_consequence='insertion')

