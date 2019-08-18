Input
-----

The CPSR workflow accepts a single input file:

-  An unannotated, single-sample **VCF file** (>= v4.2) with germline
   calls (SNVs/InDels)

VCF
~~~

-  We **strongly** recommend that the input VCF is compressed and
   indexed using `bgzip <http://www.htslib.org/doc/tabix.html>`__ and
   `tabix <http://www.htslib.org/doc/tabix.html>`__
-  If the input VCF contains multi-allelic sites, these will be subject
   to `decomposition <http://genome.sph.umich.edu/wiki/Vt#Decompose>`__
-  Variants used for reporting should be designated as ‘PASS’ in the VCF
   FILTER column

**IMPORTANT NOTE**: CPSR generates a number of VCF INFO annotation tags
that is appended to the query VCF. We will therefore encourage the users
to submit query VCF files that have not been subject to annotations by
other means, but rather a VCF file that comes directly from variant
calling. If not, there are likely to be INFO tags in the query VCF file
that coincide with those produced by CPSR.

CPSR configuration file
~~~~~~~~~~~~~~~~~~~~~~~

The cancer predisposition sequencing report can be flexibly configured
in a TOML-formatted configuration file. The default TOML configuration
file, with descriptive comments wrt. usage are shown below:

::

   # CPSR configuration options (TOML).

   [cancer_predisposition_genes]
   ABCB11 = true
   ABRAXAS1 = true
   ACD = true
   AIP = true
   AKT1 = true
   ALK = true
   APC = true
   APOBEC3B = true
   AR = true
   ATM = true
   ATR = true
   AXIN1 = true
   AXIN2 = true
   BAP1 = true
   BARD1 = true
   BLM = true
   BMPR1A = true
   BRAF = true
   BRCA1 = true
   BRCA2 = true
   BRIP1 = true
   BUB1B = true
   CASR = true
   CBL = true
   CDC73 = true
   CDH1 = true
   CDH10 = true
   CDK4 = true
   CDKN1B = true
   CDKN1C = true
   CDKN2A = true
   CEBPA = true
   CEP57 = true
   CHEK2 = true
   COL7A1 = true
   CTNNA1 = true
   CTNNB1 = true
   CTR9 = true
   CTRC = true
   CXCR4 = true
   CYLD = true
   DDB2 = true
   DICER1 = true
   DIRAS3 = true
   DIS3L2 = true
   DKC1 = true
   DOCK8 = true
   DROSHA = true
   DTX3L = true
   EGFR = true
   ELANE = true
   ENG = true
   EPCAM = true
   ERBB4 = true
   ERCC1 = true
   ERCC2 = true
   ERCC3 = true
   ERCC4 = true
   ERCC5 = true
   ETV6 = true
   EXT1 = true
   EXT2 = true
   EZH2 = true
   FAH = true
   FANCA = true
   FANCB = true
   FANCC = true
   FANCD2 = true
   FANCE = true
   FANCF = true
   FANCG = true
   FANCI = true
   FANCL = true
   FANCM = true
   FAT1 = true
   FEN1 = true
   FH = true
   FLCN = true
   GALNT12 = true
   GATA2 = true
   GBA = true
   GJB2 = true
   GPC3 = true
   GREM1 = true
   HABP2 = true
   HFE = true
   HMBS = true
   HNF1A = true
   HNF1B = true
   HOXB13 = true
   HRAS = true
   ITK = true
   JMJD1C = true
   KDR = true
   KIF1B = true
   KIT = true
   KRAS = true
   LMO1 = true
   LZTR1 = true
   MAP2K1 = true
   MAP2K2 = true
   MAX = true
   MEN1 = true
   MET = true
   MITF = true
   MLH1 = true
   MLH3 = true
   MPL = true
   MRE11 = true
   MSH2 = true
   MSH3 = true
   MSH6 = true
   MTAP = true
   MUTYH = true
   NBN = true
   NF1 = true
   NF2 = true
   NHP2 = true
   NOP10 = true
   NRAS = true
   NSD1 = true
   NTHL1 = true
   OGG1 = true
   PALB2 = true
   PAX5 = true
   PDGFRA = true
   PHOX2B = true
   PIK3CA = true
   PINK1 = true
   PMS1 = true
   PMS2 = true
   POLD1 = true
   POLE = true
   POLH = true
   POLQ = true
   POT1 = true
   PPM1D = true
   PRDM9 = true
   PRF1 = true
   PRKAR1A = true
   PRSS1 = true
   PTCH1 = true
   PTEN = true
   PTPN11 = true
   PTPN13 = true
   RAD50 = true
   RAD51 = true
   RAD51B = true
   RAD51C = true
   RAD51D = true
   RAF1 = true
   RB1 = true
   RCC2 = true
   RECQL = true
   RECQL4 = true
   RET = true
   RFWD3 = true
   RHBDF2 = true
   RING1 = true
   RINT1 = true
   RMRP = true
   RUNX1 = true
   SBDS = true
   SCG5 = true
   SDHA = true
   SDHAF2 = true
   SDHB = true
   SDHC = true
   SDHD = true
   SERPINA1 = true
   SETBP1 = true
   SH2B3 = true
   SH2D1A = true
   SHOC2 = true
   SLC25A13 = true
   SLX4 = true
   SMAD4 = true
   SMARCA4 = true
   SMARCB1 = true
   SMARCE1 = true
   SOS1 = true
   SPINK1 = true
   SPOP = true
   SPRED1 = true
   SPRTN = true
   SRY = true
   STAT3 = true
   STK11 = true
   SUFU = true
   TERF2IP = true
   TERT = true
   TGFBR1 = true
   TGFBR2 = true
   TMEM127 = true
   TNFRSF6 = true
   TP53 = true
   TP63 = true
   TRIM37  = true
   TSC1 = true
   TSC2 = true
   TSHR = true
   UROD = true
   VHL = true
   WAS = true
   WRN = true
   WT1 = true
   XPA = true
   XPC = true
   XRCC2 = true

   [maf_limits]
   ## choose upper MAF thresholds for report of unclassified variants
   maf_tgp = 0.001
   maf_gnomad = 0.001

   [popgen]
   ## choose population source in gnomAD and 1000 Genomes Project, defaults to the global set

   ## For gnomaAD, this can by any of the following values (three-letter codes):
   ## "afr" - African/American (12,020 individuals (7,652 WES / 4,368 WGS))
   ## "amr" - Admixed American (17,210 individuals (16,791 WES / 419 WGS))
   ## "eas" - East Asian (9,435 individuals (8,624 WES / 811 WGS))
   ## "sas" - Sout Asian (15,391 individuals (15,391 WES / 0 WGS))
   ## "asj" - Ashkenazi Jewish (5,076 individuals (4,925 WES / 151 WGS))
   ## "nfe" - Non-Finnish European (63,369 individuals (55,860 WES / 7,509 WGS))
   ## "fin" - Finnish (12,897 individuals (11,150 WES / 1,747 WGS))
   ## "oth" - Other (3,234 individuals (2,743 WES / 491 WGS))
   ## "global" - All populations (138,632 individuals (123,136 WES / 15,496 WGS))
   pop_gnomad = "global"

   ## For 1000 Genomes Project this can be any of "afr","amr","sas","eas","eur","global"
   pop_tgp = "global"

   [visual]
   # Choose visual theme of report, any of: "default", "cerulean", "journal", "flatly", "readable", "spacelab", "united", "cosmo", "lumen", "paper", "sandstone", "simplex", or "yeti" (https://bootswatch.com/)
   report_theme = "default"

   [custom_tags]
   ## list VCF info tags that should be present in JSON output
   ## tags should be comma separated, i.e. custom_tags = "GATK_FILTER,VARSCAN_FILTER"
   custom_tags = ""

   [dbnsfp]
   ## CPSR performs a ranking of unclassified/novel variants according to its likelihood of being pathogenic, adopting the
   ## same implementations of ACMG-AMP guidelines as outlined by Huang et al., 2018 (Cell)
   ##
   ## One criteria relates to whether insilico predictions can support a pathogenic versus benign nature
   ## The user can here configure the minimum number of algorithms that must have called a variant
   ## as 'damaging', and the maximum number of algorithms that called the variant as 'tolerated'
   ## in order for it to have a consensus call as 'damaging' (and vice versa for consensus calls
   ## as 'tolerated')
   ## min_majority should not exceed 8 and should not be less than 5, max_minority should not exceed 2,
   ## min_majority + max_minority should not exceed 8

   min_majority = 5
   max_minority = 1

   [gwas]
   gwas_hits = false
   ## Required p-value for reporting of GWAS hits
   p_value_min = 5e-8

   [other]
   vcf_validation = true
   n_vcfanno_proc = 4
   n_vep_forks = 4
   vep_skip_intergenic = false
