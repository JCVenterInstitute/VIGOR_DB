Viral Databases for VIGOR4
====================

This repository contains viral reference databases and configurations for use with [VIGOR4](https://github.com/JCVenterInstitute/VIGOR4)


# Organization

The viral database contains three kinds of files: the viral databases
themselves, virus specific configuration files and mature peptide
databases.

##  Viral Database Files

The viral databases are simple fasta files.

## Configuration Files

Each viral database may optionally have its own VIGOR4 .INI
configuration file which overrides the default settings for that virus
and matched genes. By default, VIGOR4 looks for a virus specific
configuration file by appending the .ini suffix. For example, the Flu
A database is named "flua_db" and the flua configuration files is
"flua_db.ini"".

### Configuration Sections

Any configuration parameters set outside of a section are defaults for
the virus. The configuration file may also contain one or more gene
specific sections. These gene specific sections should be named
gene:GENE_NAME.

### Configuration Parameters

Below is a list of the configuration parameters that may be set in the
deflines of the viral databases or in the virus specific configuration
files.  The list was produced by using the --list-config-parameters
option to VIGOR4 which is the recommended way to get the most update
list of settable parameters.


alignment_score_factor

    Environment variable:          VIGOR_ALIGNMENT_SCORE_FACTOR
    System.property:               vigor.alignment_score_factor
    Settable levels                virus config,gene config/variant
    Description:                   Weight of alignment score in evaluating alignment

alternate_startcodon

    Environment variable:          VIGOR_ALTERNATE_STARTCODON
    System.property:               vigor.alternate_startcodon
    Settable levels                gene config/variant
    Description:                   Alternate start codons for gene. Format is CODON[,CODON,..]

circular_genome

    Environment variable:          VIGOR_CIRCULAR_GENOME
    System.property:               vigor.circular_genome
    Settable levels                virus config,gene config/variant
    Note:                          This parameters is currently unimplemented
    Description:                   When this parameter is set to TRUE, VIGOR consider the genome as circular, enabling annotating genes spanning both ends of the sequence (which would be continuous when circularized).

complete_genome

    Environment variable:          VIGOR_COMPLETE_GENOME
    System.property:               vigor.complete_genome
    Settable levels                virus config,gene config/variant
    Note:                          This parameters is currently unimplemented
    Description:                   When this parameter is set to TRUE, VIGOR assumes that the sequence represents the whole genome and all the genes must be present (unless a gene has the is_optional attribute)

excludes_gene

    Environment variable:          VIGOR_EXCLUDES_GENE
    System.property:               vigor.excludes_gene
    Settable levels                gene config/variant
    Description:                   Gene/sequence-level attribute describing the incompatibility of presence of two or more alternate genes (generally used only in large DNA viruses).

frameshift_sensitivity

    Environment variable:          VIGOR_FRAMESHIFT_SENSITIVITY
    System.property:               vigor.frameshift_sensitivity
    Settable levels                virus config,gene config/variant
    Note:                          This parameters is currently unimplemented
    Description:                   Dictates the sensitivity VIGOR should use in handling frame-shifts. Accepted values: 0, 1, 2, with 2 being the strictest, forcing VIGOR to create pseudogenes and raising error messages whenever a perfect model for a given gene cannot be created (this is used most for validating assemblies).

is_optional

    Environment variable:          VIGOR_IS_OPTIONAL
    System.property:               vigor.is_optional
    Settable levels                gene config/variant
    Description:                   This parameter works in combination with complete_genome and frameshift_sensitivity: if it is set to TRUE, the absence of that particular gene for which is set suppresses the stricter behavior set by the use of the other parameters.

is_required

    Environment variable:          VIGOR_IS_REQUIRED
    System.property:               vigor.is_required
    Settable levels                gene config/variant
    Description:                   By setting this parameter to TRUE, VIGOR will produce an error message in the case it cannot create a valid gene model for that given gene.

leakystop_notFound_score

    Environment variable:          VIGOR_LEAKYSTOP_NOTFOUND_SCORE
    System.property:               vigor.leakystop_notFound_score
    Settable levels                virus config,gene config/variant
    Description:                   Weight to be assigned to the gene model candidates scoring when an expected leaky stop codon is not found (note: in some viruses the leaky stop codon is not always present).

leakystop_score_factor

    Environment variable:          VIGOR_LEAKYSTOP_SCORE_FACTOR
    System.property:               vigor.leakystop_score_factor
    Settable levels                virus config,gene config/variant
    Description:                   Weight for scoring leaky stops for gene model selection.

matpepdb

    Environment variable:          VIGOR_MATPEPDB
    System.property:               vigor.matpepdb
    Settable levels                gene config/variant
    Description:                   Location and name of the mature peptides reference database (e.g. matpepdb=<vigordata>/flua_ha_mp)

mature_pep_mincoverage

    Environment variable:          VIGOR_MATURE_PEP_MINCOVERAGE
    System.property:               vigor.mature_pep_mincoverage
    Settable levels                virus config,gene config/variant
    Description:                   Minimum percent coverage of the sequence of the reference mature peptide to consider the prediction valid.

mature_pep_minidentity

    Environment variable:          VIGOR_MATURE_PEP_MINIDENTITY
    System.property:               vigor.mature_pep_minidentity
    Settable levels                virus config,gene config/variant
    Description:                   Minimum percent identity between the sequence of the reference mature peptide and the sequence of the candidate mature peptide.

mature_pep_minsimilarity

    Environment variable:          VIGOR_MATURE_PEP_MINSIMILARITY
    System.property:               vigor.mature_pep_minsimilarity
    Settable levels                virus config,gene config/variant
    Description:                   Minimum percent similarity between the sequence of the reference mature peptide and the sequence of the candidate mature peptide. Note: similarity is calculated using compatibility values found in Blosum40 distance matrix

max_aa_overlap

    Environment variable:          VIGOR_MAX_AA_OVERLAP
    System.property:               vigor.max_aa_overlap
    Settable levels                virus config,gene config/variant
    Description:                   Maximum number of proteins that may overlap for alignment fragments to be considered compatible when generating a gene model

max_align_merge_aa_gap

    Environment variable:          VIGOR_MAX_ALIGN_MERGE_AA_GAP
    System.property:               vigor.max_align_merge_aa_gap
    Settable levels                virus config,gene config/variant
    Description:                   Maximum number of proteins in a gap between two alignments to consider them for merging.

max_gene_overlap

    Environment variable:          VIGOR_MAX_GENE_OVERLAP
    System.property:               vigor.max_gene_overlap
    Settable levels                virus config,gene config/variant
    Description:                    In reporting gene models, maximum overlap of genes allowed.

max_intron_size

    Environment variable:          VIGOR_MAX_INTRON_SIZE
    System.property:               vigor.max_intron_size
    Settable levels                virus config,gene config/variant
    Description:                   Maximum sequence length of an intron

max_nt_overlap

    Environment variable:          VIGOR_MAX_NT_OVERLAP
    System.property:               vigor.max_nt_overlap
    Settable levels                virus config,gene config/variant
    Description:                   Maximum number of nucleotides that may overlap for alignment fragments to be considered compatible when generating a gene model

min_functional_len

    Environment variable:          VIGOR_MIN_FUNCTIONAL_LEN
    System.property:               vigor.min_functional_len
    Settable levels                gene config/variant
    Description:                   Minimum functional length for a protein (expressed in aa) to be functional: if a premature stop codon makes it shorter than that, it should be annotated as pseudogene.

min_gene_coverage

    Environment variable:          VIGOR_MIN_GENE_COVERAGE
    System.property:               vigor.min_gene_coverage
    Settable levels                virus config,gene config/variant
    Description:                   Minimum coverage of genes

min_intron_size

    Environment variable:          VIGOR_MIN_INTRON_SIZE
    System.property:               vigor.min_intron_size
    Settable levels                virus config,gene config/variant
    Description:                   Minimum sequence length of an intron

min_missing_AA_size

    Environment variable:          VIGOR_MIN_MISSING_AA_SIZE
    System.property:               vigor.min_missing_AA_size
    Settable levels                virus config,gene config/variant
    Description:                   Minimum number of proteins missing in a given alignment to search for missing exons.

min_pseudogene_coverage

    Environment variable:          VIGOR_MIN_PSEUDOGENE_COVERAGE
    System.property:               vigor.min_pseudogene_coverage
    Settable levels                virus config,gene config/variant
    Note:                          This parameters is currently unimplemented
    Description:                   Minimum percentage of coverage in the alignment between candidate pseudogene and reference protein, based on the longest of the two.

min_pseudogene_identity

    Environment variable:          VIGOR_MIN_PSEUDOGENE_IDENTITY
    System.property:               vigor.min_pseudogene_identity
    Settable levels                virus config,gene config/variant
    Note:                          This parameters is currently unimplemented
    Description:                   Minimum percentage of identity in the alignment between candidate pseudogene and reference protein, based on the longest of the two.

min_pseudogene_similarity

    Environment variable:          VIGOR_MIN_PSEUDOGENE_SIMILARITY
    System.property:               vigor.min_pseudogene_similarity
    Settable levels                virus config,gene config/variant
    Note:                          This parameters is currently unimplemented
    Description:                   Minimum percentage of similarity in the alignment between candidate pseudogene and reference protein, based on the longest of the two.

min_seq_gap_length

    Environment variable:          VIGOR_MIN_SEQ_GAP_LENGTH
    System.property:               vigor.min_seq_gap_length
    Settable levels                virus config,gene config/variant
    Description:                   Minimum number of undefined nucleotides (i.e. Ns) to consider it a sequencing gap.

noncanonical_splicing

    Environment variable:          VIGOR_NONCANONICAL_SPLICING
    System.property:               vigor.noncanonical_splicing
    Settable levels                virus config,gene config/variant
    Description:                   List of alternative splicing donor and acceptor sequence pairs. Format: noncanonical_splicing=donor+acceptor,donor+acceptor,... (e.g. noncanonical_splicing=AA+GT)

note

    Environment variable:          VIGOR_NOTE
    System.property:               vigor.note
    Settable levels                virus config,gene config/variant
    Description:                   Information associated to a gene or protein variant (to be reported in the .tbl and GFF 3 outputs).

ribosomal_slippage

    Environment variable:          VIGOR_RIBOSOMAL_SLIPPAGE
    System.property:               vigor.ribosomal_slippage
    Settable levels                gene config/variant
    Description:                   Ribosomal slippage. Format is ribosomal_slippage=offset/frameshift/regex (e.g. ribosomal_slippage=-7/+1/[BDHKNTWY][BCHMNSVY][BCHMNSVY][BDHKNTWY][BDHKNTWY][BDHKNTWY][BCHMNSVY][BDGKNRSV][BCDHKMNSTVWY][BCHMNSVY])

rna_editing

    Environment variable:          VIGOR_RNA_EDITING
    System.property:               vigor.rna_editing
    Settable levels                gene config/variant
    Description:                   RNA editing. Format is rna_editing=offset/regex/insertion string/note (e.g. rna_editing=0/GGGG/[ARWMDHVN][ARWMDHVN][ARWMDHVN][ARWMDHVN][GRSKBDVN][GRSKBDVN][GRSKBDVN]/four non-templated G's inserted during transcription)

shared_cds

    Environment variable:          VIGOR_SHARED_CDS
    System.property:               vigor.shared_cds
    Settable levels                gene config/variant
    Description:                   List of other genes (gene symbols) sharing the same region of the viral genome. Format: shared_cds=Gene1_ID,Gene_2ID (e.g. shared_cds=NSP1-2; shared_cds=NSP1-1,NSP1-3)

splice_form

    Environment variable:          VIGOR_SPLICE_FORM
    System.property:               vigor.splice_form
    Settable levels                virus config,gene config/variant
    Description:                   Describes the expected gene structure with actual lengths of exons and introns of the reference protein aligned back to its originating genome (or the closest possible genome, if the original one is not available). Format: e\d+(i-?\d+e\d+)* The string always starts and ends with exons. (e.g. splice_form="e26i686e268"; splice_form="e1656") Note: given that this parameter is specific for each single isoform, it is allowed only on the sequence header.

splicing_score_factor

    Environment variable:          VIGOR_SPLICING_SCORE_FACTOR
    System.property:               vigor.splicing_score_factor
    Settable levels                virus config,gene config/variant
    Description:                   Weight to apply on splice sites (e.g. canonical v.s. non-canonical, how far from expected location, etc.) in the scoring of gene models.

start_codon_search_window

    Environment variable:          VIGOR_START_CODON_SEARCH_WINDOW
    System.property:               vigor.start_codon_search_window
    Settable levels                virus config,gene config/variant
    Description:                   Number of nucleotides before and after a candidate site to check for a start codon

start_codons

    Environment variable:          VIGOR_START_CODONS
    System.property:               vigor.start_codons
    Settable levels                virus config,gene config/variant
    Description:                   Comma separated list of expected start codons

start_score_factor

    Environment variable:          VIGOR_START_SCORE_FACTOR
    System.property:               vigor.start_score_factor
    Settable levels                virus config,gene config/variant
    Description:                   Weight to apply on start codons (e.g. canonical v.s. non-canonical, how far from expected location, etc.) in the scoring of gene models.

stop_codon_readthrough

    Environment variable:          VIGOR_STOP_CODON_READTHROUGH
    System.property:               vigor.stop_codon_readthrough
    Settable levels                gene config/variant
    Description:                   Format is amino acid/offset/regex (e.g.stop_codon_readthrough=-11/R/[NATC][NATC][NT][NG][NA][NC][NTG][NAG][NATG][NATCG][NTC][NAG][NAG])

stop_codon_search_window

    Environment variable:          VIGOR_STOP_CODON_SEARCH_WINDOW
    System.property:               vigor.stop_codon_search_window
    Settable levels                virus config,gene config/variant
    Description:                   Number of nucleotides before and after a candidate site to check for a stop codon

stop_score_factor

    Environment variable:          VIGOR_STOP_SCORE_FACTOR
    System.property:               vigor.stop_score_factor
    Settable levels                virus config,gene config/variant
    Description:                   Weight to apply on stop codons (e.g. how far from expected location, etc.) in the scoring of gene models.

tiny_exon3

    Environment variable:          VIGOR_TINY_EXON3
    System.property:               vigor.tiny_exon3
    Settable levels                gene config/variant
    Description:                   Tiny exon 3. Format is regex:[offset]

tiny_exon5

    Environment variable:          VIGOR_TINY_EXON5
    System.property:               vigor.tiny_exon5
    Settable levels                gene config/variant
    Description:                   Tiny exon 5. Format is regex:[offset]

### Mature Peptide Files

## Current status
VIGOR database supports the following viruses:

* Influenza (A & B for human, avian, and swine, and C for human)
* SARS-CoV-2
* West Nile Virus (I and II)
* Zika Virus
* Chikungunya Virus
* Eastern Equine Encephalitis Virus
* Respiratory Syncytial Virus
* Rotavirus (A, B, C, F, and G)
* Enterovirus
* Lassa Mammarenavirus
* Alphaviruses (VEEV and EEEV)
* Antennavirus
* Bandavirus
* Beidivirus
* Cicadellivirus
* Coguvirus
* Dengue Virus
* Embecovirus
* Feravirus
* Goukovirus
* Hantaviridae
* Hartmanivirus
* Hibecovirus
* Horwuvirus
* Hudivirus
* Hudovirus
* Inshuvirus
* Ixovirus
* Jonvirus
* Laulavirus
* Lentinuvirus
* Measles Virus
* Merbecovirus
* Mobuvirus
* Monkeypox Virus
* Nairoviridae
* Nobecovirus
* Orthophasmavirus
* Peribunyaviridae
* Phasivirus
* Phlebovirus
* Pidchovirus
* Reptarenavirus
* Rubodvirus
* Sarbecovirus
* Sapovirus
* Sawastrivirus
* Tenuivirus
* Uukuvirus
* Wenrivirus
* Wuhivirus


