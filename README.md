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

AAOverlap_offset

        Environment variable:          VIGOR_AAOVERLAP_OFFSET
        System.property:               vigor.AAOverlap_offset
        Settable levels                gene,virus

alignment_score_factor

        Environment variable:          VIGOR_ALIGNMENT_SCORE_FACTOR
        System.property:               vigor.alignment_score_factor
        Settable levels                gene,virus
        Description:                   weight of alignment score in evaluating alignment

alternate_startcodon

        Environment variable:          VIGOR_ALTERNATE_STARTCODON
        System.property:               vigor.alternate_startcodon
        Settable levels                gene
        Description:                   Alternate start codons for gene. Format is CODON[,CODON,..]

candidate_selection

        Environment variable:          VIGOR_CANDIDATE_SELECTION
        System.property:               vigor.candidate_selection
        Settable levels                gene,virus
        Description:                   TODO

circular_gene

        Environment variable:          VIGOR_CIRCULAR_GENE
        System.property:               vigor.circular_gene
        Settable levels                gene,virus
        Description:                   Gene is circular

complete_gene

        Environment variable:          VIGOR_COMPLETE_GENE
        System.property:               vigor.complete_gene
        Settable levels                gene,virus
        Description:                   Gene is complete

excludes_gene

        Environment variable:          VIGOR_EXCLUDES_GENE
        System.property:               vigor.excludes_gene
        Settable levels                gene
        Description:                   Excludes gene. TODO

frameshift_sensitivity

        Environment variable:          VIGOR_FRAMESHIFT_SENSITIVITY
        System.property:               vigor.frameshift_sensitivity
        Settable levels                gene,virus
        Description:                   How to handle frameshifts

is_optional

        Environment variable:          VIGOR_IS_OPTIONAL
        System.property:               vigor.is_optional
        Settable levels                gene
        Description:                   Gene is optional for valid model

is_required

        Environment variable:          VIGOR_IS_REQUIRED
        System.property:               vigor.is_required
        Settable levels                gene
        Description:                   Gene is required for valid model

jcvi_rules

        Environment variable:          VIGOR_JCVI_RULES
        System.property:               vigor.jcvi_rules
        Settable levels                gene,virus

leakystop_notFound_score

        Environment variable:          VIGOR_LEAKYSTOP_NOTFOUND_SCORE
        System.property:               vigor.leakystop_notFound_score
        Settable levels                gene,virus

leakystop_score_factor

        Environment variable:          VIGOR_LEAKYSTOP_SCORE_FACTOR
        System.property:               vigor.leakystop_score_factor
        Settable levels                gene,virus

locus_tag

        Environment variable:          VIGOR_LOCUS_TAG
        System.property:               vigor.locus_tag
        Settable levels                gene,virus
        Description:                   Locus tag prefix to use in output

matpepdb

        Environment variable:          VIGOR_MATPEPDB
        System.property:               vigor.matpepdb
        Settable levels                gene
        Description:                   Mature peptide database to use for gene

mature_pep_mincoverage

        Environment variable:          VIGOR_MATURE_PEP_MINCOVERAGE
        System.property:               vigor.mature_pep_mincoverage
        Settable levels                gene,virus

mature_pep_minidentity

        Environment variable:          VIGOR_MATURE_PEP_MINIDENTITY
        System.property:               vigor.mature_pep_minidentity
        Settable levels                gene,virus

mature_pep_minsimilarity

        Environment variable:          VIGOR_MATURE_PEP_MINSIMILARITY
        System.property:               vigor.mature_pep_minsimilarity
        Settable levels                gene,virus

max_align_merge_aa_gap

        Environment variable:          VIGOR_MAX_ALIGN_MERGE_AA_GAP
        System.property:               vigor.max_align_merge_aa_gap
        Settable levels                gene,virus

max_exon_size

        Environment variable:          VIGOR_MAX_EXON_SIZE
        System.property:               vigor.max_exon_size
        Settable levels                gene,virus
        Description:                   Maximum sequence length of an exon

max_gene_overlap

        Environment variable:          VIGOR_MAX_GENE_OVERLAP
        System.property:               vigor.max_gene_overlap
        Settable levels                gene,virus
        Description:                   In reporting gene models, maximum overlap of genes allowed.

max_intron_size

        Environment variable:          VIGOR_MAX_INTRON_SIZE
        System.property:               vigor.max_intron_size
        Settable levels                gene,virus
        Description:                   Maximum sequence length of an intron

min_candidate_pctsimilarity

        Environment variable:          VIGOR_MIN_CANDIDATE_PCTSIMILARITY
        System.property:               vigor.min_candidate_pctsimilarity
        Settable levels                gene,virus

min_candidate_sbjcoverage

        Environment variable:          VIGOR_MIN_CANDIDATE_SBJCOVERAGE
        System.property:               vigor.min_candidate_sbjcoverage
        Settable levels                gene,virus

min_exon_size

        Environment variable:          VIGOR_MIN_EXON_SIZE
        System.property:               vigor.min_exon_size
        Settable levels                gene,virus
        Description:                   Minimum sequence length of an exon

min_functional_length

        Environment variable:          VIGOR_MIN_FUNCTIONAL_LENGTH
        System.property:               vigor.min_functional_length
        Settable levels                gene
        Description:                   Minimum functional length

min_gene_coverage

        Environment variable:          VIGOR_MIN_GENE_COVERAGE
        System.property:               vigor.min_gene_coverage
        Settable levels                gene,virus
        Description:                   Minimum coverage of genes

min_gene_size

        Environment variable:          VIGOR_MIN_GENE_SIZE
        System.property:               vigor.min_gene_size
        Settable levels                gene,virus
        Description:                   Minimum sequence length to be considered as a gene

min_intron_size

        Environment variable:          VIGOR_MIN_INTRON_SIZE
        System.property:               vigor.min_intron_size
        Settable levels                gene,virus
        Description:                   Minimum sequence length of an intron

min_missing_AA_size

        Environment variable:          VIGOR_MIN_MISSING_AA_SIZE
        System.property:               vigor.min_missing_AA_size
        Settable levels                gene,virus
        Description:                   TODO

min_pseudogene_coverage

        Environment variable:          VIGOR_MIN_PSEUDOGENE_COVERAGE
        System.property:               vigor.min_pseudogene_coverage
        Settable levels                gene,virus

min_pseudogene_identity

        Environment variable:          VIGOR_MIN_PSEUDOGENE_IDENTITY
        System.property:               vigor.min_pseudogene_identity
        Settable levels                gene,virus

min_pseudogene_similarity

        Environment variable:          VIGOR_MIN_PSEUDOGENE_SIMILARITY
        System.property:               vigor.min_pseudogene_similarity
        Settable levels                gene,virus

min_seq_gap_length

        Environment variable:          VIGOR_MIN_SEQ_GAP_LENGTH
        System.property:               vigor.min_seq_gap_length
        Settable levels                gene,virus
        Description:                   Minimum length of a sequence gap to ....

nancanonical_splicing

        Environment variable:          VIGOR_NANCANONICAL_SPLICING
        System.property:               vigor.nancanonical_splicing
        Settable levels                gene,virus
        Description:                   Alternate splice sites

NTOverlap_offset

        Environment variable:          VIGOR_NTOVERLAP_OFFSET
        System.property:               vigor.NTOverlap_offset
        Settable levels                gene,virus


relax_align_merge_aa_gap

        Environment variable:          VIGOR_RELAX_ALIGN_MERGE_AA_GAP
        System.property:               vigor.relax_align_merge_aa_gap
        Settable levels                gene,virus

ribosomal_slippage

        Environment variable:          VIGOR_RIBOSOMAL_SLIPPAGE
        System.property:               vigor.ribosomal_slippage
        Settable levels                gene
        Description:                   Ribosomal slippage. Format is offset/frameshift/regex

rna_editing

        Environment variable:          VIGOR_RNA_EDITING
        System.property:               vigor.rna_editing
        Settable levels                gene
        Description:                   RNA editing. Format is offset/regex/insertion string/note

shared_cds

        Environment variable:          VIGOR_SHARED_CDS
        System.property:               vigor.shared_cds
        Settable levels                gene
        Description:                   Shared CDS. Format is CDS[,CDS]

splice_form

        Environment variable:          VIGOR_SPLICE_FORM
        System.property:               vigor.splice_form
        Settable levels                gene,virus
        Description:                   Splice form. Format is ([ie]\\d+)+

splicing_score_factor

        Environment variable:          VIGOR_SPLICING_SCORE_FACTOR
        System.property:               vigor.splicing_score_factor
        Settable levels                gene,virus

start_codon_search_window

        Environment variable:          VIGOR_START_CODON_SEARCH_WINDOW
        System.property:               vigor.start_codon_search_window
        Settable levels                gene,virus
        Description:                   Number of nucleotides before and after a candidate site to check for a start codon

start_score_factor

        Environment variable:          VIGOR_START_SCORE_FACTOR
        System.property:               vigor.start_score_factor
        Settable levels                gene,virus

StartCodons

        Environment variable:          VIGOR_STARTCODONS
        System.property:               vigor.StartCodons
        Settable levels                gene,virus
        Description:                   Comma separated list of expected start codons

stop_codon_readthrough

        Environment variable:          VIGOR_STOP_CODON_READTHROUGH
        System.property:               vigor.stop_codon_readthrough
        Settable levels                gene
        Description:                   format is amino acid/offset/regex

stop_codon_search_window

        Environment variable:          VIGOR_STOP_CODON_SEARCH_WINDOW
        System.property:               vigor.stop_codon_search_window
        Settable levels                gene,virus
        Description:                   Number of nucleotides before and after a candidate site to check for a stop codon

stop_score_factor

        Environment variable:          VIGOR_STOP_SCORE_FACTOR
        System.property:               vigor.stop_score_factor
        Settable levels                gene,virus

tiny_exon3

    Environment variable:          VIGOR_TINY_EXON3
    System.property:               vigor.tiny_exon3
    Settable levels                gene
    Description:                   Tiny exon 3. Format is regex:[offset]

tiny_exon5

        Environment variable:          VIGOR_TINY_EXON5
        System.property:               vigor.tiny_exon5
        Settable levels                gene
        Description:                   Tiny exon 5. Format is regex:[offset]

### Mature Peptide Files

## Current status

The databases are in the process of migration from VIGOR3 to VIGOR4, but the following virus are considered production ready:

* Influenza (A & B for human, avian, and swine)
* West Nile Virus
* Zika Virus
* Chikungunya Virus
* Eastern Equine Encephalitis Virus
* Respiratory Syncytial Virus
* Rotavirus
* Enterovirus

