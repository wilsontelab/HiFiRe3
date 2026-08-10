# HiFiRe3 Tool Suite

**HiFiRe3** ("high fire") carries pipelines and apps for analyzing
sequencing data from libraries that achieve **Hi**gh **Fi**delity 
error-corrected scoring of single nucleotide variants (SNVs) and 
structural variants (SVs) by controlling input fragment sizes and/or
**Re**striction enzyme-mediated **Re**duced **Re**presentation.

**IMPORTANT**: HiFiRe3 assumes that any insert size selection during 
library preparation was performed _before adapter ligation_, targeting 
molecules from 1N to 2N bp, and that input libararies are either 
PCR-free or used pooled PCR with unique dual indices to minimize 
template switching.

The steps to using HiFiRe3 are to:
- install the codebase in this repository
- obtain or build the required runtime environment
- if needed, basecall reads using `basecall PacBio` or `basecall ONT`
- align and analyze basecalled reads using `analyze fragments`
- characterize variants:
    - in single samples using `analyze SVs` and/or `analyze SNVs`
    - comparing across multiple samples using `compare SVs` and/or `compare SNVs`
    - comparing across muliple library types using `merge SVs`
- visualize results in the R Shiny apps

## Single-suite installation (recommended)

HiFiRe3 is implemented in the
[Michigan Data Interface](https://midataint.github.io/) (MDI)
for developing, installing and running Stage 1 HPC **pipelines** 
and Stage 2 interactive web **apps**. For HiFiRe3, we recommend 
a single-suite installation, which is accomplished by:
- cloning this tool suite repository
- running _install.sh_
- optionally running _alias.pl_ to create an `hf3` alias to the command line interface (CLI)

See the bottom of this README for information on an alternative 
container-only way of using of HiFiRe3 pipelines.

### Install this tool suite

Run the following to install the tool suite and command line interface (CLI).

```bash
git clone https://github.com/wilsontelab/HiFiRe3.git
cd HiFiRe3
./install.sh 1
```

### Create an alias to the command line interface (optional)

```bash
# you can use a different alias if you'd like, e.g., replace hf3 with HiFiRe3
perl alias.pl hf3 
```

Answer 'y' to add the alias to your bash profile, then 
reload a new shell to activate the alias for use.

If you prefer not to use an alias, 
add the installation directory to your PATH variable
or `cd` into the directory prior to calling `./hf3`.

### Test the command line interface (CLI)

For help, call the CLI with no arguments, which describes the format for pipeline calls. 

```bash
# use the alias, if you created it as described above
hf3 
hf3 --help

# or call the CLI directly without an alias
cd HiFiRe3
./hf3
./hf3 --help
```

## Obtain or build the required runtime environment

HiFiRe3 pipelines use version-controlled 3rd-party software built into a 
[conda](https://docs.conda.io/)
runtime environment. There are two ways to obtain or create that environment.

### Use Singularity containers, i.e., Apptainers (recommended)

HiFiRe3 supports Singularity containers that have the required
conda environments pre-installed. No action is needed to support their use
if either the `singularity` or `module load singularity` command is available 
on your server - the containers will download and be used automatically.

If you wish to use containers but need to run a different command
to make `singularity` available, follow the instructions in
`.../mdi/config/singularity.yml` to communicate that command to the CLI.

### Build the Conda environments locally

If you don't have Singularity or Apptainer available on your system,
or prefer or need to build the environments yourself, you can
build them in your HiFiRe3 installation using
[micromamba](https://mamba.readthedocs.io/)
as follows:

```sh
hf3 analyze conda --create
```

In a shared server environment, the environment build command may get killed by the host.
If that happens, run the command on a cluster worker node with sufficient resources, e.g.:

```sh
# example for a Slurm-based cluster server
salloc --account <your_slurm_account> --cpus-per-task 4 --mem-per-cpu 4G 
hf3 analyze conda --create
exit
```

### Communicate your environment choice using option `--runtime`

Option `--runtime` defaults to value 'auto', which will prefer to use
Singularity containers and fall back to locally built environments only if 
that fails. To force the use of a runtime environment, set the 
option to either `--runtime singularity` or `--runtime conda`.

## Execute a pipeline from the command line

### Job files

HiFiRe3 pipelines can be called entirely using the CLI introduced above. However, you 
are encouraged to create YAML-format job configuration files that define the
parameters for your job and execution steps.

See [the templates folder](https://github.com/wilsontelab/HiFiRe3/tree/main/templates)
for job file templates for all HiFiRe3 pipelines
and actions, and <https://midataint.github.io/mdi/docs/job_config_files.html>
for extended help on using job files. Job file templates can also be generated with 
command `hf3 <pipeline> template`, e.g., `hf3 basecall template`.

The HiFiRe3 CLI and job files can run pipeline actions either 
inline in the calling shell or by submitting jobs to your server job scheduler.
The latter is recommended for most use cases. Thus, our most common usage pattern is:

```sh
hf3 inspect myJob.yml          # check the formatting of your job file
hf3 mkdir myJob.yml            # create any missing output directories
hf3 submit --dry-run myJob.yml # test the job file to see what will happen
hf3 submit myJob.yml           # submit the job to Slurm or your scheduler
hf3 myJob.yml status           # show the state of all submitted jobs
hf3 myJob.yml top              # monitor a running job
hf3 myJob.yml report           # show a job log report
hf3 myJob.yml ls               # show the contents of a job's output diretory
```

### Workflow sequence

HiFiRe3 has _pipelines_ each with associated _actions_, 
listed here in execution order of the most common actions:
- `basecall ONT` or `basecall PacBio` (where `basecall` is a _pipeline_ and `ONT` is an _action_)
- `analyze fragments`
- `analyze SVs`  or `compare SVs`
- `analyze SNVs` or `compare SNVs`
- `merge SVs`

Required/common options are described below; use 
`hf3 <pipeline> <action> --help` or `hf3 <pipeline> template` 
for complete option information, or see the action help 
[here](https://github.com/wilsontelab/HiFiRe3/tree/main/options).

### Universally required options

Options `--output-dir/-O` and `--data-name/-N` are required by all pipeline actions.
Output files are placed into directory `<--output-dir>/<--data-name>`.

Option `--genome` (default hg38) is required by nearly all actions. 
This is usually all that is needed, but see below for special use cases. 

## Input read data

HiFiRe3 accepts read input files from the following sequencing platforms,
read configurations, and library types, communicated via options 
`--sequencing-platform` and `--library-type`:

| --sequencing-platform | --library-type |
|---|---|
| Illumina_2x150 | Nextera |
| Illumina_2x150 | TruSeq |
| Aviti_2x150 | Nextera |
| Aviti_2x150 | TruSeq |
| Aviti_2x150 | Elevate |
| Aviti_1x300 | Nextera |
| Aviti_1x300 | TruSeq |
| Aviti_1x300 | Elevate |
| Ultima | Solaris |
| ONT | Rapid |
| ONT | Ligation |
| PacBioMolecule | HiFi |
| PacBioStrand | HiFi |

All tools assume any insert size selection was performed _before adapter ligation_. 
However, the pipeline is flexible with regard to which error correction methods you 
exploit, e.g., you can forego size selection and/or RE-mediated DNA fragmentation as 
suits your needs. 

Different library types use restriction enyzmes in different ways.
Most commonly, adapters are ligated onto DNA fragments generated by a blunt RE (ligFree).
Alternatively, tagmentation can be applied to initial genomic DNA fragments generated 
with a 5' overhanging RE and blocked by ddNTP filling (tagFree).

## Recommended HiFiRe3 workspace organization

The following is an optional but time-tested strategy for organizing input,
output, job, and resource files in your HiFiRe3 workspace 
(create *** folders manually as needed).

```sh
HiFiRe3                    # root folder you created above
├── input                  # *** folder for your input data files
│   └── project1           # *** subfolder for a related set of samples
├── jobFiles               # *** folder for your job configuration files
│   └── project1
├── mdi                    # MDI codebase folder created by HiFiRe3 installation
│   ├── config             # folder with configuration files you may need
│   └── resources/genomes  # folder where genome files are placed by default
└── output                 # *** folder for your output data files
    └── project1
```

## Basecalling long reads

The `basecall` pipeline finishes processing PacBioStrand or ONT
reads into the unaligned BAM format required for later pipeline actions. 

```sh
hf3 basecall --help
hf3 basecall PacBio --help
hf3 basecall ONT --help
hf3 basecall ...
```

For short-read libraries, input reads will have been basecalled by 
the sequencing platform and you should proceed directly to fragment analysis.

### PacBioStrand, from unaligned single-strand consensus bam files

For PacBio HiFiRe3 with sufficiently small inserts, initial basecalling is 
performed during your Revio or Vega run, which must be configured with the 
following settings, identified in HiFiRe3 as the **PacBioStrand** platform:
- Library Type = Standard (e.g., WGS)
- Insert Size = the middle of your selected size range, i.e., (N + 2N) / 2 = 1.5N
- Movie Acquisition Time = at least 24 hours (or the maximum allowed)
- Data Options:
    - Include Base Kinetics = YES (or NO, if you do not intend to look for damaged bases)
    - Consensus Mode = Strand **<<<< IMPORTANT!**

These settings ensure you will have:
- sufficient time to get around shorter HiFiRe3 inserts _more times than a normal HiFi library_ 
- consensus output files generated for each strand independently

Point `basecall PacBio` at your strand consensus unaligned BAM files using
option `--pacbio-dir`. HiFiRe3 compares the strand consensuses to each other 
and to the reference genome to assign a final base call and quality score. 

When the two strands agree on a base, that base is reported with the highest 
quality score of the two strands regardless of the reference alignment.

When heteroduplex DNA is detected, it might result from:
- a consensus basecalling error on one strand but not the other (usually an indel)
- a base mismatch characterized by two high-quality but non-Watson-Crick paired bases
- a damaged base on one of the two strands opposite an undamaged base

If one heteroduplex strand matches the reference, that base is reported with 
its quality score. If neither base matches the reference, an N base and zero
quality are reported. In either case, the values of the mismatched bases and 
kinetic information are reported in 
[SAM/BAM tags enumerated here](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/formats/hf3_tags.rs).

PacBioStrand reads processed in this way can be used for error-corrected calling of 
both SVs and SNVs. Alternatively, larger PacBioMolecule reads can be subjected to 
normal on-device basecalling with molecule-level consensuses built from both strands
if only error-corrected SV calling is needed. The `basecall PacBio` action is not 
needed for PacBioMolecule reads.

### ONT, from POD5 files

HiFiRe3 must perform basecalling of ONT reads from primary trace data 
to make best use of the fixed bases at RE-cleaved DNA ends. ONT reads cannot 
support high-fidelity SNV calling, but fully support error-corrected calling 
of rare mosaic SVs.

The required input for ONT is one or more POD5 format files obtained from your
nanopore sequencing run. Use option `--pod5-dir` to point at your files in your 
call to `basecall ONT`. HiFiRe3 uses 
[Dorado](https://github.com/nanoporetech/dorado)
for high accuracy basecalling and custom scripts for RE-aware adapter trimming 
that preserves ONT channel information. See the options help for information 
on how to tune ONT basecalling and/or to invoke modified basecalling.

### HiFiRe3 tags in unaligned bam files

The final output of either HiFiRe3 `basecall` action is one or more unaligned
BAM files in folder `<--output-dir>/<--data-name>/ubam`. Important metatadata 
specific to HiFiRe3 or derived from vendor files is propagated to output files as
[SAM/BAM tags enumerated here](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/formats/hf3_tags.rs).

## Shared fragment analysis upstream of SV and SNV calling

The `analyze` pipeline aligns and processes reads with the goal of finding
rare mosaic variants, especially singleton variants, with extremely low baseline 
artifact rates for SVs and, for PacBioStrand, SNVs/indels.

The `analyze fragments` action filters, aligns and indexes reads to prepare for 
variant analysis. Various filters ensure that only high-quality conforming reads 
and bases are used for variant calling downstream. Indexing allows for efficient 
comparison of reads to RE sites and read recovery during visualization.

```sh
hf3 analyze --help
hf3 analyze fragments --help
hf3 analyze fragments ...
```

Custom metadata generated by fragment analysis is propagated to output BAM files as
[SAM/BAM tags enumerated here](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/formats/hf3_tags.rs).
These tags codify information on RE site matching, size distributions, alignment 
and junction quality filterting, and more.

### Read alignment using minimap2

The first step in `analyze fragments` aligns reads to the reference genome using 
[minimap2](https://github.com/lh3/minimap2), 
for both short and long reads. Use option --`read-file-dir` and `--read-file-prefix` to communicate 
the path to your input read files if they were not generated with `basecall PacBio` or `basecall ONT`.

### Inferring RE site locations in reduced representation reads

When option `--enzyme-name` is set to something other than NA, HiFiRe3 expects reads
to arise from (a subset of) all possible RE fragments in a genome.

For adapter ligation (ligFree) libraries, sequencing reads end at RE sites.
Each fragment is sequenced enough times to establish its consensus genotype without 
external information. This is known as a _reduced representation_ libary. While 
most fragments can be predicted from the _in silico_ RE site analysis, others 
will be genotype-specific due to RFLPs. Therefore, `analyze fragments` uses the 
ends of sequenced inserts to locate recurring sample-specific RE sites. 

Clonal SNV-containing fragments will be among the set of located fragments if 
they are within the (selected) insert size range, where reference and 
SNV/indel-containing reads will be ~the same size.

In contrast, rare mosaic SV-containing fragments may not end at located RE sites 
depending on the relative sizes of clonal and SV fragments, i.e.,  an SV fragment 
may be in a library's insert size range even though the parental fragment(s) it 
bridges are not. For this reason, both _in silico_ and located RE sites are used 
when matching read outer ends to RE sites. Any fragment ending at expected RE sites 
of either type are used for SV discovery, which, unlike clonal fragments, can be 
located ~anywhere in the genome or target regions.

These steps do not apply and are skipped for randomly sheared library inputs, 
even when tagmentation is applied to inital RE-cleaved fragments (tagFree). 
Reference RE sites at junctions are still used to filter against chimeric SVs, 
but read outer endpoints are not expected to match RE sites. 

You can also skip RFLP assessment by setting option `--skip-rflp-detection`, 
in which case only reference RE sites will be tracked during SV error correction.
Alternatively, you can set option `--site-override-file` to a file containing
a properly constructed HiFiRe3 RE site list instead of performing _de novo_ 
RE site discovery from the current sample's reads.

## Clonal and mosaic variant calling

Depending on your application and sequencing platform, you may wish to call 
SVs, SNVs/indels, or both from your processed reads, which is done in distinct 
actions of the `analyze` pipeline:

```sh
hf3 analyze SVs --help
hf3 analyze SVs ...
hf3 analyze SNVs --help
hf3 analyze SNVs ...
```

Error-corrected SNV calling requires sufficiently small PacBioStrand reads 
processed using `basecall PacBio` as described above. The process uses
the tags added during basecalling to determine which bases of a read are allowed
to call SNVs and indels based on homoduplex strand validation. SNV calling
occurs in two passes. The first pass includes all reads, even non-duplex reads,
for maximal sensitivity for clonal variant calling. The second pass includes 
only duplex reads for maximal specificity for rare mosaic variant calling.

SV error correction enforces filters against end-to-end chimeras that may variably 
reject junctions:
- with adapters present in inserted non-reference bases (always applied)
- flanked by low quality bases or alignments (always applied)
- that match a known RE fragment endpoint (if a RE was used to fragment the genomic DNA)
- more than 1N bp away from either insert end (if size selection was performed before adapter ligation)

Use options `--sequencing-platform` and `--library-type` to tell the pipeline how
to perform adapter and alignment quality checks based on the nature of the reads.

Use option `--enzyme-name` to communicate the RE used to prepare a HiFiRe3 library
upstream of adapter ligation or tagmentation. Leave `--enzyme-name`as the 
default value 'NA' if your library is a sheared DNA library.

Use options `--min-selected-size`, `--selected-size-cv`, and `--min-allowed-size` 
to communicate the selected DNA insert size range. Leave both `--min-selected-size` and 
`--min-allowed-size` at the default of 0 if size selection was not performed. Otherwise,
you will most commonly set `--min-selected-size` to communicate the lower insert size cutoff,
which should be strictly established, e.g., using a BluePippin device. If needed, the values
calculated from `--min-selected-size` and `--selected-size-cv` can be overridden by setting
`--min-allowed-size` to a non-zero value. In either usage, these options establish the
1N insert size value used to enforce size-based SV error correction.

RE-based and size-based SV error correction methods are substantially redundant but 
independent, so both can be used to achieve maximal SV error correction. The pipelines
will also call SVs when neither error correction method was in use.

Importantly, each SV error correction method only allows the pipeline to reliably reject 
end-to-end chimeras. Many middle-to-middle chimeras arising from DNA fragmentation after 
size-selection cannot be error corrected, so care must be taken to aggressively limit 
DNA fragmentation that occurs between size selection and adapter ligation.

The final output of both variant calling actions is a data package suitable for
loading into the interactive visualization app, in additional variant-specific files.
Many output files have the `.bgz` extension to indicate that they are bgzipped
and tabix-indexed for random access data retrieval. They can be decompressed with either
`bgzip` or `gzip`.

### SV-specific output files

The primary output of interest from SV variant calling is the "final junctions"
file named `XXXX.analysis.final_junctions_1.txt.bgz`, a custom tab-delimited format. 

Addtional SV output files are mainly intended for use in the R Shiny app and 
include alignment coverage maps and read-level metadata for 
all SV junction-containing reads, including SEQ, QUAL, CIGAR fields.

Details on these file formats, including column definitions and data types,
can be found in the Rust structures that define them:
- [final_junctions table format](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/junctions/junction.rs#L277)
- [read_path table format](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/junctions/read_path.rs#L22)
- [alignment map table format](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/junctions/alignment.rs#L24)

Notably, final junction files retain ALL junctions detected in all (on-target) reads 
that passed read-level quality filtering, including junctions identified as chimeric
artifacts, single-molecule junctions, etc. This supports comparisons of the properties 
of real vs. artifact junctions. Filtering the table yields lists of SV junctions of 
greatest interest for different applications, e.g., clonal vs. mosaic SVs. 

Some important filtering information is found in bit-encoded alignment-level and 
junction-level "failure flags", where a bit that is set means that quality check failed. 
Please see this summary of the 
[failure flag bit encoding](https://github.com/wilsontelab/HiFiRe3/blob/main/shiny/shared/session/utilities/jxn_filters.R#L16)
for details.

### SNV-specific output files

The primary outputs of SNV/indel calling are:
- two variant list files named `XXXX.<READ_LEVEL>.snv_indel.txt.bgz`
- two genome pileup files named `XXXX.<READ_LEVEL>.pileup.bed.bgz`.

Each file type has two variants where READ_LEVEL is replaced with either:
- `all_reads`, where non-duplex reads were included during variant calling
- `error_corrected`, where only duplex error-corrected reads were used

The all_reads variant list provides the most sensitive and accurate detection 
of clonal variants, where the error_corrected list provides the most
specific asssessment of rare variants. Pileup files are used 
for visualization in the app TrackBrowser.

As with SV files, ALL variants from the indicated reads are included, even 
variants that will eventually be filtered away as not passing homoduplex 
strand validation or having too few or too many observations.

Details on these file formats can be found in the Rust structures that define them:
- [snv_indel variant list format](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/snvs/variant.rs#L35)
- [pileup file format](https://github.com/wilsontelab/HiFiRe3/blob/main/shared/crates/hf3_tools/src/snvs/pileup.rs#L34)

**PENDING: conversion of HiFiRe3 file formats to VCF format

## Comparing variants across samples

HiFiRe3 applications often seek to distinguish single-molecule variants,
i.e., SVs and SNVs detected only once in a sample, from those detected
multiple times. Multiply-sequenced variants must have arisen from different 
input DNA molecules in PCR-free libraries, whereas single-molecule variants 
are consistent with mutations that occurred during an experiment or recently 
in a tissue or cell lineage.

Sensitivity for detecting multiply-sequenced variants, including variants
that arose in a clonal ancestor of the samples under study, increases if data  
from many samples from the same source genotype can be combined, which is the
job of `compare SVs` and `compare SNVs`. 

The processes and file formats for comparing variants across samples are the 
same as those applied to a single sample above. Aligned reads from multiple samples 
are analyzed together with tracking of the samples that contributed to each called 
variant. This approach requires that all input samples are of the same library 
type, including sequencing platform and restriction enzyme, and that `analyze fragments`
has been run on all samples individually. However, `analyze SVs` and `analyze SNVs` 
are not required to run `compare SVs` or `compare SNVs`, as variant calling is 
repeated anew.

To compare SVs across different types of libaries, e.g., to compare the final junction 
tables of short and long reads from the sample source, use `merge SVs` instead. 
Variant merging does not apply to SNVs since only PacBioStrand libraries support 
SNV calling.

## Targeted genomic analysis

Many error-corrected sequencing approaches benefit from focusing reads to specific
genomic regions, e.g., by hybridization capture or ONT adaptive sampling.
If applicable, use options `--targets-bed` and `--region-padding` to communicate 
the genomic regions to which your reads were targeted. 

During SV calling, HiFiRe3 pipelines only consider reads where the 5' end 
aligned to a target region. This approach is critical for ONT adaptive sampling
where off-target 5' ends will be present in the library in a truncated form that
should not be used for equivalent SV calling to on-target reads.

Note that region targeting is compounded on top of RE-mediated reduced representation
sampling to allow considerable focusing of reads onto specific genomic fragments. 

## Launch the interactive apps server

Once all pipelines have finished running, you will want to view and interact
with your data.

To install and launch the HiFiRe3 apps server, we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.

After following the instructions to run the Desktop on your local machine
or server, load a HiFiRe3 data package file ending in `.mdi.package.zip`,
into the app interface. If possible, the relevant container
will automatically be downloaded and used to run the apps server.

## Download and digest a reference genome (special cases only)

HiFiRe3 requires a reference genome assembly and supporting annotation files,
including tables of _in silico_ restriction enzyme (RE) sites.
Usually, these files are automatically prepared on first use of a given `--genome`.
However, a few special use cases require you to prepare genome files manually.

### Manually prepare a genome before running jobs

You may want or need to prepare the genome in advance as follows:

```sh
hf3 prepare genome --help
# substitute hs1 with your desired genome
hf3 prepare genome --genome hs1 --output-dir $PWD/_prepare_hs1 --data-name hs1
```

Even better, use the 
[templates/prepare_genome.yml](https://github.com/wilsontelab/HiFiRe3/blob/main/templates/prepare_genome.yml)
job file template. 

### Optionally create a composite genome for mixed species analysis

Quality assessment of SV artifact rates can be enhanced by mixing samples from 
two species prior to libary preparation and sequencing. The `combine genomes` 
pipeline action creates the required composite reference, which should then be 
provided as option `--genome` in later analysis steps while also setting 
option `--is-composite-genome`.

```sh
hf3 combine genomes --help
# substitute hs1 and dm6 with your desired genomes
hf3 combine genomes --genome1 hs1 --genome2 dm6 --output-dir $PWD/_prepare_hs1_dm6 --data-name hs1_dm6
```

Combining genomes requires that each individual genome was previously prepared as above.
Even better, use the 
[templates/combine_genomes.yml](https://github.com/wilsontelab/HiFiRe3/blob/main/templates/combine_genomes.yml)
job file template.

The resulting composite reference assembly carries the original chromosome names
suffixed with the respective source genome, e.g., `chr1_hs1` and `chr4_dm6`.

## Obtain the required `hf3_tools` utility (special cases only)

Many actions in HiFiRe3 pipelines are executed by the `hf3_tools` utility, 
an executable binary compiled from 
[source Rust code in this repository](https://github.com/wilsontelab/HiFiRe3/tree/main/shared/crates/hf3_tools).

For most use cases, `hf3_tools` will be automatically downloaded from
GitHub on first use. However, two special cases require your action
to obtain or create the binary.

### Manually download the binary

You may want or need to pre-download the binary for the matching version of 
the tool suite code.

First, go to <https://github.com/wilsontelab/HiFiRe3/releases/latest>
to get the latest version number of the code you installed above, 
e.g., `v1.2.3`. Use that version number to complete the following:

```sh
VERSION=<version> # replace <version> with your desired version number, e.g., v1.2.3
URL=https://github.com/wilsontelab/HiFiRe3/releases/download/$VERSION/hf3_tools-x86_64-unknown-linux-gnu.tar.gz
DIR=bin/HiFiRe3/$VERSION
cd /path/to/HiFiRe3/mdi # the mdi folder within your HiFiRe3 installation
mkdir -p $DIR
curl -sLf ${URL} | tar -xz -C $DIR
```

### Compile the binary using Rust (developers only)

Developers most often need to compile the Rust code for themselves 
using the support features provided by the mdi-pipelines-framework.

```bash
cd /path/to/HiFiRe3/mdi/suites/developer-forks/HiFiRe3/shared/crates/hf3_tools # must compile from within the crate directory
hf3 analyze rust --help
hf3 analyze rust --create  1.92 # create a versioned Rust development environment
hf3 analyze rust --compile 1.92 # compile hf3_tools using the created environment
hf3 analyze rust --gcc "module load gcc/15.1.0" --compile 1.92 # if a command is required to make C compilers available 
```

It is also possible to compile the Rust code from first principles
if all prerequisites are met. The compiled executable binary must be 
copied into file `/path/to/HiFiRe3/mdi/bin/HiFiRe3/dev/hf3_tools`.

Alternatively, if you are not developing Rust code, you can download
or copy a valid `hf3_tools` binary into `/path/to/HiFiRe3/mdi/bin/HiFiRe3/dev`
without compiling it yourself. 

Use of the developer binary in the `dev` folder is activated using
the `hf3 -d` developer option on all HiFiRe3 calls.

## Alternative container-only method of using HiFiRe3 pipelines

Some users may only be interested in using HiFiRe3 pipelines
as a standalone program, e.g., if HiFiRe3 actions are to be
incorporated into a pre-existing workflow management system.
HiFiRe3 pipeline containers support this installation-free usage.

Importantly, running pipelines via a direct call to a container,
rather than via the CLI installed on your server, disables
all of the worflow/job manager helpers in favor of your 
pre-existing solution. All data processing tasks and configuration 
remain exactly the same since the container has the same installation
in it as described above.

### Download the relevant container

You do not need to clone or install this repository, simply download the relevant 
container image from:
- <https://github.com/wilsontelab/HiFiRe3/pkgs/container/hifire3> 

e.g., using command (update vX.X to an actual version, e.g., v0.5):
- `singularity pull oras://ghcr.io/wilsontelab/hifire3:vX.X`

### Use the container with directory bind mounts and direct action calls

As always, you must bind mount all relevant server paths to 
the container so that running pipeline jobs have access to your files. The 
recommended directory organization above makes this easy by requiring only one bind
mount of the HiFiRe3 root directory (even that is unnecessary
if you work from that directory as the current working directory is implicitly
mounted).

Then simply follow `singularity run <image>.sif` with  your pipeline, action, 
and options as you would for any typical data analysis program, e.g.:

```sh
# every relevant directory option should be prefixed with a bind mount directory
ROOT_DIR=/path/to/HiFiRe3
singularity run \
    --bind $ROOT_DIR \
    <image>.sif \
        <pipeline> <action> \
            --tmp-dir        $ROOT_DIR/tmp \
            --output-dir     $ROOT_DIR/output \
            --data-name      my-data-name \
            --option-name    123
            # etc.
```

The job report log is printed to STDOUT for you to consume as needed in your workflow.
