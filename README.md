[![codecov](https://codecov.io/gh/alnfedorov/reat/branch/main/graph/badge.svg?token=9TF8UQD0PD)](https://codecov.io/gh/alnfedorov/reat)
![license](https://img.shields.io/github/license/alnfedorov/reat)

### REAT

**REAT** is an RNA editing analysis toolkit designed with a focus on performance, low memory footprint, and ease of use.

Refer to the [Wikipedia](https://en.wikipedia.org/wiki/RNA_editing) for classification and basic overview of known
editing events. **REAT** can handle all types of edits (eg, A->I, C->U, etc.), except for editing by insertions or
deletions.

The target platform for **REAT** is an x86-64 computer running Linux; working on Windows or ARM is not guaranteed.

Please, feel free to open issues regarding bugs, installation issues, feature requests, or unintuitive behavior.
Quality, correctness, and ease of use are vital to us.

### Not released yet

**REAT** is under active development and is not yet been officially released. Despite good tests coverage (both unit and
integration), the tool still lacks user feedback. I.e., use at your own risk.

### Features

* Summarizing editing for provided regions or separate loci
* Efficient multithreading
* Strand prediction for unstranded libraries
* Autoref: simple yet useful inference of single nucleotide polymorphisms (SNP)
* Editing Index (EI) for the given set of ROI
* Flexible filtering options with reasonable default settings

See [details](#details) section for more in-depth explanation of some features.

### Limitations

Here is a list of known limitations; feel free to open an issue if one of them is critical for you:

* Potential RNA editing by insertions/deletions is ignored
* No Python or R interface
* Lack of support for ARM and Windows/macOS builds

### Installation

#### Pre-built binaries

For a quick start, you can try the pre-built **REAT** binary posted
via [GitHub releases](https://github.com/alnfedorov/REAT/releases). Besides the tagged versions, we also provide the
latest tested build from the [main branch](https://github.com/alnfedorov/REAT/releases/tag/latest).

Note that you will have to add executable permissions to the file: `sudo chmod +x reat`.

#### Cargo

Rust package manager _cargo_ is the recommended way to build **REAT** from sources:

```shell
cargo install --git https://github.com/alnfedorov/reat
reat --help
```

Follow [this](https://www.rust-lang.org/tools/install) page to install rust toolchain if you don't have one.

In addition, one will need CMake (for zlib-ng), which should be available in most package managers (
e.g, `apt install cmake`).

### Basic usage

**REAT** supports two modes: ROI-based and site-based.

#### ROI mode

ROI stands for Region of Interest, a **REAT** mode in which edits are summarized for each provided genomic region. For
instance, it is useful in the context of studying Alu repeats editing.

Example:

* _BED-like file with ROIs_(repeats.bed):

```text
chrY	76783396	76783529	RSINE   .    -   ... # The rest of the columns are ignored
... # Remaining regions omitted
```

* _Command:_

```shell
reat rois --input techrep1.bam techrep2.bam --rois repeats.bed \
          --reference GRCm38.fa --stranding "f/s" --threads $(nproc)
```

Output columns:

* **contig, start, end, strand, name** - ROI coordinate and name from the input file
* **trstrand** - transcription strand; predicted for unstranded libraries and deducted from the design for stranded
  experiments
* **coverage** - number of unique reads covering ROI (after applying all filters)
* **#X** - number of _X_ nucleotides in the sequence of a given ROI (always forward strand sequence)
* **X->Y** - the total number of events observed in a given ROI where a reference nucleotide _X_ was replaced by _Y_.
  That is, A->A is a number of A matches, and A->G denotes the total number of observed A->I edits

Note that **X->Y** notation always denotes matches/mismatches relative to the forward strand. For example, T->C
mismatches for reverse strand ROI are, in fact, A->G _RNA_ mismatches.

If autoref feature is enabled, edits are summarised after correcting for any potential SNPs.

#### ROI editing index

For a given set of regions, one can always calculate an editing index (EI) for all possible matches and mismatches.

For example, formally EI for mismatches A->G is defined as follows: EI(A->G) = P(transcribed ROIs are A->G edited) = (
forward(A->G) + reverse(T->C)) / (forward(A->A + A->C + A->G + A->T) + reverse(T->A + T->C + T->G + T->T)).

In other words, it is the probability that a given A in the provided ROI will be A->G edited after the transcription.

Note that the EI is strand independent as it is defined for the _transcribed_ regions of RNA. For example, to study A->I
editing, you need to refer exclusively to the A->G column and completely ignore T->C.

If a given set of ROIs represent all Alu repeats in the genome, then A->G EI is the so-called Alu Editing Index (AEI).

Example:

```shell
for bamfile in *.bam
do
  reat rois --input $bamfile  --rois alu-repeats.bed \
            --reference hg19.fasta --threads $(nproc) \
            --stranding "s/f" --annotation gencode.gff3.gz \
            --ei=AEI.csv --saveto /dev/null # append Alu editing index values and skip ROI results completely
done
```

Output columns:

* **name** - name of the experiment
* **ROI-file** - path to the file with rois
* **unstranded** - number of ROIs for which no transcription strand was deduced/predicted
* **X->Y** - editing index for _X->Y_ pair

One can call **REAT** multiple times with the same CSV file to append rows to the EI table.

#### Site mode

The **REAT** site-based mode is a classic scenario for estimating RNA editing for each genomic locus.

Example:

```shell
reat site --input rnaseq.bam \
          --reference hg38.fa --threads $(nproc) \
          --stranding "u" --annotation gencode.gff3.gz \
          --saveto /dev/stdout 
```

Output columns:

* **contig,pos** - coordinate of the locus
* **trstrand** - transcription strand; predicted for unstranded libraries and deducted from the design for stranded
  experiments
* **refnuc** - reference nucleotide from the FASTA assembly
* **prednuc** - predicted reference nucleotide(assembly nucleotide if autoref feature is disabled)
* **X** - the total number of sequenced nucleotides X; X is one of \[A, C, G, T\].

Similarly to the ROI mode, the reference and sequenced nucleotides **X** are always reported with respect to the forward
strand. That is, a minus strand locus with ten A's corresponds to ten sequenced T's from RNA fragments.

### Details

#### Strand prediction

To predict transcription strand for ROI/loci in unstranded experiments, **REAT** uses two strategies.

First, ROI/loci strand will be derived from the overlapping genes/exons strand(only works if genome annotation is
provided). This approach can be summarized in the following table:

| overlapping genes | overlapping exons | predicted strand |
|:-----------------:|:-----------------:|:----------------:|
|         +         |         +         |         +        |
|         -         |         -         |         -        |
|        +/-        |         +         |         +        |
|        +/-        |         -         |         -        |
|        +/-        |        +/-        |         .        |

That is, **REAT** checks overlapping genes first. If they are genes on the + and the - strand, exons are considered. In
the worst-case scenario, an unknown(`.`) strand is returned.

Second, for ROIs / loci for which **REAT** could not predict the strand from the annotation, **REAT** attempts to derive
the strand based on the observed A->I editing.

For the + strand transcripts, A->I edits are A->G mismatches, and for the - strand, T->C mismatches. Note that in many
cases, this heuristic fails (no A->I editing at all), and such ROIs / loci will be left unstranded in the final table.

#### Autoref

With sufficient coverage, we can automatically adjust the reference sequence for observed SNVs based on RNA-seq data.
**REAT** uses a simple heuristic for this: if coverage is sufficient and the frequency of the most abundant nucleotide
exceeds the threshold, then the reference nucleotide is the most abundant nucleotide at that locus.

Note that hyper-editing flag allows one to skip A->G and T->C corrections to explore potential hyperedited ROI/loci.

#### How `N`s are handled?

`N` is routinely used to indicate unknown nucleotides in assemblies and sequencing data. Here are a few notes on how `N`
s are handled by **REAT**:

* Ignored in sequencing data (reads). That is, _X->N_ mismatches are always ignored
* For unknown nucleotide `N`, all canonical nucleotides (A, C, G, T) are treated as mismatches
* In **rois** mode, `N` reference positions are skipped

Note that the above notes apply to `N`s after _Autoref_ (if enabled). That is, in most cases, `N`s will be replaced by
an appropriate nucleotide during the _Autoref_ pass.

#### What are include/exclude lists?

In short, these lists specify DNA regions that will be included or excluded from the analysis completely. I.e. counting 
only positions inside provided regions(include) or completely skipping them(exclude).

Regardless of the overlap with excluded / included regions, ROIs will be printed with their original coordinates and
names to make them distinguishable in the subsequent analysis. This is what makes usage of include/exclude regions
different from simply subtracting/intersting ROIs with them - original ROIs won't be splitted in the output.

[//]: # (### CLI arguments)

[//]: # ()
[//]: # (Here is a list of arguments&#40;`--help`&#41; supported by the Command Line Interface &#40;CLI&#41;:)

[//]: # ()
[//]: # (#### Shared args)

[//]: # ()
[//]: # (```)

[//]: # (Core:)

[//]: # (        --binsize <binsize>)

[//]: # (            Summarize mismatches per locus/ROI, providing each thread with genome bins&#40;job share&#41; of at most X base)

[//]: # (            pairs[default: 64000])

[//]: # ()
[//]: # (    -i, --input <input>...)

[//]: # (            Path to the input BAM file&#40;s&#41;. May contain a space-separated list of files, in which case they are treated)

[//]: # (            as technical replicates and pulled together)

[//]: # ()
[//]: # (    -n, --name <name>)

[//]: # (            Name of the run.[default: NA])

[//]: # ()
[//]: # (    -o, --saveto <saveto>)

[//]: # (            Path to the output tsv file. By default, the results are printed to stdout[default: /dev/stdout])

[//]: # ()
[//]: # (    -r, --reference <reference>)

[//]: # (            Indexed fasta file with a reference genome assembly. Contig / chromosome names must match the entries in the)

[//]: # (            BAM header &#40;s&#41;)

[//]: # ()
[//]: # (    -s, --stranding <stranding>)

[//]: # (            Strand-specificity of the experiment, i.e. matching between the read strand and the gene strand. Use "u" for)

[//]: # (            unstranded experiments; other available options based on the RSeQC nomenclature&#40;see infer_experiment.py)

[//]: # (            docs&#41;: same:"s" &#40;++,--&#41;, flip:"f" &#40;+-,-+&#41;, same read1/flip read2:"s/f" &#40;1++,1--/2+-,2-+&#41;, flip read1/same)

[//]: # (            read2:"f/s" &#40;1+-,1-+/2++,2--&#41;[possible values: u, s, f, s/f, f/s])

[//]: # ()
[//]: # (    -t, --threads <threads>)

[//]: # (            Maximum number of threads to spawn at once[default: 1])

[//]: # ()
[//]: # ()
[//]: # (Reads filtering:)

[//]: # (    -3, --trim3 <trim3>)

[//]: # (            Trim bases from the 3’ &#40;right&#41; end of each read before processing. Can be used to hard skip low-quality)

[//]: # (            bases at the end of reads if no trimming was done before / during the alignment.[default: 0])

[//]: # ()
[//]: # (    -5, --trim5 <trim5>)

[//]: # (            Trim bases from the 5’ &#40;left&#41; end of each read before processing. In particular, one can skip the first ~12)

[//]: # (            bases with non-random composition due to priming biases, what is a common anomaly in short-read RNA-seq)

[//]: # (            experiments.[default: 0])

[//]: # ()
[//]: # (        --ex-flags <ex-flags>)

[//]: # (            Exclude reads for which any of the specified BAM flags are set. For example, a value of 2820 will result in)

[//]: # (            skipping unmapped reads, supplementary and secondary alignments, reads that fail platform/vendor quality)

[//]: # (            checks. Use zero&#40;0&#41; to disable this filter[default: 2820])

[//]: # ()
[//]: # (        --in-flags <in-flags>)

[//]: # (            Include only reads for which all the specified BAM flags are set. For example, a value of 3 will result in)

[//]: # (            keeping only reads that were mapped in proper pairs. Use zero&#40;0&#41; to disable this filter[default: 0])

[//]: # ()
[//]: # (        --mapq <mapq>)

[//]: # (            Count only reads with mapq ≥ threshold. Note that reads with mapq = 255 are skipped by default&#40;mapq 255)

[//]: # (            means "not available" according to the SAM spec&#41;[default: 1])

[//]: # ()
[//]: # (        --mapq-255)

[//]: # (            Count reads with mapq=255. Useful for aligners that do not fully conform to the SAM specification &#40;e.g. STAR)

[//]: # (            with default parameters&#41;)

[//]: # ()
[//]: # (        --phread <phread>)

[//]: # (            Count only bases with phread ≥ threshold. For reference, Phread is defined as -10 log₁₀[error probability],)

[//]: # (            so Phread = 20 means 1 error in 100 base calls[default: 20])

[//]: # ()
[//]: # ()
[//]: # (Stranding:)

[//]: # (        --annotation <annotation>)

[//]: # (            Genome annotation in the GFF3 format. Genomic features &#40;exons and genes&#41; are used only to inference loci/ROI)

[//]: # (            strand based on the most likely direction of transcription &#40;see the GitHub documentation for details&#41;. It is)

[//]: # (            recommended to provide genome annotation for unstranded libraries, otherwise stranding will be highly)

[//]: # (            inaccurate.)

[//]: # ()
[//]: # (        --str-min-freq <str-min-freq>)

[//]: # (            Automatically predict strand based on the observed A->I editing for locus/ROI with A->G freq >= threshold)

[//]: # (            &#40;freq = ∑ A->G / &#40;∑ A->G + ∑ A->A&#41;&#41;[default: 0.05])

[//]: # ()
[//]: # (        --str-min-mismatches <str-min-mismatches>)

[//]: # (            Automatically predict strand based on the observed A->I editing for locus/ROI with A->G mismatches >=)

[//]: # (            threshold. It is a fallback strand prediction heuristic, used only for the unstranded libraries. Not)

[//]: # (            relevant for organisms without active ADAR-like enzymes.[default: 50])

[//]: # ()
[//]: # ()
[//]: # (Autoref:)

[//]: # (        --hyperedit)

[//]: # (            Turn on the "hyperediting" mode, i.e. do not correct&#40;replace&#41; A with G and T with C. This will ensure that)

[//]: # (            potentially hyper-editable sites are not accidentally lost)

[//]: # ()
[//]: # (        --ref-min-cov <ref-min-cov>)

[//]: # (            Automatically correct reference sequence for loci with coverage ≥ the threshold. In short, there is no)

[//]: # (            reason to use the assembly nucleotide "T " if we have sequenced 100% "A ". This heuristic is especially)

[//]: # (            useful in regions of low complexity&#40;or simple repeats&#41;, where such SNPs can affect the editing)

[//]: # (            estimation.[default: 20])

[//]: # (```)

[//]: # ()
[//]: # (#### ROI mode specific)

[//]: # ()
[//]: # (```)

[//]: # (Stats:)

[//]: # (        --ei <ei>)

[//]: # (            File for saving Editing Indexes &#40;EI&#41;. If the file already exists, EI for the current experiments will be)

[//]: # (            appended to it)

[//]: # ()
[//]: # ()
[//]: # (Special information:)

[//]: # (        --rois <rois>)

[//]: # (            Path to a BED file with regions of interest&#40;ROIS&#41; with at least 4 first BED columns&#40;chr, start, end, name&#41;)

[//]: # ()
[//]: # ()
[//]: # (Output filtering:)

[//]: # (        --out-min-cov <out-min-cov>)

[//]: # (            Output only ROIs covered by at least X unique reads&#40;after reads/bases filtering&#41;[default: 10])

[//]: # ()
[//]: # (        --out-min-freq <out-min-freq>)

[//]: # (            Output only ROI having total mismatches frequency ≥ threshold &#40;freq = ∑ mismatches / coverage&#41;[default:)

[//]: # (            0.01])

[//]: # ()
[//]: # (        --out-min-mismatches <out-min-mismatches>)

[//]: # (            Output only ROI having total number of mismatches ≥ threshold. Mismatches are counted jointly, i.e. for the)

[//]: # (            "A" reference we have "C" + "G" + "T". For "N" reference all nucleotides are considered as mismatches. This)

[//]: # (            is a deliberate choice to allow a subsequent user to work through / filter such records[default: 5])

[//]: # (```)

[//]: # ()
[//]: # (#### Loci mode specific)

[//]: # ()
[//]: # (```)

[//]: # (Output filtering:)

[//]: # (        --out-min-cov <out-min-cov>)

[//]: # (            Output only loci covered by at least X unique reads&#40;after reads/bases filtering&#41;[default: 10])

[//]: # ()
[//]: # (        --out-min-freq <out-min-freq>)

[//]: # (            Output only loci having total mismatches frequency ≥ threshold &#40;freq = ∑ mismatches / coverage&#41;[default:)

[//]: # (            0.01])

[//]: # ()
[//]: # (        --out-min-mismatches <out-min-mismatches>)

[//]: # (            Output only loci having total number of mismatches ≥ threshold. Mismatches are counted jointly, i.e. for the)

[//]: # (            "A" reference we have "C" + "G" + "T". For "N" reference all nucleotides are considered as mismatches. This)

[//]: # (            is a deliberate choice to allow a subsequent user to work through / filter such records.[default: 3])

[//]: # (```)

## TODO:

### Must:

1. Check if the integration tests are correct.
2. Test the include/exclude functionality.
3. Command to create a BED track with edits based on the output for each site.
4. More unit-tests.
5. Provide CLI args description in the README

### Maybe:

1. Auto-inference for an exclude list based on the provided fasta (homopolymers, special simple repeats). A module to do
   this on the fly for autoref corrected positions?
2. Handle heterozygous sites and known SNPs (excluding RNA based, i.e. cDNA)
3. Robust error handling (remove unwrap/expect)
4. Lower autoref/coverage thresholds
5. Remove PCR duplicates in PE libraries by default?
6. Count each mate as one half in PE experiments?
7. Better stranding prediction algorithm (include downstream regions)

### Want:

1. Joint analysis of several samples
2. Advanced EI for a set of sites / particular repeats (currently can be done using Python/R)
3. [Base Quality Score Recalibration](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)
4. [Indels Realignment](https://qcb.ucla.edu/wp-content/uploads/sites/14/2016/03/GATKwr12-3-IndelRealignment.pdf)
5. Statistical inference for replicated editing sites / group of sites
