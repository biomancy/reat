[![codecov](https://codecov.io/gh/alnfedorov/rada/branch/main/graph/badge.svg?token=9TF8UQD0PD)](https://codecov.io/gh/alnfedorov/rada)
![license](https://img.shields.io/github/license/alnfedorov/rada)

### RADA

**RADA** is an RNA editing analysis package designed with a focus on performance, low memory footprint, and ease of use.

Refer to the [Wikipedia](https://en.wikipedia.org/wiki/RNA_editing) for classification and basic overview of known
editing events. **RADA** can handle all types of edits (eg, A->I, C->U, etc.), except for editing by insertions or
deletions.

The target platform for **RADA** is an x86-64 computer running Linux; working on Windows or ARM is not guaranteed.

Please, feel free to open issues regarding bugs, installation issues, feature requests, or unintuitive behavior.
Quality, correctness, and ease of use are vital to us.

### Not released yet

**RADA** is under active development and is not yet been officially released. While having unit tests whenever possible,
the tool still lacks integration tests and users feedback. We will release in a few months; until then, use at your own
risk.

### Features

* Summarizing editing for provided regions or separate loci
* Efficient multithreading
* Strand prediction for unstranded libraries
* Autoref: automatic detection and accounting for single nucleotide variants (SNV)
* Editing Index (EI) for the given set of ROI

See [details](#details) section for more in-depth explanation of some features.

### Limitations

Here is a list of known limitations; feel free to open an issue if one of them is critical for you:

* Potential RNA editing by insertions/deletions is ignored
* No Python or R interface
* Lack of support for ARM and Windows/macOS builds

### Installation

#### Pre-built binaries

For a quick start, you can try the pre-built **RADA** binary posted
via [GitHub releases](https://github.com/alnfedorov/rada/releases). Besides the tagged versions, we also provide the
latest tested build from the [main branch](https://github.com/alnfedorov/rada/releases/tag/latest).

Note that you will have to add executable permissions to the file: `sudo chmod +x rada`.

#### Cargo

Rust package manager _cargo_ is the recommended way to build **RADA** from sources:

```shell
cargo install --git https://github.com/alnfedorov/rada
rada --help
```

Follow [this](https://www.rust-lang.org/tools/install) page to install rust toolchain if you don't have one.

### Basic usage

**RADA** supports two modes: ROI-based and loci-based.

#### ROI mode

ROI stands for Region of Interest, a **RADA** mode in which edits are summarized for each provided genomic region. For
instance, It is useful in the context of studying Alu repeats editing.

Example:

* _BED-like file with ROIs_(repeats.bed):

```text
chrY	76783396	76783529	RSINE   .    -   ... # The rest of the columns are ignored
... # Remaining regions omitted
```

* _Command:_

```shell
rada --input techrep1.bam techrep2.bam --roi repeats.bed \
     --reference GRCm38.fa --stranding "f/s" --threads $(nproc) \
     --saveto /dev/stdout 
```

Output columns:

* **chr, start, end, name** - ROI coordinate and name from the input file
* **strand** - transcription strand; predicted for unstranded libraries and deducted from the design for stranded
  experiments
* **X->Y** - Total number of events observed in a given ROI where a reference nucleotide _X_ was replaced by _Y_. That
  is, A->A is a number of A matches, and A->G denotes the total number of observed A->I edits.

Note that **X->Y** notation always denotes matches/mismatches relative to the forward strand. For example, T->C
mismatches for REVERSE strand ROI are, in fact, A->G edits within ROI in a pool of sequenced RNAs.

If autoref feature is enabled, edits are summarised after correcting for any potential SNVs.

#### ROI editing index

For a given set of regions, one can always calculate an editing index (EI) for all possible matches and mismatches.

For example, formally EI for mismatches A->G is defined as follows: EI(A->G) = P(transcribed ROI is A->G edited) = (
forward(A->G) + REVERSE(T->C)) / (forward(A->A + A->C + A->G + A->T) + REVERSE(T->A + T->C + T->G + T->T)). In other
words, it is the probability that a given A in the provided ROI will be A->G edited after the transcription.

Note that the EI is strand independent as it is defined for the transcribed regions of RNA. For example, to edit A-> I,
you need to refer exclusively to the A->G column and completely ignore T->C.

If a given set of ROIs represent all Alu repeats in the genome, then A->G EI is the so-called Alu Editing Index (AEI).

Example:

```shell
for bamfile in *.bam
do
  rada --input bamfile  --roi alu-repeats.bed \
       --reference hg19.fasta --threads $(nproc) \
       --stranding "s/f" --annotation gencode.gff3.gz \
       --ei=AEI.tsv --saveto /dev/null # append Alu editing index values and skip ROI results completely
done
```

Output columns:

* **files** - comma-separated input BAM files
* **X->Y** - editing index for _X->Y_ pair

One can call **RADA** multiple times with the same TSV file to append rows to the EI table.

#### Loci mode

The **RADA** loci-based mode is a classic scenario for estimating RNA editing for each genomic locus.

Example:

```shell
rada --input rnaseq.bam  --binsize 50000 \
     --reference hg38.fa --threads $(nproc) \
     --stranding "u" --annotation gencode.gff3.gz \
     --saveto /dev/stdout 
```

Output columns:

* **chr,position** - coordinate of the locus
* **strand** - transcription strand; predicted for unstranded libraries and deducted from the design for stranded
  experiments
* **X** - the total number of sequenced nucleotides X; X is one of \[A, C, G, T\].
* **reference** - predicted reference nucleotide; just a reference assembly if autoref feature is disabled

Similarly to the ROI mode, the reference and sequenced nucleotides **X** are always reported with respect to the forward
strand. That is, a minus strand locus with ten A's corresponds to ten sequenced T's from RNA fragments originated from
the minus strand.

### CLI arguments

Here is a list of arguments(`rada --help`) supported by the Command Line Interface (CLI):

```
Core:
    -i, --input <input>...
            Path to the input BAM file(s). May contain a space-separated list of files, in which case they are
            treated as technical replicates and pulled together.

    -r, --reference <reference>
            Indexed fasta file with a reference genome assembly. Contig/chromosome names must match the names in
            the input BAM headers.

        --roi <roi>
            Path to a BED-like file with 4 columns (chr, start, end, name) in which the target regions of
            interest (ROI) are declared. If specified, total mismatches will be reported for each region of
            interest rather than loci.

        --binsize <binsize>
            Summarize the mismatches for each locus, providing each worker thread with genome bins(job share) of
            approximately the specified size (in base pairs).

    -s, --stranding <stranding>
            Strand-specificity of the experiment, i.e. matching between the read strand and the transcript
            strand. Use "u" for unstranded experiments; other available options based on the RSeQC
            nomenclature(see infer_experiment.py docs): same:"s" (++,--), flip:"f" (+-,-+), same read1/flip
            read2:"sf" (1++,1--/2+-,2-+), flip read1/same read2:"fs" (1+-,1-+/2++,2--).[possible values: u, s,
            f, s/f, f/s]

    -o, --saveto <saveto>
            Path to the output tsv file. Use /dev/stdout to print result to the standard output.

    -t, --threads <threads>
            Maximum number of threads to spawn at once.[default: 1]


Filtering:
        --mapq <mapq>
            Count only reads with mapq ≥ threshold. Note that reads with mapq = 255 are always skipped(mapq is
            not available according to the SAM spec).[default: 1]

        --phread <phread>
            Count only bases with phread ≥ threshold. For reference, Phread is defined as -10 log₁₀[error
            probability], so Phread = 20 means 1 error in 100 base calls.[default: 20]

        --out-min-mismatches <out-min-mismatches>
            Output only loci/ROI having total number of mismatches ≥ threshold. Mismatches are counted jointly,
            i.e. for the "A" reference we have "C" + "G" + "T". For "N" reference all nucleotides are considered
            as mismatches. This is a deliberate choice to allow a subsequent user to work through / filter such
            records.[default: 5]

        --out-min-freq <out-min-freq>
            Output only loci/ROI having total mismatches frequency ≥ threshold (freq = ∑ mismatches /
            coverage)[default: 0.01]


Autoref:
        --ref-min-cov <ref-min-cov>
            Automatically correct reference sequence for loci with coverage ≥ the threshold. In short, there is
            no reason to use the assembly nucleotide "T " if we have sequenced 100% "A ". This heuristic is
            especially useful in regions of low complexity(or simple repeats), where such SNPs can affect the
            editing estimation.[default: 20]

        --ref-min-freq <ref-min-freq>
            Automatically correct reference sequence for loci with the most common nucleotide frequency ≥
            cutoff[default: 0.95]

        --hyperedit
            Turn on the "hyperediting" mode, i.e. do not correct(replace) A with G and T with C. This will
            ensure that potentially hyper-editable sites are not accidentally lost.


Stranding:
        --annotation <annotation>
            Genome annotation in the GFF3 format. Genomic features (exons and genes) are used to inference
            loci/ROI strand based on the most likely direction of transcription (see the GitHub documentation
            for details). It is recommended to provide genome annotation for unstranded libraries, otherwise
            stranding will be highly inaccurate.

        --str-min-mismatches <str-min-mismatches>
            Automatically predict strand based on the observed A->I editing for locus/ROI with A->G mismatches
            >= threshold. It is a fallback strand prediction heuristic, used only for the unstranded libraries.
            Not relevant for organisms without active ADAR-like enzymes.[default: 50]

        --str-min-freq <str-min-freq>
            Automatically predict strand based on the observed A->I editing for locus/ROI with A->G freq >=
            threshold (freq = ∑ A->G / (∑ A->G + ∑ A->A))[default: 0.05]


Statistics:
        --ei <ei>
            File for saving Editing Indexes (EI). If the file already exists, the EI results for the current
            experiments will be appended to it.
```

### Details

#### Strand prediction

To predict transcription strand for ROI/loci in unstranded experiments, **RADA** uses two strategies.

First, ROI/loci strand will be derived from the overlapping genes/exons strand. This approach only works if the GFF3
genome annotation is provided and can be summarized in the following table:

| overlapping genes | overlapping exons | predicted strand |
|:-----------------:|:-----------------:|:----------------:|
|         +         |         +         |         +        |
|         -         |         -         |         -        |
|        +/-        |         +         |         +        |
|        +/-        |         -         |         -        |
|        +/-        |        +/-        |         ?        |

That is, **RADA** checks overlapping genes first. If they are genes on the + and the - strand, exons are considered. In
the worst-case scenario, an unknown(`.`) strand is returned.

Second, for ROIs / loci for which **RADA** could not predict the strand from the annotation, **RADA** attempts to derive
the strand based on the observed A->I editing.

For the + strand transcripts, A-> I edits are A-> G mismatches, and for the - strand, T-> C mismatches. Note that in
many cases, this heuristic fails (no A->I editing at all), and such ROIs / loci will be left unstranded in the final
table.

#### Autoref

With sufficient coverage, we can automatically adjust the reference sequence for observed SNVs based on RNA-seq data.
**RADA** uses a simple heuristic for this: if coverage is sufficient and the frequency of the most abundant nucleotide
exceeds the threshold, then the reference nucleotide is the most abundant nucleotide at that locus.

Note that hyper-editing flag allows one to skip A->G and T->C corrections to explore potential hyperedited ROI/loci.

# TODO:

- Check grammar
- Integration tests
- Editing indices (alu editing index, Non-synonimous editing index)
- Graph: (RAM usage per thread) vs (#threads)
- Graph: (Speedup against REDItools) vs (#threads) + note about run time
- Handle errors properly!

- Optimization: group nearby regions together and fetch reads for all the regions inside the batch
- Сoverage details for each region(cov per nucleotide, total overlapped fragments)
