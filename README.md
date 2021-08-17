[![codecov](https://codecov.io/gh/alnfedorov/rada/branch/main/graph/badge.svg?token=9TF8UQD0PD)](https://codecov.io/gh/alnfedorov/rada)
![license](https://img.shields.io/github/license/alnfedorov/rada)

### RADA

**RADA** is an RNA editing analysis package, designed with emphasis on performance, low memory footprint, and ease of
use.

Refer to the wikipedia [article](https://en.wikipedia.org/wiki/RNA_editing) for a classification and basic overview of
known editing events. **RADA** can handle all types of editing (i.e. A->I, C->U, etc), except for editing by insertions
or deletions.

**RADA** officially supports only Linux distributions running on an x86-64 machine. Working on ARM or Windows systems is
not guaranteed.

Please, feel free to open issues concerning bugs, installation problems, feature requests or unintuitive behaviour.
Quality, correctness, and ease of use matter a lot for us.

### Not release yet

**RADA** is under an active development and is not yet formally released. While having unit-tests whenever possible, the
tool still lacks integration tests and users feedback. We will release in a few months, until then use at your own risk.

### Features

See details section for more in-depth explanation of some options.

### Limitations

This a list of known limitations, feel free to open an issue if one of them is critical for you:

* RNA editing by insertions and deletions is not supported
* No Python or R interface
* Lack of official support for ARM and Windows builds

### Installation

Usage of the rust _cargo_ package manager is the recommended way of installing **RADA**:

```shell
cargo install --git https://github.com/alnfedorov/rada
rada --help
```

Follow [this](https://www.rust-lang.org/tools/install) page to install rust toolchain if you don't have one.

TODO: Docker container and prebuilt binaries.

### Basic usage

**RADA** supports running in two modes: ROI based and loci based.

#### ROI mode

ROI stands for region of interest, a **RADA** mode when editing is summarised for each provided region. It is useful in
a context of studying, for example, Alu repeats editing. Example:

```shell

```

Output example:

#### Loci mode

Loci-based **RADA** mode is a typical scenario of quantifying editing per-locus for the whole genome. Example:

```shell

```

Output example

### CLI arguments

Running `rada --help` will give you the following list of supported arguments:

```
Core:
    -i, --input <input>...
            Path to the input BAM file(s). May contain a comma-separated list of files, in which case they are treated
            as technical replicates and pulled together.

    -r, --reference <reference>
            Indexed fasta file with a reference genome assembly. Contig/chromosome names must match the names in the
            input BAM headers.

        --roi <roi>
            Path to a BED-like file with 4 columns (chr, start, end, name) in which the target regions of interest (ROI)
            are declared. If specified, cumulative mismatches will be reported for each region of interest rather than
            loci.

        --binsize <binsize>
            Summarize the mismatches for each locus, providing each worker thread with genome bins(job share) of
            approximately the specified size (in base pairs).

    -s, --stranding <stranding>
            Strand-specificity of the experiment, i.e. matching the read strand and the transcript strand. Use "u" for
            unstranded experiments; other available options based on the RSeQC nomenclature(see infer_experiment.py
            docs): "s" (++,--), "f" (+-,-+), "sf" (1++,1--,2+-,2-+), "fs" (1+-,1-+,2++,2--).[possible values: u, s, f,
            sf, fs]

    -o, --saveto <saveto>
            Path to the output tsv file. May end with ".gz", in which case the stream is automatically gzipped.

    -t, --threads <threads>
            Maximum number of threads to spawn at once.[default: 1]


Filtering:
        --mapq <mapq>
            Count only reads with mapq ≥ threshold. Note that reads with mapq = 255 are always skipped(mapq is not
            available according to the SAM spec).[default: 1]

        --phread <phread>
            Count only bases with phread ≥ threshold. For reference, Phread is defined as -10 log₁₀[error probability],
            so Phread = 20 means 1 error in 100 base calls.[default: 20]

        --out-min-mismatches <out-min-mismatches>
            Output only loci/ROI having total number of mismatches ≥ threshold. Mismatches are counted jointly, i.e. for
            the "A" reference we have "C" + "G" + "T". For "N" reference all nucleotides are considered as mismatches.
            This is a deliberate choice to allow a subsequent user to work through / filter such records.[default: 5]

        --out-min-freq <out-min-freq>
            Output only loci/ROI having total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)[default:
            0.01]


Autoref:
        --ref-min-cov <ref-min-cov>
            Automatically correct reference sequence for loci with coverage ≥ the threshold. In short, there is no
            reason to use the assembly nucleotide "T " if we have sequenced 100% "A ". This heuristic is especially
            useful in regions of low complexity(or simple repeats), where such SNPs can affect the editing
            estimation.[default: 20]

        --ref-min-freq <ref-min-freq>
            Automatically correct reference sequence for loci with the most common nucleotide frequency ≥ cutoff
            [default: 0.95]


Stranding:
        --annotation <annotation>
            Genome annotation in the GFF3 format. Genomic features (exons and genes) are used to inference loci/ROI
            strand based on the most likely direction of transcription (see the GitHub documentation for details). It is
            highly recommended to provide genome annotation for unstranded libraries, otherwise stranding will be highly
            inaccurate.

        --str-min-mismatches <str-min-mismatches>
            Automatically predict strand based on the observed A->I editing for locus/ROI with A->G mismatches >=
            threshold. It is a fallback strand prediction heuristic, used only for the unstranded libraries. Not
            relevant for organisms without active ADAR-like enzymes.[default: 50]

        --str-min-freq <str-min-freq>
            Automatically predict strand based on the observed A->I editing for locus/ROI with A->G frequency >=
            threshold (freq = ∑ A->G / (∑ A->G + ∑ A->A)) [default: 0.05]
```

### Details

#### Stranding

TODO

#### Autoref

TODO

# TODO:

- Check grammar
- Finish Docs
- Hyperediting mode (keep G's in the refnucpred module)
- Editing indices (alu editing index, Non-synonimous editing index)
- Graph: (RAM usage per thread) vs (#threads)
- Graph: (Speedup against REDItools) vs (#threads) + note about running time
