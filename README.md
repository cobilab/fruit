# FRUIT &nbsp;&nbsp;&nbsp;&nbsp; [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
Fast Relative Uniqueness fInder for proTein sequences

## RUN
FRUIT includes three tools: `fruit-map`, `fruit-filter` and `fruit-visual`.

### FRUIT-MAP
The `fruit-map` finds uniqueness in a group of target proteins, relatively to a group of reference proteins. The output of this tool is a group of files (*.sng) including zeros and ones associated with positions of the characters in the target sequences. Note that each output file and its corresponding target file are of equal size. This output can then be fed to `fruit-filter` tool.

```bash
./fruit-map  [OPTIONS]...  -r [REF1:REF2:...]  -t [TAR1:TAR2:...]
```
The followings are two examples:
```bash
./fruit-map -k 7 -r ref -t tar
```
or
```bash
./fruit-map -k 10 -r ref1:ref2 -t tar1:tar2:tar3
```
It is recommended to choose short names for reference and target sequences.

To see the possible options, type:
```bash
./fruit-map
```
which provides the following:
```text
DESCRIPTION
  The output: *.sng

  Mandatory arguments:
  -r   FILES           reference file(s) (SEQ/FASTA/FASTQ)
  -t   FILES           target file(s)    (SEQ/FASTA)

  Options:
  -k   INT             k-mer size
  -bs  INT             size of bloom filter, in log2 format
                       e.g. 10 corresponds to 2^10
  -v                   more information
  -h                   usage guide
```
Reference proteins can have SEQ (M, A, R, D, etc), FASTA or FASTQ formats, and target proteins can be of SEQ or FASTA formats. Note, the maximum allowed k-mer size depends on the cardinality (number of different characters) of the input files. Also note that bloom filter is the structure used to save the relatively unique regions, and the size of this structure can be customized by the user; although, the optimal values are automatically selected, by default.

### FRUIT-FILTER
The `fruit-filter` takes \*.sng files, generated by `fruit-map` tool, and filters runs of zeros, which represent relatively unique regions. The output of this tool is a group of position files (*.pos) including the beginning and the end positions of relatively unique regions. This output can be later fed to `fruit-visual` tool.

```bash
./fruit-filter  [OPTIONS]...  -i [FILE1:FILE2:...]
```

The following shows an example:
```bash
./fruit-filter -w 2 -i file1:file2
```

To see the possible options, type:
```bash
./fruit-filter
```
which provides the following:
```text
DESCRIPTION
  The output: *.pos

  Mandatory arguments:
  -i   FILES           map file(s), generated by
                       fruit-map tool (*.sng)
  Options:
  -w   INT             window size
  -v                   more information
  -h                   usage guide
```
`-i` stands for input(s). The window size (`-w`) determines the minumum size of any unique region. It is 10, by default.

### FRUIT-VISUAL
The `fruit-visual` takes \*.pos files, generated by `fruit-filter` tool, and visualizes the relative singular regions. The output of this tool is an "svg" file.

```bash
./fruit-visual  [OPTIONS]...  -i [POS_FILE1:POS_FILE2:...]  -o SVG_FILE
```

The following shows an example:
```bash
./fruit-visual -i ab.pos:cd.pos:ef.pos -o out.svg
```

To see the possible options, type:
```bash
./fruit-visual
```
which provides the following:
```text
DESCRIPTION
  Mandatory arguments:
  -i  POS-FILES        position file(s), generated by
                       fruit-filter tool (*.pos)
  Options:
  -o  SVG-FILE         output image name (*.svg). Default: map.svg
  -sl STRINGS          sequence labels. If any space in the label,
                       surround it by "s, e.g. "prot seq label".
                       Default: names of position files
  -t  INT              tick length: [1, 4294967295]
  -am INT              axis max: [1, 4294967295]
  -th 0|1              tick label human readable: 0=false, 1=true.
                       if true: e.g. 1200 -> 1.2 K. Default: true
  -v                   more information
  -h                   usage guide
```
`-i` stands for input(s). `-o` stands for output, and if is not determined by user, it will be named "map.svg".

## Experiments
To replicate the results presented in the paper, run the bash script *run.sh*, provided in the "script/" directory. Note that to generate the synthetic datasets, we have used the "mutate" tool, available in "prog/" directory, in order to mutate an input sequence with a specific rate. For example, the following command takes 50 as the seed, mutates the "input_seq" with the rate of 10%, and writes the result into "mutated_seq":
```bash
./mutate -s 50 -r 10 -i input_seq -o mutated_seq
```
`-s`, `-r`, `-i` and `-o` stand for seed, rate, input and output, respectively. Consider the input_seq as ACG**T**ACGTAC, the mutated_seq might be ACG**C**ACGTAC, that is 10% of the symbols (1 out of 10) has been mutated.

<!-- ## Example -->

<!-- ### Compare fruit with other methods
In order for comparison, you might set the parameters in 
"run.sh" bash script, then run it:
```bash
./run.sh
```
With this script, you can download the datasets, install the dependencies, 
install the other tools, run all the tools, and finally, visualize the results. -->

## CITE
Please cite the following, if you use **FRUIT**:
* M. Hosseini, D. Pratas and A.J. Pinho, ``A Probabilistic Method to Find and Visualize Distinct Regions in Protein Sequences,'' *27th European Signal Processing Conference* (*EUSIPCO*), pp. 1-5, IEEE, 2019.

<!-- ## RELEASES
* [Release](https://github.com/smortezah/fruit/releases) 1: . -->

## ISSUES
Please let us know if there is any 
[issues](https://github.com/smortezah/fruit/issues).
