Cpg_Islands
===========

Search for cpg islands in genome using HMM

Installation
===========

make

Usage
=====

Prediction mode:
    HMM_CpG [-p <parameters_file>] -cn <chromosome name> < <data> > <output_file>

    -p(--parameters) - allows to specify model parameters file. You can use train mode to obtain parameters file.
                        By default 'parameters' file will be used.
    -cn(--chromosome-name) - allows to specify chromosome name. It will be used as prefix to each cpg (required for bed format).

Training mode:
    HMM_CpG [--cpg <file_with_cpg_examples>] [--non-cpg <file_with_non_cpg_examples>] [-acl <positive integer>] [-ancl <positive integer>]
    Note that you should specify at least one of --cpg or --non-cpg keys to activate training mode.

    --cpg - allows to specify file with training cpg examples. Examples must be separated with 'N'.
    --non-cpg - allows to specify file with training non-cpg examples. Examples must be separated with 'N'.
    -acl(--average-cpg-length) - allows to specify average cpg length. If it doesn't specified, average training cpg length will be taken.
        If training cpg file doesn't specified, default value (765) will be taken.
    -ancl(--average-non-cpg-length) - allows to specify average non-cpg length. If it doesn't specified, average training non-cpg length will be taken.
    If training nno-cpg file doesn't specified, default value (100415) will be taken.

Examples:
./HMM_CpG --cpg training_cpg_set -acl 900 > my_parameters
./HMM_CpG --non-cpg training_non_cpg_set > my_parameters
./HMM_CpG --cpg training_cpg_set --non-cpg training_non_cpg_set -acl 900 -ancl 100500 > my_parameters

cat genome.fa | ./HMM_CpG
./HMM_CpG -p my_parameters < chr1.fa > out.bed
