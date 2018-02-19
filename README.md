[![Build Status](https://travis-ci.org/seqan/vaquita.svg?branch=master)](https://travis-ci.org/seqan/vaquita)
<p align="center"><img height="200" src="http://jongkyu.kim/images/vaquita_420_340.png"></p>

[Vaquita](http://www.worldwildlife.org/species/vaquita) accurately identifies __structural variations__ using split-reads, discordant read-pairs, soft-clipped reads, and read-depth information. Vaquita does not depend on external tools and very fast. You can analyze __50x WGS sample within an hour__.


Download & Compile
-----------------
    git clone https://github.com/seqan/vaquita.git
    mkdir vaquita-build && cd vaquita-build
    cmake ../vaquita && make vaquita -j 4

Vaquita supports `GCC≥4.9` and `Clang≥3.8`.

Usage
-----------------
    vaquita call -r [reference.fa] [input.bam] > [output.vcf]

The `.bam` file must be sorted by coordinate (eg. `samtools sort`).
You can find more options using `vaquita call --help`.

Citation
-----------------
Jongkyu Kim and Knut Reinert, Vaquita: Fast and Accurate Identification of Structural Variations using Combined Evidence. _Workshop on Algorithmic Bioinformatics (WABI) 2017_

DOI: [10.4230/LIPIcs.WABI.2017.13](http://drops.dagstuhl.de/opus/volltexte/2017/7635/)

* You can find all the scripts and information about raw datasets that I used for benchmarking at [this repository](https://github.com/xenigmax/vaquita_WABI2017).

Contact
-----------------
Jongkyu Kim (vaquita@jongkyu.kim)
