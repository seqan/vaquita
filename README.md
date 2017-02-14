<p align="center"><img height="200" src="http://jongkyu.kim/images/vaquita_420_340.png" align="middle"></p>
Vaquita identifies __structural variations__ using split-reads, abnormal read-pairs, soft-clipped reads, and read-depth information.

Installation
-----------------
    git clone --recursive https://github.com/xenigmax/vaquita.git
    mkdir vaquita-build && cd vaquita-build
    cmake ../vaquita && make vaquita
 

Usage
-----------------
    vaquita -cg [reference.fa] [input.bam] > [output.vcf]

* the `.bam` file must be sorted by coordinates. eg. `samtools sort`

You can find more options using `vaquita --help.`
