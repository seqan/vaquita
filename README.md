<p align="center"><img height="200" src="http://jongkyu.kim/images/vaquita_420_340.png"></p>

[Vaquita](http://www.worldwildlife.org/species/vaquita) accurately identifies __structural variations__ using split-reads, abnormal read-pairs, soft-clipped reads, and read-depth information. Vaquita does not depend on external tools and very fast. You can analyze 50x WGS sample within 50 minutes. 

The current version only supports __short-reads__. For long-reads, you can try Vaquita with `-v 1  --no-pe --no-ce --no-re`. However, it's performance is not rigorously tested.

Download & Complie
-----------------
    git clone --recursive https://github.com/xenigmax/vaquita.git
    mkdir vaquita-build && cd vaquita-build
    cmake ../vaquita && make vaquita
 

Usage
-----------------
    vaquita -cg [reference.fa] [input.bam] > [output.vcf]

* the `.bam` file must be sorted by coordinates. eg. `samtools sort`

You can find more options using `vaquita --help.`

