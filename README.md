<p align="center"><img height="200" src="http://jongkyu.kim/images/vaquita_420_340.png"></p>

[Vaquita](http://www.worldwildlife.org/species/vaquita) accurately identifies __structural variations__ using split-reads, discordant read-pairs, soft-clipped reads, and read-depth information. Vaquita does not depend on external tools and very fast. You can analyze __50x WGS sample within an hour__.

The current version is developed for __short-reads__ datasets. For long and noisy reads, it's performance is not rigorously tested. The next version will support long-reads officially.

Download & Complie
-----------------
    git clone https://github.com/seqan/vaquita.git
    mkdir vaquita-build && cd vaquita-build
    cmake ../vaquita && make vaquita -j 4
    
Vaquita supports GCC ≥ 4.9 and Clang ≥ 3.6.

Usage
-----------------
    vaquita call -r [reference.fa] [input.bam] > [output.vcf]

* the `.bam` file must be sorted by coordinates. eg. `samtools sort`
You can find more options using `vaquita --help.`

Contact
-----------------
Jongkyu Kim (xenigmax@gmail.com)
