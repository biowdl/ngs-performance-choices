ngs-performance-choices
=======================
Which tools to choose based on various performance metrics.

Performance testing
+++++++++++++++++++

For fair comparison between tools it is best to use containerized environments.
These eliminate any user specific installs or compilations. 

For running the tools `Singularity <https://github.com/hpcng/singularity>`_ was
used. The containers were provided by `biocontainers
<https://quay.io/biocontainers>`_.

Tools were run with ``singularity exec -eip
docker://quay.io/biocontainers/<tool>:<version-tag> command``. The ``-e``,
``-i`` and ``-p`` flags isolate environment, IPC namespace and PID namespace
respectively. By default singularity mounts ``$PWD`` and user ``HOME``
directories, which was convenient when testing.

For benchmarking [hyperfine](https://github.com/sharkdp/hyperfine) was used.
This command line tool is useful as it can perform warmup runs to ensure 
images are downloaded and disk caches are full.

Multithreading
++++++++++++++

Multithreading is a programming technique that allows a program to split up
its execution in multiple independent threads that allow for parallel execution
on the CPU.

Advantages:

+ More efficient use of resources (less CPU's in idle state)
+ Significant decrease in wall clock time if used correctly.

Disadvantages:

+ Performance overhead compared to single-threaded algorithms
+ More complex code required

When threads share resources some conditions
must be satisfied to ensure the threads don't alter the shared resource at
the same time (race condition) or wait on each other to change the resource
and therefore never changing it (deadlock). Taking care these problems is
what makes the code more complex and comes with a performance penalty.

Some algorithms come with a very small performance overhead with regards to
multithreading some with a much bigger one.

How much threads should be used?
--------------------------------
Given the overhead of using more than one thread, the default should be one
thread. Situations in which more threads can be considered are:

+ A reduction in wall clock time is required or desirable
+ A program is bottlenecking other applications
+ Using a program in single-threaded mode will leave resources idle.

Examples:

+ When time is no object, it may be cheaper to use one thread per task on cloud
  backends. This will decrease the amount of CPU time needed.
+ On clinical pipelines, getting the diagnose one day earlier may be worth more
  than computational overhead.
+ The program is used in a pipe and other programs spend a lot of time waiting
  on it.
+ mapping/alignment tools require a lot of memory for the index, but
  not much additional memory per thread. On a node that has 16 GB and 4 threads
  it is more efficient to use all 4 threads and 16 GB of ram, than to use 12
  GB of ram on one thread.

Compression and compression levels
++++++++++++++++++++++++++++++++++

Compression tools
-----------------
The bioinformatics world seems to have standardised on `zlib
<https://www.zlib.net/>`_. Most bioinformatics tools can open ``.gz`` files and
output to them. Therefore tools that implement zlib should be used for
compressing files.

Not all zlib implementations are created equal.
The most well-known zlib implementation is `gzip <http://www.gzip.org/>`_. Gzip
is very old and is written as a single-threaded algorithm. Another
implementation is `pigz <http://zlib.net/pigz/>`_, which was written to make use
of multiple threads. Pigz is **always** faster than gzip. Even on a single core
the algorithm is simply more efficient. ``pigz -p 1`` (no multithreading) beats
``gzip`` in both compression and decompression. [Cutadapt](
https://github.com/marcelm/cutadapt) therefore does its compression and
decompression with pigz.

Why pigz:

.. code-block::

    $ hyperfine -w 3 -r 10 "gzip -1 -c big.fastq > /dev/null"
    Benchmark #1: gzip -1 -c big.fastq > /dev/null
      Time (mean ± σ):     22.335 s ±  0.284 s    [User: 22.123 s, System: 0.202 s]
      Range (min … max):   21.963 s … 22.890 s    10 runs

    $ hyperfine -w 3 -r 10 "pigz -p 1 -1 -c big.fastq > /dev/null"
    Benchmark #1: pigz -p 1 -1 -c big.fastq > /dev/null
      Time (mean ± σ):     19.343 s ±  0.225 s    [User: 19.085 s, System: 0.243 s]
      Range (min … max):   19.093 s … 19.737 s    10 runs

    $ hyperfine -w 3 -r 10 "gzip -cd big.fastq.gz > /dev/null"
    Benchmark #1: gzip -cd big.fastq.gz > /dev/null
      Time (mean ± σ):     10.363 s ±  0.086 s    [User: 10.274 s, System: 0.059 s]
      Range (min … max):   10.217 s … 10.471 s    10 runs

    $ hyperfine -w 3 -r 10 "pigz -p 1  -cd big.fastq.gz > /dev/null"
    Benchmark #1: pigz -p 1  -cd big.fastq.gz > /dev/null
      Time (mean ± σ):      6.196 s ±  0.037 s    [User: 6.134 s, System: 0.053 s]
      Range (min … max):    6.132 s …  6.250 s    10 runs

Pigz is more than 10% faster when compressing and more than 40 % (!!!) faster
when decompressing.

which compression level to use?
-------------------------------
There are several zlib compression levels, from 1 to 9. Alternatively a file
can not be compressed. Should a file be compressed, and if so, at which level?

Compression was done using ``pigz -p 1``

============ ================ ========= ============= =============
level        time (seconds)   size      relative time relative size
============ ================ ========= ============= =============
uncompressed 0.914            2,3G      0.05          4.44
1            19.770           513M      1.00          1.00
2            21.713           496M      1.10          0.97
3            30.370           467M      1.54          0.91
4            28.196           465M      1.43          0.91
5            47.708           446M      2.41          0.87
6            108.467          402M      5.49          0.78
7            174.239          386M      8.81          0.75
8            220.316          383M      11.14         0.75
9            223.419          383M      11.30         0.75
============ ================ ========= ============= =============

The used data was quite repetitive so might not have been a best benchmark for
the highest compression levels. What we see is that anything above compression
level 1 uses disproportionally more compute time for the benefits it gives.
Compression level 4 might be worth it but takes 43% compute time for a 9%
smaller filesize.

The default compression should be level 1. What we see in for example ``Picard
Markduplicates`` is that the execution time is halved when compression level
is set from 5 (default) to 1.

Alignment tools
+++++++++++++++

DNA Alignment
-------------
For short read DNA alignment, `BWA <https://github.com/lh3/bwa>`_ has been the
tool of choice. It is used by many insitutions as the go to standard.

When evaluating the performance of `BIOWDL's germline-DNA pipeline
<https://github.com/biowdl/germline-DNA>`_, it was found that the alignment step
using BWA takes ~50% of the wall-clock time and 70% of the compute time. This
was *after* all the remaing steps in the pipeline had optimized compression
levels. Alignment is by far the most expensive step. The reason that the wall
clock time is less than the compute time is that the task can be parallelized
very well with virtually no loss of efficiency.

To solve the long alignment times a new tool has come along: `bwa-mem2
<https://github.com/bwa-mem2/bwa-mem2>`_. It promises to deliver increased
speeds at the cost of an increase in memory. This increase in memory was quite
substantial (10x) but since bwa-mem2 2.1 the increase has been reduced to 3x.
It has this been more viable for use generally, and can now be tested locally
on a workstation which makes it easier to have apples-to-apples performance
comparisons.

Benchmarks
----------
To compare bwa and bwa-mem2 the reference genome for pipelines: `hg38 including
alt contigs and hs38d1 decoy sequences
<ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz>`_
was used. It is provided by the `NCBI's ftp site
<ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/>`_
. The `BWA index can also be downloaded from that page
<ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz>`_.
The index for bwa-mem2 was generated from the fasta file, with the ``.alt``
index file from the bwa index (it's compatible).

The input FASTQ file was a 5 million record FASTQ file from the Genome in a
bottle sample with read length 150. For a total of 1,6 GB. The file was
uncompressed and stored on a NVME SSD drive. The indexes were stored on the
same drive. The output was written as uncompressed SAM to a in-memory disk of
4GB.

Bwa-mem2 comes with several executables, each optimized for a particular
`instruction set <https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions>`_.
By default it chooses the fastest instruction set allowed by the CPU.
These instruction sets allow fasts operations on vectors. Meaning that larger
amounts of numbers can be manipulated at the same time for much less cost. The
different instruction sets were benchmarked. SSE4.1, SSE4.2, AVX and AVX2 were
tested as these were supported by the used CPU for testing (AMD Ryzen 5 3600).
AVX512 is only available on modern Intel CPUs and was not tested.

Bwa-mem2 is provided as both source and as binaries on its github page.
The binaries are compiled using the Intel C compiler. `This compiler has
been known to deliberately reduce performance on non-Intel CPU's
<https://en.wikipedia.org/wiki/Intel_C%2B%2B_Compiler#Reception>`_. It is also
proprietary. Therefore the `bioconda distribution
<http://bioconda.github.io/recipes/bwa-mem2/README.html>_` uses the GNU
compiler collection. Both Intel and GNU compiled binaries were tested. Also
the source was compiled locally using GCC and AOCC (`AMD Optimizing C Compiler
<https://developer.amd.com/amd-aocc/>`_) to see how much value was added by
local compilation.

Versions used were BWA 0.7.17 and bwa-mem2 2.1

Test results
............

The baseline is BWA. This uses slightly over 6 GB in memory (6.3 GB).

.. code-block::

    $ hyperfine -w 2 -r 5 'bwa mem -t 8 -o ramdisk/bwa-mem.sam GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa mem -t 8 -o ramdisk/bwa-mem.sam GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     179.499 s ±  0.700 s    [User: 1425.645 s, System: 3.179 s]
      Range (min … max):   179.042 s … 180.719 s    5 runs

bwa-mem2 uses almost 19GB (18.9 GB) maximum which is approximately 3 times the
memory used in BWA. Below follow the GNU Compiler compiled binary results.
Note how the instruction sets have a big impact on performance. They are
ordered from older to newer.

.. code-block::

    $ hyperfine -w 2 -r 5 'bwa-mem2.sse41 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2.sse41 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     139.950 s ±  4.528 s    [User: 974.764 s, System: 9.624 s]
      Range (min … max):   136.166 s … 145.029 s    5 runs

    $ hyperfine -w 2 -r 5 'bwa-mem2.sse42 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2.sse42 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     124.358 s ±  0.194 s    [User: 881.606 s, System: 9.801 s]
      Range (min … max):   124.187 s … 124.636 s    5 runs

    $ hyperfine -w 2 -r 5 'bwa-mem2.avx mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2.avx mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     115.176 s ±  0.594 s    [User: 807.904 s, System: 9.872 s]
      Range (min … max):   114.439 s … 116.068 s    5 runs

    $ hyperfine -w 2 -r 5 'bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     96.902 s ±  0.184 s    [User: 662.478 s, System: 9.813 s]
      Range (min … max):   96.756 s … 97.118 s    5 runs

Below are the results using the binaries found on the bwa-mem2 github page.

.. code-block::

    $ hyperfine -w 2 -r 5 'bwa-mem2-2.1_x64-linux/bwa-mem2.sse41 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2-2.1_x64-linux/bwa-mem2.sse41 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     124.173 s ±  0.089 s    [User: 874.309 s, System: 9.747 s]
      Range (min … max):   124.056 s … 124.292 s    5 runs

    $ hyperfine -w 2 -r 5 'bwa-mem2-2.1_x64-linux/bwa-mem2.sse42 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2-2.1_x64-linux/bwa-mem2.sse42 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     121.848 s ±  0.966 s    [User: 856.020 s, System: 9.917 s]
      Range (min … max):   121.142 s … 123.442 s    5 runs

    $ hyperfine -w 2 -r 5 'bwa-mem2-2.1_x64-linux/bwa-mem2.avx mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2-2.1_x64-linux/bwa-mem2.avx mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     111.642 s ±  0.145 s    [User: 774.780 s, System: 9.681 s]
      Range (min … max):   111.502 s … 111.861 s    5 runs

    $ hyperfine -w 2 -r 5 'bwa-mem2-2.1_x64-linux/bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2-2.1_x64-linux/bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     93.252 s ±  0.477 s    [User: 626.916 s, System: 10.021 s]
      Range (min … max):   92.757 s … 93.881 s    5 runs

Compiled on the test machine using GCC:

.. code-block::

    $ hyperfine -w 2 -r 5 'bwa-mem2-source/bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2-source/bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     95.301 s ±  0.446 s    [User: 645.498 s, System: 9.800 s]
      Range (min … max):   94.923 s … 95.942 s    5 runs

Compiled on the test machine using AOCC:

.. code-block::

    $ hyperfine -w 2 -r 5 'bwa-mem2-source/bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq'
    Benchmark #1: bwa-mem2-source/bwa-mem2.avx2 mem -t 8 -o ramdisk/bwa-mem2.sam bwa-mem2/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna ../big2.fastq
      Time (mean ± σ):     91.842 s ±  0.105 s    [User: 618.004 s, System: 9.729 s]
      Range (min … max):   91.685 s … 91.961 s    5 runs

Relative performances
.....................
Relative time to completion of executables distributed on bioconda.
(Compiled using GCC). Lower is better

=============== ==============
program         relative time
=============== ==============
BWA             100.00%
bwa-mem2.sse41  77.97%
bwa-mem2.sse42  69.22%
bwa-mem2.avx    64.16%
bwa-mem2.avx2   53.98%
=============== ==============

Relative time to completion of bwa-mem2 given compilation methods.
This only compares bwa-mem2.avx2 as it was the fastest on the test system.

=================================== =====================
method                              relative time
=================================== =====================
Conda distributed (GCC-compiled)    100.00%
Github distributed (ICPC-compiled)  96.23%
Compiled locally (GCC)              98.35%
Compiled locally (AOCC)             94.78%
=================================== =====================

Conclusion
----------
Bwa-mem2 2.1 delivers a substantial decrease in time to completion for a
three-fold increase in memory use. The three-fold increase in memory seems big
but for the human genome this means going from 6.3 to 19 GB of memory usage.
Given the circumstances in which bwa-mem2 is most practical (8 threads or
higher for WGS samples) the core-to-thread ratio of 2.5G or lower per thread
is very reasonable.

Bwa-mem2 successfully utilizes SIMD instructions which means that it performs
much better on more modern architectures. However this should also be in reach
for most institutions as AVX2 was introduced with Intel's Haswell architecture
in 2013. AVX512 could not be benchmarked on the test machine. It was introduced
in Skylake architecture server processors in 2017.

The bwa-mem2 README's claim that the Intel-compiled binaries perform better is
verified even on this AMD test machine. However, the improvement is about 4%.
Compiling locally using an open-source compiler gives a smaller improvement of
less than 2%. Using a proprietary compiler that was optimized by the CPU vendor
(AMD) led to a performance increase of 5%. These percentages are very small.

The bioconda distribution comes with numerous advantages such as an open
toolchain and automatic update of a biocontainer for reproducibility purposes.
It is debatable whether small performance increases due to using better
compilers are worth losing these perks. This may differ on other types of
processors used.

Sorting tools
+++++++++++++

Which sorting tool should be used?
----------------------------------
There are three contenders: ``samtools sort``, ``sambamba sort`` and
``picard SortSam``. Testing was performed on a 760 MB unsorted BAM file. Memory
settings and MAX_RECORDS_IN_RAM were set to low levels to simulate a 
WGS scenario where the resulting files are so big they cannot possible kept in
memory and need to be split in multiple files.

Test results
............

Samtools sort and samtools index were run in the same run because sambamba and
picard can index on the fly. So indexing needs to be done as well for fair 
comparison.

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'
      Time (mean ± σ):     41.780 s ±  0.339 s    [User: 37.899 s, System: 1.833 s]
      Range (min … max):   41.460 s … 42.287 s    5 runs

For picard we use ``VALIDATION_STRINGENCY=SILENT``. We can assume the aligner
outputs correct BAM formatted records. The tool only needs coordinates or names
depending on the sort order.

.. code-block::

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate COMPRESSION_LEVEL=1'
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate COMPRESSION_LEVEL=1
      Time (mean ± σ):     33.420 s ±  0.369 s    [User: 46.486 s, System: 1.750 s]
      Range (min … max):   33.152 s … 34.067 s    5 runs

For comparison here is ``picard SortSam`` with the default validation:

.. code-block::

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 SORT_ORDER=coordinate COMPRESSION_LEVEL=1'
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 SORT_ORDER=coordinate COMPRESSION_LEVEL=1
      Time (mean ± σ):     39.163 s ±  0.372 s    [User: 55.183 s, System: 1.956 s]
      Range (min … max):   38.688 s … 39.728 s    5 runs

Sambamba has an interface similar to samtools and was run with the same settings.

.. code-block::

    hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba sort -t0 -m 128M -l 1 -o test.bam unsorted.bam'
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba sort -t0 -m 128M -l 1 -o test.bam unsorted.bam
      Time (mean ± σ):     89.088 s ±  0.770 s    [User: 85.938 s, System: 1.622 s]
      Range (min … max):   87.847 s … 89.940 s    5 runs

``picard SortSam`` seems to be the fastest but a quick look at user time shows
that it is slower than ``samtools sort``. ``picard`` seems to implement some
multithreading. Also the resulting file size for 
``samtools`` (744 MB) is smaller than that of ``picard`` (960 MB).

Since the sorting tool should be used in a pipe behind an aligner, a sort 
tool that uses the least CPU time is preferred, as more time can go towards the 
alignment. Also the output filesize is preferably small. ``samtools`` is the
best tool to use here.

Sambamba is not in the same ball park as the other two tools with regards to
sorting and should therefore not be considered.

How much threads should be used?
--------------------------------

1

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'
      Time (mean ± σ):     42.022 s ±  0.319 s    [User: 38.012 s, System: 1.860 s]
      Range (min … max):   41.720 s … 42.539 s    5 runs


2

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@2 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@2 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'
      Time (mean ± σ):     23.124 s ±  0.244 s    [User: 41.531 s, System: 3.432 s]
      Range (min … max):   22.764 s … 23.355 s    5 runs

3

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@3 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@3 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'
      Time (mean ± σ):     19.528 s ±  1.319 s    [User: 43.423 s, System: 3.657 s]
      Range (min … max):   18.475 s … 21.207 s    5 runs


Using additional threads decreases the wall clock time but increases the total
CPU time.

How much memory should be used?
-------------------------------
All sorting tools work in the following way:

- A file is read in. The reads are sorted in a in-memory buffer.
- Once the buffer is full, it is written to disk to a tmp file.
- Once the entire file is read all tmp files and the memory buffer are merged.

If the sorting tool can hold the entire BAM into memory then no disk I/O is 
needed, giving significantly better performance.

When the BAM file is bigger than the in-memory buffer, part of it will be
written to disk. In WGS the BAM file to be sorted is usually very big.
160GB for a level 1 compressed BAM file is not uncommon. Sorting this file
with a 4GB in-memory buffer will create ~150 ~1GB temporary files (these are
compressed 4GB bam files).

Increasing memory does only affect the number of temporary files written. The
number of temporary files does not have a significant impact on the time as
most of the time is spent sorting. The number of temporary files written is
important, is the maximum number of open file handles may be reached.
Using a default in the 2-4GB range seems reasonable for sorting WGS BAM data.

When should the bam be sorted?
------------------------------

The BAM should be sorted directly after alignment using a unix pipe. 
Writing the BAM to a file and then using sort afterwards is a waste. The sort 
algorithm will chunk up the bam file in sorted small bam files before merging
these in the resulting bam file. Therefore a sort algorithm will write the 
entire bam file to disk twice. To not use a pipe from the aligner will increase
that to three times. Also additional time will be needed to compress and 
decompress the file from disk.

Marking duplicates
++++++++++++++++++

Which program should be used?
-----------------------------

Samtools, picard and sambamba can all mark duplicate reads. Samtools
requires a more complex pipeline. GATK may have some requirements on how the
duplicates are marked. Therefore ``picard MarkDuplicates`` is a good candidate
as well as ``sambamba markdup`` which promises to work according to the
picard criteria.

Test results
............

For testing the samtools sorted test.bam was used.

.. code-block::

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics'
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     51.616 s ±  0.325 s    [User: 59.094 s, System: 1.549 s]
      Range (min … max):   51.168 s … 51.959 s    5 runs

.. code-block::

    hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam
      Time (mean ± σ):     86.467 s ±  0.567 s    [User: 83.899 s, System: 1.494 s]
      Range (min … max):   86.023 s … 87.431 s    5 runs


Sambamba requires more CPU seconds 84 vs 59 for picard. But, the Picard file
is significantly bigger 960M vs 766M. That's a big difference, especially when
handling big WGS files. This can be multiple gigabytes.

Picard uses the intel deflater by default which gives very large files for
compression level 1. We can also use the jdk deflater which should yield the
same filesize as sambamba.

.. code-block::

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics USE_JDK_INFLATER=true USE_JDK_DEFLATER=true' && du -h markdup.bam
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics USE_JDK_INFLATER=true USE_JDK_DEFLATER=true
      Time (mean ± σ):     65.913 s ±  0.449 s    [User: 73.844 s, System: 1.458 s]
      Range (min … max):   65.480 s … 66.503 s    5 runs


This generates files of 765M which is virtually the same as the 766M
by sambamba. But compute time  for picard (74 vs 89) is better.

If you feel better using the intel deflater and inflater compression level 3
yields similar results.

.. code-block::

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=2 METRICS_FILE=markdup.metrics' && du -h markdup.bam
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=2 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     47.503 s ±  0.661 s    [User: 54.885 s, System: 1.463 s]
      Range (min … max):   46.647 s … 48.283 s    5 runs

    960M	markdup.bam

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=3 METRICS_FILE=markdup.metrics' && du -h markdup.bam
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=3 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     76.495 s ±  0.268 s    [User: 85.189 s, System: 1.381 s]
      Range (min … max):   76.171 s … 76.848 s    5 runs

    742M	markdup.bam

It is slightly slower 77 vs 74 seconds with a slightly smaller bam file 753M vs 765M.

For comparison here is picard's execution time with default settings.

.. code-block::

    $ hyperfine -r 2 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=5 METRICS_FILE=markdup.metrics'
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=5 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     103.899 s ±  0.180 s    [User: 111.700 s, System: 1.445 s]
      Range (min … max):   103.772 s … 104.026 s    2 runs

Which generates a bam file of 711M.

But sambamba has a multithreaded advantage. How does it scale?

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam
      Time (mean ± σ):     87.504 s ±  0.489 s    [User: 84.816 s, System: 1.547 s]
      Range (min … max):   86.847 s … 88.011 s    5 runs

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 2 -l 1 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 2 -l 1 test.bam markdup.bam
      Time (mean ± σ):     51.343 s ±  0.312 s    [User: 87.543 s, System: 1.984 s]
      Range (min … max):   50.995 s … 51.812 s    5 runs

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 3 -l 1 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 3 -l 1 test.bam markdup.bam
      Time (mean ± σ):     42.527 s ±  0.287 s    [User: 95.170 s, System: 2.402 s]
      Range (min … max):   42.220 s … 42.943 s    5 runs

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 4 -l 1 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 4 -l 1 test.bam markdup.bam
      Time (mean ± σ):     37.447 s ±  0.494 s    [User: 99.267 s, System: 2.478 s]
      Range (min … max):   37.000 s … 38.188 s    5 runs

Sambamba uses some threads effectively to its advantage, however diminishing
returns hit quickly. At 4 threads CPU time/ wall clock time is less than 3. Meaning
that 1 thread is not sufficiently utilized.

Maybe the scaling gets better when using better compression?

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 4 -l 4 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 4 -l 4 test.bam markdup.bam
      Time (mean ± σ):     43.028 s ±  0.428 s    [User: 118.324 s, System: 2.617 s]
      Range (min … max):   42.612 s … 43.727 s    5 runs

CPU time/ wall clock time is still less than 3. So it is not possible to fully
utilize the 4 threads. Requiring heavier compression does not increase
thread utilization.

Since x threads do not utilize x cores, we can see what the effect would be
if we run sambamba with x threads on a 2-core machine by using taskset.

.. code-block::

    $ hyperfine -w 2 -r 5 "taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam"
    Benchmark #1: taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam
      Time (mean ± σ):     81.394 s ±  0.252 s    [User: 78.575 s, System: 1.279 s]
      Range (min … max):   81.172 s … 81.828 s    5 runs

    $ hyperfine -w 2 -r 5 "taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 2 -l 1 test.bam markdup.bam"
    Benchmark #1: taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 2 -l 1 test.bam markdup.bam
      Time (mean ± σ):     48.222 s ±  0.125 s    [User: 80.145 s, System: 1.734 s]
      Range (min … max):   48.091 s … 48.380 s    5 runs

    $ hyperfine -w 2 -r 5 "taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 3 -l 1 test.bam markdup.bam"
    Benchmark #1: taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 3 -l 1 test.bam markdup.bam
      Time (mean ± σ):     45.894 s ±  0.230 s    [User: 79.484 s, System: 1.940 s]
      Range (min … max):   45.651 s … 46.204 s    5 runs

    $ hyperfine -w 2 -r 5 "taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 4 -l 1 test.bam markdup.bam"
    Benchmark #1: taskset -c 0,1 singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 4 -l 1 test.bam markdup.bam
      Time (mean ± σ):     47.583 s ±  0.085 s    [User: 80.978 s, System: 2.030 s]
      Range (min … max):   47.455 s … 47.666 s    5 runs

We can see that using 3 threads on 2 cores leads to te lowest wall clock time.
Being able to better utilize the CPU resources.

Xeon server results
...................

Picard
.. code-block::

    Intel compression level 5
    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=5 METRICS_FILE=markdup.metrics' && stat -c %s markdup.bam | numfmt --to=iec-i
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=5 METRICS_FILE=markdup.metrics
    Time (mean ± σ):     202.168 s ±  6.038 s    [User: 241.791 s, System: 8.639 s]
      Range (min … max):   191.935 s … 207.536 s    5 runs
    711Mi

    Intel compression level 1
    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics' && stat -c %s markdup.bam | numfmt --to=iec-i
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     112.341 s ±  4.293 s    [User: 159.056 s, System: 7.862 s]
      Range (min … max):   106.777 s … 117.376 s    5 runs
    960Mi

    Intel compression level 2
    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=2 METRICS_FILE=markdup.metrics' && stat -c %s markdup.bam | numfmt --to=iec-i
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=2 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     124.673 s ±  6.990 s    [User: 167.482 s, System: 11.114 s]
      Range (min … max):   112.443 s … 129.351 s    5 runs
    960Mi

    Intel compression level 3
    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=3 METRICS_FILE=markdup.metrics' && stat -c %s markdup.bam | numfmt --to=iec-i
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=3 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     163.478 s ±  5.353 s    [User: 204.590 s, System: 8.764 s]
      Range (min … max):   160.655 s … 173.039 s    5 runs
      Warning: Statistical outliers were detected. Consider re-running this benchmark on a quiet PC without any interferences from other programs. It might help to use the '--warmup' or '--prepare' options.
    742Mi

    JDK level 1
     hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics USE_JDK_INFLATER=true USE_JDK_DEFLATER=true' && stat -c %s markdup.bam | numfmt --to=iec-i
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=1 METRICS_FILE=markdup.metrics USE_JDK_INFLATER=true USE_JDK_DEFLATER=true
      Time (mean ± σ):     147.671 s ±  5.456 s    [User: 194.152 s, System: 8.623 s]
      Range (min … max):   141.741 s … 156.001 s    5 runs
    765Mi

sambamba

.. code-block::

    hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam"  && stat -c %s markdup.bam | numfmt --to=iec-i
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 1 -l 1 test.bam markdup.bam
      Time (mean ± σ):     143.724 s ±  1.610 s    [User: 136.320 s, System: 4.303 s]
      Range (min … max):   141.083 s … 145.420 s    5 runs

    766Mi

Conclusion
..........
Using picard with compression level 1 and defaults will yield a very big file
which is undesirable. When the jdk deflater is utilized the file size is
similar to sambamba and we can make an apples-to-apples comparison.

When utilizing a single core picard is most efficient with 74 seconds vs 84 for
sambamba. Sambamba can use more cores however and finish in a wall clock time
of 51 seconds when utilizing 2 cores. Beyond that it does not scale well.

A thing to consider is the memory usage. Sambamba will not use more than
3 gigabyte of memory on its default settings. Broad institute runs Picard
MarkDuplicates on VMs with 7 GB.

Another thing to consider is that Picard outputs a metrics file stating how
many reads were marked as duplicate.

It is not a decisive win for any application.

======= ====== =============== ================
.       Picard sambamba 1 core sambamba 2 cores
======= ====== =============== ================
cores   1      1               2
memory  7G     3G              3G
time    74s    84s             51s
metrics yes    no              no
======= ====== =============== ================

Merging BAM files
+++++++++++++++++

Which program should be used
-----------------------------
Samtools, sambamba and picard all provide tools for merging several sorted
bam files together.

Two test BAMs with each 5 million reads were generated. The test BAMs were
distinct from each other. Both BAM files were sorted using samtools.

Test results
............

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools merge -f -l1 -@0 merged.bam test.bam test2.bam && samtools index merged.bam'"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools merge -f -l1 -@0 merged.bam test.bam test2.bam && samtools index merged.bam'
      Time (mean ± σ):     45.746 s ±  0.223 s    [User: 42.968 s, System: 1.143 s]
      Range (min … max):   45.547 s … 46.012 s    5 runs

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba merge -l1 -t0 merged.bam test.bam test2.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba merge -l1 -t0 merged.bam test.bam test2.bam
      Time (mean ± σ):     91.850 s ±  1.192 s    [User: 87.480 s, System: 2.778 s]
      Range (min … max):   90.252 s … 93.187 s    5 runs


Alignment with BWA
+++++++++++++++++++

To pipe or not to pipe?
-----------------------

Pipes can be used to couple BWA with the postalt.js script and samtools sort
in order to let the alignment process output a sorted BAM file. This massively
limits the amount of files written to disk which is a great advantage.

A disadvantage of pipes is that they are blocking when read while empty or
while written and full. The default pipe size on Linux is 64K. It is possible
for a pipe to get clogged up repeatedly and thus lengthen computation time.

Scripts
-------
These scripts were used:

+ 8 cores BWA, postalt script, samtools sort 3 cores. Sequentially.
+ 8 cores BWA, postalt script, samtools sort 3 cores. Pipes.
+ 8 cores BWA, postalt script, samtools sort 3 cores. Named pipes.
+ 8 cores BWA, postalt script, samtools sort 3 cores. Named pipes altered with
  [mkbigfifo](https://github.com/biowdl/mkbigfifo).
+ 12 cores BWA, postalt script, samtools sort 4 cores. Sequentially.
+ 12 cores BWA, postalt script, samtools sort 4 cores. Pipes.
+ 12 cores BWA, postalt script, samtools sort 4 cores. Named pipes.
+ 12 cores BWA, postalt script, samtools sort 4 cores. Named pipes altered with
  [mkbigfifo](https://github.com/biowdl/mkbigfifo).
Test results
------------

.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_no_pipes_8cores.sh'
    Benchmark #1: bash bwa_no_pipes.sh
      Time (mean ± σ):     205.159 s ±  1.112 s    [User: 1186.185 s, System: 8.004 s]
      Range (min … max):   203.323 s … 206.748 s    10 runs

    $ du -h ramdisk/*
    435M	ramdisk/no_pipes.aln.bam
    2,4G	ramdisk/no_pipes.postalt.sam
    2,2G	ramdisk/no_pipes.sam

.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_with_pipes_8cores.sh'
    Benchmark #1: bash bwa_with_pipes.sh
      Time (mean ± σ):     171.717 s ±  0.520 s    [User: 1240.633 s, System: 9.746 s]
      Range (min … max):   170.844 s … 172.695 s    10 runs

    $ du -h ramdisk/with_pipes.aln.bam
    435M	ramdisk/with_pipes.aln.bam

.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_fifo_8cores.sh'
    Benchmark #1: bash bwa_fifo_8cores.sh
      Time (mean ± σ):     171.638 s ±  1.279 s    [User: 1247.522 s, System: 10.922 s]
      Range (min … max):   170.079 s … 174.891 s    10 runs

.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_mkbigfifo_8cores.sh'
    Benchmark #1: bash bwa_mkbigfifo_8cores.sh
      Time (mean ± σ):     170.907 s ±  2.347 s    [User: 1240.188 s, System: 6.786 s]
      Range (min … max):   168.868 s … 175.535 s    10 runs


.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_no_pipes_12cores.sh'
    Benchmark #1: bash bwa_no_pipes_12cores.sh
      Time (mean ± σ):     155.212 s ±  0.614 s    [User: 1271.680 s, System: 8.780 s]
      Range (min … max):   154.406 s … 156.669 s    10 runs

    $ du -h ramdisk/*
    435M	ramdisk/no_pipes.aln.bam
    2,4G	ramdisk/no_pipes.postalt.sam
    2,2G	ramdisk/no_pipes.sam


.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_with_pipes_12cores.sh'
    Benchmark #1: bash bwa_with_pipes_12cores.sh
      Time (mean ± σ):     123.853 s ±  0.928 s    [User: 1295.435 s, System: 10.535 s]
      Range (min … max):   122.943 s … 125.636 s    10 runs

    $ du -h ramdisk/with_pipes.aln.bam
    435M	ramdisk/with_pipes.aln.bam

.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_fifo_12cores.sh'
    Benchmark #1: bash bwa_fifo_12cores.sh
      Time (mean ± σ):     122.921 s ±  0.541 s    [User: 1288.907 s, System: 11.465 s]
      Range (min … max):   122.244 s … 123.700 s    10 runs

.. code-block::

    $ hyperfine -w 2 -r 10 'bash bwa_mkbigfifo_12cores.sh'
    Benchmark #1: bash bwa_mkbigfifo_12cores.sh
      Time (mean ± σ):     125.581 s ±  0.955 s    [User: 1303.329 s, System: 8.499 s]
      Range (min … max):   123.360 s … 126.660 s    10 runs
