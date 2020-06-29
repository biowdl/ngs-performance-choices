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

which compression level to use?
-------------------------------
There are several zlib compression levels, from 1 to 9. Alternatively a file
can not be compressed. Should a file be compressed, and if so, at which level?

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
Samtools works with additional threads.

0

.. code-block::

    $ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'
      Time (mean ± σ):     42.022 s ±  0.319 s    [User: 38.012 s, System: 1.860 s]
      Range (min … max):   41.720 s … 42.539 s    5 runs

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
CPU time. Since the samtools program will be used in a pipe with an alignment
program it is best to use ``-@0`` to prevent overhead.

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

    hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 0 -l 1 test.bam markdup.bam"
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba markdup -t 0 -l 1 test.bam markdup.bam
      Time (mean ± σ):     86.467 s ±  0.567 s    [User: 83.899 s, System: 1.494 s]
      Range (min … max):   86.023 s … 87.431 s    5 runs


Sambamba requires more CPU seconds 84 vs 59 for picard. But, the Picard file
is significantly bigger 960M vs 766M. That's a big difference, especially when
handling big WGS files. This can be multiple gigabytes.

Further investigation is needed. How much compute time is needed for picard
to get smaller files?

.. code-block::

    $ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=3 METRICS_FILE=markdup.metrics'
    Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx4G -XX:ParallelGCThreads=1 MarkDuplicates INPUT=test.bam OUTPUT=markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=3 METRICS_FILE=markdup.metrics
      Time (mean ± σ):     81.179 s ±  0.395 s    [User: 89.123 s, System: 1.453 s]
      Range (min … max):   80.819 s … 81.840 s    5 runs

Compression level 3 generates files of 742M which is close enough to the 766
by sambamba. Also in compute time (84 vs 89) seconds both are similar.

But sambamba has a multithreaded advantage. How does it scale?

