# ngs-performance-choices
Which tools to choose based on various performance metrics.

## Performance testing.

For fair comparison between tools it is best to use containerized environments.
These eliminate any user specific installs or compilations. 

For running the tools [Singularity](https://github.com/hpcng/singularity) was 
used. The containers were provided by [biocontainers](
https://quay.io/biocontainers).

Tools were run with `singularity exec -eip 
docker://quay.io/biocontainers/<tool>:<version-tag> command`. The `-e`, `-i` 
and `-p` flags isolate environment, IPC namespace and PID namespace 
respectively. By default singularity mounts `$PWD` and user `HOME` directories,
which was convenient when testing.

For benchmarking [hyperfine](https://github.com/sharkdp/hyperfine) was used.
This command line tool is useful as it can perform warmup runs to ensure 
images are downloaded and disk caches are full.

## Compression and compression levels

### Compression tools
The bioinformatics world seems to have standardised on [zlib](
https://www.zlib.net/). Most bioinformatics tools can open `.gz` files and 
output to them. Therefore tools that implement zlib should be used for compressing files.

Not all zlib implementations are created equal.
The most well-known zlib implementation is [gzip](http://www.gzip.org/). Gzip
is very old and is written as a single-threaded algorithm. Another 
implementation is [pigz](http://zlib.net/pigz/), which was written to make use 
of multiple threads. Pigz is **always** faster than gzip. Even on a single core
the algorithm is simply more efficient. `pigz -p 1` (no multithreading) beats
`gzip` in both compression and decompression. [Cutadapt](
https://github.com/marcelm/cutadapt) therefore does its compression and 
decompression with pigz.

### Which compression level to use?
There are several zlib compression levels, from 1 to 9. Alternatively a file
can not be compressed. Should a file be compressed, and if so, at which level?

 

## Sorting tools.

### Which sorting tool should be used?
There are three contenders: `samtools sort`, `sambamba sort` and 
`picard SortSam`. Testing was performed on a 760 MB unsorted BAM file. Memory 
settings and MAX_RECORDS_IN_RAM were set to low levels to simulate a 
WGS scenario where the resulting files are so big they cannot possible kept in
memory and need to be split in multiple files.

Testing results:

Samtools sort and samtools index were run in the same run because sambamba and
picard can index on the fly. So indexing needs to be done as well for fair 
comparison.
```
$ hyperfine -w 2 -r 5 "singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'"
Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 bash -c 'samtools sort -@0 -m 128M -l 1 -o test.bam unsorted.bam && samtools index test.bam test.bai'
  Time (mean ± σ):     41.780 s ±  0.339 s    [User: 37.899 s, System: 1.833 s]
  Range (min … max):   41.460 s … 42.287 s    5 runs
```

For picard we use `VALIDATION_STRINGENCY=SILENT`. We can assume the aligner 
outputs correct BAM formatted records. The tool only needs coordinates or names
depending on the sort order.
```
$ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate COMPRESSION_LEVEL=1'
Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate COMPRESSION_LEVEL=1
  Time (mean ± σ):     33.420 s ±  0.369 s    [User: 46.486 s, System: 1.750 s]
  Range (min … max):   33.152 s … 34.067 s    5 runs
```
For comparison here is `picard SortSam` with the default validation:
```
$ hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 SORT_ORDER=coordinate COMPRESSION_LEVEL=1'
Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/picard:2.23.1--h37ae868_0 picard -Xmx3G -XX:ParallelGCThreads=1 SortSam INPUT=unsorted.bam OUTPUT=test.bam CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000 SORT_ORDER=coordinate COMPRESSION_LEVEL=1
  Time (mean ± σ):     39.163 s ±  0.372 s    [User: 55.183 s, System: 1.956 s]
  Range (min … max):   38.688 s … 39.728 s    5 runs
```

Sambamba has an interface similar to samtools and was run with the same settings.
```
hyperfine -w 2 -r 5 'singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba sort -t0 -m 128M -l 1 -o test.bam unsorted.bam'
Benchmark #1: singularity exec -eip docker://quay.io/biocontainers/sambamba:0.7.1--h148d290_2 sambamba sort -t0 -m 128M -l 1 -o test.bam unsorted.bam
  Time (mean ± σ):     89.088 s ±  0.770 s    [User: 85.938 s, System: 1.622 s]
  Range (min … max):   87.847 s … 89.940 s    5 runs
```

`picard SortSam` seems to be the fastest but a quick look at user time shows 
that it is slower than `samtools sort`. `picard` seems to implement some 
multithreading. Also the resulting file size for 
`samtools` (744 MB) is smaller than that of `picard` (960 MB). 

Since the sorting tool should be used in a pipe behind an aligner, a sort 
tool that uses the least CPU time is preferred, as more time can go towards the 
alignment. Also the output filesize is preferably small. `samtools` is the best
tool to use here. 

Sambamba is not in the same ball park as the other two tools with regards to
sorting and should therefore not be considered.

### How much threads should be used?
???

### How much memory should be used?
All sorting tools work in the following way:

- A file is read in. The reads are sorted in a in-memory buffer.
- Once the buffer is full, it is written to disk to a tmp file.
- Once the entire file is read all tmp files and the memory buffer are merged.

If the sorting tool can hold the entire BAM into memory then no disk I/O is 
needed, giving significantly better performance.

The percentage 

### When should the bam be sorted?

The BAM should be sorted directly after alignment using a unix pipe. 
Writing the BAM to a file and then using sort afterwards is a waste. The sort 
algorithm will chunk up the bam file in sorted small bam files before merging
these in the resulting bam file. Therefore a sort algorithm will write the 
entire bam file to disk twice. To not use a pipe from the aligner will increase
that to three times. Also additional time will be needed to compress and 
decompress the file from disk.


