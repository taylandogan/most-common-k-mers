# most-common-k-mers

#####Name: MostCommonK-mers
#####Author: Taylan Dogan
#####License: The MIT License (MIT)
#####Requirements: C++11

###Description:
Given a .fastq file, MostCommonK-mers is a small piece of code that attempts to extract most common N K-mers.
Due to complexity issues, it eliminates less common K-mers much more harshly as the input size gets bigger and bigger,
hence its accuracy reduces. This accuracy problem could be solved by dividing the input file into equal sized parts(~100MB),
then merging the results of these parts.

The algorithm gives you the full results as a file. But it prints the top N to the command line, along with the counts of K-mers.
These counts of K-mers are not exact, the program misses some K-mers here and there.

###Attention: 
**This algorithm may find inaccurate results!** It may miss some common K-mers, especially when a K-mer is distributed over the file uniformly.
Hopefully, this issue will be fixed when the items in todo section are completed. 
Also, note that it just uses the genomic sequences in the given .fastq file. (It does not consider the quality information or other info.)

###Information Concerning Execution Times:
(CPU: Intel i7-2670QM)

Input File              | # of lines| Size of the file  | Execution Time

ERR068396.filt.fastq    |   ~270K   |       17.2 MB     |   ~3.2 sec

ERR059924.filt.fastq    |   ~283K   |       18.1 MB     |   ~3 sec

ERR059932.filt.fastq    |   ~636K   |       40.6 MB     |   ~9 sec

ERR047698.filt.fastq    |   ~800K   |       47.4 MB     |   ~11 sec

ERR055763_1.filt.fastq  |   ~70M    |       4.5 GB      |   ~28.5 min

*Note that you may decrease the execution time by reducing the fairness constant.
Though, a higher fairness constant would give us more accurate results.
The value of the fairness constant should be between 1.0 and 2.0. Default value is 1.25.

*Increasing K will most likely shorten the execution time too.

###How to Use:
An example is given below.
Arguments: [input .fastq file] [K] [N] [fairness const] [output file]

- make
- ./MostCommonK-mers ERR047698.filt.fastq 30 25 1.25 'ERR047698_output.txt'

###To Do:
- Divide the input file into segments(~100MB), execute the algorithm for each segment, then merge the results.

