MMBIRFinder
===========

MMBIRFinder is a bioinformatics tool to detect microhomology-mediated break-induced replication(MMBIR) events. 

Specifically, MMBIRFinder detects template-switching events associated with MMBIR.

The tool writes to a Log file for easy recording of important parameters, run time settings, and results.

Installation
--------------

The included executable should work in most Linux environments. Otherwise, the source code can be compiled using the included Make file.
```sh
make -f Makefile
```

Usage
--------------
Run the program by invoking the executable followed by the configuration file.

```sh
./mmbirFinder configFile [options]
    -?                  Command line options
    configFile          Program configuration file name
    -j jobIdString      Identifier to be used for this job
    -OconfOption=Value  Override options in the configuration file
```


Configuration File
--------------

    --- General parameters ---
    sLogDirName                The directory where log files should be stored
    projectDirectory           The directory to store and retrive all created files. 
                               Leave empty for current working directory.
    referenceFile              The FASTA reference file.
    chromosome                 Limit the search to a single chromosome. Enter 0 to search them all.
    readsFile                  The FASTA or FASTQ reads file.
    alignedFile                The output alignment file from BWA.
    unalignedFile              The SAM and FASTA file name for the unaligned reads.
    outputFile                 The name of the SAM file for the BWA alignment output.
    locationsFile              The name of the temporary file for the possible MMBIR locations.
    clusterFile                The name of the file to output the cluster information.
    pairedEnd                  Is the dataset paired-end (yes/no)?
    
    --- Alignment steps parameters ---
    - All are in the format true/false -
    runAlignment               A testing method to save time.
    onlyAlign                  Only run the alignment, not the mmbirFinder. 
                               ** WARNING *** If true, runAlignment must also be true
    indexGenome                Index the genome
    fullAlign                  Run the full alignment on the full-read dataset.
                               Usefule for when you already have an aligned BAM/SAM file.
    bamFile                    Is the alignment in a BAM format?
    extractUnalignedReads      Extract unaligned reads? Should set to TRUE, unless error occurred in earlier step.
    halfAlign                  Perform half read alignment? Should set to TRUE, unless error occurred in earlier step.
    extractHalfReads           Extract half reads? Should set to TRUE, unless error occurred in earlier step.
    filterOut                  Filter out reads? Should set to TRUE, unless error occurred in earlier step.
    
    --- Clustering parameters ---
    performClustering          Should the program perform clustering? 
    minConsolidate             The minimum number of reads (evidence) to consolidate (IMPORTANT).
    mysql                      Should the pipeline use MySQL to find reads? If true, the program will stop and the 
                               included MySQL script is required to continue. The MySQL is faster, but requires more 
                               manual steps.
    mysqlFile                  The mysql_results.txt file from the script output
    
    --- Search parameters ---
    searchLength               The length in b.p. to search for a possible template strand upstream from possible MMBIR.
    minSeqLength               The minimum length of the split-read DNA sequence.
    minBirLength               The minimum length of a MMBIR region (IMPORTANT).
    minAlignedLength           The minimum length of the aligned region to be considered a successful MMBIR.
    missCount                  For the FSM, the number of CONSECUTIVE misses to open a MMBIR region.
    hitCount                   For the FSM, the number of CONSECUTIVE hits to close a MMBIR region.
    tolerance                  The distance from the beginning of the read to consider a MMBIR region.




License
----

GPLv2


Version
----

0.1
