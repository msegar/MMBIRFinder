MMBIRFinder
=========

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



MySQL Script
----
The included MySQL script greatly decreases the computation time to run the program. If set in the configuration file,
the program will stop after clustering the anchored reads. The program now needs to find the other half of the read (i.e.
the unaligned read). Since there are millions of unaligned reads that need to be paired with the anchored/aligned read,
it is a computationally expensive problem. 

After cluserering, the program will output a file called 'half_read_clusters.txt' where the individual columns are as followed:
1. Read name
2. Sequence
3. Location
4. Chromosome
5. Cluster

Similarly, the unaligned files from the second BWA alignment are outputted. While the filename will changed depending on the
configuration settings, the default name is 'unaligned2_2.sam'. A file named 'columns.txt' should be created using any method
that has unaligned read name in column 1 and the DNA sequence in column 2 extracted from 'unaligned2_2.sam'.
```linux
less unaligned2_2.sam | sed '1,25d' | awk '{print $1 "\t" $10}' > columns.txt
```

The following steps need to then be executed at the command prompt and MySQL:

```bash
mysql --local-infile -u root -p
```
```mysql
create database test;
create table unaligned ( `read` varchar(100) NOT NULL PRIMARY KEY, `sequence` varchar(100) NOT NULL ) ENGINE = MYISAM;
LOAD DATA LOCAL '.../path/to/file/columns.txt' INTO TABLE unaligned;

mysql> show tables;
+----------------+
| Tables_in_test |
+----------------+
| unaligned      |
+----------------+

mysql> show columns from unaligned;
+-------+--------------+------+-----+---------+-------+
| Field | Type         | Null | Key | Default | Extra |
+-------+--------------+------+-----+---------+-------+
| read  | varchar(100) | NO   | PRI | NULL    |       |
| seq   | varchar(100) | NO   |     | NULL    |       |
+-------+--------------+------+-----+---------+-------+
```

Back at the command prompt, since the engine is MYISAM, the database needs to be indexed (this may take up to 20 hours!):
```linux
myisamchk -o /path/to/mysql/test/unaligned.MYI --key_bufer_size=2G
```

Finally, run the Perl script db.pl changing the username and password as appropriate:
```perl
perl db.pl
```
The outputted file 'mysql_results.txt' is then used as an input for the second half of the MMBIRFinder tool. 

**Don't forget to changed the configuration file to not run the initial alignment and clustering again!!**  


License
----

GPLv2


Version
----

0.1


    
