# Cluster Extraction for the DeepClust Database:

This script extract clusters from the DeepClust Database efficiently and fast with [Apache Parquet](https://parquet.apache.org/) and [DuckDB](https://duckdb.org/).

## General

The DeepClust Database is formatted into an [Apache Parquet](https://parquet.apache.org/) file.  
Parquet Files are structured into RowGroups where the data is stored column wise.
This allows extremely fast querying without the need to load the whole file into memory.
By limiting the RowGroup size of 169 125 rows, the groups are small enough to fit into ~ 1GB of memory, allowing the data extraction on laptops.

The first step in the script is to determine where in the parquet file the searched Cluster is situated.
For this, an inedxed [DuckDB](https://duckdb.org/) Database is created. 
In the second step, the clusters are then extracted from the DeepClust Database.  

Here is a benchmark done with an external SanDisk Extreme portable SSD with 4TB and 1GB I/O and an Intel Core i5-6500 CPU @ 3.2 GHz.

| Number of Clusters | Time in [s] | Number of sequences |
|--------------------|----------|---------------------|
| 1                  | 1        | 119                 |
| 100                | 5        | 7 391               |
| 1000               | 124      | 39 546              |
| 32000              | 2 608    | 1 649 332           |
| 100000             | 5 887    | 5 196 357           |

Unfortunately the Parallelization only works on Linux systems and not on macOS or windows.
The script can still be used on those systems by leaving the thread option as "1".  

The DuckDB version should be 0.10.0 or newer.

### References and helpfull links:
- [DuckDB](https://duckdb.org/) 
- [Apache Parquet](https://parquet.apache.org/)
- [Apache Arrow](https://arrow.apache.org/docs/index.html)
- [PyArrow](https://arrow.apache.org/docs/python/index.html) 
## Options:
[--centroids CENTROIDS]  
*Comma separated list of centroids to extract*  
Clusters are represented by centroid sequences. The script takes a comma separated list of centroids and extracts the corresponding member sequences.  
[--centroid_file PATH/TO/CENTROID_LIST]  
*File with centroids to extract, one centroid per line*  
The script also accepts a file where each line contains a centroid ID, both input modes can be combined.  
[path_to_DCD PATH/TO/DeepClustDatabase]  
*Location of the DeepClust database in parquet format*  
[path_to_output PATH/TO/OUTPUT/DIRECTORY]  
*Path to Output Directory*  
[path_to_index /PATH/TO/INDEX]  
*Path to index in parquet format, DuckDB persistent Database will be created here if not already existent.*   
[--per-clust-output INT]   
*0: All Sequences are written to a single FASTA file; 1: For each cluster a Fasta file is written.*   
[--threads INT]  
*Number of threads to use*   
[--max_num_of_cluster_at_once INT]
If the Memory usage is to high, it can be reduced by limiting the number of Clusters extracted at once. This refers to the number of clusters per CPU. With 8GB a limit of 2000 seems to work well.
*Maximum number of cluster to extract at once, 0 means all; Default is 0*   
[--verbose INT]   
*Print verbose outut*

## Example: 
To extract the biggest cluster with 1 960 408 Members call:  
~~~  
python3 cluster_index.py --centroids NR_UGC79169.1 --per-clust-output 1 --threads 1 --verbose 1 -path_to_DCD /PATH/TO/DeepClustParquet/joined_with_index_RowGroupFinal.parquet -path_to_output /PATH/TO/OUTPUT -path_to_index /PATH/TO/DeepClustParquet/clust_index_RowGroup.parquet 
~~~

Extraction of the cluster will be written to /PATH/TO/OUTPUT/ into a FASTA file named "NR_UGC79169.1.fa" and done on one CPU. 

## Datafile Description
The database file ("joined_with_index_RowGroupFinal.parquet) contains 4 columns and 19 388 011 922 rows: 

 | f0 (Sequence ID)                                                          | f1 (Sequence)                                         | f2 (Cluster ID)                                                           | f3 (Row Number) |
 |---------------------------------------------------------------------------|-------------------------------------------------------|---------------------------------------------------------------------------|-----------------|
 | AGNOSTOS_GB_GCA_000007185.1_AE009439.1_1000_-_906357_906515_orf-1000      | METGGTCSVRPTLIEAGDTPPYGTCPGAERVSDVRSWERDLLSKGTKLEDAV* | AGNOSTOS_GB_GCA_000007185.1_AE009439.1_1000_-_906357_906515_orf-1000      | 0               |
 | .<br/>.<br/>.                                                             | .<br/>.<br/>.                                         | .<br/>.<br/>.                                                             | .<br/>.<br/>.   |
 | metaclust_all_tagenome__1003787_1003787.scaffolds.fasta_scaffold9999973_2 | MSAALVAGFVTVLLWGSAFVGIR                               | metaclust_all_tagenome__1003787_1003787.scaffolds.fasta_scaffold9999973_2 | 19388011921     |

Since all members of a cluster are grouped and the start of a cluster is known, the row-number column can be used to search for a cluster efficiently with predicate pushdown.
Additionally, more columns could be added in the future containing more information for each sequence without influencing the cluster extraction process.



 The index file ("clust_index_RowGroup.parquet") contains 4 Columns and 335 434 803 Rows:

|RowStart|NMEMBER|CLUSTER| RowGroup      |
|-----|----|----|---------------|
|1        | 5       | AGNOSTOS_GB_GCA_000007185.1_AE009439.1_1018_+_924582_925004_orf-1018 | 0             |               
| .<br/>.<br/>.                                                             | .<br/>.<br/>.                                         | .<br/>.<br/>.                                                             | .<br/>.<br/>.<br/> |
|19388011913 | 3       | metaclust_all_tagenome__1003787_1003787.scaffolds.fasta_scaffold9999794_1 | 116209        |

Only clusters with at least 3 members are listed in the index file, since this reduces the row number extremely.
The first column indicates where the cluster starts in the database file, the second the number of members and the third one contains the cluster ID. 
In the fourth column, the RowGroup number refers to the RowGroups in the parquet file, which contains the clusters. This makes it easy to only extract the needed part of the Database file where the cluster is situated.
If a cluster extends over several RowGroups, the numbers are seperated by a ",".

Similar to the Database File, columns could be added in the future containg more information about each of the clusters.

This file is provided as a parquet file or as a persistent DuckDB database.
The script extracts indices from the DuckDB database and if only the parquet file is downloaded, the script will create the database at the same place with the name "persistent".

Inspired by MMseqs2 internal cluster format and the [FFindex](https://github.com/ahcm/ffindex) format.

## Mapping any Sequence ID onto the corresponding cluster
The file called "SeqIdMapClustId.parquet" contains two columns: SEQID and CLUSTERID.
First contains the sequence ID, and the file is sorted by this column.
The second contains the corresponding centroid ID or cluster ID.  
With [DuckDB](https://duckdb.org/), the file can be queried fast:  
~~~ 
duckdb -c "SELECT * FROM read_parquet('SeqIdMapClustId.parquet') WHERE SEQID = 'SEQUENCE ID';"
~~~ 
Replace SEQUENCE ID by the desired sequence.
