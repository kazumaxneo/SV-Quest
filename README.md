######################################################################################## 
SV-Quest 1.0

A Perl scripts to call indel position from mapped.bam.   

SV Quest: Sensitive Structure Variation detection tool

Kazuma Uesaka, Hiroshi Kouno, Kazuki Terauchi, Yuichi Fujita, Tatsuo Omata, and Kunio Ihara  



Input: 
  bam file and reference.fasta for the mapping   
Outnput:	
  indel and deletion position printed to STDOUT  

Usage:  
  perl SV-Quest.pl    


 The mapped.bam and it's reference.fasta should be included in the same folder,  
 and copy the Indel_Hunter.pl in this folder.
########################################################################################

Structural variations (SVs) are large genomic rearrangements that vary significantly in size, making them challenging to detect with the relatively short reads from next-generation sequencing (NGS). A number of heuristic method have been developed to overcome these restrictions. They works well if SV subtype or SV size were fitted with their methodology. But they are not comprehensive nor tolerant for discovering all classs of SV, and demand high sequence throughput for accurate SV detection. For detecting SV in comprehensive, we developed a new method named SV Quest. Unique feature of SV Quest is the use of mismatch accumulation for detecting SV. The performance of SV Quest was compared with five well known SV detection tools. The results show that among these tools SV Quest is highly sensitive at the relatively low false discovery rates in many cases.  
The software is implement in Perl. It is freely available  at https://github.com/kazumaxneo/SV-Quest under the GPLv3 license.

<p align="center"><img src="Figure.png" alt="workflow" width="800"></p>

<p align="center"><img src="Table.png" alt="workflow" width="800"></p>
    
## Requirements  
- SAMTools  (version >= 1.3.1)  
- BWA (version >= 0.7.17)  
- sambamba  (version >= 0.8.0)  
- circos (v0.67. only required if drawing indel map))  


Install Anaconda (Mac OS X, Linux).  

```
conda install -c bioconda bwa 
conda install -c bioconda samtools 
conda install -c bioconda sambamba 
conda install -c bioconda circos 
```
    


## Source
```
cd $HOME 
git clone https://github.com/kazumaxneo/SV-Quest.git
cd SV-Quest/
echo export PATH=\$PATH:`pwd`\ >> ~/.bash_profile && source ~/.bash_profile
SV-Quest.pl
```
    


## Calling SVs and other event （e.g. mis-assembly junction, small indels）
1、call from fastq
```
./SV-Quest.pl -f refererence.fa -1 forward.fastq -2 reverse.fastq
```

2、call from bam
```
./SV-Quest.pl -f reference.fa -i error_permitted_aligned.bam
```
    
## Test run
```
chmod u+x SV-Quest.pl
mkdir test
tar zxvf sample.tar.gz -C test/
cd test/
../SV-Quest.pl -f chromosome.fasta -1 forward.fq -2 reverse.fq
```  
    


## Licence ##

GPL v3.


