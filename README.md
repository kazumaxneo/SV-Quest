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

    
## Requirements  
- SAMTools  (version >= 1.3.1)   
- BWA (version >= 0.7.17)  
- circos (v0.67. only required if drawing indel map))  


Install [HomeBrew](http://brew.sh/) (Mac OS X) or [LinuxBrew](http://brew.sh/linuxbrew/) (Linux).  
```
brew tap brewsci/science
brew install bwa
brew install samtools
brew install circos
```
    


## Source
```
cd $HOME 
git clone git@github.com:kazumaxneo/SV-Quest.git
echo': echo export PATH=\$PATH:`pwd`\ >> ~/.bash_profile && source ~/.bash_profile
SV-Quest.pl
```
    


## Calling SVs
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


