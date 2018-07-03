

## Module 1
Hard-coded file names were used for pacbio and fastq data, so I updated them to use the filenames in config, at least initially as they are copied into the $TMPDIR.  
In these snippets, commented code was the original.
```
# cp ${pacDIR}/raw_pac.fasta $TMPDIR
cp ${pacM4} ${TMPDIR}/raw_pac.fasta
```
and
```
# cp ${illDIR}/raw_f.fastq XXXXX/raw_f.fastq
cp $illF XXXXX/raw_f.fastq
```
and
```
$cp ${illDIR}/raw_m.fastq XXXXX/raw_m.fastq
cp $illM XXXXX/raw_m.fastq
```

The qsub/sbatch command generated for aligning male fastq to the MtDNA has as its last step the deletion of the tmp directory. This is problematic because it is very possible that the female alignment is still running, or hasn't even started yet, causing that process to fail. For the timebeing, I'm just commenting it out, and allowing the tmpdir to persist.

The original scripts called for 32G memory and 20 cores to run mitochondrial alignments. This seems like overkill, and greatly restricts where it will run on the ITTC resources, which don't have some many big cores.  I've reduced it to 10G and 8 cores, so it gets onto a node more quickly. 


## Module 2
modified initial resource request to 6 cores and 10G.  Original request seems unnecessarily high since this is just a wrapper for submitting other jobs.  


Module 2 starts out with a strange step, copying the original raw fasta into the $TMPDIR. This seems to have no purpose, because the code proceeds to operate on the filtered data file, which is hard-coded as m_pac.fasta. I think this is an error, because m_pac.fasta wouldn't exist in a new tmp directory. It just happens that the default configuration file has the "raw" pacbio data file defined as m_pac.fas in the configuration, and hard-coded "raw_pac.fasta" in module 1 to get the original data into play.
  
So the simple solution (probably the original intention) is for module 2 to copy "m_pac.fasta" from the pacbio directory. This way it is operating on the filtered data as well.
```
# cp $pacM $TMPDIR
cp ${pacDIR}/m_pac.fasta $TMPDIR/m_pac.fasta  
```

Reduced resource request for node array to 8 cores and 10G mem, since very few nodes will have the resources as initially configured.
  
For bowtie alignmments of filtered reads in arrays, the original fastq files are copied, but names are not updated.  Should be using filtered files produced by first module.  
```
 # cp $illM XXXXXTMPDIR
   cp ${illDIR}/m.fastq  XXXXXTMPDIR/m.fastq  # correcting so mt-filtered reads are used.                 

 # cp $illF XXXXXTMPDIR
   cp ${illDIR}/f.fastq XXXXXTMPDIR   # correcting so mt-filtered reads are used.
```


## Module 3

In the merging of separate unique read count files, the original code explicitly reverses ordering of columns, from Scaffold,Count to Count,Scaffold. However, subsequent steps in the pipeline seemingly assumed the original ordering of columns.  As is, no useful output is generated.  But if I swap the awk statement to keep the initial ordering of columns it all seems to work. For instance, in the sorting and merging of the files, the (Count,Sort) ordering sorts and joins  by Count, which really makes no sense; it should be by scaffold.  
  
Here is the updated code:  

```
cat $CWD/pacBio_illmapping/mapping_rawdata/*_female_uniq | awk '{print $1, $2}'> $CWD/pacBio_illmapping/mapping_rawdata/female_unsort    
cat $CWD/pacBio_illmapping/mapping_rawdata/*_male_uniq | awk '{print $1, $2}'> $CWD/pacBio_illmapping/mapping_rawdata/male_unsort 
```


For calculating library sizes, should be using mt-filtered fastq files for getting library size, not original raw data. This probably a minor difference, though depends on what portion of Illumina data correspond to MtDNA

```
# illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBMsize=$(wc -l ${illDIR}/m.fastq | awk '{print ($1/4)}')
# illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illLIBFsize=$(wc -l ${illDIR}/f.fastq | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))
```


