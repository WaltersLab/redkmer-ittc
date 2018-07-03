

# Module 1
Hard-coded file names were used for pacbio and fastq data, so I updated them to use the filenames in config, at least initially as they are copied into the $TMPDIR
In these snippets, commeted code was the original.
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

# Module 2
Module 2 starts out with a strange step, copying the original raw fasta into the $TMPDIR. This seems to have no purpose, because the code proceeds to operate on the filtered data file, which is hard-coded as m_pac.fasta. I think this is an error, because m_pac.fasta wouldn't exist in a new tmp directory. It just happens that the default configuration file has the "raw" pacbio data file defined as m_pac.fas in the configuration, and hard-coded "raw_pac.fasta" in module 1 to get the original data into play.
  
So the simple solution (probably the original intention) is for module 2 to copy "m_pac.fasta" from the pacbio directory. This way it is operating on the filtered data as well.
```
# cp $pacM $TMPDIR
cp ${pacDIR}/m_pac.fasta $TMPDIr/m_pac.fasta  
```

 


