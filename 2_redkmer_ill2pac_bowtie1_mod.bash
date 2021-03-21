#!/bin/bash
#SBATCH -J redkmer2
#SBATCH -t 10:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --tmp=500G
#SBATCH -e redkmer2_%j_error.log
#SBATCH -o redkmer2_%j_output.log

source $SLURM_SUBMIT_DIR/redkmer_mod.cfg
# module load slurm-torque  

mkdir -p $CWD/pacBio_illmapping
mkdir -p $CWD/pacBio_illmapping/logs
mkdir -p $CWD/pacBio_illmapping/mapping_rawdata
mkdir -p $CWD/pacBio_illmapping/index
mkdir -p $CWD/pacBio_bins
mkdir -p $CWD/pacBio_bins/fasta
#rm -f $CWD/pacBio_illmapping/mapping_rawdata/*

echo "==================================== Generating pacBio data chunks ======================================="

# cp $pacM $TMPDIR
# cp ${pacDIR}/m_pac.fasta $TMPDIR/m_pac.fasta  # JRW: no need to copy, just read m_pac.fasta in siute
grep -n ">" ${pacDIR}/m_pac.fasta |cut -f1 -d: > ${pacDIR}/pacMsplitter
READNpacM=$(cat ${pacDIR}/pacMsplitter | echo $((`wc -l`)))
echo "Total number of reads $READNpacM !"

READNUNIT=$(((($READNpacM))/$NODES))
READSTART=1
READEND=$READNUNIT
	
for i in $(eval echo "{1..$NODES}")
	do
   	echo "Align chunk $i (out of $NODES) from read $READSTART to read $READEND !"
	
	ACTUALSTART=$(sed -n "$READSTART"p ${pacDIR}/pacMsplitter)
	ACTUALEND=$(sed -n "$READEND"p ${pacDIR}/pacMsplitter)
	
	if [ "$i" -eq "$NODES" ];
		then
		ACTUALEND=$(wc -l ${pacDIR}/m_pac.fasta | awk '{print $1}')
		ACTUALEND=$(($ACTUALEND+1))
		echo $ACTUALEND
	else
		echo "next.."
	fi
	sed -n "$ACTUALSTART,$(($ACTUALEND-1))"p ${pacDIR}/m_pac.fasta > ${pacDIR}/${i}_m_pac.fasta

	READSTART=$(($READSTART + $READNUNIT))
	READEND=$(($READEND + $READNUNIT))

done

echo "==================================== Done step 2A! ======================================="

cat > ${CWD}/qsubscripts/pacbins.bashX <<EOF
#!/bin/bash
#SBATCH -J redkmer2B
#SBATCH -t 80:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --tmp=500G
#SBATCH -e ${CWD}/reports/redkarray_%j_error.log
#SBATCH -o ${CWD}/reports/redkarray_%j_output.log
#SBATCH --array=1-${NODES}

source $SLURM_SUBMIT_DIR/redkmer_mod.cfg

module load Bowtie/1.2.2
module load GCC/6.2.0-2.27

	echo "==================================== Indexing chunk XXXXX{SLURM_ARRAY_TASK_ID} ======================================="
		cp ${pacDIR}/XXXXX{SLURM_ARRAY_TASK_ID}_m_pac.fasta XXXXXTMPDIR  # copy chunk to tmpdir
		$BOWTIEB XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_m_pac.fasta XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_m_pac

	#echo "==================================== make counting tool ======================================="	
		#cp ${BASEDIR}/Cscripts/* XXXXXTMPDIR
		#make

	echo "==================================== Working on male chunk XXXXX{SLURM_ARRAY_TASK_ID} ======================================="

		cp ${illDIR}/m.fastq.gz  XXXXXTMPDIR/m.fastq.gz   
		$BOWTIE -a -t -5 ${TRIMM5} -3 ${TRIMM3} -p $ARRAYCORES -v $VVAL XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_m_pac --suppress 1,2,4,5,6,7,8,9 XXXXXTMPDIR/m.fastq.gz 1> XXXXXTMPDIR/male.txt 2> $CWD/pacBio_illmapping/logs/XXXXX{SLURM_ARRAY_TASK_ID}_male_log.txt
		rm XXXXXTMPDIR/m.fastq.gz  # clean up

	echo "==================================== Counting, sorting for male chunck XXXXX{SLURM_ARRAY_TASK_ID} ===================================="

		${BASEDIR}/Cscripts/count XXXXXTMPDIR/male.txt > XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_male_uniq
		cp XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_male_uniq $CWD/pacBio_illmapping/mapping_rawdata/
		rm XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_male_uniq
	
	echo "==================================== Done male chunk XXXXX{SLURM_ARRAY_TASK_ID} ! ===================================="

	echo "==================================== Working on female chunk XXXXX{SLURM_ARRAY_TASK_ID} ======================================="
		# cp $illF XXXXXTMPDIR
		cp ${illDIR}/f.fastq.gz XXXXXTMPDIR/f.fastq.gz   # correcting so mt-filtered reads are used. # JRW: don't copy, use single source
		$BOWTIE -a -t -5 ${TRIMM5} -3 ${TRIMM3} -p $ARRAYCORES -v $VVAL XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_m_pac --suppress 1,2,4,5,6,7,8,9 XXXXXTMPDIR/f.fastq.gz  1> XXXXXTMPDIR/female.txt 2> $CWD/pacBio_illmapping/logs/XXXXX{SLURM_ARRAY_TASK_ID}_female_log.txt
		rm XXXXXTMPDIR/f.fastq.gz # clean up

	echo "==================================== Counting, sorting for male chunck XXXXX{SLURM_ARRAY_TASK_ID} ===================================="
		${BASEDIR}/Cscripts/count XXXXXTMPDIR/female.txt > XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_female_uniq
		cp XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_female_uniq $CWD/pacBio_illmapping/mapping_rawdata/
		rm XXXXXTMPDIR/XXXXX{SLURM_ARRAY_TASK_ID}_female_uniq
		rm -rf XXXXXTMPDIR 
	echo "==================================== Done female chunk XXXXX{SLURM_ARRAY_TASK_ID} ! ===================================="


echo "==================================== Done step 2B chunk XXXXX{SLURM_ARRAY_TASK_ID} ! ======================================="

EOF
sed 's/XXXXX/$/g' ${CWD}/qsubscripts/pacbins.bashX > ${CWD}/qsubscripts/pacbins.bash

sbatch ${CWD}/qsubscripts/pacbins.bash

#rm -rf $TMPDIR

echo "==================================== Done step 2B! ======================================="
