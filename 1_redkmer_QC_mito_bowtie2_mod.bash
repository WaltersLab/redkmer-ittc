#!/bin/bash
#SBATCH -J redkmer1
#SBATCH -t 70:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --tmp=200G
#SBATCH -e redkmer1_%j_error.log
#SBATCH -o redkmer1_%j_output.log

echo "========== starting up step 1 =========="

source $SLURM_SUBMIT_DIR/redkmer_mod.cfg
module load Bowtie2/2.2.9
module load SAMtools
# module load slurm-torque

echo "========== setting up directories =========="

mkdir -p $CWD/qsubscripts
mkdir -p $CWD/QualityReports
mkdir -p $CWD/plots
mkdir -p $CWD/MitoIndex
mkdir -p $CWD/reports

echo "========== filtering pacBio libary by read length =========="

# cp ${pacDIR}/raw_pac.fasta $TMPDIR
# cp ${pacM} ${TMPDIR}/raw_pac.fasta
pigz -p $CORES -dc ${pacM} > ${TMPDIR}/raw_pac.fasta
$SAMTOOLS faidx $TMPDIR/raw_pac.fasta
awk -v pl="$pac_length" -v plm="$pac_length_max" '{if($2>=pl && $2<=plm)print $1}' $TMPDIR/raw_pac.fasta.fai | xargs samtools faidx $TMPDIR/raw_pac.fasta > $TMPDIR/m_pac.fasta
cp $TMPDIR/m_pac.fasta ${pacDIR}/m_pac.fasta
rm -rf $TMPDIR   # clean up here. don't need tmpdir sub-scripts below, as edited, because reading/writing in place to panassis

echo "========== building mitochondiral index =========="

$BOWTIE2B --threads $CORES $MtREF ${CWD}/MitoIndex/MtRef_bowtie2

echo "========== filtering for mitochondiral reads =========="

cat > ${CWD}/qsubscripts/femalemito.bashX <<EOF
#!/bin/bash
#SBATCH -J redk_f_mito
#SBATCH -t 20:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH -e ${CWD}/reports/redk_f_mito_%j_error.log
#SBATCH -o ${CWD}/reports/redk_f_mito_%j_output.log


module load Bowtie2/2.2.9
#module load FastQC

# cp ${illDIR}/raw_f.fastq XXXXX/raw_f.fastq
# cp $illF XXXXX/raw_f.fastq  #JRW: no need to copy this to temp, just read directly from original
echo "========== NOT producing quality report for female illumina library =========="
# $FASTQC XXXXX/raw_f.fastq -o ${CWD}/QualityReports # JRW: Skipping fastqc. who cares...?
echo "========== removing female illumina reads mapping to mitochondrial DNA =========="
#JRW: align directly from original gzipped fastq, and write directly to a gzipped file.
$BOWTIE2 -p $CORES -x ${CWD}/MitoIndex/MtRef_bowtie2 -U $illF --un-gz ${illDIR}/f.fastq.gz 1>/dev/null 2> ${illDIR}/f_bowtie2.log
#cp XXXXX/f.fastq ${illDIR}

EOF
sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/femalemito.bashX > ${CWD}/qsubscripts/femalemito.bash
sbatch  ${CWD}/qsubscripts/femalemito.bash

cat > ${CWD}/qsubscripts/malemito.bashX <<EOF
#!/bin/bash
#SBATCH -J redk_m_mito
#SBATCH -t 20:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH -e ${CWD}/reports/redk_m_mito_%j_error.log
#SBATCH -o ${CWD}/reports/redk_m_mito_%j_output.log

module load Bowtie2/2.2.9
#module load FastQC

# cp ${illDIR}/raw_m.fastq XXXXX/raw_m.fastq
# cp $illM XXXXX/raw_m.fastq
echo "========== NOT producing quality report for male illumina library =========="
# $FASTQC XXXXX/raw_m.fastq -o ${CWD}/QualityReports
echo "========== removing male illumina reads mapping to mitochondrial DNA =========="
$BOWTIE2 -p $CORES -x $CWD/MitoIndex/MtRef_bowtie2 -U $illM --un-gz ${illDIR}/m.fastq.gz 1>/dev/null 2> ${illDIR}/m_bowtie2.log
# cp XXXXX/m.fastq ${illDIR}
#rm -rf XXXXX 
EOF
sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/malemito.bashX > ${CWD}/qsubscripts/malemito.bash
sbatch ${CWD}/qsubscripts/malemito.bash

printf "======= Done step 1 =======\n"

