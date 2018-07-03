#!/bin/bash
#SBATCH -J redkmer1
#SBATCH -t 70:00:00
#SBATCH -c 20
#SBATCH --mem=32G
#SBATCH -e redkmer1_%j_error.log
#SBATCH -o redkmer1_%j_output.log

echo "========== starting up step 1 =========="

source $SLURM_SUBMIT_DIR/redkmer.cfg
module load Bowtie2/2.2.9
module load SAMtools
module load slurm-torque

echo "========== setting up directories =========="

mkdir -p $CWD/qsubscripts
mkdir -p $CWD/QualityReports
mkdir -p $CWD/plots
mkdir -p $CWD/MitoIndex
mkdir -p $CWD/reports

echo "========== filtering pacBio libary by read length =========="

# cp ${pacDIR}/raw_pac.fasta $TMPDIR
cp ${pacM} ${TMPDIR}/raw_pac.fasta
$SAMTOOLS faidx $TMPDIR/raw_pac.fasta
awk -v pl="$pac_length" -v plm="$pac_length_max" '{if($2>=pl && $2<=plm)print $1}' $TMPDIR/raw_pac.fasta.fai | xargs samtools faidx $TMPDIR/raw_pac.fasta > $TMPDIR/m_pac.fasta
cp $TMPDIR/m_pac.fasta ${pacDIR}/m_pac.fasta

echo "========== building mitochondiral index =========="

$BOWTIE2B --threads $CORES $MtREF ${CWD}/MitoIndex/MtRef_bowtie2

echo "========== filtering for mitochondiral reads =========="

cat > ${CWD}/qsubscripts/femalemito.bashX <<EOF
#!/bin/bash
#SBATCH -J redkmer_f_mito
#SBATCH -t 20:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH -e ${CWD}/reports/%j_error.log
#SBATCH -o ${CWD}/reports/%j_output.log

module load Bowtie2/2.2.9
module load FastQC

# cp ${illDIR}/raw_f.fastq XXXXX/raw_f.fastq
cp $illF XXXXX/raw_f.fastq
echo "========== producing quality report for female illumina library =========="
$FASTQC XXXXX/raw_f.fastq -o ${CWD}/QualityReports
echo "========== removing female illumina reads mapping to mitochondrial DNA =========="
$BOWTIE2 -p $CORES -x ${CWD}/MitoIndex/MtRef_bowtie2 -U XXXXX/raw_f.fastq --un XXXXX/f.fastq 1>/dev/null 2> ${illDIR}/f_bowtie2.log
cp XXXXX/f.fastq ${illDIR}
EOF
sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/femalemito.bashX > ${CWD}/qsubscripts/femalemito.bash
qsub ${CWD}/qsubscripts/femalemito.bash

cat > ${CWD}/qsubscripts/malemito.bashX <<EOF
#!/bin/bash
#SBATCH -J redkmer_m_mito
#SBATCH -t 20:00:00
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH -e ${CWD}/reports/%j_error.log
#SBATCH -o ${CWD}/reports/%j_output.log

module load Bowtie2/2.2.9
module load FastQC

$cp ${illDIR}/raw_m.fastq XXXXX/raw_m.fastq
cp $illM XXXXX/raw_m.fastq
echo "========== producing quality report for male illumina library =========="
$FASTQC XXXXX/raw_m.fastq -o ${CWD}/QualityReports
echo "========== removing male illumina reads mapping to mitochondrial DNA =========="
$BOWTIE2 -p $CORES -x $CWD/MitoIndex/MtRef_bowtie2 -U XXXXX/raw_m.fastq --un XXXXX/m.fastq 1>/dev/null 2> ${illDIR}/m_bowtie2.log
cp XXXXX/m.fastq ${illDIR}
#rm -rf XXXXX 
EOF
sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/malemito.bashX > ${CWD}/qsubscripts/malemito.bash
qsub ${CWD}/qsubscripts/malemito.bash

printf "======= Done step 1 =======\n"

