set -euo pipefail

REPO="/mnt/c/Users/현종/OneDrive - 인하대학교/문서/Qoka/IBS6503-2026-02"
DATA="$REPO/data/week4"
BIO="$HOME/qoka_envs/bio/bin"
export PATH="$BIO:$PATH"

hisat2 --version | head -1
samtools --version | head -1

set -x

echo "### 1. 결과 폴더 준비 ###"
mkdir -p index bam fastq

echo "### 1b. FASTQ를 공백 없는 로컬 경로로 복사 (hisat2 래퍼가 경로 내 공백을 잘못 파싱하는 문제 회피) ###"
for sample in 1_control 1_treated 2_control 2_treated; do
  cp "$DATA/${sample}.R1.fastq.gz" "fastq/${sample}.R1.fastq.gz"
  cp "$DATA/${sample}.R2.fastq.gz" "fastq/${sample}.R2.fastq.gz"
done

echo "### 2. 참조 서열 압축 해제 ###"
gunzip -k -c "$DATA/hg19_chr1_11.4Mb.fa.gz" > index/ref.fa

echo "### 3. hisat2 인덱스 생성 ###"
hisat2-build index/ref.fa index/ref < /dev/null

echo "### 4. 네 샘플 정렬 -> 위치순 정렬(sort) -> 인덱싱(index) ###"
for sample in 1_control 1_treated 2_control 2_treated; do
  hisat2 -p "$(nproc)" -x index/ref \
    -1 "fastq/${sample}.R1.fastq.gz" \
    -2 "fastq/${sample}.R2.fastq.gz" \
    -S "bam/${sample}.sam" < /dev/null

  samtools sort -@ "$(nproc)" -o "bam/${sample}.sorted.bam" "bam/${sample}.sam"
  samtools index "bam/${sample}.sorted.bam"
  rm "bam/${sample}.sam"
done

rm -rf fastq

set +x
echo "### 결과 파일 ###"
ls -la bam/
echo "### flagstat ###"
for sample in 1_control 1_treated 2_control 2_treated; do
  echo "--- ${sample} ---"
  samtools flagstat "bam/${sample}.sorted.bam"
done
