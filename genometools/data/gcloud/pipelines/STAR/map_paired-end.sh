#!/bin/bash
# align paired-end reads using STAR

# variables
SCRIPT_DIR='{{ script_dir }}'
INDEX_DIR='{{ index_dir }}'
FASTQ_FILE1='{{ fastq_file1 }}'
FASTQ_FILE2='{{ fastq_file2 }}'
OUTPUT_DIR='{{ output_dir }}'
NUM_THREADS='{{ num_threads }}'
MAX_MATE_DIST='{{ max_mate_dist }}'
CHIM_SEGMENT_MIN='{{ chim_segment_min }}'
SEED_START_LMAX='{{ seed_start_lmax }}'
SCORE_MIN_ABS='{{ score_min_abs }}'
SCORE_MIN_REL='{{ score_min_rel }}'
MATCH_NMIN_ABS='{{ match_nmin_abs }}'
MATCH_NMIN_REL='{{ match_nmin_rel }}'

# store everything in /root directory
cd /root

# copy startup scripts
echo "Downloading startup scripts..."
gsutil -q cp -r "${SCRIPT_DIR}/*" ./

# install docker
echo "Installing docker"
chmod a+x install_docker.sh
./install_docker.sh

# copy index
echo "Downloading the STAR index"
gsutil -q cp -r "${INDEX_DIR}" ./
INDEX_NAME=${INDEX_DIR##*/}

# copy FASTQ file (.fa.gz)
echo "Downloading the reads (.fastq.gz)"
gsutil -q cp "${FASTQ_FILE1}" "${FASTQ_FILE2}" ./
FASTQ_FILE_NAME1=${FASTQ_FILE1##*/}
FASTQ_FILE_NAME2=${FASTQ_FILE2##*/}

# generate output folder
mkdir output

# run STAR
echo "Performing alignment with STAR..."
docker run -v "/root/:/data/" quay.io/biocontainers/star:2.5.3a--0 STAR \
    --runThreadN $NUM_THREADS \
    --genomeDir "/data/${INDEX_NAME}" \
    --readFilesIn "/data/${FASTQ_FILE_NAME1}" "/data/${FASTQ_FILE_NAME2}" \
    --readFilesCommand "zcat" \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "/data/output/" \
    --outSAMstrandField intronMotif \
    --outSAMmultNmax 1 \
    --alignMatesGapMax $MAX_MATE_DIST \
    --chimSegmentMin $CHIM_SEGMENT_MIN \
    --seedSearchStartLmax $SEED_START_LMAX \
    --outFilterScoreMin $SCORE_MIN_ABS \
    --outFilterScoreMinOverLread $SCORE_MIN_REL \
    --outFilterMatchNmin $MATCH_NMIN_ABS \
    --outFilterMatchNminOverLread $MATCH_NMIN_REL

#--alignEndsType Extend5pOfRead1 \

echo "Uploading the output files to the bucket..."
gsutil -q cp output/* "${OUTPUT_DIR}/"

# delete instance
echo "Deleting instance..."
chmod a+x delete_instance.sh
./delete_instance.sh
