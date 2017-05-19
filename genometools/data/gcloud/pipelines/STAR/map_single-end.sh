#!/bin/bash
# align single-end reads with STAR

# variables
SCRIPT_DIR='{{ script_dir }}'
INDEX_DIR='{{ index_dir }}'
FASTQ_FILE='{{ fastq_file }}'
OUTPUT_DIR='{{ output_dir }}'
NUM_THREADS='{{ num_threads }}'
SEED_START_LMAX='{{ seed_start_lmax }}'
MISMATCH_NMAX='{{ mismatch_nmax }}'
MULTIMAP_NMAX='{{ multimap_nmax }}'
OUT_MULT_NMAX='{{ out_mult_nmax }}'
SPLICE_MIN_OVERHANG='{{ splice_min_overhang }}'
OUT_SAM_TYPE='{{ out_sam_type }}'

# store everything in /root directory
cd /root

# copy startup scripts
mkdir scripts
echo "Downloading startup scripts..."
gsutil -q cp -r "${SCRIPT_DIR}/*" scripts/
chmod -R a+x scripts/*

# install compiled crcmod module
echo "Installing compiled version of crcmod module..."
./scripts/install_crcmod.sh

# install docker
echo "Installing docker"
./scripts/install_docker.sh

# copy index
echo "Downloading the STAR index"
gsutil -q cp -r "${INDEX_DIR}" ./
INDEX_NAME=${INDEX_DIR##*/}

# copy FASTQ file (.fa.gz)
echo "Downloading the reads (.fastq.gz)"
gsutil -q cp "${FASTQ_FILE}" ./
FASTQ_FILE_NAME=${FASTQ_FILE##*/}

# generate output folder
mkdir output

# run STAR
echo "Performing alignment with STAR..."
docker run -v "/root/:/data/" quay.io/biocontainers/star:2.5.3a--0 STAR \
    --runThreadN $NUM_THREADS \
    --genomeDir "/data/${INDEX_NAME}" \
    --readFilesIn "/data/${FASTQ_FILE_NAME}" \
    {% if compressed %} --readFilesCommand "zcat" {% endif %} \
    --outFileNamePrefix "/data/output/" \
    --outSAMstrandField intronMotif \
    --outMultimapperOrder "Random" \
    --seedSearchStartLmax $SEED_START_LMAX \
    --outFilterMismatchNmax $MISMATCH_NMAX \
    --outFilterMultimapNmax $MULTIMAP_NMAX \
    --outSAMmultNmax $OUT_MULT_NMAX \
    --alignSJDBoverhangMin $SPLICE_MIN_OVERHANG \
    --outSAMtype $OUT_SAM_TYPE

echo "Uploading the output files to the bucket..."
gsutil -q cp output/* "${OUTPUT_DIR}/"

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
./scripts/delete_instance.sh
{% endif %}
