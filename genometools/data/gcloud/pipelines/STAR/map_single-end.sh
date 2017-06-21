#!/bin/bash
# align single-end reads with STAR

set -x  # echo on

# variables
SCRIPT_DIR='{{ script_dir }}'
INDEX_DIR='{{ index_dir }}'
OUTPUT_DIR='{{ output_dir }}'
#NUM_THREADS='{{ num_threads }}'
#SEED_START_LMAX='{{ seed_start_lmax }}'
#MISMATCH_NMAX='{{ mismatch_nmax }}'
#MULTIMAP_NMAX='{{ multimap_nmax }}'
#OUT_MULT_NMAX='{{ out_mult_nmax }}'
#SPLICE_MIN_OVERHANG='{{ splice_min_overhang }}'
#OUT_SAM_TYPE='{{ out_sam_type }}'

# store everything in /root directory
cd /root

# create data directory
mkdir data

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
gsutil -q cp -r "${INDEX_DIR}" ./data/
INDEX_NAME=${INDEX_DIR##*/}

# copy all FASTQ files (.fa.gz)
{% for fastq_file in fastq_files %}

FASTQ_FILE='{{ fastq_file }}'
FASTQ_FILE_NAME=${FASTQ_FILE##*/}
echo "Downloading the reads for '${FASTQ_FILE_NAME}' (.fastq.gz)"
gsutil -q cp "${FASTQ_FILE}" ./data/

{% endfor %}

# generate output folder
mkdir output

### perform alignment for all FASTQ files
{% for fastq_file in fastq_files %}

FASTQ_FILE='{{ fastq_file }}'
FASTQ_FILE_NAME="${FASTQ_FILE##*/}"
FASTQ_NAME="${FASTQ_FILE_NAME%%.*}"

echo "Performing alignment of '${FASTQ_FILE_NAME}' with STAR..."

# create output directory for this FASTQ file
mkdir "output/${FASTQ_NAME}"

docker run -v "/root/:/host/" quay.io/biocontainers/star:2.5.3a--0 STAR \
    --genomeDir "/host/data/${INDEX_NAME}" \
    --readFilesIn "/host/data/${FASTQ_FILE_NAME}" \
    --outFileNamePrefix "/host/output/${FASTQ_NAME}/" \
    --outMultimapperOrder "Random" \
    {% if num_threads is defined and num_threads is not none 
        %}--runThreadN {{ num_threads }} \
    {% endif
    %}{% if compressed is defined and compressed
        %}--readFilesCommand "zcat" \
    {% endif
    %}{% if seed_start_lmax is defined and seed_start_lmax is not none
        %}--seedSearchStartLmax {{ seed_start_lmax }} \
    {% endif
    %}{% if mismatch_nmax is defined and mismatch_nmax is not none
        %}--outFilterMismatchNmax {{ mismatch_nmax }} \
    {% endif
    %}{% if multimap_nmax is defined and multimap_nmax is not none
        %}--outFilterMultimapNmax {{ multimap_nmax }} \
    {% endif
    %}{% if out_mult_nmax is defined and out_mult_nmax is not none
        %}--outSAMmultNmax {{ out_mult_nmax }} \
    {% endif
    %}{% if splice_min_overhang is defined and splice_min_overhang is not none
        %}--alignSJDBoverhangMin {{ splice_min_overhang }} \
    {% endif
    %}{% if out_sam_type is defined and out_sam_type is not none
        %}--outSAMtype {{ out_sam_type }} \
    {% endif
    %}{% if keep_unmapped is defined and keep_unmapped
        %}--outSAMunmapped Within \
    {% endif %}

{# --outSAMstrandField intronMotif \ #}

{% endfor %}

echo "Uploading the output files to the bucket..."
gsutil -q cp -r output/* "${OUTPUT_DIR}/"

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
./scripts/delete_instance.sh
{% endif %}
