#!/bin/bash
# create a STAR index

# variables
SCRIPT_DIR='{{ script_dir }}'
GENOME_FILE='{{ genome_file }}'
ANNOTATION_FILE='{{ annotation_file if annotation_file is not none else '' }}'
SPLICE_OVERHANG={{ splice_overhang }}
OUTPUT_DIR='{{ output_dir }}'
NUM_THREADS={{ num_threads }}
GENOME_MEMORY_LIMIT={{ genome_memory_limit }}
CHROMOSOME_BIN_BITS={{ chromosome_bin_bits }}

# store everything in /root directory
cd /root

# copy startup scripts
echo "Downloading startup scripts..."
gsutil -q cp -r "${SCRIPT_DIR}/*" ./

# install docker
chmod a+x install_docker.sh
./install_docker.sh

# download the genome file (.fa.gz)
echo "Downloading the genome file..."
gsutil -q cp "$GENOME_FILE" ./

# decompress the genome file
echo "Decompressing the genome file..."
GENOME_FILE_NAME=${GENOME_FILE##*/}
DECOMPRESSED_GENOME_FILE='genome.fa'
zcat "$GENOME_FILE_NAME" > "$DECOMPRESSED_GENOME_FILE"

{% if annotation_file is not none %}
# download the annotation file (.gtf.gz)
echo "Downloading the annotation file..."
gsutil -q cp "$ANNOTATION_FILE" ./

# decompress the annotation file
echo "Decompressing the annotation file..."
ANNOTATION_FILE_NAME=${ANNOTATION_FILE##*/}
DECOMPRESSED_ANNOTATION_FILE='annotations.gtf'
zcat "$ANNOTATION_FILE_NAME" > "$DECOMPRESSED_ANNOTATION_FILE"
{% endif %}

echo "Generating the STAR index..."
mkdir index
docker run -v "/root/:/data/" quay.io/biocontainers/star:2.5.3a--0 STAR \
    --runThreadN $NUM_THREADS \
    --limitGenomeGenerateRAM $GENOME_MEMORY_LIMIT \
    --runMode genomeGenerate \
    --genomeDir "/data/index" \
    --genomeFastaFiles "/data/${DECOMPRESSED_GENOME_FILE}" \
    --genomeChrBinNbits $CHROMOSOME_BIN_BITS \
    {% if annotation_file is not none %} \
        --sjdbGTFfile "/data/${DECOMPRESSED_ANNOTATION_FILE}" \
        --sjdbOverhang $SPLICE_OVERHANG \
    {% endif %}

echo "Copying the index to the bucket..."
gsutil -q cp -r index/* "${OUTPUT_DIR}/"

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
chmod a+x delete_instance.sh
./delete_instance.sh
{% endif %}
