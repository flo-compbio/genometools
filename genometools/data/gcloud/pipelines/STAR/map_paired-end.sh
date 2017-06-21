#!/bin/bash
# align paired-end reads using STAR

set -x  # echo on

# store everything in /root directory
cd /root

# copy startup scripts
mkdir scripts
echo "Downloading startup scripts..."
SCRIPT_DIR='{{ script_dir }}'
gsutil -q cp -r "${SCRIPT_DIR}/*" scripts/
chmod -R a+x scripts/*

# install compiled crcmod module
echo "Installing compiled version of crcmod module..."
./scripts/install_crcmod.sh

# install docker
echo "Installing docker"
./scripts/install_docker.sh

# copy index
mkdir index
echo "Downloading the STAR index"
INDEX_DIR='{{ index_dir }}'
gsutil -q cp -r "${INDEX_DIR}" ./index/
INDEX_NAME=${INDEX_DIR##*/}

# create data directory
mkdir reads

# create output directory
mkdir output

# process all pairs, one at a time
{% for r1_fastq_file, r2_fastq_file in fastq_file_pairs %}

R1_FASTQ_FILE='{{ r1_fastq_file }}'
R1_FASTQ_FILE_NAME=${R1_FASTQ_FILE##*/}
R2_FASTQ_FILE='{{ r2_fastq_file }}'
R2_FASTQ_FILE_NAME=${R2_FASTQ_FILE##*/}
FASTQ_NAME=${R1_FASTQ_FILE_NAME%_*}

echo "Now processing '${FASTQ_NAME}'..."
mkdir output/${FASTQ_NAME}

echo "Downloading the reads (FASTQ files)..."
gsutil -q cp "$R1_FASTQ_FILE" "$R2_FASTQ_FILE" ./reads/

# run STAR
echo "Performing alignment with STAR..."
docker run -v "/root:/host" quay.io/biocontainers/star:2.5.3a--0 STAR \
    --genomeDir "/host/index/${INDEX_NAME}" \
    --readFilesIn "/host/reads/${R1_FASTQ_FILE_NAME}" "/host/reads/${R2_FASTQ_FILE_NAME}" \
    --outFileNamePrefix "/host/output/${FASTQ_NAME}/" \
    {% if num_threads is defined and num_threads is not none 
        %}--runThreadN {{ num_threads }} \
    {% endif
    %}{% if compressed is defined and compressed
        %}--readFilesCommand "zcat" \
    {% endif
    %}{% if out_sam_type is defined and out_sam_type is not none
        %}--outSAMtype {{ out_sam_type }} \
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
    %}{% if keep_unmapped is defined and keep_unmapped
        %}--outSAMunmapped Within \
    {% endif
    %}{% if max_intron_size is defined and max_intron_size is not none
        %}--alignIntronMax {{ max_intron_size }} \
    {% endif
    %}{% if max_mate_dist is defined and max_mate_dist is not none
        %}--alignMatesGapMax {{ max_mate_dist }} \
    {% endif %}

    {#

    # --alignIntronMax = (2ˆwinBinNbits)*winAnchorDistNbins
    # --alignMatesGapMax = (2ˆwinBinNbits)*winAnchorDistNbins
    # --winBinNbits 16
    # --winAnchorDistNbins 9
    # => 589815

    #--outSAMstrandField intronMotif \
    #--outSAMmultNmax 1 \
    #--alignMatesGapMax $MAX_MATE_DIST \
    #--chimSegmentMin $CHIM_SEGMENT_MIN \
    #--seedSearchStartLmax $SEED_START_LMAX \
    #--outFilterScoreMin $SCORE_MIN_ABS \
    #--outFilterScoreMinOverLread $SCORE_MIN_REL \
    #--outFilterMatchNmin $MATCH_NMIN_ABS \
    #--outFilterMatchNminOverLread $MATCH_NMIN_REL
    #--alignEndsType Extend5pOfRead1 \

    #}

echo "Uploading the output files to the bucket..."
OUTPUT_DIR='{{ output_dir }}'
gsutil -q cp output/${FASTQ_NAME}/* "${OUTPUT_DIR}/${FASTQ_NAME}/"

{% endfor %}

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
./scripts/delete_instance.sh
{% endif %}
