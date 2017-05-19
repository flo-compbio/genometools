#!/bin/bash

SCRIPT_DIR='{{ script_dir }}'
BARCODE_READ_FILE='{{ barcode_read_file }}'
MRNA_READ_FILE='{{ mrna_read_file }}'
BARCODE1_FILE='{{ barcode1_file }}'
BARCODE2_FILE='{{ barcode2_file }}'
OUTPUT_DIR='{{ output_dir }}'
MAX_READS={{ max_reads }}

# store everything in /root directory
cd /root

# download startup scripts
echo "Downloading startup scripts from Google Storage..."
mkdir scripts
gsutil -q cp -r "${SCRIPT_DIR}/*" scripts/
chmod -R a+x scripts/*

# install compiled crcmod module
echo "Installing compiled version of crcmod module..."
./scripts/install_crcmod.sh

# download input data
echo "Downloading input data from Google Storage..."
mkdir data
gsutil -q cp "$BARCODE1_FILE" data/barcode1_list.txt
gsutil -q cp "$BARCODE2_FILE" data/barcode2_list.txt
gsutil -q cp "$BARCODE_READ_FILE" data/barcode_reads.fastq.gz
gsutil -q cp "$MRNA_READ_FILE" data/mrna_reads.fastq.gz

# install GCC
./scripts/install_gcc.sh

# install Python 3
./scripts/install_python.sh

# install numpy
pip3 install numpy

# install cython
pip3 install cython

# install pandas
pip3 install pandas

# generate C file with cython
mkdir build
echo "Cythoning the barcode processing script..."
cython --embed -o build/process_barcodes.c \
        scripts/inDrop/process_barcodes.pyx

# compile and link C file to generate executable
echo "Compiling and linking the barcode processing script..."
PYTHON_INCLUDES=`python3-config --includes`
PYTHON_LIB_DIRS=`python3-config --ldflags`
PYTHON_LIBS=`python3-config --libs`
NUMPY_INCLUDE=`python3 -c "import numpy as np; print(np.get_include())"`
#gcc script.c $PYTHON_INCLUDES $PYTHON_LIB_DIRS $PYTHON_LIBS -o script
gcc build/process_barcodes.c $PYTHON_INCLUDES -I${NUMPY_INCLUDE} \
        $PYTHON_LIB_DIRS $PYTHON_LIBS \
        -o build/process_barcodes

# run the executable
echo "Running the barcode processing script..."
mkdir output
./build/process_barcodes \
        -rb "data/barcode_reads.fastq.gz" -rm "data/mrna_reads.fastq.gz" \
        -b1 "data/barcode1_list.txt" -b2 "data/barcode2_list.txt" \
        -o "output" \
        --max-reads $MAX_READS

# copy the results to the bucket
echo "Uploading the results to Google Storage..."
gsutil -q cp -r "output/*" "${OUTPUT_DIR}/"

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
chmod a+x delete_instance.sh
./delete_instance.sh
{% endif %}
