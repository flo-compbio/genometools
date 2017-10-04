#!/bin/bash

# install the compiled version of the crcmod Python package
# this dramatically increases performance of gsutil operations

apt-get update
apt-get install -y python-pip python-dev build-essential
pip install --upgrade pip

/usr/bin/yes | pip uninstall crcmod
pip install --no-cache crcmod