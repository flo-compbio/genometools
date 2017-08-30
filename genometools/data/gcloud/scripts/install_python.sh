#!/bin/bash

apt-get update
apt-get install -y python3-pip python3-dev build-essential
pip3 install --upgrade pip  # required to enable pip cache (under /root/.cache/pip/http)
pip3 install --upgrade pyasn1

# install gcloud, also installs oauth2client
#pip3 install google-cloud
pip3 install google-cloud-storage oauth2client