#!/bin/bash
# deletes the instance

# get instance name
INSTANCE_NAME=`curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google"`

# get instance zone
zone=`curl -s "http://metadata.google.internal/computeMetadata/v1/instance/zone" -H "Metadata-Flavor: Google"`
arrzone=(${zone//// })
INSTANCE_ZONE=${arrzone[-1]}

# echo "Self-deletion..."
gcloud compute instances delete "$INSTANCE_NAME" --zone="$INSTANCE_ZONE" -q