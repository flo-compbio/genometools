# using google.auth, see https://google-auth.readthedocs.io/en/latest/

"""
Notes on authentication and permissions:

We use a service account to access Google Cloud Storage.
A service account is an account that belongs to an application or a VM,
instead of an individual. For more details, see here:
https://cloud.google.com/iam/docs/service-accounts

There are two ways of granting permissions to the service account for
performing operations on the bucket. (I prefer the second way.)

The first way uses roles - you can add these roles to the service
account using the Cloud Console. Under "IAM & Admin"/IAM, click on the
"Role(s)" drop-down menu, and select the roles from the "Storage"
submenu. The "Storage Object Viewer" role allows the service account
to download files, and the "Storage Object Creator" role allows the
service account to upload files. However, these roles do not allow the
service account to list files, or to delete or overwrite exisitng files.
Surprisingly, this requires the "Storage Admin" role (which seems overly
general --- the "Storage Object Admin" role does not work).

I therefore prefer to instead grant permissions to the service account
by making the account a READER (for read-only access) or a WRITER (for
read and write access) on the bucket. You can do this using the Cloud
Console. Under "Storage"/Browser, click on the three dots to the right
of the bucket name, and add a row with
  "ENTITY"=User,
  "NAME"=(email address of the service account), and
  "ACCESS"=Reader (or Writer).

Finally, you have to create a JSON key for the storage account (if you
haven't already) You can download a JSON key file using the Cloud Console.
Under "IAM & Admin"/Service accounts, click on the three dots to the right of
the service account name, and select "Create key".
"""
import os
import logging
import json

#from oauth2client.service_account import ServiceAccountCredentials
from google.cloud import storage
from google.cloud.exceptions import Conflict
# from google.oauth2 import service_account

from . import get_storage_credentials

#from genometools import misc

#GCLOUD_PROJECT = 'bo-project'
#GCLOUD_BUCKET = 'bo-scripts'

LOGGER = logging.getLogger(__name__)

def get_client(key, project):
    """Gets a `Client` object (required by the other functions).
    
    TODO: docstring"""
    cred = get_storage_credentials(key)
    return storage.Client(project=project, credentials=cred)


def get_files(client, bucket, prefix=''):
    """Lists files/objects on a bucket.
    
    TODO: docstring"""
    bucket = client.get_bucket(bucket)
    files = list(bucket.list_blobs(prefix=prefix))    
    return files


def download_file(client, bucket, remote_path, local_path, overwrite=False):
    """Downloads a file from a bucket.
    
    TODO: docstring"""
    bucket = client.get_bucket(bucket)
    if (not overwrite) and os.path.isfile(local_path):
        raise OSError('File already exists!')
    with open(local_path, 'wb') as ofh:
        bucket.get_blob(remote_path).download_to_file(ofh)

        
def delete_file(client, bucket, remote_path):
    """Deletes a file from a bucket.
    
    TODO: docstring"""
    bucket = client.get_bucket(bucket)
    bucket.delete_blob(remote_path)


def upload_file(client, bucket, local_path, remote_path, overwrite=False):
    """Uploads a file to a bucket.
    
    TODO: docstring"""
    bucket = client.get_bucket(bucket)
    blob = storage.Blob(remote_path, bucket)
    if (not overwrite) and blob.exists():
        raise Conflict('File/object already exists on the bucket!')
    blob.upload_from_filename(local_path)
