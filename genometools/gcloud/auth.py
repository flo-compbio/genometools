"""Authentication functions."""

#from google.cloud import storage
#from google.cloud.exceptions import Conflict

# for authenticating for storage
from google.oauth2 import service_account


# for authenticating for compute
from oauth2client.service_account import ServiceAccountCredentials


def authenticate_storage(client_secrets, read_only=False):
    """Authenticates a service account for reading and/or writing on a bucket.
    
    TODO: docstring"""
    if read_only:
        scopes = ['https://www.googleapis.com/auth/devstorage.read_only']
    else:
        scopes = ['https://www.googleapis.com/auth/devstorage.read_write']

    credentials = service_account.Credentials.from_service_account_info(
            client_secrets)
    
    scoped_credentials = credentials.with_scopes(scopes)
    
    return scoped_credentials


def authenticate_compute(client_secrets):
    """Authenticates a service account for the compute engine.
    
    TODO: docstring"""
    scopes = ['https://www.googleapis.com/auth/compute']

    credentials = ServiceAccountCredentials.from_json_keyfile_dict(
            client_secrets, scopes=scopes)
    
    return credentials