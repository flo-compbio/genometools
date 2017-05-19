import logging

from google.oauth2 import service_account
from oauth2client.service_account import ServiceAccountCredentials

LOGGER = logging.getLogger(__name__)


def get_storage_credentials(key, read_only=False):
    """Authenticates a service account for reading and/or writing on a bucket.
    
    This uses the `google.oauth2.service_account` module to obtain "scoped
    credentials". These can be used with the `google.storage` module.

    TODO: docstring"""
    if read_only:
        scopes = ['https://www.googleapis.com/auth/devstorage.read_only']
    else:
        scopes = ['https://www.googleapis.com/auth/devstorage.read_write']

    credentials = service_account.Credentials.from_service_account_info(key)
    
    scoped_credentials = credentials.with_scopes(scopes)
    
    return scoped_credentials


def get_compute_credentials(key):
    """Authenticates a service account for the compute engine.

    This uses the `oauth2client.service_account` module. Since the `google`
    Python package does not support the compute engine (yet?), we need to make
    direct HTTP requests. For that we need authentication tokens. Obtaining
    these based on the credentials provided by the `google.auth2` module is
    much more cumbersome than using the `oauth2client` module.

    See:
    - https://cloud.google.com/iap/docs/authentication-howto
    - https://developers.google.com/identity/protocols/OAuth2ServiceAccount
    
    TODO: docstring"""
    scopes = ['https://www.googleapis.com/auth/compute']

    credentials = ServiceAccountCredentials.from_json_keyfile_dict(
        key, scopes=scopes)
    
    return credentials
