import json
import time
import logging

import requests

# from .. import service_account
from .. import wait_for_zone_op

_LOGGER = logging.getLogger(__name__)


def create_instance(credentials, project, zone, name,
                    startup_script=None, startup_script_url=None,
                    metadata=None,
                    machine_type='f1-micro', tags=None,
                    disk_size_gb=10, wait_until_done=False):
    """Create instance with startup script.

    TODO: docstring"""

    if startup_script is not None and startup_script_url is not None:
        raise ValueError('Cannot specify a startup script string and URL '
                         'at the same time!')

    access_token = credentials.get_access_token()

    if metadata is None:
        metadata = {}
    meta_items = [{'key': k, 'value': v} for k, v in metadata.items()]

    if tags is None:
        tags = []

    if startup_script is not None:
        meta_items.insert(
            0, {'key': 'startup-script', 'value': startup_script}
        )
    elif startup_script_url is not None:
        meta_items.insert(
            0, {'key': 'startup-script-url', 'value': startup_script_url})

    payload = {
        "name": name,
        "zone": "projects/%s/zones/%s" % (project, zone),
        "machineType": "projects/%s/zones/%s/machineTypes/%s"
                       % (project, zone, machine_type),
        "metadata": {
            "items": meta_items
        },
        "tags": {
            "items": tags
        },
        "disks": [
            {
                "type": "PERSISTENT",
                "boot": True,
                "mode": "READ_WRITE",
                "autoDelete": True,
                "deviceName": "test-instance",
                "initializeParams": {
                    "sourceImage": "https://www.googleapis.com/compute/v1/projects/ubuntu-os-cloud/global/images/ubuntu-1604-xenial-v20161115",
                    "diskType": "projects/%s/zones/%s/diskTypes/pd-standard" % (project, zone),
                    "diskSizeGb": str(disk_size_gb)
                }
            }
        ],
        "canIpForward": False,
        "networkInterfaces": [
            {
                "network": "projects/%s/global/networks/default" % project,
                "subnetwork": "projects/%s/regions/%s/subnetworks/default" % (project, zone[:-2]),
                "accessConfigs": [
                    {
                        "name": "External NAT", "type": "ONE_TO_ONE_NAT"
                    }
                ]
            }
        ],
        "description": "",
        "scheduling": {
            "preemptible": False,
            "onHostMaintenance": "MIGRATE",
            "automaticRestart": True
        },
        "serviceAccounts": [
            {
                "email": "default",
                "scopes": [
                    'https://www.googleapis.com/auth/compute',
                    "https://www.googleapis.com/auth/devstorage.read_write",
                    "https://www.googleapis.com/auth/logging.write",
                    "https://www.googleapis.com/auth/monitoring.write",
                    "https://www.googleapis.com/auth/servicecontrol",
                    "https://www.googleapis.com/auth/service.management.readonly",
                    "https://www.googleapis.com/auth/trace.append"
                ]
            }
        ]
    }
    #header = 'Authorization: Bearer 1/fFBGRNJru1FQd44AzqT3Zg'
    headers = {
        'Authorization': 'Bearer %s' % access_token.access_token
    }
    
    _LOGGER.debug('Access token: %s' % access_token.access_token)

    r = requests.post('https://www.googleapis.com/compute/v1/'
                      'projects/%s/zones/%s/instances' % (project, zone),
                      headers=headers, json=payload)

    r.raise_for_status()

    op_name = r.json()['name']

    _LOGGER.info('Submitted request to create intsance '
                '(HTTP code: %d).',
                r.status_code)

    if wait_until_done:
        _LOGGER.info('Waiting until operation is done...')
        wait_for_zone_op(access_token, project, zone, op_name)

    return op_name


def delete_instance(credentials, project, zone, name, wait_until_done=False):
    """Delete an instance.

    TODO: docstring
    """

    access_token = credentials.get_access_token()

    headers = {
        'Authorization': 'Bearer %s' % access_token.access_token
    }

    r = requests.delete('https://www.googleapis.com/compute/v1/'
                        'projects/%s/zones/%s/instances/%s'
                        % (project, zone, name),
                        headers=headers)

    r.raise_for_status()

    op_name = r.json()['name']

    _LOGGER.info('Submitted request to create intsance '
                '(HTTP code: %d).',
                r.status_code)

    if wait_until_done:
        _LOGGER.info('Waiting until operation is done...')
        wait_for_zone_op(access_token, project, zone, op_name)

    return op_name


def wait_for_instance_deletion(credentials, project, zone, instance_name,
                               interval_seconds=5):
    """Wait until an instance is deleted.
    
    We require that initially, the specified instance exists.
    TODO: docstring
    """
    
    t0 = time.time()
    access_token = credentials.get_access_token()
    headers = {
        'Authorization': 'Bearer %s' % access_token.access_token
    }

    r = requests.get('https://www.googleapis.com/compute/v1/'
                     'projects/%s/zones/%s/instances/%s'
                     % (project, zone, instance_name),
                     headers=headers)
    
    if r.status_code == 404:
        raise AssertionError('Instance "%s" does not exist!' % instance_name)
        
    r.raise_for_status()
    _LOGGER.debug('Instance "%s" exists.', instance_name)

    while True:
        time.sleep(interval_seconds)
        r = requests.get('https://www.googleapis.com/compute/v1/'
                         'projects/%s/zones/%s/instances/%s'
                         % (project, zone, instance_name),
                         headers=headers)
        if r.status_code == 404:
            break
        r.raise_for_status()
        _LOGGER.debug('Instance "%s" still exists.', instance_name)
        
    t1 = time.time()
    t = t1-t0
    t_min = t/60.0
    _LOGGER.info('Instance was deleted after %.1f s (%.1f m).', t, t_min)
