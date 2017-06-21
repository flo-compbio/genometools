"""Operation functions."""

import logging
import json
import time

import requests

LOGGER = logging.getLogger(__name__)


def wait_for_zone_op(access_token, project, zone, name, interval=1.0):
    """Wait until a zone operation is finished.

    TODO: docstring"""

    assert isinstance(interval, (int, float))
    assert interval >= 0.1
    
    status = 'RUNNING'
    progress = 0

    LOGGER.info('Waiting for zone operation "%s" to finish...', name)

    while status != 'DONE':

        # print('\rRunning! Progress: %d %%' % progress, end='')
        r = requests.get('https://www.googleapis.com/compute/v1/'
                         'projects/%s/zones/%s/operations/%s?access_token=%s'
                         % (project, zone, name, access_token.access_token))
        r.raise_for_status()

        # print(r.status_code)
        #result = json.loads(r.content.decode('utf-8'))
        result = r.json()
        status = result['status']
        progress = result['progress']

        time.sleep(interval)

    LOGGER.info('Zone operation "%s" done!', name)


def wait_for_global_op(access_token, project, name, interval=1.0):
    """Wait until a global operation is finished.

    TODO: docstring"""
    
    assert isinstance(interval, (int, float))
    assert interval >= 0.1
    
    status = 'RUNNING'
    progress = 0

    LOGGER.info('Waiting for global operation "%s" to finish...', name)

    while status != 'DONE':

        # print('\rRunning! Progress: %d %%' % progress, end='')
        r = requests.get('https://www.googleapis.com/compute/v1/'
                         'projects/%s/global/operations/%s?access_token=%s'
                         % (project, name, access_token.access_token))
        r.raise_for_status()

        # print(r.status_code)
        #result = json.loads(r.content.decode('utf-8'))
        result = r.json()
        status = result['status']
        progress = result['progress']

        time.sleep(interval)

    LOGGER.info('Operation %s done!', name)