import logging

import requests

LOGGER = logging.getLogger(__name__)


def add_tcp_firewall_rule(project, access_token, name, tag, port):
    """Adds a TCP firewall rule.

    TODO: docstring"""
    
    headers = {
        'Authorization': 'Bearer %s' % access_token.access_token
    }
    
    payload = {
        "name": name,
        "kind": "compute#firewall",
        "sourceRanges": [ "0.0.0.0/0" ],
        "targetTags": [ tag ],
        "allowed": [ { "IPProtocol": "tcp", "ports": [ str(port), ] } ],
        "network": "projects/%s/global/networks/default" % project
    }
    
    r = requests.post('https://www.googleapis.com/compute/v1/'
                      'projects/%s/global/firewalls'
                      % project, headers=headers, json=payload)
    r.raise_for_status()
    op = r.json()['name']

    LOGGER.info('Requested to add firewall rule for tag "%s" '
                '(operation "%s")...', tag, op)

    return op
