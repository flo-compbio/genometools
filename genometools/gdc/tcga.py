# Copyright (c) 2016 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Functions for downloading TCGA data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging
import json
import io
from collections import Iterable, OrderedDict

import requests
import pandas as pd

logger = logging.getLogger(__name__)


def get_clinical_data(tcga_id):
    """Get clinical data for a TCGA project.
    
    Parameters
    ----------
    tcga_id : str
        The TCGA project ID.
    
    Returns
    -------
    `pandas.DataFrame`
        The clinical data.abs

    Notes
    -----
    Clinical data is associated with individual cases (patients). These
    correspond to rows in the returned data frame, and are identified by
    12-character TCGA barcodes. 
    """

    payload = {
        'attachment': 'true',
        "filters": json.dumps({
            "op": "and",
            "content": [
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.program.name",
                        "value":["TCGA"]}},
                {
                    "op": "in",
                    "content": {
                        "field": "project.project_id",
                        "value": [tcga_id]}}]
        }),
        'fields': 'case_id',
        'expand': 'demographic,diagnoses,family_histories,exposures',
        'format': 'JSON',
        'pretty': 'true',
        'size': 10000,
        'filename': 'clinical.project-%s' % tcga_id,
    }
    r = requests.get('https://gdc-api.nci.nih.gov/cases', params=payload)
    
    j = json.loads(r.content.decode())
    clinical = {}
    valid = 0
    for s in j:
        if 'diagnoses' not in s:
            continue
        valid += 1
        assert len(s['diagnoses']) == 1
        diag = s['diagnoses'][0]
        tcga_id = diag['submitter_id'][:12]
        clinical[tcga_id] = diag

    logger.info('Found clinical data for %d cases.', valid)
    df = pd.DataFrame.from_dict(clinical).T
    df.sort_index(inplace=True)
    return df
    #clin.to_csv(download_file, sep='\t')


def get_biospecimen_data(tcga_id):
    """Get biospecimen data for a TCGA project.
    
    Parameters
    ----------
    tcga_id : str
        The TCGA project ID.
    
    Returns
    -------
    `pandas.DataFrame`
        The biospecmin data.

    Notes
    -----
    Biospecimen data is associated with individual vials. TCGA vials correspond
    to portions of a sample, and are uniquely identified by a 16-character
    barcode. For example, one vial can contain FFPE material and the other
    fresh-frozen material from the same sample. Each row in the returned data
    frame corresponds to a vial.
    """

    payload = {
        'attachment': 'true',
        "filters": json.dumps({
            "op": "and",
            "content": [
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.program.name",
                        "value":["TCGA"]}},
                {
                    "op": "in",
                    "content": {
                        "field": "project.project_id",
                        "value": [tcga_id]}}]
        }),
        'fields': 'case_id',
        'expand': ('samples,samples.portions,'
                   'samples.portions.analytes,'
                   'samples.portions.analytes.aliquots,'
                   'samples.portions.analytes.aliquots.annotations,'
                   'samples.portions.analytes.annotations,'
                   'samples.portions.submitter_id,'
                   'samples.portions.slides,'
                   'samples.portions.annotations,'
                   'samples.portions.center'),
        'format': 'JSON',
        'pretty': 'true',
        'size': 10000,
        'filename': 'biospecimen.project-%s' % tcga_id,
    }
    r = requests.get('https://gdc-api.nci.nih.gov/cases', params=payload)
    
    j = json.loads(r.content.decode())
    biospec = {}
    valid = 0
    for case in j:
        if 'samples' not in case:
            continue
        valid += 1
        for s in case['samples']:
            tcga_id = s['submitter_id'][:16]
            del s['portions']
            del s['submitter_id']
            biospec[tcga_id] = s

    logger.info('Found biospecimen data for %d cases.', valid)
    df = pd.DataFrame.from_dict(biospec).T
    df.sort_index(inplace=True)
    return df
    #df.to_csv(download_file, sep='\t')


def get_masked_cnv_manifest(tcga_id):
    """Get manifest for masked TCGA copy-number variation data.
    
    Params
    ------
    tcga_id : str
        The TCGA project ID.
    download_file : str
        The path of the download file.
        
    Returns
    -------
    `pandas.DataFrame`
        The manifest.
    """
    payload = {
        "filters": json.dumps({
            "op": "and",
            "content" : [
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.program.name",
                        "value":["TCGA"]}},
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.project_id",
                        "value":[tcga_id]}},
                {
                    "op":"in",
                    "content":{
                        "field":"files.data_category",
                        "value":["Copy Number Variation"]}},
                {
                    "op":"in",
                    "content":{
                        "field":"files.data_type",
                        "value":["Masked Copy Number Segment"]}}]
        }),
        "return_type":"manifest",
        "size":10000,
    }

    r = requests.get('https://gdc-api.nci.nih.gov/files', params=payload)
    df = pd.read_csv(io.StringIO(r.text), sep='\t', header=0)
    logger.info('Obtained manifest with %d files.', df.shape[0])
    return df


def get_expression_manifest(tcga_id):
    """Get manifest for masked TCGA gene expression data.
    
    Params
    ------
    tcga_id : str
        The TCGA project ID.
        
    Returns
    -------
    `pandas.DataFrame`
        The manifest.
    """
    payload = {
        "filters":json.dumps({
            "op":"and",
            "content":[
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.program.name",
                        "value":["TCGA"]}},
                {
                    "op":"in",
                    "content":{
                        "field":"cases.project.project_id",
                        "value":[tcga_id]}},
                {
                    "op":"in",
                    "content":{
                        "field":"experimental_strategy",
                        "value":["RNA-Seq"]}},
                {
                    "op":"in",
                    "content":{
                        "field":"files.analysis.workflow_type",
                        "value":["HTSeq - FPKM"]}}
            ]}),
        "return_type":"manifest",
        "size":10000,
    }

    r = requests.get('https://gdc-api.nci.nih.gov/files', params=payload)
    df = pd.read_csv(io.StringIO(r.text), sep='\t', header=0)
    logger.info('Obtained manifest with %d files.', df.shape[0])
    return df


def get_file_samples(file_ids):
    """Get TCGA associated sample barcodes for a list of file IDs. 
    
    Params
    ------
    file_ids : Iterable
        The file IDs.
        
    Returns
    -------
    `pandas.Series`
        Series containing file IDs as index and corresponding sample barcodes.
    """
    assert isinstance(file_ids, Iterable)

    # query TCGA API to get sample barcodes associated with file IDs
    payload = {
        "filters":json.dumps({
            "op":"in",
            "content":{
                "field":"files.file_id",
                "value": list(file_ids),
            }
        }),
        "fields":"file_id,cases.samples.submitter_id",
        "size":10000
    }
    r = requests.post('https://gdc-api.nci.nih.gov/files', data=payload)
    j = json.loads(r.content.decode('utf-8'))
    file_samples = OrderedDict()
    for hit in j['data']['hits']:
        file_id = hit['file_id']
        assert len(hit['cases']) == 1
        case = hit['cases'][0]
        assert len(case['samples']) == 1
        sample = case['samples'][0]
        sample_barcode = sample['submitter_id']
        file_samples[file_id] = sample_barcode

    df = pd.DataFrame.from_dict(file_samples, orient='index')
    df = df.reset_index()
    df.columns = ['file_id', 'sample_barcode']
    return df


def get_unique_sample_files(file_samples):
    """Filter file_sample data frame to only keep one file per sample.
    
    Params
    ------
    file_samples : `pandas.DataFrame`
        A data frame containing a mapping between file IDs and sample barcodes.
        This type of data frame is returned by :meth:`get_file_samples`.
        
    Returns
    -------
    `pandas.DataFrame`
        The filtered data frame.

    Notes
    -----
    In order to remove redundant files in a consistent fashion, the samples are
    sorted by file ID, and then the first file for each sample is kept.
    """
    assert isinstance(file_samples, pd.DataFrame)

    df = file_samples  # use shorter variable name

    # sort data frame by file ID
    df = df.sort_values('file_id')

    # - some samples have multiple files with the same barcode,
    #   corresponding to different aliquots
    # get rid of those duplicates
    logger.info('Original number of files: %d', len(df.index))
    df.drop_duplicates('sample_barcode', keep='first', inplace=True)
    logger.info('Number of files after removing duplicates from different '
                'aliquots: %d', len(df.index))

    # - some samples also have multiple files corresponding to different vials
    # add an auxilliary column that contains the sample barcode without the
    # vial tag (first 15 characters)
    df['sample_barcode15'] = df['sample_barcode'].apply(lambda x: x[:15])

    # use auxilliary column to get rid of duplicate files
    df.drop_duplicates('sample_barcode15', keep='first', inplace=True)
    logger.info('Number of files after removing duplicates from different '
                'vials: %d', len(df.index))

    # get rid of auxilliary column 
    df.drop('sample_barcode15', axis=1, inplace=True)

    # restore original order using the (numerical) index
    df.sort_index(inplace=True)
    
    # return the filtered data frame
    return df


def get_fresh_primary_tumors(biospecimen):
    """Filter biospecimen data to only keep non-FFPE primary tumor samples.
    
    Parameters
    ----------
    biospecimen : `pandas.DataFrame`
        The biospecimen data frame. This type of data frame is returned by
        :meth:`get_biospecimen_data`.
    
    Returns
    -------
    `pandas.DataFrame`
        The filtered data frame.
    """
    df = biospecimen  # use shorter variable name

    # get rid of FFPE samples
    num_before = len(df.index)
    df = df.loc[~df['is_ffpe']]
    logger.info('Excluded %d files associated with FFPE samples '
                '(out of %d files in total).',
                num_before - len(df.index), num_before)

    # only keep primary tumors
    num_before = len(df.index)
    df = df.loc[df['sample_type'] == 'Primary Tumor']
    logger.info('Excluded %d files not corresponding to primary tumor '
                'samples (out of %d files in total).',
                num_before - len(df.index), num_before)

    return df
