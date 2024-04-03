import logging

import requests

from requests import RequestException
from retry import retry


logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

eutils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
esearch_url = eutils_url + 'esearch.fcgi'
esummary_url = eutils_url + 'esummary.fcgi'
efetch_url = eutils_url + 'efetch.fcgi'


@retry(exceptions=(ConnectionError, RequestException), tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def get_rsid_summaries(rsids, api_key=None):
    payload = {'db': 'snp', 'id': ','.join(rsids), 'retmode': 'JSON'}
    if api_key:
        payload['api_key'] = api_key
    req = requests.get(esummary_url, params=payload)
    req.raise_for_status()
    data = req.json()
    if 'result' in data:
        results = data['result']
        results.pop('uids')
        return results
    return {}


def parse_spdi_from_ncbi_result(ncbi_result):
    """
    Parse SPDI expressions from NCBI result.

    :param ncbi_result: dict response from NCBI summary endpoint for a single rsID
    :return: SPDI coords (chrom, pos, ref, list of alts)
    """
    spdis = ncbi_result['spdi'].split(',')
    # Parse each SPDI into coordinates
    chrom = pos = ref = None
    alts = []
    for spdi in spdis:
        # TODO warn or something if chrom/pos/ref differs among these
        chrom, pos, ref, alt = spdi.split(':')
        alts.append(alt)
    return chrom, pos, ref, alts


def get_spdi_coords_for_rsid(rsid):
    """
    Return SPDI coordinates for a single rsIDs by querying NCBI.

    :param rsid:
    :return: dict mapping rsid to SPDI coords as tuple (chrom, pos, ref, list of alts)
    """
    if rsid.startswith('rs'):
        rsid = rsid[2:]
    summary_result = get_rsid_summaries([rsid])
    if rsid not in summary_result:
        logger.warning(f'Could not get SPDI from NCBI for rs{rsid}')
        return None, None, None, []
    return parse_spdi_from_ncbi_result(summary_result[rsid])


def get_spdi_coords_for_rsids(rsids):
    """
    Return SPDI coordinates for a list of rsIDs, using single batch query to NCBI.

    :param rsids: list of rsIDs
    :return: dict mapping rsid to SPDI coords as tuple (chrom, pos, ref, list of alts)
    """
    rsids = [rs[2:] if rs.startswith('rs') else rs for rs in rsids]
    summary_results = get_rsid_summaries(rsids)
    rsid_to_spdis = {}
    for k in summary_results:
        rsid = f'rs{k}'
        rsid_to_spdis[rsid] = parse_spdi_from_ncbi_result(summary_results[k])
    if set(rsids) != set(rsid_to_spdis.keys()):
        logger.warning('Did not get SPDI from NCBI for all rsids')
    return rsid_to_spdis
