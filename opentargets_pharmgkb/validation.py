import json
import logging

import jsonschema
import requests

logging.basicConfig()
logger = logging.getLogger(__package__)
logger.setLevel(level=logging.DEBUG)


def get_ot_json_schema():
    # Temporary schema url for testing
    schema_url = 'https://raw.githubusercontent.com/opentargets/json_schema/ds_3110_pharmacogenomics_schema/schemas/pharmacogenomics.json'
    response = requests.get(schema_url)
    if response.ok:
        return response.json()
    raise ValueError('Problem getting JSON schema')


def validate_evidence_string(ev_string):
    try:
        ot_schema_contents = get_ot_json_schema()
        jsonschema.validate(ev_string, ot_schema_contents, format_checker=jsonschema.FormatChecker())
        return True
    except jsonschema.exceptions.ValidationError as err:
        logger.error('Error: evidence string does not validate against schema.')
        logger.error(f'Error message: {err}')
        logger.error(f'Complete evidence string: {json.dumps(ev_string)}')
        return False
    except jsonschema.exceptions.SchemaError:
        logger.error('Error: OpenTargets schema file is invalid')
        return False
