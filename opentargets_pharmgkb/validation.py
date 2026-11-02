import json
import logging
import os

import jsonschema

logging.basicConfig()
logger = logging.getLogger(__package__)
logger.setLevel(level=logging.DEBUG)

OT_JSON_SCHEMA = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'etc/pharmacogenomics.json')

def get_ot_json_schema():
    with open(OT_JSON_SCHEMA) as f:
        return json.load(f)


def validate_evidence_string(ev_string, ot_schema_contents=None):
    try:
        if not ot_schema_contents:
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
