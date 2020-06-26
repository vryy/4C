#
# Maintainer: Amadeus Gebauer
#

import numpy as np
import re

CURRENT_SECTION_FORMAT = 0
CURRENT_SECTION_GEOMETRY = 1
CURRENT_SECTION_VARIABLE = 2
CURRENT_SECTION_TIME = 3
CURRENT_SECTION_FILE = 4


def read_case(file):
    case = {}
    case['timesets'] = {}
    case['filesets'] = {}
    case['geometries'] = []
    case['variables'] = {}
    with open(file, 'r') as f:
        current_section = None
        current_index = None

        lines = f.readlines()
        # first read file- and timesets
        for line in lines:
            case_line = _read_case_line(line)

            if case_line == None:
                continue
            elif line.strip() == 'TIME':
                current_section = CURRENT_SECTION_TIME
            elif line.strip() == 'FILE':
                current_section = CURRENT_SECTION_FILE
            elif line.strip() == 'FORMAT':
                current_section = CURRENT_SECTION_FORMAT
            elif line.strip() == 'GEOMETRY':
                current_section = CURRENT_SECTION_GEOMETRY
            elif line.strip() == 'VARIABLE':
                current_section = CURRENT_SECTION_VARIABLE
            elif current_section == CURRENT_SECTION_TIME:

                if case_line['name'] == 'time set':
                    current_index = int(case_line['values'][0])
                    case['timesets'][current_index] = None
                elif case_line['name'] == 'time values':
                    case['timesets'][current_index] = np.asarray(
                        [float(i) for i in case_line['values']])
                elif case_line['name'] == None:
                    case['timesets'][current_index] = np.append(case['timesets'][current_index], np.asarray([
                        float(i) for i in case_line['values']]))

            elif current_section == CURRENT_SECTION_FILE:

                if case_line['name'] == 'file set':
                    current_index = int(case_line['values'][0])
                    case['filesets'][current_index] = 0
                elif case_line['name'] == 'number of steps':
                    case['filesets'][current_index] = int(
                        case_line['values'][0])

        # now read geometry and variable
        for line in lines:
            case_line = _read_case_line(line)

            if case_line == None:
                continue
            elif line.strip() == 'TIME':
                current_section = CURRENT_SECTION_TIME
            elif line.strip() == 'FILE':
                current_section = CURRENT_SECTION_FILE
            elif line.strip() == 'FORMAT':
                current_section = CURRENT_SECTION_FORMAT
            elif line.strip() == 'GEOMETRY':
                current_section = CURRENT_SECTION_GEOMETRY
            elif line.strip() == 'VARIABLE':
                current_section = CURRENT_SECTION_VARIABLE
            elif current_section == CURRENT_SECTION_GEOMETRY:
                if case_line['name'] == 'model':
                    case['geometries'].append({
                        'timeset': int(case_line['values'][0]),
                        'fileset': int(case_line['values'][0]),
                        'geofile': case_line['values'][2]
                    })
            elif current_section == CURRENT_SECTION_VARIABLE:
                case['variables'][case_line['values'][2]] = {
                    'timeset': int(case_line['values'][0]),
                    'fileset': int(case_line['values'][1]),
                    'type': case_line['name'],
                    'file': case_line['values'][3]
                }
        return case

    return case


def _read_case_line(line):
    """Reads a line of a case file and tries to interpret it"""
    stripped_line = line.strip()

    result = {}

    if stripped_line == '':
        return None

    if ':' in stripped_line:
        # this is a name with values
        stripped_line_part = stripped_line.partition(':')
        result['name'] = stripped_line_part[0]
        result['values'] = _read_case_values(stripped_line_part[2])
    else:
        # that are just values
        result['name'] = None
        result['values'] = _read_case_values(stripped_line)

    return result


def _read_case_values(line_part):
    return re.sub('( |\t)+', ' ', line_part.strip()).split(' ')
