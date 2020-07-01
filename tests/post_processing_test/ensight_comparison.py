import argparse
from datfile_reader import read_result_description
from casefile_reader import read_case
from geometry_reader import read_geometry
from variable_reader import read_variable_per_node
import os
import sys


def get_corresponding_quantity(dis, quantity):
    if dis == 'structure':

        # structure quantitities
        if quantity == 'dispx':
            return 'displacement', 0
        elif quantity == 'dispy':
            return 'displacement', 1
        elif quantity == 'dispz':
            return 'displacement', 2

    elif dis == 'cardiovascular0d':
        # cardiovascular quantities
        return None, None

    # quantity not implemented
    raise NotImplementedError()


if __name__ == '__main__':

    # create parser for input arguments
    parser = argparse.ArgumentParser(
        description='Python scripts that checks the results from the result description of the dat file after postprocessing with the ensight filter')

    # add input arguments to parser
    parser.add_argument('datfile', help='Path to the .dat input file')
    parser.add_argument('casefile', help='Path to the ensight output file')

    # read input arguments
    args = parser.parse_args()

    result_description = read_result_description(args.datfile)
    case_file = read_case(args.casefile)

    if (len(case_file['geometries']) != 1):
        raise ValueError('Currently only 1 geometry is supported. Found {0} geometries'.format(
            len(case_file['geometries'])))

    basepath = os.path.dirname(args.casefile)
    geometry = read_geometry(os.path.join(
        basepath, case_file['geometries'][0]['geofile']))

    errors = 0
    success = 0
    for r in result_description:
        try:
            variable_name, dim = get_corresponding_quantity(
                r['dis'], r['quantity'])
        except NotImplementedError:
            print (
                '{0}: {1} is NOT IMPLEMENTED'.format(r['dis'], r['quantity']))
            errors += 1
            continue

        # check whether variable is implemented for postprocessing tests
        if variable_name == None:
            print (
                '{0}: {1} is NOT TESTABLE. SKIPPING!'.format(r['dis'], r['quantity']))
            continue

        # check whether variable exists in case file
        if variable_name not in case_file['variables']:
            print (
                '{0}: {1} at node {2}: variable could NOT BE FOUND in ensight file'.format(r['dis'], r['quantity'], r['node']))
            errors += 1
            continue

        variable = case_file['variables'][variable_name]

        # check, whether node can be found in the case file
        node_index = 0
        found_nid = False
        for nid in geometry['node_ids']:
            if nid == r['node']:
                found_nid = True
                break
            node_index += 1

        # check whether node could be found
        if not found_nid:
            print (
                '{0}: {1} at node {2}: node could NOT BE FOUND in ensight file'.format(r['dis'], r['quantity'], r['node']))
            errors += 1
            continue

        if variable['type'].endswith('per node'):
            timeset = case_file['timesets'][variable['timeset']]
            timestep = len(timeset) - 1

            myvalue = read_variable_per_node(
                variable, timestep, node_index, dim, len(geometry['node_ids']), basepath)

            # check correctness
            tolerance = max(1e-7, r['tolerance'])
            if abs(myvalue - r['value']) > tolerance:
                print (
                    '{0}: {1} at node {2} is WRONG --> actresult={3}, givenresult={4}, abs(diff)={5} > {6}'.format(r['dis'], r['quantity'], r['node'], myvalue, r['value'], abs(myvalue - r['value']), tolerance))
                errors += 1
                continue

            # finally we are here: Result is correct
            print (
                '{0}: {1} at node {2} is CORRECT, abs(diff)={3} < {4}'.format(r['dis'], r['quantity'], r['node'], abs(myvalue - r['value']), tolerance))

            success += 1
        else:
            print(variable['type'])

    if success == 0 and errors == 0:
        print('Nothing has been tested. Ensight comparison is not possible for your kind of results! You need to adapt it to test your quantities as well.')
        sys.exit(1)

    if errors > 0:
        sys.exit(1)
