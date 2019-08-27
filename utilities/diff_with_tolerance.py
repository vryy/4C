# -*- coding: utf-8 -*-
"""
This script compares two files with a tolerance.
"""


# Import python modules.
import os
import sys
import numpy as np


def convert_line_to_array(line):
    """
    Try to convert a string to an array with floats.
    """

    # Split up the line.
    if ',' in line:
        line = line.split(',')
    else:
        line = line.split()
    line = [item.strip() for item in line]

    # Try to parse all parts of the line to a float.
    line_array = []
    try:
        for item in line:
            line_array.append(float(item))
    except ValueError:
        return False, None

    # All floats could be converted. We still have to check if the array is empty.
    if len(line_array) == 0:
        return False, None
    else:
        return True, line_array


def read_csv(path):
    """
    Load a csv file as a numpy array.
    """
    if not os.path.isfile(path):
        raise ValueError('The file {} does not exist.'.format(path))

    # Load the lines in the file.
    with open(path, 'r') as f:
        lines = f.readlines()

    # Go through each line an try to convert it to float.
    data = []
    for line in lines:
        is_float_line, line_array = convert_line_to_array(line)
        if is_float_line:
            data.append(line_array)

    return np.array(data)


if __name__ == '__main__':
    """
    Execution part of script.
    """

    # Read arguments.
    if not len(sys.argv) == 4:
        raise ValueError(('Wrong number of input arguments. Got {} instead '
            + 'of 4!').format(len(sys.argv)))
    eps = float(sys.argv[1])
    file_ref = sys.argv[2]
    file_comp = sys.argv[3]

    print('\n\nCompare files with tolerance{}:\n{}\n{}'.format(eps, file_ref, file_comp))

    # Load each file as a real array.
    data_ref = read_csv(file_ref)
    data_comp = read_csv(file_comp)

    # Calculate the absolute error in the data entries.
    abs_diff = np.abs(data_comp - data_ref)
    if eps > np.max(abs_diff):
        print('CSV comparison successful!')
    else:
        print('Largest error is {}'.format(np.max(abs_diff)))
        raise ValueError('CSV comparison failed!')
