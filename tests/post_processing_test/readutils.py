#
# Maintainer: Amadeus Gebauer
#

import numpy as np
import struct


def read_string(file_handle, length):
    """Reads a string of specified length in ASCII"""
    result_str = [None] * length
    pos = 0
    while pos < length:
        result_str[pos] = file_handle.read(1)
        pos += 1

    return b''.join(result_str).decode("ascii")


def read_int(file_handle):
    """Reads integer"""
    return struct.unpack('<i', file_handle.read(4))[0]


def read_ints(file_handle, length):
    """Reads integer"""
    return np.fromfile(file_handle, dtype='<i', count=length)


def read_floats(file_handle, length):
    """reads array of length floats"""
    return np.fromfile(file_handle, dtype='<f4', count=length)
