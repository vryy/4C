#
# Maintainer: Amadeus Gebauer
#

import numpy as np
from readutils import read_string, read_int, read_ints, read_floats
import os

DIM = {
    'scalar': 1,
    'vector': 3,
    'tensor sym': 6,
}


def read_variable_per_node(variable, timestep, node, dim, numnodes, basepath):
    with VariableStream(variable, numnodes, basepath) as var_stream:

        last = None
        for i in range(0, timestep+1):

            # var_stream.seek(timestep)

            last = float(var_stream.read_step()[node, dim])
        return last


class VariableStream:

    def __init__(self, variable, numnodes, basepath):
        self.variable = variable
        self.number_items = numnodes
        self.current_timestep = 0

        self.dim = DIM[variable['type'].split(' ')[0]]

        self.path = os.path.join(basepath, variable['file'])

    def seek(self, number_timesteps):
        # skip certain amounts of timesteps
        self.current_timestep += number_timesteps

        # compute size of each timestep
        timestep_byte_size = 80+80+80+4+80+self.dim*self.number_items*4+80

        # seek reader to the timestep
        self.f.seek(number_timesteps*timestep_byte_size, 1)

    def read_step(self):
        if self.f == None or self.f.closed:
            raise Exception('You have to open the VariableStream before')

        try:
            var = np.zeros((self.number_items, self.dim))

            # read next timestep
            while True:
                filetest = read_string(self.f, 80)
                if filetest.startswith('BEGIN TIME STEP'):
                    break

            # next bytes are just uninteresting chars
            self.f.seek(80+80+4, 1)

            # read ele type if variable per element
            # self.f.seek(80)
            read_string(self.f, 80)

            # interesting floats
            var = read_floats(self.f, self.number_items *
                              self.dim).reshape(self.dim, self.number_items).transpose()

            # uninteresting chars
            self.f.seek(80, 1)

            # increase current timestep
            self.current_timestep += 1
            return var
        except EOFError:
            return None

    def __next__(self):
        var = self.read_step()
        if var is None:
            raise StopIteration

        return var

    next = __next__  # Python 2

    def __enter__(self):
        # open file
        self.f = open(self.path, 'rb')
        return self

    def __exit__(self, type, value, traceback):
        # close file
        self.f.close()

    def __iter__(self):
        return self
