from readutils import read_string, read_int, read_ints, read_floats

KNOWN_ELEMENTS = {
    'point': 1,
    'bar2': 2,
    'bar3': 3,
    'tria3': 3,
    'tria6': 6,
    'quad4': 4,
    'quad8': 8,
    'tetra4': 4,
    'tetra10': 10,
    'pyramid5': 5,
    'pyramid13': 13,
    'penta6': 6,
    'penta15': 15,
    'hexa8': 8,
    'hexa20': 20
}


def read_geometry(path):
    with open(path, 'rb') as f:

        geometry = {}
        geometry['node_ids'] = None
        geometry['coordinates'] = None
        geometry['elements'] = {}

        try:
            # ignore first 3 lines
            f.seek(80+80+80+80, 1)

            node_id_given = 'given' in read_string(f, 80)
            element_id_given = 'given' in read_string(f, 80)
            f.seek(80+4+80+80, 1)

            # read node ids if available
            number_nodes = read_int(f)
            if node_id_given:
                geometry['node_ids'] = read_ints(
                    f, number_nodes)

            # read coordinates
            geometry['coordinates'] = read_floats(
                f, number_nodes*3).reshape(3, number_nodes).transpose()

            # read different element types
            element_type = read_string(f, 80).rstrip('\x00')
            while not element_type.startswith('part') and not element_type.startswith('END'):

                # read element ids if available
                element_ids = None
                number_elements = read_int(f)
                if element_id_given:
                    element_ids = read_ints(f, number_elements)

                # read node ids

                if element_type not in KNOWN_ELEMENTS:
                    raise ValueError(
                        'Unknown element type {0}'.format(element_type))

                number_nodes_per_element = KNOWN_ELEMENTS[element_type]

                element_nodes = read_ints(f, number_elements*number_nodes_per_element) \
                    .reshape((number_elements, number_nodes_per_element))

                # create element list
                geometry['elements'][element_type] = {
                    'ids': element_ids,
                    'nodes': element_nodes
                }

                # read next element type
                element_type = read_string(f, 80).rstrip('\x00')
        except EOFError:
            pass

    return geometry
