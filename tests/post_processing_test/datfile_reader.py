#
# Maintainer: Amadeus Gebauer
#

import re


def read_result_description(datfile):

    # read all sections
    sections = read_sections(datfile)

    # extract result description
    if 'RESULT DESCRIPTION' not in sections:
        return []

    result_descriptions = sections['RESULT DESCRIPTION']

    results = []

    for r in result_descriptions:
        if r.strip().startswith('//'):
            continue
        try:
            this_result = {}
            this_result['dis'] = read_option_item(r, 'DIS', 1)[0]
            node = read_option_item(r, 'NODE', 1)[0]
            if node is not None:
                this_result['node'] = int(node)
            element = read_option_item(r, 'ELEMENT', 1)[0]
            if element is not None:
                this_result['element'] = int(element)
            this_result['quantity'] = read_option_item(r, 'QUANTITY', 1)[0]
            value = read_option_item(r, 'VALUE', 1)[0]
            if value is not None:
                this_result['value'] = float(value)
            tol = read_option_item(r, 'TOLERANCE', 1)[0]
            if tol is not None:
                this_result['tolerance'] = float(tol)

            this_result['special_quantity'] = read_option_item(
                r, 'SPECIAL QUANTITY', 1)[0]

            results.append(this_result)
        except:
            pass

    return results


def read_option_item(line, option, num):
    regex = re.compile(
        r"(^| ){0}{1}($|\s)".format(re.escape(option), num * "[ ]+([\\S]+)")
    )

    # split comment
    line = line.split("//", 1)[0]

    # read option
    match = regex.search(line)

    if not match:
        return None, None

    if num == 1:
        return match.group(2), match.span(0)
    else:
        return [match.group(i) for i in range(2, num + 2)], match.span(0)


def read_sections(filename):
    """
    Reads a dat file format and returns a dictionary containing of the sections with their lines

    Args:
        filename: Path to the dat file

    Returns:
        dict: Dictionary with section names as key and lines as value
    """
    content = {}
    with open(filename, "r") as origin:
        re_title = re.compile(r"^-{3,}(.*)")

        current_section = ""
        content[current_section] = []
        for line in origin:
            line_no_comment = line.split("//", 1)[0].strip()
            match_title = re_title.match(line_no_comment)
            if match_title:
                # this is a section title
                current_section = match_title.group(1)
                if current_section in content:
                    raise ValueError(
                        "{0} is dublicate!".format(current_section))

                content[current_section] = []
            else:
                content[current_section].append(line)

    return content
