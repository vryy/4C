# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Create json schema from 4C yaml file."""
import yaml
import warnings
import json
import argparse

FOURC_TO_JSON_SCHEMA_DICT = {
    "double": "number",
    "int": "integer",
    "bool": "boolean",
    "string": "string",
}


def map_property_to_json_schema(
    properties_original: dict,
    key_original: str,
    properties_schema: dict,
    key_schema: str,
    filter=lambda x: x,
    optional: bool = False,
):
    value = properties_original.pop(key_original, None)
    if value is None:
        if not optional:
            raise ValueError(f"Value {key_original} was not provided!")
    else:
        properties_schema[key_schema] = filter(value)


def create_properties_dict(properties_entry: dict) -> dict:
    """Generate properties dict for one option.

    Args:
        properties_entry (dict): 4C exported dict

    Returns:
        dict: properties dict for the JSON schema
    """
    properties = {}

    # Loop over the definitions of the properties
    for key, value in properties_entry.items():

        # Check if the optains contains the type
        if "type" not in value.keys():
            raise Exception(f"Ups could not load {key}: {value}")

        description = value.pop("description", "No description yet")

        # Add default to the description if provided
        if "default" in value:
            description += f"\nDefault: {value['default']}"

        properties[key] = {"description": description}
        properties[key]["title"] = f"{key} ({value['type']})"
        map_property_to_json_schema(
            value, "type", properties[key], "type", filter=FOURC_TO_JSON_SCHEMA_DICT.get
        )
        map_property_to_json_schema(
            value, "valid options", properties[key], "enum", optional=True
        )

        # Defaults are ready to go, but currently not properly provided!
        # For now we do not export them!
        # map_property_to_json_schema(value, "default", properties[key], "default", optional=True)

    return properties


def create_json_schema(fourc_metadata: dict) -> dict:
    """Create a JSON schema dict.

    Args:
        fourc_metada (dict): Metadata exported by 4C

    Returns:
        dict: data for the JSON schema
    """
    commit_hash = fourc_metadata["metadata"]["commit_hash"]
    json_schema_dict = {
        "$id": commit_hash,  # not the best solution but for now it works
        "$schema": "https://json-schema.org/draft/2020-12/schema",  # Schema standard we use
        "title": f"4C Schema for commit hash: {commit_hash}",
        "description": "A schema for 4C yaml input files.",
        "type": "object",
    }

    properties_dict = {}

    # Title section
    properties_dict["TITLE"] = {
        "title": "Section 'Title'",
        "description": "You can add a description to your input file here.",
        "type": "string",
        "default": f"Using commit hash {commit_hash}",
    }

    # Includes section
    properties_dict["INCLUDES"] = {
        "title": "Section 'INCLUDES'",
        "description": "Paths to additional input files.",
        "type": "array",
        "items": {"type": "string"},
    }

    # Loop over sections
    for section_name, section_options in fourc_metadata["parameters"].items():

        # Some sections are empty
        if not section_options:
            warnings.warn(f"For section {section_name} no metadata is available!")

        description = f"Section '{section_name}'\n"
        description += section_options.pop("description", "No description yet")
        properties_dict[section_name] = {
            "title": section_name,
            "description": description,
            "type": "object",  # this type is required by the JSON schema template
            "properties": create_properties_dict(section_options),
        }

    json_schema_dict["properties"] = properties_dict
    return json_schema_dict


def load_fourc_yaml_and_escape_bools(fourc_yaml_path: str) -> dict:
    """Load 4C yaml file.

    Since in 4C currently bools are not supported in the inputs, they need some special treatment to avoid converting them to bool while loading. Once this is supported in the code, this functionally can be dropped.

    Args:
        fourc_yaml_path (str): Path to the 4C metadata file

    Returns:
        dict: loaded options
    """
    from yaml.constructor import SafeConstructor
    from yaml.loader import Reader, Scanner, Parser, Composer, Resolver

    class CustomSafeConstructor(SafeConstructor):
        def add_bool(self, node):
            return self.construct_scalar(node)

    CustomSafeConstructor.add_constructor(
        "tag:yaml.org,2002:bool", CustomSafeConstructor.add_bool
    )

    class KeepStringSafeLoader(
        Reader, Scanner, Parser, Composer, CustomSafeConstructor, Resolver
    ):
        def __init__(self, stream):
            Reader.__init__(self, stream)
            Scanner.__init__(self)
            Parser.__init__(self)
            Composer.__init__(self)
            CustomSafeConstructor.__init__(self)
            Resolver.__init__(self)

    with open(fourc_yaml_path, "r", encoding="UTF-8") as f:
        fourc_metadata = yaml.load(f.read(), KeepStringSafeLoader)
    return fourc_metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create JSON schema from 4C metadata yaml"
    )
    parser.add_argument(
        "fourc_metadata_yaml_path", help="Path to the yaml file generated by 4C."
    )
    parser.add_argument("json_schema_path", help="Path for the JSON schema.")
    args = parser.parse_args()

    # Load the yaml file exported by 4C
    fourc_metadata = load_fourc_yaml_and_escape_bools(args.fourc_metadata_yaml_path)

    # create the json schema dict
    json_schema_dict = create_json_schema(fourc_metadata)

    # Export the JSON schema
    with open(args.json_schema_path, "w", encoding="UTF-8") as f:
        json.dump(json_schema_dict, f, indent=1)
