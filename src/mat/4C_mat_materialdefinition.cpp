// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_materialdefinition.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_mat_par_bundle.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN



/*======================================================================*/
/*======================================================================*/
Mat::MaterialDefinition::MaterialDefinition(
    std::string materialname, std::string description, Core::Materials::MaterialType mattype)
    : materialname_(std::move(materialname)),
      description_(std::move(description)),
      mattype_(mattype)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::MaterialDefinition::add_component(Core::IO::InputSpec&& c)
{
  components_.emplace_back(std::move(c));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<std::pair<int, Core::IO::InputParameterContainer>> Mat::MaterialDefinition::read(
    Core::IO::InputFile& input)
{
  std::string name = "MATERIALS";

  std::vector<std::pair<int, Core::IO::InputParameterContainer>> found_materials;
  auto input_line = Core::IO::InputSpecBuilders::group(components_);
  for (const auto& line : input.lines_in_section(name))
  {
    Core::IO::ValueParser parser(
        line, {.user_scope_message = "While reading 'MATERIALS' section: ",
                  .base_path = input.file_for_section(name).parent_path()});

    parser.consume("MAT");
    const int matid = parser.read<int>();
    const std::string name = parser.read<std::string>();

    if (name == materialname_)
    {
      if (matid <= -1) FOUR_C_THROW("Illegal negative ID provided");

      Core::IO::InputParameterContainer input_data;
      Core::IO::fully_parse(parser, input_line, input_data);

      found_materials.emplace_back(matid, std::move(input_data));
    }
  }

  return found_materials;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& Mat::MaterialDefinition::print(
    std::ostream& stream, const Core::FE::Discretization* dis)
{
  // a string holding the comment indicating symbols for DAT input file
  const std::string comment = "//";

  // the descriptive lines (comments)
  stream << comment << '\n';
  stream << comment << " " << description_ << '\n';
  for (const auto& c : components_)
  {
    stream << comment << " " << c.name() << (c.required() ? " " : " (optional) ") << c.description()
           << '\n';
  }

  // the default line
  stream << comment << "MAT 0   " << materialname_ << "   ";
  auto input_line = Core::IO::InputSpecBuilders::group(components_);

  Core::IO::InputParameterContainer container;
  input_line.set_default_value(container);
  input_line.print(stream, container);

  stream << '\n';

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::append_material_definition(std::vector<std::shared_ptr<MaterialDefinition>>& matlist,
    const std::shared_ptr<MaterialDefinition>& mat)
{
  // test if material was defined with same name or type
  std::vector<std::shared_ptr<Mat::MaterialDefinition>>::const_iterator m;
  for (m = matlist.begin(); m != matlist.end(); ++m)
  {
    std::shared_ptr<Mat::MaterialDefinition> mmd = *m;

    if (mmd->type() == mat->type())
      FOUR_C_THROW(
          "Trying to define two materials with the same type '%d'\n"
          "Please revise your definitions of valid materials",
          mmd->type());

    if (mmd->name() == mat->name())
      FOUR_C_THROW(
          "Trying to define two materials with the same name '%s'\n"
          "Please revise your definitions of valid materials",
          mmd->name().c_str());
  }

  // no coincidence found
  if (m == matlist.end())
    matlist.push_back(mat);
  else
    FOUR_C_THROW("Trouble in determining coincidences of material definitions");
}

FOUR_C_NAMESPACE_CLOSE
