/*----------------------------------------------------------------------*/
/*! \file

\brief Central storage of element input line definitions

\level 1


*/
/*----------------------------------------------------------------------*/



#include "4C_fem_general_element_definition.hpp"

#include "4C_comm_parobjectfactory.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void PrintElementDatHeader()
{
  Core::Elements::ElementDefinition ed;
  ed.print_element_dat_header_to_stream(std::cout);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Elements::ElementDefinition::print_element_dat_header_to_stream(std::ostream& stream)
{
  setup_valid_element_lines();

  print_section_header(stream, "STRUCTURE ELEMENTS");

  //  PrintElementLines(stream,"ART");
  print_element_lines(stream, "BEAM3");
  print_element_lines(stream, "BEAM3R");
  print_element_lines(stream, "BEAM3EB");
  print_element_lines(stream, "BEAM3K");
  print_element_lines(stream, "BELE3");
  print_element_lines(stream, "RIGIDSPHERE");
  print_element_lines(stream, "NSTET5");
  print_element_lines(stream, "SHELL7P");
  print_element_lines(stream, "SHELL7PSCATRA");
  print_element_lines(stream, "SOLID");
  print_element_lines(stream, "SOLIDPORO");
  print_element_lines(stream, "SOLIDSCATRA");
  print_element_lines(stream, "SOLIDH18_DEPRECATED");
  print_element_lines(stream, "SOLIDH20_DEPRECATED");
  print_element_lines(stream, "SOLIDH27_DEPRECATED");
  print_element_lines(stream, "SOLIDH27PORO");
  print_element_lines(stream, "SOLIDH27PLAST");
  print_element_lines(stream, "SOLIDH27THERMO");
  print_element_lines(stream, "SONURBS27THERMO");
  print_element_lines(stream, "SOLIDH20THERMO");
  print_element_lines(stream, "SONURBS27");
  print_element_lines(stream, "SOLIDH8_DEPRECATED");
  print_element_lines(stream, "MEMBRANE");
  print_element_lines(stream, "SOLIDH8P1J1");
  print_element_lines(stream, "SOLIDH8FBAR_DEPRECATED");
  print_element_lines(stream, "SOLIDH8FBARSCATRA_DEPRECATED");
  print_element_lines(stream, "SOLIDH8FBARTHERMO");
  print_element_lines(stream, "SOLIDH8PORO");
  print_element_lines(stream, "SOLIDH8POROSCATRA");
  print_element_lines(stream, "SOLIDH8POROP1");
  print_element_lines(stream, "SOLIDH8POROP1SCATRA");
  print_element_lines(stream, "SOLIDH8THERMO");
  print_element_lines(stream, "SOLIDH8PLAST");
  print_element_lines(stream, "SOLIDH8SCATRA");
  print_element_lines(stream, "SOLIDSH18");
  print_element_lines(stream, "SOLIDSH18PLAST");
  print_element_lines(stream, "SOLIDSH8");
  print_element_lines(stream, "SOLIDSH8PLAST");
  print_element_lines(stream, "SOLIDSH8P8");
  print_element_lines(stream, "SOLIDSHW6");
  print_element_lines(stream, "SOLIDT10_DEPRECATED");
  print_element_lines(stream, "SOLIDT4_DEPRECATED");
  print_element_lines(stream, "SOLIDT4PLAST");
  print_element_lines(stream, "SOLIDT4PORO");
  print_element_lines(stream, "SOLIDT4THERMO");
  print_element_lines(stream, "SOLIDT10THERMO");
  print_element_lines(stream, "SOLIDT4PLAST");
  print_element_lines(stream, "SOLIDT4SCATRA_DEPRECATED");
  print_element_lines(stream, "SOLIDT10SCATRA_DEPRECATED");
  print_element_lines(stream, "SOLIDW6_DEPRECATED");
  print_element_lines(stream, "SOLIDP5_DEPRECATED");
  print_element_lines(stream, "SOLIDP5FBAR_DEPRECATED");
  print_element_lines(stream, "SOLIDW6SCATRA_DEPRECATED");
  print_element_lines(stream, "TORSION3");
  print_element_lines(stream, "TRUSS3");
  print_element_lines(stream, "TRUSS3SCATRA");
  print_element_lines(stream, "WALL");
  print_element_lines(stream, "WALLNURBS");
  print_element_lines(stream, "WALLSCATRA");
  print_element_lines(stream, "WALLQ4PORO");
  print_element_lines(stream, "WALLQ4POROSCATRA");
  print_element_lines(stream, "WALLQ4POROP1");
  print_element_lines(stream, "WALLQ4POROP1SCATRA");
  print_element_lines(stream, "WALLQ9PORO");


  print_section_header(stream, "FLUID ELEMENTS");
  print_element_lines(stream, "FLUID");
  print_element_lines(stream, "FLUIDXW");
  print_element_lines(stream, "FLUIDHDG");
  print_element_lines(stream, "FLUIDHDGWEAKCOMP");
  print_element_lines(stream, "FLUIDIMMERSED");
  print_element_lines(stream, "FLUIDPOROIMMERSED");

  print_section_header(stream, "LUBRICATION ELEMENTS");
  print_element_lines(stream, "LUBRICATION");

  print_section_header(stream, "TRANSPORT ELEMENTS");
  print_element_lines(stream, "TRANSP");

  print_section_header(stream, "TRANSPORT2 ELEMENTS");
  print_element_lines(stream, "TRANSP");

  print_section_header(stream, "ALE ELEMENTS");
  print_element_lines(stream, "ALE2");
  print_element_lines(stream, "ALE3");

  // PrintElementLines(stream,"BELE3_3");
  // PrintElementLines(stream,"VELE3");

  print_section_header(stream, "THERMO ELEMENTS");
  print_element_lines(stream, "THERMO");

  print_section_header(stream, "ARTERY ELEMENTS");
  print_element_lines(stream, "ART");

  print_section_header(stream, "REDUCED D AIRWAYS ELEMENTS");
  print_element_lines(stream, "RED_AIRWAY");
  print_element_lines(stream, "RED_ACINUS");
  print_element_lines(stream, "RED_ACINAR_INTER_DEP");

  print_section_header(stream, "ELECTROMAGNETIC ELEMENTS");
  print_element_lines(stream, "ELECTROMAGNETIC");
  print_element_lines(stream, "ELECTROMAGNETICDIFF");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Elements::ElementDefinition::print_section_header(std::ostream& stream, std::string name)
{
  unsigned l = name.length();
  stream << "--";
  for (int i = 0; i < std::max<int>(65 - l, 0); ++i) stream << '-';
  stream << name << '\n';
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Elements::ElementDefinition::print_element_lines(std::ostream& stream, std::string name)
{
  if (definitions_.find(name) != definitions_.end())
  {
    std::map<std::string, Input::LineDefinition>& defs = definitions_[name];
    for (std::map<std::string, Input::LineDefinition>::iterator i = defs.begin(); i != defs.end();
         ++i)
    {
      stream << "// 0 " << name << " " << i->first << " ";
      i->second.print(stream);
      stream << '\n';
    }
  }
  else
    stream << "no element type '" << name << "' defined\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Elements::ElementDefinition::setup_valid_element_lines()
{
  Core::Communication::ParObjectFactory::instance().setup_element_definition(definitions_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Input::LineDefinition* Core::Elements::ElementDefinition::element_lines(
    std::string name, std::string distype)
{
  // This is ugly. But we want to access both maps just once.
  std::map<std::string, std::map<std::string, Input::LineDefinition>>::iterator j =
      definitions_.find(name);
  if (j != definitions_.end())
  {
    std::map<std::string, Input::LineDefinition>& defs = j->second;
    std::map<std::string, Input::LineDefinition>::iterator i = defs.find(distype);
    if (i != defs.end())
    {
      return &i->second;
    }
  }
  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
