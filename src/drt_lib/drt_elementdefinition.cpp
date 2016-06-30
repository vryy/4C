/*----------------------------------------------------------------------*/
/*!
\file drt_elementdefinition.cpp

\brief Central storage of element input line definitions

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#include "drt_elementdefinition.H"
#include "drt_parobjectfactory.H"


/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void PrintElementDatHeader()
{
  DRT::INPUT::ElementDefinition ed;
  ed.PrintElementDatHeaderToStream(std::cout);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::PrintElementDatHeaderToStream(std::ostream& stream)
{
  SetupValidElementLines();

  PrintSectionHeader(stream,"STRUCTURE ELEMENTS");

  //  PrintElementLines(stream,"ART");
  PrintElementLines(stream,"BEAM3");
  PrintElementLines(stream,"BEAM3R");
  PrintElementLines(stream,"BEAM3CL");
  PrintElementLines(stream,"BEAM3EB");
  PrintElementLines(stream,"BEAM3EBTOR");
  PrintElementLines(stream,"BEAM3EBANISOTROP");
  PrintElementLines(stream,"BEAM3K");
  PrintElementLines(stream,"BELE2");
  PrintElementLines(stream,"BELE3");
  PrintElementLines(stream,"DISCSH3");
  PrintElementLines(stream,"RIGIDSPHERE");
  //PrintElementLines(stream,"CONSTRELE2");
  //PrintElementLines(stream,"CONSTRELE3");
  //PrintElementLines(stream,"PTET4");
  PrintElementLines(stream,"NSTET4");
  PrintElementLines(stream,"NSTET5");
  //PrintElementLines(stream,"NSTET5SCATRA");
  PrintElementLines(stream,"SHELL8");
  PrintElementLines(stream,"SOLIDH18");
  PrintElementLines(stream,"SOLIDH20");
  PrintElementLines(stream,"SOLIDH27");
  PrintElementLines(stream,"SOLIDH27PORO");
  PrintElementLines(stream,"SOLIDH27PLAST");
  PrintElementLines(stream,"SOLIDH27THERMO");
  PrintElementLines(stream,"SONURBS27THERMO");
  PrintElementLines(stream,"SOLIDH20THERMO");
  PrintElementLines(stream,"SONURBS27");
  PrintElementLines(stream,"SOLIDH8");
  PrintElementLines(stream,"MEMBRANE");
  PrintElementLines(stream,"SOLIDH8P1J1");
  PrintElementLines(stream,"SOLIDH8FBAR");
  PrintElementLines(stream,"SOLIDH8FBARSCATRA");
  PrintElementLines(stream,"SOLIDH8FBARTHERMO");
  PrintElementLines(stream,"SOLIDH8PORO");
  PrintElementLines(stream,"SOLIDH8THERMO");
  PrintElementLines(stream,"SOLIDH8PLAST");
  PrintElementLines(stream,"SOLIDH8SCATRA");
  PrintElementLines(stream,"SOLIDSH18");
  PrintElementLines(stream,"SOLIDSH18PLAST");
  PrintElementLines(stream,"SOLIDSH8");
  PrintElementLines(stream,"SOLIDSH8PLAST");
  PrintElementLines(stream,"SOLIDSH8P8");
  PrintElementLines(stream,"SOLIDSHW6");
  PrintElementLines(stream,"SOLIDT10");
  PrintElementLines(stream,"SOLIDT4");
  PrintElementLines(stream,"SOLIDT4PORO");
  PrintElementLines(stream,"SOLIDT4THERMO");
  PrintElementLines(stream,"SOLIDT10THERMO");
  PrintElementLines(stream,"SOLIDT4PLAST");
  PrintElementLines(stream,"SOLIDT4SCATRA");
  PrintElementLines(stream,"SOLIDT10SCATRA");
  PrintElementLines(stream,"SOLIDW6");
  PrintElementLines(stream,"SOLIDP5");
  PrintElementLines(stream,"SOLIDP5FBAR");
  PrintElementLines(stream,"SOLIDW6SCATRA");
  PrintElementLines(stream,"TORSION2");
  PrintElementLines(stream,"TORSION3");
  PrintElementLines(stream,"TRUSS2");
  PrintElementLines(stream,"TRUSS3");
  PrintElementLines(stream,"WALL");
  PrintElementLines(stream,"WALLQ4PORO");
  PrintElementLines(stream,"WALLQ9PORO");


  PrintSectionHeader(stream,"FLUID ELEMENTS");
  PrintElementLines(stream,"COMBUST3");
  PrintElementLines(stream,"FLUID");
  PrintElementLines(stream,"FLUIDXW");
  PrintElementLines(stream,"FLUIDHDG");
  PrintElementLines(stream,"FLUIDIMMERSED");
  PrintElementLines(stream,"FLUIDPOROIMMERSED");
  PrintElementLines(stream,"MEFLUID");

  PrintSectionHeader(stream,"LUBRICATION ELEMENTS");
  PrintElementLines(stream,"LUBRICATION");

  PrintSectionHeader(stream,"TRANSPORT ELEMENTS");
  PrintElementLines(stream,"TRANSP");
  PrintElementLines(stream,"METRANSP");

  PrintSectionHeader(stream,"ALE ELEMENTS");
  PrintElementLines(stream,"ALE2");
  PrintElementLines(stream,"ALE3");

  //PrintElementLines(stream,"BELE3_3");
  //PrintElementLines(stream,"VELE3");

  PrintSectionHeader(stream,"THERMO ELEMENTS");
  PrintElementLines(stream,"THERMO");

  PrintSectionHeader(stream,"ARTERY ELEMENTS");
  PrintElementLines(stream,"ART");

  PrintSectionHeader(stream,"REDUCED D AIRWAYS ELEMENTS");
  PrintElementLines(stream,"RED_AIRWAY");
  PrintElementLines(stream,"RED_ACINUS");
  PrintElementLines(stream,"RED_ACINAR_INTER_DEP");

  PrintElementLines(stream,"DCOHESIVE");

  PrintSectionHeader(stream,"ACOUSTIC ELEMENTS");
  PrintElementLines(stream,"ACOUSTIC");
  PrintElementLines(stream,"ACOUSTICSOL");

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::PrintSectionHeader(std::ostream& stream, std::string name, bool color)
{
  std::string blue2light = "";
  std::string bluelight = "";
  std::string redlight = "";
  std::string yellowlight = "";
  std::string greenlight = "";
  std::string magentalight = "";
  std::string endcolor = "";

  if (color)
  {
  }

  unsigned l = name.length();
  stream << redlight << "--";
  for (int i=0; i<std::max<int>(65-l,0); ++i) stream << '-';
  stream << greenlight << name << endcolor << '\n';
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::PrintElementLines(std::ostream& stream, std::string name)
{
  if (definitions_.find(name)!=definitions_.end())
  {
    std::map<std::string,LineDefinition>& defs = definitions_[name];
    for (std::map<std::string,LineDefinition>::iterator i=defs.begin(); i!=defs.end(); ++i)
    {
      stream << "// 0 " << name << " " << i->first << " ";
      i->second.Print(stream);
      stream << '\n';
    }
  }
  else
    stream << "no element type '" << name << "' defined\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupValidElementLines()
{
  DRT::ParObjectFactory::Instance().SetupElementDefinition( definitions_ );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::LineDefinition* DRT::INPUT::ElementDefinition::ElementLines(std::string name, std::string distype)
{
  // This is ugly. But we want to access both maps just once.
  std::map<std::string,std::map<std::string,LineDefinition> >::iterator j = definitions_.find(name);
  if (j!=definitions_.end())
  {
    std::map<std::string,LineDefinition>& defs = j->second;
    std::map<std::string,LineDefinition>::iterator i = defs.find(distype);
    if (i!=defs.end())
    {
      return &i->second;
    }
  }
  return NULL;
}
