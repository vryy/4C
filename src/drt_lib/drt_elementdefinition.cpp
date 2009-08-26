#ifdef CCADISCRET

#include "drt_elementdefinition.H"


/*----------------------------------------------------------------------*/
//! Print function to be called from C
/*----------------------------------------------------------------------*/
extern "C"
void PrintElementDatHeader()
{
  DRT::INPUT::ElementDefinition ed;

  ed.SetupValidElementLines();

  ed.PrintSectionHeader(std::cout,"STRUCTURE ELEMENTS");

  ed.PrintElementLines(std::cout,"ART");
  ed.PrintElementLines(std::cout,"BEAM2");
  ed.PrintElementLines(std::cout,"BEAM2R");
  ed.PrintElementLines(std::cout,"BEAM3");
  ed.PrintElementLines(std::cout,"CONSTRELE2");
  ed.PrintElementLines(std::cout,"CONSTRELE3");
  ed.PrintElementLines(std::cout,"PTET4");
  ed.PrintElementLines(std::cout,"SOLID3");
  ed.PrintElementLines(std::cout,"SOLIDH20");
  ed.PrintElementLines(std::cout,"SOLIDH27");
  ed.PrintElementLines(std::cout,"SOLIDH8");
  ed.PrintElementLines(std::cout,"SOLIDH8P1J1");
  ed.PrintElementLines(std::cout,"SOLIDSH8");
  ed.PrintElementLines(std::cout,"SOLIDSH8P8");
  ed.PrintElementLines(std::cout,"SOLIDSHW6");
  ed.PrintElementLines(std::cout,"SOLIDT10");
  ed.PrintElementLines(std::cout,"SOLIDT4");
  ed.PrintElementLines(std::cout,"SOLIDW6");
  ed.PrintElementLines(std::cout,"TORSION2");
  ed.PrintElementLines(std::cout,"TORSION3");
  ed.PrintElementLines(std::cout,"TRUSS2");
  ed.PrintElementLines(std::cout,"TRUSS3");
  ed.PrintElementLines(std::cout,"WALL");

  ed.PrintSectionHeader(std::cout,"FLUID ELEMENTS");
  ed.PrintElementLines(std::cout,"COMBUST3");
  ed.PrintElementLines(std::cout,"CONDIF2");
  ed.PrintElementLines(std::cout,"CONDIF3");
  ed.PrintElementLines(std::cout,"FLUID2");
  ed.PrintElementLines(std::cout,"FLUID3");
  ed.PrintElementLines(std::cout,"TRANSP");
  ed.PrintElementLines(std::cout,"XDIFF3");
  ed.PrintElementLines(std::cout,"XFLUID3");

  ed.PrintSectionHeader(std::cout,"ALE ELEMENTS");
  ed.PrintElementLines(std::cout,"ALE2");
  ed.PrintElementLines(std::cout,"ALE3");

  //ed.PrintElementLines(std::cout,"BELE3");
  //ed.PrintElementLines(std::cout,"VELE3");

  //ed.PrintSectionHeader(std::cout,"THERMAL ELEMENTS");
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
  SetupFluid3Lines();
  //SetupSolidH8Lines();
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupFluid3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["FLUID3"];

  defs["HEX8"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedIntVector("GP",3)
    .AddNamedString("CA")
    ;

  defs["HEX20"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedIntVector("GP",3)
    .AddNamedString("CA")
    ;

  defs["HEX27"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedIntVector("GP",3)
    .AddNamedString("CA")
    ;

  defs["TET4"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    ;

  defs["TET10"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    ;

  defs["WEDGE6"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE15"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["PYRAMID5"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS8"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("NURBS8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS27"]
//     .AddInt("ID")
//     .AddTag("FLUID3")
    .AddDoubleVector("NURBS27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidH8Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDH8"];

#if 0
  defs["HEX8"]
//     .AddInt("ID")
//     .AddTag("SOLIDH8")
    .AddDoubleVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TET")
    .AddNamedInt("HYB")
    .AddNamedInt("FORM")
    .AddNamedString("STRESSES")
    .AddNamedString("TSI_COUPTYP")
    ;
#endif

  defs["HEX8"]
//     .AddInt("ID")
//     .AddTag("SOLIDH8")
    .AddDoubleVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("EAS")
    ;

}


#endif
