#ifdef CCADISCRET

#include "drt_elementdefinition.H"
#include "drt_linedefinition.H"


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
    definitions_[name]->Print(stream,false);
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
Teuchos::RCP<DRT::INPUT::Lines> DRT::INPUT::ElementDefinition::ElementLines(std::string name)
{
  if (definitions_.find(name)!=definitions_.end())
    return definitions_[name];
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupFluid3Lines()
{
  DRT::INPUT::LineDefinition fluid3_hex8;
  fluid3_hex8
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedIntVector("GP",3)
    .AddNamedString("CA")
    ;

  DRT::INPUT::LineDefinition fluid3_hex20;
  fluid3_hex20
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedIntVector("GP",3)
    .AddNamedString("CA")
    ;

  DRT::INPUT::LineDefinition fluid3_hex27;
  fluid3_hex27
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedIntVector("GP",3)
    .AddNamedString("CA")
    ;

  DRT::INPUT::LineDefinition fluid3_tet4;
  fluid3_tet4
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    ;

  DRT::INPUT::LineDefinition fluid3_tet10;
  fluid3_tet10
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    ;

  DRT::INPUT::LineDefinition fluid3_wedge6;
  fluid3_wedge6
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  DRT::INPUT::LineDefinition fluid3_wedge15;
  fluid3_wedge15
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  DRT::INPUT::LineDefinition fluid3_pyramid5;
  fluid3_pyramid5
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  DRT::INPUT::LineDefinition fluid3_nurbs8;
  fluid3_nurbs8
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("NURBS8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  DRT::INPUT::LineDefinition fluid3_nurbs27;
  fluid3_nurbs27
    .AddInt("ID")
    .AddTag("FLUID3")
    .AddNamedDoubleVector("NURBS27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("FLUID3"));
  lines->Add(fluid3_hex8);
  lines->Add(fluid3_hex20);
  lines->Add(fluid3_hex27);
  lines->Add(fluid3_tet4);
  lines->Add(fluid3_tet10);
  lines->Add(fluid3_wedge6);
  lines->Add(fluid3_wedge15);
  lines->Add(fluid3_pyramid5);
  lines->Add(fluid3_nurbs8);
  lines->Add(fluid3_nurbs27);
  definitions_["FLUID3"] = lines;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidH8Lines()
{
  DRT::INPUT::LineDefinition solid8_a;
  solid8_a
    .AddInt("ID")
    .AddTag("SOLIDH8")
    .AddNamedDoubleVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TET")
    .AddNamedInt("HYB")
    .AddNamedInt("FORM")
    .AddNamedString("STRESSES")
    .AddNamedString("TSI_COUPTYP")
    ;

  DRT::INPUT::LineDefinition solid8_b;
  solid8_b
    .AddInt("ID")
    .AddTag("SOLIDH8")
    .AddNamedDoubleVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("EAS")
    ;

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("SOLIDH8"));
  lines->Add(solid8_a);
  lines->Add(solid8_b);
  definitions_["SOLIDH8"] = lines;
}


#endif
