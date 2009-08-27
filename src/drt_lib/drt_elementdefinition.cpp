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
  ed.PrintElementLines(std::cout,"SHELL8");
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
//   SetupArtLines();
//   SetupBeam2Lines();
//   SetupBeam2rLines();
//   SetupBeam3Lines();
//   SetupConstrele2Lines();
//   SetupConstrele3Lines();
//   SetupPtet4Lines();
//   SetupShell8Lines();
//   SetupSolid3Lines();
//   SetupSolidh20Lines();
//   SetupSolidh27Lines();
//   SetupSolidh8Lines();
//   SetupSolidh8p1j1Lines();
//   SetupSolidsh8Lines();
//   SetupSolidsh8p8Lines();
//   SetupSolidshw6Lines();
//   SetupSolidt10Lines();
//   SetupSolidt4Lines();
//   SetupSolidw6Lines();
//   SetupTorsion2Lines();
//   SetupTorsion3Lines();
//   SetupTruss2Lines();
//   SetupTruss3Lines();
//   SetupWallLines();

//   SetupCombust3Lines();
//   SetupFluid2Lines();
  SetupFluid3Lines();
//   SetupTranspLines();
//   SetupXdiff3Lines();
//   SetupXfluid3Lines();

//   SetupAle2Lines();
//   SetupAle3Lines();

  // backward compatibility
  // still needed?
  //SetupCondif2Lines();
  //SetupCondif3Lines();
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
void DRT::INPUT::ElementDefinition::SetupArtLines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["ART"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupBeam2Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["BEAM2"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupBeam2rLines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["BEAM2R"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN2"]
    .AddDoubleVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LINE3"]
    .AddDoubleVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN3"]
    .AddDoubleVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LINE4"]
    .AddDoubleVector("LINE4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN4"]
    .AddDoubleVector("LIN4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LINE5"]
    .AddDoubleVector("LINE5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN5"]
    .AddDoubleVector("LIN5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupBeam3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["BEAM3"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN2"]
    .AddDoubleVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE3"]
    .AddDoubleVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN3"]
    .AddDoubleVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE4"]
    .AddDoubleVector("LINE4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN4"]
    .AddDoubleVector("LIN4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LINE5"]
    .AddDoubleVector("LINE5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;

  defs["LIN5"]
    .AddDoubleVector("LIN5",5)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("MOMIN")
    //.AddNamedDouble("MOMIN")
    .AddNamedDouble("MOMINPOL")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupConstrele2Lines()
{
  // No reading for this element! Will be created on the fly, not from a .dat file.
  //std::map<std::string,LineDefinition>& defs = definitions_["CONSTRELE2"];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupConstrele3Lines()
{
  // No reading for this element! Will be created on the fly, not from a .dat file.
  //std::map<std::string,LineDefinition>& defs = definitions_["CONSTRELE3"];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupPtet4Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["PTET4"];

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupShell8Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SHELL8"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedString("SDC")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedString("SDC")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedString("SDC")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedString("SDC")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedString("SDC")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolid3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLID3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidh20Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDH20"];

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedDouble("STRENGTH")
    .AddNamedIntVector("GP",3)
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidh27Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDH27"];

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedDouble("STRENGTH")
    .AddNamedIntVector("GP",3)
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidh8Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDH8"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedString("KINEM")
    .AddNamedString("EAS")
    .AddNamedDoubleVector("RAD",3)
    .AddNamedDoubleVector("AXI",3)
    .AddNamedDoubleVector("CIR",3)
    .AddNamedDouble("STRENGTH")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidh8p1j1Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDH8P1J1"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidsh8Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDSH8"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedString("KINEM")
    .AddNamedString("EAS")
    .AddNamedString("THICKDIR")
    .AddNamedDoubleVector("RAD",3)
    .AddNamedDoubleVector("AXI",3)
    .AddNamedDoubleVector("CIR",3)
    .AddNamedDouble("STRENGTH")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidsh8p8Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDSH8P8"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedString("KINEM")
    .AddNamedString("THICKDIR")
    .AddNamedString("STAB")
    .AddNamedString("ANS")
    .AddNamedString("EAS")
    .AddNamedString("LIN")
    .AddNamedString("ISO")
    .AddNamedDoubleVector("RAD",3)
    .AddNamedDoubleVector("AXI",3)
    .AddNamedDoubleVector("CIR",3)
    .AddNamedDouble("STRENGTH")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidshw6Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDSHW6"];

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddNamedString("EAS")
    .AddNamedString("OPTORDER")
    .AddNamedDoubleVector("RAD",3)
    .AddNamedDoubleVector("AXI",3)
    .AddNamedDoubleVector("CIR",3)
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidt10Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDT10"];

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidt4Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDT4"];

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupSolidw6Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["SOLIDW6"];

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddNamedDoubleVector("RAD",3)
    .AddNamedDoubleVector("AXI",3)
    .AddNamedDoubleVector("CIR",3)
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupTorsion2Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["TORSION2"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;

  defs["LIN2"]
    .AddDoubleVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupTorsion3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["TORSION3"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;

  defs["LIN2"]
    .AddDoubleVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupTruss2Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["TRUSS2"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;

  defs["LIN2"]
    .AddDoubleVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupTruss3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["TRUSS3"];

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;

  defs["LIN2"]
    .AddDoubleVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupWallLines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["WALL"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;

  defs["NURBS4"]
    .AddIntVector("NURBS4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;

  defs["NURBS9"]
    .AddIntVector("NURBS9",9)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",2)
    .AddString("STRESS_STRAIN")
    .AddNamedString("EAS")
    .AddNamedString("STRESSES")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupCombust3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["COMBUST3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupCondif2Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["CONDIF2"];

  // backward compatibility
  defs = definitions_["TRANSP"];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupCondif3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["CONDIF3"];

  // backward compatibility
  defs = definitions_["TRANSP"];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupFluid2Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["FLUID2"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["NURBS4"]
    .AddIntVector("NURBS4",4)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["NURBS9"]
    .AddIntVector("NURBS9",9)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;

  defs["THQ9"]
    .AddDoubleVector("THQ9",9)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",2)
    .AddNamedInt("GP_TRI")
    .AddNamedString("GP_ALT")
    .AddNamedString("NA")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupFluid3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["FLUID3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS8"]
    .AddIntVector("NURBS8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupTranspLines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["TRANSP"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    ;

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    ;

  defs["NURBS4"]
    .AddIntVector("NURBS4",4)
    .AddNamedInt("MAT")
    ;

  defs["NURBS9"]
    .AddIntVector("NURBS9",9)
    .AddNamedInt("MAT")
    ;

  defs["LINE2"]
    .AddDoubleVector("LINE2",2)
    .AddNamedInt("MAT")
    ;

  defs["LINE3"]
    .AddDoubleVector("LINE3",3)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupXdiff3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["XDIFF3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    ;

  defs["NURBS8"]
    .AddIntVector("NURBS8",8)
    .AddNamedInt("MAT")
    ;

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupXfluid3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["XFLUID3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    ;

  defs["NURBS8"]
    .AddIntVector("NURBS8",8)
    .AddNamedInt("MAT")
    ;

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupAle2Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["ALE2"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::ElementDefinition::SetupAle3Lines()
{
  std::map<std::string,LineDefinition>& defs = definitions_["ALE3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    ;
}


#endif
