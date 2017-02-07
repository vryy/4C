/*!----------------------------------------------------------------------
\file artery_input.cpp
\brief

\level 3

\maintainer Lena Yoshihara

*----------------------------------------------------------------------*/
#ifdef D_ARTNET

#include "artery.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Artery::ReadElement(const std::string& eletype,
                                        const std::string& distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  int ngp;
  linedef->ExtractInt("GP",ngp);

  switch (ngp)
  {
  case 1:
    gaussrule_ = DRT::UTILS::intrule_line_1point;
    break;
  case 2:
    gaussrule_ = DRT::UTILS::intrule_line_2point;
    break;
  case 3:
    gaussrule_ = DRT::UTILS::intrule_line_3point;
    break;
  case 4:
    gaussrule_ = DRT::UTILS::intrule_line_4point;
    break;
  case 5:
    gaussrule_ = DRT::UTILS::intrule_line_5point;
    break;
  case 6:
    gaussrule_ = DRT::UTILS::intrule_line_6point;
    break;
  case 7:
    gaussrule_ = DRT::UTILS::intrule_line_7point;
    break;
  case 8:
    gaussrule_ = DRT::UTILS::intrule_line_8point;
    break;
  case 9:
    gaussrule_ = DRT::UTILS::intrule_line_9point;
    break;
  case 10:
    gaussrule_ = DRT::UTILS::intrule_line_10point;
    break;
  default:
    dserror("Reading of ART element failed: Gaussrule for line not supported!\n");
  }

  return true;
}

#endif  // #ifdef D_ARTNET
