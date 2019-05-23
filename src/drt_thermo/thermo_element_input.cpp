/*----------------------------------------------------------------------*/
/*!
\file thermo_element_input.cpp
\brief element input routines
\level 1
\maintainer Christoph Meier
*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifdef D_THERMO

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "thermo_element.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/fouriervar.H"


/*----------------------------------------------------------------------*
 | read element and set required information                  gjb 01/08 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Thermo::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SetDisType(DRT::StringToDistype(distype));

  if (Shape() == DRT::Element::nurbs27) SetNurbsElement() = true;

  if (Material()->MaterialType() == INPAR::MAT::m_th_fourier_var)
  {
    if (Shape() != DRT::Element::hex8)
      dserror("Setup call for FourierVar only implemented for hex8.");
    Teuchos::RCP<MAT::FourierVar> mat = Teuchos::rcp_dynamic_cast<MAT::FourierVar>(Material());
    mat->Setup(8);
  }

  return true;
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef D_THERMO
