/*!----------------------------------------------------------------------
\file so_hex27_input.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#include "so_hex27.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex27::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(NUMGPT_SOH27, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  if (buffer=="linear")
  {
    kintype_ = soh27_linear;
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = soh27_nonlinear;
  }
  else dserror ("Reading SO_HEX27 element failed KINEM unknown");

  // only for linear SVK materials and small strain plastic materials
  bool admissibl_mat = false;
  if ( (mat->MaterialType() == INPAR::MAT::m_stvenant)
       or (mat->MaterialType() == INPAR::MAT::m_thermostvenant)
       or (mat->MaterialType() == INPAR::MAT::m_pllinelast)
       or (mat->MaterialType() == INPAR::MAT::m_thermopllinelast)
       or (mat->MaterialType() == INPAR::MAT::m_elpldamage) )
    admissibl_mat = true;

  // check for SVK material if geometrically linear
  if ( (kintype_ == soh27_linear) and (admissibl_mat == false) )
    dserror("ERROR: Only linear elasticity (SVK) for geometrically linear hex27 element");

  return true;
}


