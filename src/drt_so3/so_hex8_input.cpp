/*!----------------------------------------------------------------------
\file so_hex8_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "../drt_mat/robinson.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);

  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)

  Teuchos::RCP<MAT::Material> mat = Material();

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
   std::string buffer;

  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
   kintype_ = INPAR::STR::kinem_linear;
  }
  else if (buffer=="nonlinear")
  {
   kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else dserror ("Reading SO_HEX8 element failed");

  // check if material kinematics is compatible to element kinematics
  so3mat->ValidKinematics(kintype_);

  // read EAS technology flag
  linedef->ExtractString("EAS",buffer);

  // full EAS technology
  if (buffer=="full")
  {
    eastype_ = soh8_easfull;
    neas_ = 21;               // number of eas parameters for full EAS
    soh8_easinit();
  }
  // mild EAS technology
  else if (buffer=="mild")
  {
    eastype_ = soh8_easmild;
    neas_ = 9;               // number of eas parameters for mild EAS
    soh8_easinit();
  }
  // no EAS technology
  else if (buffer=="none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else
    dserror("Reading of SO_HEX8 EAS technology failed");

  return true;
}

