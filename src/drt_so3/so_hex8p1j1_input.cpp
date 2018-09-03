/*!----------------------------------------------------------------------
\file so_q1p0hex8_input.cpp
\brief

<pre>
Maintainer: Lena Wiechert
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/


#include "so_hex8p1j1.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_Hex8P1J1::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->Setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
  std::string buffer;
  // read kinematic flag
  linedef->ExtractString("KINEM", buffer);
  // no linear case implemented so far, hence just a dummy check
  if (buffer == "linear")
  {
    // kintype_ = soh8_linear;
    dserror("Only nonlinear kinematics for SO_HEX8p1j1 implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else
    dserror("Reading SO_HEX8p1j1 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}
