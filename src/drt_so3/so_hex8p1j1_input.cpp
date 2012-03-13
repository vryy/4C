/*!----------------------------------------------------------------------
\file so_q1p0hex8_input.cpp
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET


#include "so_hex8p1j1.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_Hex8P1J1::ReadElement(const std::string& eletype,
                                             const std::string& distype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_plneohooke)
  {
    MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(Material().get());
    plastic->Setup(NUMGPT_SOH8);
  }
  // temporary variable for read-in
    std::string buffer;
  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  // no linear case implemented so far, hence just a dummy check
  if (buffer=="linear")
  {
    //kintype_ = soh8_linear;
    dserror("Only nonlinear kinematics for SO_HEX8p1j1 implemented!");
  }
  else if (buffer=="nonlinear")
  {
    // kintype_ = soh8_nonlinear;
  }
   else dserror ("Reading SO_HEX8p1j1 element failed KINEM unknown");
  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

