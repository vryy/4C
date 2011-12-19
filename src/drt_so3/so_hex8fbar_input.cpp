/*!----------------------------------------------------------------------
\file so_hex8fbar_input.cpp
\brief

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8fbar.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8fbar::ReadElement(const std::string& eletype,
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
  else if (Material()->MaterialType() == INPAR::MAT::m_growth)
  {
    MAT::Growth* grow = static_cast <MAT::Growth*>(Material().get());
    grow->Setup(NUMGPT_SOH8, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture){
    MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(Material().get());
    comix->Setup(NUMGPT_SOH8, linedef);
  }

  // temporary variable for read-in
  std::string buffer;

	// read kinematic flag
	linedef->ExtractString("KINTYP",buffer);
	if (buffer=="lin")
	{
		kintype_ = soh8_geolin;
		dserror("Only Total Lagrange for SO_HEX8FBAR implemented!");
	}
	else if (buffer=="nln")
	{
		kintype_ = soh8_totlag;
	}
	else dserror ("Reading SO_HEX8FBAR element failed");

  return true;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
