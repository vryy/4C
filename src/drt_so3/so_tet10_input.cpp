/*!----------------------------------------------------------------------**###
\file so_tet10_input.cpp
\brief

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_tet10.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet10::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);
  
  if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular){
    MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(Material().get());
    holzcard->Setup(NUMGPT_SOTET10, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
    MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(Material().get());
    humcard->Setup(NUMGPT_SOTET10, linedef);
  }

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if(buffer=="linear")
  {
    kintype_ = so_tet10_linear;
    dserror("Reading of SO_TET10 element failed only nonlinear kinematics implemented");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer=="nonlinear")    kintype_ = so_tet10_nonlinear;
  // geometrically non-linear with Updated Lagrangean approach
  else dserror("Reading of SO_TET10 element failed KINEM unknown");

  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
