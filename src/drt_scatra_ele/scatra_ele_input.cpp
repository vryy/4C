/*----------------------------------------------------------------------*/
/*!
\file scatra_element_input.cpp
\brief

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/myocard.H"


/*----------------------------------------------------------------------*
 | read element input (public)                                          |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  SetDisType(DRT::StringToDistype(distype));

  if (Material()->MaterialType() == INPAR::MAT::m_myocard)
      {
      Teuchos::RCP<MAT::Myocard> myocard = Teuchos::rcp_dynamic_cast<MAT::Myocard>(Material());
      myocard->Setup(linedef);
      }
  return true;
}

