//-----------------------------------------------------------------------
/*!
\file ale3_input.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE


#include "ale3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Ale3::ReadElement(const std::string& eletype,
                                      const std::string& distype,
                                      DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  return true;
}


#endif
