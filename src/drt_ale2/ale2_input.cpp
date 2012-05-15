//-----------------------------------------------------------------------
/*!
\file ale2_input.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE


#include "ale2.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Ale2::ReadElement(const std::string& eletype,
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
