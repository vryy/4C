/*!----------------------------------------------------------------------**###
\file so_nstet_input.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nstet.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::NStet::ReadElement(const std::string& eletype,
                                      const std::string& distype,
                                      DRT::INPUT::LineDefinition* linedef)
{
  string buffer;
  
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractString("STAB",buffer);
  if (buffer=="none")         stabtype_ = DRT::ELEMENTS::NStet::so_tet4_stab_none;
  else if (buffer=="voldev" || buffer=="devvol")  
                              stabtype_ = DRT::ELEMENTS::NStet::so_tet4_voldev;
  else if (buffer=="vol")     stabtype_ = DRT::ELEMENTS::NStet::so_tet4_vol;
  else if (buffer=="dev")     stabtype_ = DRT::ELEMENTS::NStet::so_tet4_dev;
  else if (buffer=="puso")    stabtype_ = DRT::ELEMENTS::NStet::so_tet4_puso;
  else dserror("Unknown type of stabilization for NStet: {VolDev,Vol,Dev,Puso}");

  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
