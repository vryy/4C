/*!----------------------------------------------------------------------
\file so_nurbs27_input.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "so_nurbs27.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::NURBS::So_nurbs27::ReadElement(const std::string& eletype,
                                                   const std::string& distype,
                                                   DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // read possible gaussian points, obsolete for computation
  std::vector<int> ngp;
  linedef->ExtractIntVector("GP",ngp);
  for (int i=0; i<3; ++i)
    if (ngp[i]!=3)
      dserror("Only version with 3 GP for So_N27 implemented");

  // we expect kintype to be total lagrangian
  kintype_ = sonurbs27_totlag;

  return true;
}
