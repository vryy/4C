/*!
\file scatra_ele_impl_utils.cpp

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "scatra_ele_impl_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

namespace SCATRA
{

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool IsBinaryElectrolyte(const std::vector<double>& valence)
{
  int numions(0);
  for (int k=0; k < (int) valence.size(); k++)
  {
    if (abs(valence[k]) > EPS10)
      numions++;
  }
  return (numions == 2);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double CalResDiffCoeff(
    const std::vector<double>& valence,
    const std::vector<double>& diffus,
    const std::vector<double>& diffusvalence
    )
{
  if ((abs(valence[0])<EPS10) or (abs(valence[1])<EPS10))
    dserror("ion species with ids 0 or 1 cannot be neutral.");
  const double n = (diffusvalence[0]-diffusvalence[1]);
  if (abs(n) < EPS15)
    dserror("denominator in resulting diffusion coefficient is zero");

  return diffus[0]*diffus[1]*(valence[0]-valence[1])/n;
}


} // namespace SCATRA

#endif  // #ifdef CCADISCRET
