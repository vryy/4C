/*!----------------------------------------------------------------------
\file patspec.H

\brief A collection of methods to modify patient specific geometries

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifndef PATSPEC_H
#define PATSPEC_H

//#include "../drt_lib/drt_discret.H"

namespace DRT
{
  class Discretization;
}

namespace Teuchos
{
  class ParameterList;
}

namespace PATSPEC
{

  void PatientSpecificGeometry(DRT::Discretization& dis);

  void ComputeEleNormalizedLumenDistance(DRT::Discretization& dis);

  void GetILTDistance(const int eleid, 
                      Teuchos::ParameterList& params, 
                      DRT::Discretization& dis);

}
#endif /* PATSPEC_H */
