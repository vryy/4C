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

#include <Teuchos_RefCountPtr.hpp>

class Epetra_Vector;

namespace DRT
{
  class Discretization;
}

namespace Teuchos
{
  class ParameterList;
}

namespace LINALG
{
  class SparseOperator;
}

namespace PATSPEC
{

  void PatientSpecificGeometry(DRT::Discretization& dis);

  void ComputeEleNormalizedLumenDistance(DRT::Discretization& dis);

  void ComputeEleLocalRadius(DRT::Discretization& dis);

  void GetILTDistance(const int eleid, 
                      Teuchos::ParameterList& params, 
                      DRT::Discretization& dis);

  void GetLocalRadius(const int eleid, 
                      Teuchos::ParameterList& params, 
                      DRT::Discretization& dis);
                      
  void CheckEmbeddingTissue(DRT::Discretization& discret,
                            Teuchos::RCP<LINALG::SparseOperator> stiff,
                            Teuchos::RCP<Epetra_Vector> fint);

}
#endif /* PATSPEC_H */
