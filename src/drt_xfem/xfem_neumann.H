/*----------------------------------------------------------------------*/
/*!
\file xfem_neumann.H

\brief base XFEM Neumann boundary conditions

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning think about removing these routines!!!

*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_NEUMANN_H_
#define XFEM_NEUMANN_H_


#include <Teuchos_RCP.hpp>


class Epetra_Vector;

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
}

namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
  }
}

namespace LINALG
{
  class SparseOperator;
}
namespace Teuchos
{
  class ParameterList;
}


namespace XFEM
{

/// evaluate Neumann boundary conditions
void EvaluateNeumann(  Teuchos::RCP<GEO::CutWizard>         wizard,
                       Teuchos::ParameterList &             params,
                       Teuchos::RCP<DRT::Discretization>    discret,
                       Teuchos::RCP<Epetra_Vector>          systemvector,
                       Teuchos::RCP<LINALG::SparseOperator> systemmatrix = Teuchos::null );

/// evaluate Neumann boundary conditions
void EvaluateNeumann(  Teuchos::RCP<GEO::CutWizard>         wizard,
                       Teuchos::ParameterList &             params,
                       Teuchos::RCP<DRT::Discretization>    discret,
                       Epetra_Vector&                       systemvector,
                       LINALG::SparseOperator*              systemmatrix = NULL );

/// evaluate standard Neumann boundary conditions
void EvaluateNeumannStandard( std::multimap<std::string,DRT::Condition* >&  condition,
                              bool                                          usetime,
                              const double                                  time,
                              bool                                          assemblemat,
                              Teuchos::ParameterList &                      params,
                              Teuchos::RCP<DRT::Discretization>             discret,
                              Epetra_Vector&                                systemvector,
                              LINALG::SparseOperator*                       systemmatrix);

void CutNeumannSurf(
    Teuchos::RCP<DRT::Element> neumann_surface,
    DRT::Element* parentele, GEO::CUT::ElementHandle * parentele_handle,
    bool & eval_Neumann,
    bool & xfem_eval_Neumann);

/// evaluate XFEM Neumann boundary conditions
void EvaluateNeumannXFEM(     Teuchos::RCP<GEO::CutWizard>                  wizard,
                              std::multimap<std::string,DRT::Condition* >&  condition,
                              bool                                          usetime,
                              const double                                  time,
                              bool                                          assemblemat,
                              Teuchos::ParameterList &                      params,
                              Teuchos::RCP<DRT::Discretization>             discret,
                              Epetra_Vector&                                systemvector,
                              LINALG::SparseOperator*                       systemmatrix);

int getParentElementId(       Teuchos::RCP<DRT::Discretization>  discret,
                              Teuchos::RCP<DRT::Element>  surf_ele  );



}

#endif /* XFEM_NEUMANN_H_ */
