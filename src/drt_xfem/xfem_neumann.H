/*!-----------------------------------------------------------------------------------------------*
\file xfem_neumann.cpp

\brief base XFEM Neumann boundary conditions

<pre>
Maintainer: Dipl.-Math. Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


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
  class FluidWizard;
}

namespace XFEM
{

/// evaluate Neumann boundary conditions
void EvaluateNeumann(  Teuchos::RCP<XFEM::FluidWizard>      wizard,
                       Teuchos::ParameterList &             params,
                       DRT::Discretization &                discret,
                       Teuchos::RCP<DRT::Discretization>    cutdiscret,
                       Teuchos::RCP<Epetra_Vector>          systemvector,
                       Teuchos::RCP<LINALG::SparseOperator> systemmatrix = Teuchos::null );

/// evaluate Neumann boundary conditions
void EvaluateNeumann(  Teuchos::RCP<XFEM::FluidWizard>      wizard,
                       Teuchos::ParameterList &             params,
                       DRT::Discretization &                discret,
                       Teuchos::RCP<DRT::Discretization>    cutdiscret,
                       Epetra_Vector&                       systemvector,
                       LINALG::SparseOperator*              systemmatrix = NULL );

/// evaluate standard Neumann boundary conditions
void EvaluateNeumannStandard( std::multimap<std::string,DRT::Condition* >&  condition,
                              bool                                          usetime,
                              const double                                  time,
                              bool                                          assemblemat,
                              Teuchos::ParameterList &                      params,
                              DRT::Discretization &                         discret,
                              Epetra_Vector&                                systemvector,
                              LINALG::SparseOperator*                       systemmatrix);

void CutNeumannSurf(
    Teuchos::RCP<DRT::Element> neumann_surface,
    DRT::Element* parentele, GEO::CUT::ElementHandle * parentele_handle,
    bool & eval_Neumann,
    bool & xfem_eval_Neumann);

/// evaluate XFEM Neumann boundary conditions
void EvaluateNeumannXFEM(     Teuchos::RCP<XFEM::FluidWizard>               wizard,
                              std::multimap<std::string,DRT::Condition* >&  condition,
                              bool                                          usetime,
                              const double                                  time,
                              bool                                          assemblemat,
                              Teuchos::ParameterList &                      params,
                              DRT::Discretization &                         discret,
                              Teuchos::RCP<DRT::Discretization>             cutdiscret,
                              Epetra_Vector&                                systemvector,
                              LINALG::SparseOperator*                       systemmatrix);

int getParentElementId(       DRT::Discretization&        discret,
                              Teuchos::RCP<DRT::Element>  surf_ele  );



}

#endif /* XFEM_NEUMANN_H_ */
