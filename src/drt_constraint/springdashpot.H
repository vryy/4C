/*!----------------------------------------------------------------------
\file patspec.H

\brief Methods for spring and dashpot constraints / boundary conditions

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#ifndef SPRINGDASHPOT_H
#define SPRINGDASHPOT_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

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

namespace IO
{
  class DiscretizationWriter;
}

namespace SPRINGDASHPOT
{

  void SpringDashpot(Teuchos::RCP<DRT::Discretization> dis);


  void EvaluateSpringDashpot(Teuchos::RCP<DRT::Discretization> discret,
                            Teuchos::RCP<LINALG::SparseOperator> stiff,
                            Teuchos::RCP<Epetra_Vector> fint,
                            Teuchos::RCP<Epetra_Vector> disn,
                            Teuchos::RCP<Epetra_Vector> veln,
                            Teuchos::ParameterList parlist);




}
#endif /* SPRINGDASHPOT_MANAGER_H */
