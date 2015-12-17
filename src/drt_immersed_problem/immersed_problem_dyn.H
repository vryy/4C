/*!----------------------------------------------------------------------
\file immersed_problem_dyn.H

\brief global algorithm control class for all immersed algorithms

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PROBLEM_DYN_H
#define IMMERSED_PROBLEM_DYN_H

#include <Teuchos_RCP.hpp>

namespace DRT{
  class Discretization;
}


  //! global control routine for all immersed problems.
  void immersed_problem_drt();

  //! global cell migration algorithm
  void CellMigrationControlAlgorithm();

  //! Ghost the discretization handed in to this method on all procs
  void CreateGhosting(const Teuchos::RCP<DRT::Discretization> distobeghosted);

/*----------------------------------------------------------------------*/
#endif /* #ifndef IMMERSED_PROBLEM_DYN_H*/
