/*!----------------------------------------------------------------------
\file immersed_problem_dyn.H

\brief global algorithm control file for all immersed algorithms

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 2

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


  /*!
    \brief Ghost the discretization handed in to this method on all procs

    \param distobeghosted (in) : RCP to the discretisation which is ghosted after execution of this method

    \return void

   */
  void CreateGhosting(const Teuchos::RCP<DRT::Discretization> distobeghosted);


  /*!
    \brief Redistribute discretizations by binning

    \param vector_of_discretizations (in) : vector of RCPs to the discretisations which are supposed to be
                                            redistributed by binning.

    \return void

   */
  void CreateBinning(const std::vector<Teuchos::RCP<DRT::Discretization> > vector_of_discretizations);

/*----------------------------------------------------------------------*/
#endif /* #ifndef IMMERSED_PROBLEM_DYN_H*/
