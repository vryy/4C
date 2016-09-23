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

  //! global control routine for all immersed problems.
  void immersed_problem_drt();

  //! global cell migration algorithm
  void CellMigrationControlAlgorithm();

/*----------------------------------------------------------------------*/
#endif /* #ifndef IMMERSED_PROBLEM_DYN_H*/
