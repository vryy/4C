/*----------------------------------------------------------------------*/
/*!
\file elemag_timeint.H

\brief Base class functions for time integration of electromagnetics

<pre>
\level 3

\maintainer Volker Gravemeier
            gravemeier@lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ELEMAG_TIMINT_H
#define ELEMAG_TIMINT_H

/*----------------------------------------------------------------------*
 | headers                                             gravemeier 06/17 |
 *----------------------------------------------------------------------*/

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>

#include "../drt_inpar/inpar_elemag.H"

/*----------------------------------------------------------------------*
 | forward declarations                                gravemeier 06/17 |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseOperator;
  class MapExtractor;
  class Solver;
}
namespace DRT
{
  class Discretization;
  class DiscretizationHDG;
  class Node;
  class ResultTest;
}
namespace IO
{
  class DiscretizationWriter;
}

/*----------------------------------------------------------------------*
 | general time integration for electromagnetics       gravemeier 06/17 |
 *----------------------------------------------------------------------*/
namespace ELEMAG
{

class ElemagTimeInt
{
public:
  // constructor
  ElemagTimeInt(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
  );

  /// virtual destructor
  virtual ~ElemagTimeInt();

  /// initialization routine
  void Init();

  /// print information to screen
  void PrintInformationToScreen();

  /// time integration
  //virtual void Integrate();
  void Integrate();

  /// set initial field by given function
  //virtual void SetInitialField(int startfuncno);
  void SetInitialField(int startfuncno);

  /// call elements to calculate system matrix
  //virtual void AssembleMatAndRHS();
  void AssembleMatAndRHS();

  /// apply Dirichlet boudnary conditions to system
  void ApplyDirichletToSystem();

  /// solve the system for trace and interior variables
  //virtual void Solve();
  void Solve();

  /// output to screen
  void OutputToScreen();

  /// output routine
  //virtual void Output();
  void Output();

  /// write restart
  //virtual void WriteRestart();

  /// read restart information for given time step
  //virtual void ReadRestart(int step);

  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dtp_;
  }

  double Time()  { return time_; }
  int    Step()  { return step_; }
  int    UpRes() { return upres_; }
  double TimeStep() { return dtp_; }

  /// create result test
  //virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

protected:

  /// discretization, solver, parameter list and output
  Teuchos::RCP<DRT::DiscretizationHDG>   discret_;
  Teuchos::RCP<LINALG::Solver>           solver_;
  Teuchos::RCP<Teuchos::ParameterList>   params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  INPAR::ELEMAG::DynamicType elemagdyna_; /// time integration scheme

  int     myrank_;                        /// processor id

  double  time_;                          /// physical time
  int     step_;                          /// time step
  int     restart_;                       /// restart step

  double  maxtime_;                       /// maximum time
  int     stepmax_;                       /// maximum step

  int     uprestart_;                     /// write restart data every uprestart_ steps
  int     upres_;                         /// write output every upres_ steps

  int numdim_;                            /// number of spatial dimensions
  double dtp_;                            /// time step size

  double  dtele_;                         /// element evaluation time
  double  dtsolve_;                       /// solver time

  bool calcerr_;                          /// flag for error calculation


  /// system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// residual vector
  Teuchos::RCP<Epetra_Vector> residual_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// vector of zeros to be used for enforcing zero Dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector> zeros_;
};

} // namespace ELEMAG

#endif /* ELEMAG_TIMINT_H */
