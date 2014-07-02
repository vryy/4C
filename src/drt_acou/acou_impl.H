/*----------------------------------------------------------------------*/
/*!
\file acou_impltimint.H

\brief Acoustics Base Algorithm for implicit time integration

       Derived classes are implicit Euler, trapezoidal rule, BDF2,
       BDF3, BDF4 and several diagonally implicit RK-schemes

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ACOU_IMPLTIMINT_H
#define ACOU_IMPLTIMINT_H


/*----------------------------------------------------------------------*
 | headers                                               schoeder 01/14 |
 *----------------------------------------------------------------------*/

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>

#include "../drt_inpar/inpar_acou.H"

/*----------------------------------------------------------------------*
 | forward declarations                                  schoeder 01/14 |
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
 | acoustics algorithm                                   schoeder 01/14 |
 *----------------------------------------------------------------------*/
namespace ACOU
{
// acoustics solver
class AcouImplicitTimeInt
{
public:

  /// constructor
  AcouImplicitTimeInt(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
  );

  /// virtual destructor
  virtual ~AcouImplicitTimeInt();

  /// read restart information for given time step
  virtual void ReadRestart(int step);

  /// set initial field to zero (needed for adjoint problem)
  virtual void SetInitialZeroField();

  /// set initial field by given function
  virtual void SetInitialField(int startfuncno, double pulse);

  /// set initial field by given light distribution
  virtual void SetInitialPhotoAcousticField(double pulse, Teuchos::RCP<Epetra_Vector> light, Teuchos::RCP<DRT::Discretization> scatradis, bool meshconform);

  /// print out user information
  void PrintInformationToScreen();

  /// run a complete simulation
  virtual void Integrate(Teuchos::RCP<Epetra_MultiVector> history = Teuchos::null, Teuchos::RCP<LINALG::MapExtractor> splitter = Teuchos::null);

  /// solve the system for trace and interior variables
  virtual void Solve();

  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dtp_;
  }

  /// call elements to calculate system matrix
  virtual void AssembleMatAndRHS();

  /// apply Dirichlet boudnary conditions to system
  void ApplyDirichletToSystem();

  /// update: current solution becomes old solution of next timestep
  virtual void TimeUpdate();

  /// output
  void Output(Teuchos::RCP<Epetra_MultiVector> history = Teuchos::null, Teuchos::RCP<LINALG::MapExtractor> splitter = Teuchos::null);

  /// write restart
  virtual void WriteRestart();

  /// function to output user information
  void OutputToScreen();

  /// evaluate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

  /// Call elements to update interior variabels and to calculate the RHS vector (in this way, we can save some computing time)
  virtual void UpdateInteriorVariablesAndAssemebleRHS();

  /// Output nodal internal pressure field at new time step
  void NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec);

  /// Output nodal internal psi field at new time step
  void NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec);

  /// Fill Touch Count Vector (needed for inverse analysis)
  void FillTouchCountVec(Teuchos::RCP<Epetra_Vector> touchcount);

  /// For result tests
  Teuchos::RCP<DRT::Discretization> Discretization();
  Teuchos::RCP<Epetra_Vector> Velnp()   { return velnp_; }


  /// create result test for encapsulated structure algorithm
  virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  /// return the name of the time integration scheme
  virtual std::string Name() = 0;


protected:

  /// discretization, solver, parameter list and output
  Teuchos::RCP<DRT::DiscretizationHDG>   discret_;
  Teuchos::RCP<LINALG::Solver>           solver_;
  Teuchos::RCP<Teuchos::ParameterList>   params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  bool    adjoint_;                     /// flag indicating whether this is standard or an adjoint solve

  int     myrank_;                      /// processor id

  double  time_;                        /// physical time
  int     step_;                        /// time step
  int     restart_;                     /// restart step

  double  maxtime_;                     /// maximum time
  int     stepmax_;                     /// maximum step

  int     uprestart_;                   /// write restart data every uprestart_ steps
  int     upres_;                       /// write output every upres_ steps

  INPAR::ACOU::DynamicType dyna_;       /// time integration scheme
  INPAR::ACOU::PhysicalType phys_;      /// physical type

  int numdim_;                          /// number of spatial dimensions
  double dtp_;                          /// time step size

  /// cpu-time measures
  double  dtele_;                       /// element evaluation time
  double  dtsolve_;                     /// solver time

  bool invana_;                         /// flag indicating whether inverse analysis is performed
  bool errormaps_;                      /// flag indicating whether error maps are calculated
  bool padaptivity_;                    /// p-adaptivity flag
  double padapttol_;                    /// p-adaptivity error tolerance
  INPAR::ACOU::CalcError calcerr_;      /// flag for error calculation

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// needed in case of adjoint problem solve
  Teuchos::RCP<Epetra_MultiVector>    adjoint_rhs_;

  //! @name velocity at time n+1, n, and n-1
  Teuchos::RCP<Epetra_Vector>    velnp_;
  Teuchos::RCP<Epetra_Vector>    veln_;
  Teuchos::RCP<Epetra_Vector>    velnm_;

  //! @name velocity and pressure at time n+1, n, n-1
  //!  and n+alpha_F for element interior in HDG
  Teuchos::RCP<Epetra_Vector>    intvelnp_;
  Teuchos::RCP<Epetra_Vector>    intveln_;
  Teuchos::RCP<Epetra_Vector>    intvelnm_;

  /// element based error vector, only if desired by user
  Teuchos::RCP<Epetra_Vector>    error_;
};

} // namespace ACOU

#endif /* ACOU_IMPLTIMINT_H */
