
/*!----------------------------------------------------------------------
\file artnetexplicitintegration.H
\brief Associated with control routine for artery solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef ARTNETEXPLICITINTEGRATION_H
#define ARTNETEXPLICITINTEGRATION_H

#include "../drt_fluid/fluidimpedancecondition.H"

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_function.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>


using namespace std;
using namespace Teuchos;

/// artery algorithms
namespace ART
{

/*!
\brief time integration for arterial network problems

\author ismail (ismail@lnm.mw.tum.de)
*/
class ArtNetExplicitTimeInt
{
  friend class ArtNetResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  ArtNetExplicitTimeInt(RCP<DRT::Discretization>  dis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output,
                       bool                      alefluid=false);


  /*!
  \brief Destructor

  */
  virtual ~ArtNetExplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop();



  /// setup the variables to do a new time step
  void PrepareTimeStep();


  /*!
  \brief solve linearised artery and bifurcation

  */
  void LinearSolve();

  /*!
    \brief build linear system matrix and rhs


    \param vel new guess at velocity, cross-sectional area, and pressure
  */
  void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);

  /*!
  \brief Update the solution after convergence of the linear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();


  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

//MAH: I think I dont need this part!!!
  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
    int whichinitialfield,
    int startfuncno
    );

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  void EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief calculate pressure

  */
  Teuchos::RCP<Epetra_Vector> CalcPressure();

  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; } // what is this??
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<Epetra_Vector> Areanp() { return areanp_; }
  Teuchos::RCP<Epetra_Vector> Presnp() { return presnp_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<Epetra_Vector> Arean() { return arean_; }
  Teuchos::RCP<Epetra_Vector> Presn() { return presn_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  Teuchos::RCP<LINALG::SparseMatrix> MassMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(massmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> StiffnessMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(stiffmat_); }
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  double Dt() const { return dta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }



  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool splitmatrix=true);


  //@}


protected:

  // don't want = operator and cctor
//  ArtNetExplicitTimeInt operator = (const ArtNetExplicitTimeInt& old);
//  ArtNetExplicitTimeInt (const ArtNetExplicitTimeInt& old);

  //! @name general algorithm parameters
  //! arterial network discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! the processor ID from the communicator
  int                           myrank_;
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;

  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  int     stepmax_; ///< maximal number of timesteps
  double  maxtime_; ///< maximal physical computation time
  //@}


  /// constant density extracted from element material for incompressible flow
  double density_;

  //! use (or not) linearisation of reactive terms on the element
  string newton_;

  /// whether we do some extrapolation at the begin of each time step
  bool extrapolationpredictor_;

  //! @name restart variables
  int uprestart_;
  int upres_;
  int writestresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;

  //@}

  /// (standard) mass matrix
  RCP<LINALG::SparseOperator> massmat_;

  /// (standard) stiffness matrix
  RCP<LINALG::SparseOperator> stiffmat_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;

  /// rhs: right hand side vector for linearised flow solution
  RefCountPtr<Epetra_Vector>    rhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

  //! @name velocities and pressures at time n+1, n and n-1
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;
  //@}

  //! @name pressures at time n+1, n and n-1
  RCP<Epetra_Vector>    presnp_;
  RCP<Epetra_Vector>    presn_;
  RCP<Epetra_Vector>    presnm_;
  //@}

  //! @name cross-sectional area at time n+1, n and n-1
  RCP<Epetra_Vector>    arean_;
  RCP<Epetra_Vector>    areanp_;
  RCP<Epetra_Vector>    areanm_;
  //@}


  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;


  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;

  /// bio related special (outflow) boudaries
  //  RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

  Teuchos::RCP<LINALG::SparseMatrix> meshmatrix_;

  //@}

}; // class ArtNetExplicitTimeInt

}


#endif  // #ifndef ARTNETEXPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
