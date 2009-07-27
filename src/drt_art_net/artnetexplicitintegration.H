
/*!----------------------------------------------------------------------
\file artnetexplicitintegration.H
\brief Associated with control routine for artery solvers,

     including instationary solvers based on

     o two-step Taylor-Galerkin


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

//#include "../drt_fluid/fluidimpedancecondition.H"

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
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ml_common.h>
#include <ml_include.h>
#include <ml_epetra_utils.h>
#include <ml_epetra.h>
#include <ml_epetra_operator.h>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_agg_genP.h>
#include <ml_operator.h>
#include <MLAPI_Error.h>
#include <MLAPI_CompObject.h>
#include <MLAPI_TimeObject.h>
#include <MLAPI_Operator.h>
#include <MLAPI_Operator_Utils.h>
#include <MLAPI_MultiVector.h>
#include <MLAPI_InverseOperator.h>
#include <MLAPI_Expressions.h>
#include <MLAPI_BaseOperator.h>
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>
#include <MLAPI_Eig.h>

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

using namespace std;
using namespace Teuchos;
using namespace MLAPI;

/// artery algorithms
namespace ART
{

/*!
\brief time integration for arterial network problems

\author ismail (ismail@lnm.mw.tum.de)
*/
class ArtNetExplicitTimeInt
{
  //friend class ArtNetResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  ArtNetExplicitTimeInt(RCP<DRT::Discretization>  dis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output);


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
  void Solve();

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

  /*!
  \brief Assembling of the RHS Vector and the LHS Matrix 

  */
  void AssembleMatAndRHS();


  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  //void EvaluateErrorComparedToAnalyticalSol();


  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  //! @name access methods for composite algorithms

  //  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; } //This variable might be needed in future!
  Teuchos::RCP<Epetra_Vector> Qnp()    { return qnp_; }
  Teuchos::RCP<Epetra_Vector> QAnp()   { return qanp_; }
  Teuchos::RCP<Epetra_Vector> Areanp() { return areanp_; }
  //Teuchos::RCP<Epetra_Vector> Presnp() { return presnp_; }
  Teuchos::RCP<Epetra_Vector> Qn()     { return qn_; }
  Teuchos::RCP<Epetra_Vector> QAn()    { return qan_; }
  Teuchos::RCP<Epetra_Vector> Arean()  { return arean_; }
  //Teuchos::RCP<Epetra_Vector> Presn()  { return presn_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

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
  Teuchos::RCP<DRT::Discretization > Discretization() { return discret_; }

  double  Dt()   const { return dta_; }

  double  Time() const { return time_; }
  int     Step() const { return step_; }

  int     Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void    SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  //@}


protected:

  //! @name general algorithm parameters
  //! arterial network discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! the processor ID from the communicator
  int                           myrank_;

  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  int     stepmax_; ///< maximal number of timesteps
  double  maxtime_; ///< maximal physical computation time
  //@}


  /// constant density extracted from element material for incompressible flow
  double density_;

  //! @name restart variables
  int uprestart_;
  int upres_;
  int writestresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;

  //@}

  /// cpu-time measures
  double  dtele_;
  double  dtfilter_;
  double  dtsolve_;

  /// (standard) mass matrix
  RCP<LINALG::SparseOperator> massmat_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// rhs: right hand side vector
  RefCountPtr<Epetra_Vector>    rhs_;

  //! @name Volumetric Flow rate and Cross-Sectional area at time n+1, n and n-1
  RCP<Epetra_Vector>    qanp_;
  RCP<Epetra_Vector>    qan_;
  RCP<Epetra_Vector>    qanm_;
  //@}

  //! @name Volumetric Flow rate at time n+1, n and n-1
  RCP<Epetra_Vector>    qnp_;
  RCP<Epetra_Vector>    qn_;
  RCP<Epetra_Vector>    qnm_;
  //@}

  //! @name cross-sectional area at time n+1, n and n-1
  RCP<Epetra_Vector>    arean_;
  RCP<Epetra_Vector>    areanp_;
  RCP<Epetra_Vector>    areanm_;
  //@}

  //! @name Dirichlet boundary condition vectors
  RCP<Epetra_Vector>    bcval_;
  RCP<Epetra_Vector>    dbctog_;
  //@}

  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;

  /// bio related special (outflow) boudaries
  //  RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  //@}

}; // class ArtNetExplicitTimeInt

}


#endif  // #ifndef ARTNETEXPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
