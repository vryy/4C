
/*!----------------------------------------------------------------------
\file airwayimplicitintegration.H
\brief Associated with control routine for reduced dimensional airways
  solvers,

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/

#ifndef AIRWAYIMPLICITINTEGRATION_H
#define AIRWAYIMPLICITINTEGRATION_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */


#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
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

// forward declarations
namespace DRT
{
  class ResultTest;
  class RedAirwayResultTest;
}


/// airway algorithms
namespace AIRWAY
{

/*!
\brief time integration for reduced dimesnional airway network problems

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAirwayImplicitTimeInt
{

public:


  /*!
  \brief Standard Constructor

  */
  RedAirwayImplicitTimeInt(RCP<DRT::Discretization>  dis,
                           LINALG::Solver&            solver,
                           ParameterList&             params,
                           IO::DiscretizationWriter&  output);


  /*!
  \brief Destructor

  */
  virtual ~RedAirwayImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate(bool CoupledTo3D, RCP<ParameterList> CouplingParams);

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop(bool CoupledTo3D, RCP<ParameterList> CouplingParams);

  /*!
  \brief Do one step time integration (time loop)

  */
  void OneStepTimeLoop(bool CoupledTo3D = false, RCP<ParameterList> CouplingParams = Teuchos::null);


  /*!
  \brief Integrate one step

  */
  void IntegrateStep(RCP<ParameterList> CouplingTo3DParams = Teuchos::null);


  /// setup the variables to do a new time step
  void PrepareTimeStep();


  /*!
  \brief solve of airways

  */
  void Solve(Teuchos::RCP<ParameterList> CouplingTo3DParams);

  /*!
  \brief nonlinear solver of airways

  */
  void NonLin_Solve(Teuchos::RCP<ParameterList> CouplingTo3DParams);


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
  void Output  (bool               CoupledTo3D = false,
                RCP<ParameterList> CouplingParams = Teuchos::null);

  /*!
  \brief Assembling of the RHS Vector and the LHS Matrix

  */
  void AssembleMatAndRHS();

  /*!
  \brief Evaluate the error residual
  */
  void EvalResidual(Teuchos::RCP<ParameterList> CouplingTo3DParams);


  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  //void EvaluateErrorComparedToAnalyticalSol();


  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();


  //! @name access methods for composite algorithms

  Teuchos::RCP<Epetra_Vector> Pnp(){ return pnp_; }
  Teuchos::RCP<Epetra_Vector> Pn() { return pn_; }
  Teuchos::RCP<Epetra_Vector> Pnm(){ return pnm_; }

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


  //! @name methods related to coupling with 3D tissue models

  void SetupForCoupling();
  void SetAirwayFluxFromTissue(Teuchos::RCP<Epetra_Vector> coupflux);
  void ExtractPressure(Teuchos::RCP<Epetra_Vector> couppres);

  //@}

protected:

  //! @name general algorithm parameters
  //! reduced dimensional airway network discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  ParameterList                 params_;
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
  RCP<LINALG::MapExtractor>   dbcmaps_;

  /// rhs: right hand side vector
  RefCountPtr<Epetra_Vector>  rhs_;


  //! @name pressures at time n+1, n and n-1
  RCP<Epetra_Vector>          pnp_;
  RCP<Epetra_Vector>          pn_;
  RCP<Epetra_Vector>          pnm_;
  RCP<Epetra_Vector>          p_nonlin_;

  RCP<Epetra_Vector>          sysmat_iad_;
  //@}

  //! @name inlet volumetric flow rates at time n+1, n and n-1
  RCP<Epetra_Vector>          qin_np_;
  RCP<Epetra_Vector>          qin_n_;
  RCP<Epetra_Vector>          qin_nm_;

  RCP<Epetra_Vector>          qi_nl_np_;
  //@}

  //! @name outlet volumetric flow rates at time n+1, n and n-1
  RCP<Epetra_Vector>          qout_np_;
  RCP<Epetra_Vector>          qout_n_;
  RCP<Epetra_Vector>          qout_nm_;

  RCP<Epetra_Vector>          qo_nl_np_;

  RCP<Epetra_Vector>          qexp_;
  RCP<Epetra_Vector>          pexp_;
  //@}

  //! @name node Id vector
  RCP<Epetra_Vector>          nodeIds_;
  //@}

  //! @name radii vector
  RCP<Epetra_Vector>          radii_;
  //@}

  //! @name generations vector
  RCP<Epetra_Vector>          generations_;
  //@}

  //! @name Dirichlet boundary condition vectors
  RCP<Epetra_Vector>          bcval_;
  RCP<Epetra_Vector>          dbctog_;
  RCP<Epetra_Vector>          acini_bc_;

  //! @name acinar elementel values
  RCP<Epetra_Vector>          acini_e_volume0_;
  RCP<Epetra_Vector>          acini_e_volumenm_;
  RCP<Epetra_Vector>          acini_e_volumen_;
  RCP<Epetra_Vector>          acini_e_volumenp_;
  RCP<Epetra_Vector>          acini_e_volume_strain_;

  //@}

  //! Error vector that shows the convergenece of the nonlinear problem
  RCP<Epetra_Vector>          residual_;
  RCP<Epetra_Vector>          bc_residual_;
  //@}

  //! connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;
  //@}

  //!
  bool coupledTo3D_;
  //@}

  //! @name Nonlinear solution parameters
  int     maxiter_;
  double  non_lin_tol_;
  //@}


  //! @name Additional stuff related to coupling with 3D tissue models

  /// map between coupling ID and conditions on structure
  std::map<int,DRT::Condition*> coupcond_;

  /// map of coupling IDs
  Teuchos::RCP<Epetra_Map> coupmap_;

  std::map<int,double> pres_;

  //@}

}; // class RedAirwayImplicitTimeInt

}


#endif  // #ifndef AIRWAYIMPLICITINTEGRATION_H
