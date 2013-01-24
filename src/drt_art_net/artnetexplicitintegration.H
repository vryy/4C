
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

#ifndef ARTNETEXPLICITINTEGRATION_H
#define ARTNETEXPLICITINTEGRATION_H

//#include "../drt_fluid/fluidimpedancecondition.H"

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"
#include "art_junction.H"
#include "art_write_gnuplot.H"

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
                       LINALG::Solver&            solver,
                       Teuchos::ParameterList&    params,
                       IO::DiscretizationWriter&  output);


  /*!
  \brief Destructor

  */
  virtual ~ArtNetExplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate(bool CoupledTo3D, RCP<Teuchos::ParameterList> CouplingParams);

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop(bool CoupledTo3D, RCP<Teuchos::ParameterList> CouplingParams);



  /// setup the variables to do a new time step
  void PrepareTimeStep();


  /*!
  \brief solve linearised artery and bifurcation

  */
  void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

  void SolveScatra();

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
  void Output(bool               CoupledTo3D,
              RCP<Teuchos::ParameterList> CouplingParams);

  /*!
  \brief Assembling of the RHS Vector and the LHS Matrix 

  */
  void AssembleMatAndRHS();


  /*!
  \brief calculate values that could be used for postprocessing
         such as pressure and flowrate.
  */
  void CalcPostprocessingValues();


  void CalcScatraFromScatraFW(    RCP<Epetra_Vector>   scatra,
                                  RCP<Epetra_Vector>   scatran_fb);

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
  Teuchos::ParameterList&                params_;
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
  RCP<Epetra_Vector>    rhs_;

  /// (scatra) system matrix
  RCP<LINALG::SparseOperator> scatra_sysmat_;

  /// maps for scatra Dirichlet and free DOF sets
  RCP<Epetra_Vector>    scatra_bcval_;
  RCP<Epetra_Vector>    scatra_dbctog_;

  /// rhs: right hand side vector of scatra
  RCP<Epetra_Vector>    scatra_rhs_;

  //! @name Volumetric Flow rate and Cross-Sectional area at time n+1, n and n-1
  RCP<Epetra_Vector>    qanp_;
  RCP<Epetra_Vector>    qan_;
  RCP<Epetra_Vector>    qanm_;
  //@}

  //! @name Volumetric Flow rate and Cross-Sectional area at time n before solving Fluid 3D
  RCP<Epetra_Vector>    qan_3D_;
  //@}

  //! @name Volumetric Flow rate at time n+1, n and n-1
  RCP<Epetra_Vector>    qnp_;
  RCP<Epetra_Vector>    qn_;
  RCP<Epetra_Vector>    qnm_;
  //@}

  //! @name Pressure at time n
  RCP<Epetra_Vector>    pn_;
  //@}

  //! @name Area at time n
  RCP<Epetra_Vector>    an_;
  //@}

  //! @name Forward and backwar characteristic wave speeds at time n+1, n and n-1
  RCP<Epetra_Vector>    Wfo_;
  RCP<Epetra_Vector>    Wbo_;
  RCP<Epetra_Vector>    Wfnp_;
  RCP<Epetra_Vector>    Wfn_;
  RCP<Epetra_Vector>    Wfnm_;
  RCP<Epetra_Vector>    Wbnp_;
  RCP<Epetra_Vector>    Wbn_;
  RCP<Epetra_Vector>    Wbnm_;
  //@}

  //! @name scalar transport vectors at time n+1, n and n-1
  RCP<Epetra_Vector>    scatraO2n_;
  RCP<Epetra_Vector>    scatraO2np_;
  RCP<Epetra_Vector>    scatraO2wfn_;
  RCP<Epetra_Vector>    scatraO2wfnp_;
  RCP<Epetra_Vector>    scatraO2wbn_;
  RCP<Epetra_Vector>    scatraO2wbnp_;

  RCP<Epetra_Vector>    scatraCO2n_;
  RCP<Epetra_Vector>    scatraCO2np_;
  RCP<Epetra_Vector>    scatraCO2wfn_;
  RCP<Epetra_Vector>    scatraCO2wfnp_;
  RCP<Epetra_Vector>    scatraCO2wbn_;
  RCP<Epetra_Vector>    scatraCO2wbnp_;

  RCP<Epetra_Vector>    export_scatra_;
  bool solvescatra_;
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

  //! @name Junction boundary condition
  RCP<UTILS::ArtJunctionWrapper> artjun_;
  //@}

  //! @name 1D artery values at the junctions
  RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals_;
  //@}

  //! @name A condition to export 1D arteries as a gnuplot format
  RCP<UTILS::ArtWriteGnuplotWrapper> artgnu_;
  //@}


  /// connection between master and slave nodes on this proc
  RCP<map<int,std::vector<int> > >     pbcmapmastertoslave_;

  //!
  bool coupledTo3D_;
  //@}

  /// bio related special (outflow) boudaries
  //  RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  //@}

}; // class ArtNetExplicitTimeInt

}


#endif  // #ifndef ARTNETEXPLICITINTEGRATION_H
