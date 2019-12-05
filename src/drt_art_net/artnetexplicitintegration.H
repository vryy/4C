
/*! \file
\brief Associated with control routine for artery solvers,

     including instationary solvers based on

     o two-step Taylor-Galerkin

\level 2

\maintainer Johannes Kremheller



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
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_io/io.H"
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

#include "art_net_timint.H"

namespace ART
{
  /*!
  \brief time integration for arterial network problems

  \author ismail (ismail@lnm.mw.tum.de)
  */
  class ArtNetExplicitTimeInt : public TimInt
  {
    // friend class ArtNetResultTest;

   public:
    /*!
    \brief Standard Constructor

    */
    ArtNetExplicitTimeInt(Teuchos::RCP<DRT::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        FILE* errfile, IO::DiscretizationWriter& output);


    /*!
    \brief Destructor

    */
    virtual ~ArtNetExplicitTimeInt();

    /*!
    \brief Initialization

    */
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname);

    // create field test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();


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
    \brief Initialize the saving state vectors
    */
    void InitSaveState();

    /*!
    \brief Save the current vectors into the saving state vectors
    */
    void SaveState();

    /*!
    \brief Load the currently saved state vectors into the currently used vectors
    */
    void LoadState();

    /*!
    \brief update configuration and output to file/screen

    */
    void Output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

    /*!
    \brief Test results

    */
    void TestResults();

    /*!
    \brief calculate values that could be used for postprocessing
           such as pressure and flowrate.
    */
    void CalcPostprocessingValues();


    void CalcScatraFromScatraFW(
        Teuchos::RCP<Epetra_Vector> scatra, Teuchos::RCP<Epetra_Vector> scatran_fb);

    /*!
    \brief read restart data

    */
    void ReadRestart(int step, bool CoupledTo3D = false);

    //! @name access methods for composite algorithms

    //  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; } //This variable might be needed
    //  in future!
    Teuchos::RCP<Epetra_Vector> Qnp() { return qnp_; }
    Teuchos::RCP<Epetra_Vector> QAnp() { return qanp_; }
    Teuchos::RCP<Epetra_Vector> Areanp() { return areanp_; }
    // Teuchos::RCP<Epetra_Vector> Presnp() { return presnp_; }
    Teuchos::RCP<Epetra_Vector> Qn() { return qn_; }
    Teuchos::RCP<Epetra_Vector> QAn() { return qan_; }
    Teuchos::RCP<Epetra_Vector> Arean() { return arean_; }
    // Teuchos::RCP<Epetra_Vector> Presn()  { return presn_; }

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

    Teuchos::RCP<LINALG::SparseMatrix> MassMatrix()
    {
      return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(massmat_);
    }

    //@}


   protected:
    /// (standard) mass matrix
    Teuchos::RCP<LINALG::SparseOperator> massmat_;

    /// maps for scatra Dirichlet and free DOF sets
    Teuchos::RCP<Epetra_Vector> nodeIds_;
    Teuchos::RCP<Epetra_Vector> scatra_bcval_;
    Teuchos::RCP<Epetra_Vector> scatra_dbctog_;

    //! @name Volumetric Flow rate and Cross-Sectional area at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> qanp_;
    Teuchos::RCP<Epetra_Vector> qan_;
    Teuchos::RCP<Epetra_Vector> qanm_;
    //@}

    //! @name Volumetric Flow rate and Cross-Sectional area at time n before solving Fluid 3D
    Teuchos::RCP<Epetra_Vector> qan_3D_;
    //@}

    //! @name Volumetric Flow rate at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> qnp_;
    Teuchos::RCP<Epetra_Vector> qn_;
    Teuchos::RCP<Epetra_Vector> qnm_;
    //@}

    //! @name Pressure at time n
    Teuchos::RCP<Epetra_Vector> pn_;
    //@}

    //! @name Area at time n
    Teuchos::RCP<Epetra_Vector> an_;
    //@}

    //! @name Forward and backwar characteristic wave speeds at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> Wfo_;
    Teuchos::RCP<Epetra_Vector> Wbo_;
    Teuchos::RCP<Epetra_Vector> Wfnp_;
    Teuchos::RCP<Epetra_Vector> Wfn_;
    Teuchos::RCP<Epetra_Vector> Wfnm_;
    Teuchos::RCP<Epetra_Vector> Wbnp_;
    Teuchos::RCP<Epetra_Vector> Wbn_;
    Teuchos::RCP<Epetra_Vector> Wbnm_;
    //@}

    //! @name scalar transport vectors at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> scatraO2n_;
    Teuchos::RCP<Epetra_Vector> scatraO2np_;
    Teuchos::RCP<Epetra_Vector> scatraO2wfn_;
    Teuchos::RCP<Epetra_Vector> scatraO2wfnp_;
    Teuchos::RCP<Epetra_Vector> scatraO2wbn_;
    Teuchos::RCP<Epetra_Vector> scatraO2wbnp_;

    Teuchos::RCP<Epetra_Vector> scatraCO2n_;
    Teuchos::RCP<Epetra_Vector> scatraCO2np_;
    Teuchos::RCP<Epetra_Vector> scatraCO2wfn_;
    Teuchos::RCP<Epetra_Vector> scatraCO2wfnp_;
    Teuchos::RCP<Epetra_Vector> scatraCO2wbn_;
    Teuchos::RCP<Epetra_Vector> scatraCO2wbnp_;

    Teuchos::RCP<Epetra_Vector> export_scatra_;
    //@}

    //! @name saving state vectors
    Teuchos::RCP<Epetra_Vector> saved_qanp_;
    Teuchos::RCP<Epetra_Vector> saved_qan_;
    Teuchos::RCP<Epetra_Vector> saved_qanm_;

    Teuchos::RCP<Epetra_Vector> saved_Wfnp_;
    Teuchos::RCP<Epetra_Vector> saved_Wfn_;
    Teuchos::RCP<Epetra_Vector> saved_Wfnm_;

    Teuchos::RCP<Epetra_Vector> saved_Wbnp_;
    Teuchos::RCP<Epetra_Vector> saved_Wbn_;
    Teuchos::RCP<Epetra_Vector> saved_Wbnm_;

    Teuchos::RCP<Epetra_Vector> saved_scatraO2np_;
    Teuchos::RCP<Epetra_Vector> saved_scatraO2n_;
    Teuchos::RCP<Epetra_Vector> saved_scatraO2nm_;
    //@}

    //! @name cross-sectional area at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> arean_;
    Teuchos::RCP<Epetra_Vector> areanp_;
    Teuchos::RCP<Epetra_Vector> areanm_;
    //@}

    //! @name Dirichlet boundary condition vectors
    Teuchos::RCP<Epetra_Vector> bcval_;
    Teuchos::RCP<Epetra_Vector> dbctog_;
    //@}

    //! @name Junction boundary condition
    Teuchos::RCP<UTILS::ArtJunctionWrapper> artjun_;
    //@}

    //! @name 1D artery values at the junctions
    Teuchos::RCP<std::map<const int, Teuchos::RCP<ART::UTILS::JunctionNodeParams>>>
        junc_nodal_vals_;
    //@}

    //! @name A condition to export 1D arteries as a gnuplot format
    Teuchos::RCP<UTILS::ArtWriteGnuplotWrapper> artgnu_;
    //@}

    //@}

  };  // class ArtNetExplicitTimeInt

}  // namespace ART


#endif  // #ifndef ARTNETEXPLICITINTEGRATION_H
