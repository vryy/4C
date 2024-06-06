/*---------------------------------------------------------------------*/
/*! \file

\brief Associated with control routine for reduced dimensional airways
  solvers


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_RED_AIRWAYS_IMPLICITINTEGRATION_HPP
#define FOUR_C_RED_AIRWAYS_IMPLICITINTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class ResultTest;
  class RedAirwayResultTest;
}  // namespace Discret

namespace Core::LinAlg
{
  class MapExtractor;
  class Solver;
}  // namespace Core::LinAlg

namespace Airway
{
  /*!
  \brief time integration for reduced dimensional airway network problems

  */
  class RedAirwayImplicitTimeInt
  {
   public:
    /*!
    \brief Standard Constructor

    */
    RedAirwayImplicitTimeInt(Teuchos::RCP<Discret::Discretization> dis,
        std::unique_ptr<Core::LinAlg::Solver> solver, Teuchos::ParameterList& params,
        Core::IO::DiscretizationWriter& output);


    /*!
    \brief Destructor

    */
    virtual ~RedAirwayImplicitTimeInt() = default;

    /*!
    \brief start time loop for startingalgo, normal problems and restarts

    */
    void Integrate();

    /*!
    \brief start time loop for startingalgo, normal problems and restarts

    */
    void Integrate(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

    /*!
    \brief Do time integration (time loop)

    */
    void TimeLoop(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
    \brief Do one step time integration (time loop)

    */
    void TimeStep(bool CoupledTo3D = false,
        Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams = Teuchos::null);

    /*!
    \brief Integrate one step

    */
    void IntegrateStep(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams = Teuchos::null);


    /// setup the variables to do a new time step
    void prepare_time_step();


    /*!
    \brief solve of airways

    */
    void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
    \brief nonlinear solver of airways

    */
    void NonLin_Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
    \brief solve the scalar transport

    */
    void SolveScatra(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
      \brief build linear system matrix and rhs

      \param vel new guess at velocity, cross-sectional area, and pressure
    */
    void Evaluate(Teuchos::RCP<const Epetra_Vector> vel){};

    /*!
    \brief Update the solution after convergence of the linear
           iteration. Current solution becomes old solution of next
           timestep.
    */
    void TimeUpdate();

    /*!
    \brief initialize state saving vector

    */
    void InitSaveState();

    /*! Save state vectors into saving vectors
    \brief update configuration and output to file/screen

    */
    void SaveState();

    /*! load saved vectors into state vectors
    \brief update configuration and output to file/screen

    */
    void LoadState();

    /*!
    \brief update configuration and output to file/screen

    */
    void Output(bool CoupledTo3D = false,
        Teuchos::RCP<Teuchos::ParameterList> CouplingParams = Teuchos::null);

    /*!
    \brief Output for UQ problemtype

    */
    void OutputUQ(Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

    /*!
    \brief Adjust acini_volume with prestress

    */
    void compute_vol0_for_pre_stress();

    /*!
    \brief Assembling AIRWAY_ACINUS_DEP Vector for getting pext of nearest acinus

    */
    void compute_nearest_acinus(Teuchos::RCP<Discret::Discretization const> search_discret,
        std::set<int>* elecolset, std::set<int>* nodecolset,
        Teuchos::RCP<Epetra_Vector> airway_acinus_dep);

    /*!
    \brief Assembling of the RHS Vector and the LHS Matrix

    */
    void assemble_mat_and_rhs();

    /*!
    \brief Evaluate the error residual

    */
    void EvalResidual(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);


    /*!
    \brief Summing all the values in an Element ColMap
           We use this to find the total lung volume by summing up
           the acinar volumes
    */
    bool SumAllColElemVal(
        Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<Epetra_Vector> sumCond, double& sum);


    /*!
    \brief read restart data

    */
    void read_restart(int step, bool coupledTo3D = false);

    Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest();


    //! @name access methods for composite algorithms
    /// Return nodal values
    Teuchos::RCP<Epetra_Vector> Pnp() { return pnp_; }
    Teuchos::RCP<Epetra_Vector> Pn() { return pn_; }
    Teuchos::RCP<Epetra_Vector> Pnm() { return pnm_; }
    Teuchos::RCP<Epetra_Vector> Qin_np() { return qin_np_; }
    Teuchos::RCP<Epetra_Vector> Qout_np() { return qout_np_; }

    /// provide access to the Dirichlet map
    Teuchos::RCP<const Core::LinAlg::MapExtractor> DirichMaps() { return dbcmaps_; }

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

    Teuchos::RCP<Core::LinAlg::SparseMatrix> MassMatrix()
    {
      return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(massmat_);
    }
    Teuchos::RCP<Discret::Discretization> discretization() { return discret_; }

    double Dt() const { return dta_; }
    double Time() const { return time_; }
    int Step() const { return step_; }

    int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
    void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

    Teuchos::RCP<Epetra_Vector> Pext_np() { return p_extnp_; }

    /// Return elemental acini volume
    Teuchos::RCP<Epetra_Vector> AciniVolume() { return acini_e_volumenp_; }

    /// Return elemental airway volume
    Teuchos::RCP<Epetra_Vector> AirwayVolume() { return elemVolumenp_; }

    /// Return elemental airway status regarding opening
    Teuchos::RCP<Epetra_Vector> Open() { return open_; }

    /// Return elemental airway trajectory describing opening tendency
    Teuchos::RCP<Epetra_Vector> OpeningTrajectory() { return x_np_; }

    //@}


    //! @name methods related to coupling with 3D tissue models

    void SetupForCoupling();
    void set_airway_flux_from_tissue(Teuchos::RCP<Epetra_Vector> coupflux);
    void ExtractPressure(Teuchos::RCP<Epetra_Vector> couppres);

    /// Hand over outputwriter to redairway_tissue
    Core::IO::DiscretizationWriter& GetOutputWriter() { return output_; }

    /// Hand over restartreader to redairway_tissue
    Teuchos::RCP<Core::IO::DiscretizationReader> GetRestartReader(int step)
    {
      return Teuchos::rcp(new Core::IO::DiscretizationReader(
          discret_, Global::Problem::Instance()->InputControlFile(), step));
    }

    //@}

   protected:
    //! @name general algorithm parameters
    //! reduced dimensional airway network discretization
    Teuchos::RCP<Discret::Discretization> discret_;
    std::unique_ptr<Core::LinAlg::Solver> solver_;
    Teuchos::ParameterList params_;
    Core::IO::DiscretizationWriter& output_;
    //! the processor ID from the communicator
    int myrank_;

    //@}

    //! @name time stepping variables
    double time_;     ///< physical time
    int step_;        ///< timestep
    int stepmax_;     ///< maximal number of timesteps
    double maxtime_;  ///< maximal physical computation time
    //@}


    /// constant density extracted from element material for incompressible flow
    double density_;

    //! @name restart variables
    int uprestart_;
    int upres_;
    //@}

    //! @name time step sizes
    double dta_;
    double dtp_;

    //@}

    /// cpu-time measures
    double dtele_;
    double dtfilter_;
    double dtsolve_;

    /// (standard) mass matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> massmat_;

    /// (standard) system matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;

    /// rhs: right hand side vector
    Teuchos::RCP<Epetra_Vector> rhs_;


    //! @name pressures at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> pnp_;
    Teuchos::RCP<Epetra_Vector> pn_;
    Teuchos::RCP<Epetra_Vector> pnm_;
    Teuchos::RCP<Epetra_Vector> p_nonlin_;

    Teuchos::RCP<Epetra_Vector> n_intr_ac_ln_;
    //@}

    //! @name inlet volumetric flow rates at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> qin_np_;
    Teuchos::RCP<Epetra_Vector> qin_n_;
    Teuchos::RCP<Epetra_Vector> qin_nm_;

    Teuchos::RCP<Epetra_Vector> qi_nl_np_;
    //@}

    //! @trajectory vector x at time n+1 and n
    Teuchos::RCP<Epetra_Vector> x_np_;
    Teuchos::RCP<Epetra_Vector> x_n_;

    //! @state of airway open or closed
    Teuchos::RCP<Epetra_Vector> open_;

    //! @neighbouring acinus of airway
    Teuchos::RCP<Epetra_Vector> airway_acinus_dep_;

    //! @external pressure of airway
    Teuchos::RCP<Epetra_Vector> p_extnp_;
    Teuchos::RCP<Epetra_Vector> p_extn_;
    // p vector in colmap format for computing p_extnp_ and p_extn_
    Teuchos::RCP<Epetra_Vector> pnp_colmap_;
    Teuchos::RCP<Epetra_Vector> pn_colmap_;

    //! @name outlet volumetric flow rates at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> qout_np_;
    Teuchos::RCP<Epetra_Vector> qout_n_;
    Teuchos::RCP<Epetra_Vector> qout_nm_;

    Teuchos::RCP<Epetra_Vector> qo_nl_np_;

    Teuchos::RCP<Epetra_Vector> qexp_;
    Teuchos::RCP<Epetra_Vector> qexp2_;
    Teuchos::RCP<Epetra_Vector> pexp_;
    //@}

    //! @name element volume at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector> elemVolumenp_;
    Teuchos::RCP<Epetra_Vector> elemVolumen_;
    Teuchos::RCP<Epetra_Vector> elemVolumenm_;
    Teuchos::RCP<Epetra_Vector> elemVolume0_;
    Teuchos::RCP<Epetra_Vector> elemArea0_;

    //! @name element radius at time n+1
    Teuchos::RCP<Epetra_Vector> elemRadiusnp_;
    //@}

    //@}
    //! @name node Id vector
    Teuchos::RCP<Epetra_Vector> nodeIds_;
    //@}

    //! @name radii vector
    Teuchos::RCP<Epetra_Vector> radii_;
    //@}

    //! @name generations vector
    Teuchos::RCP<Epetra_Vector> generations_;
    //@}

    //! @name Dirichlet boundary condition vectors
    Teuchos::RCP<Epetra_Vector> bcval_;
    Teuchos::RCP<Epetra_Vector> dbctog_;
    Teuchos::RCP<Epetra_Vector> acini_bc_;

    //! @name acinar elementel values
    Teuchos::RCP<Epetra_Vector> acini_e_volume0_;
    Teuchos::RCP<Epetra_Vector> acini_e_volumenm_;
    Teuchos::RCP<Epetra_Vector> acini_e_volumen_;
    Teuchos::RCP<Epetra_Vector> acini_e_volumenp_;
    Teuchos::RCP<Epetra_Vector> acini_e_volume_strain_;
    Teuchos::RCP<Epetra_Vector> acini_max_strain_location_;

    bool calcV0PreStress_;
    double transpulmpress_;

    //@}

    //! @name scalar transport variables inside airways
    Teuchos::RCP<Epetra_Vector> scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> e1scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> e2scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> scatraO2n_;
    Teuchos::RCP<Epetra_Vector> e1scatraO2n_;
    Teuchos::RCP<Epetra_Vector> e2scatraO2n_;
    Teuchos::RCP<Epetra_Vector> scatraO2np_;
    Teuchos::RCP<Epetra_Vector> e1scatraO2np_;
    Teuchos::RCP<Epetra_Vector> e2scatraO2np_;
    Teuchos::RCP<Epetra_Vector> dscatraO2_;
    Teuchos::RCP<Epetra_Vector> dVolumeO2_;
    Teuchos::RCP<Epetra_Vector> acinarDO2_;

    Teuchos::RCP<Epetra_Vector> scatraCO2nm_;
    Teuchos::RCP<Epetra_Vector> scatraCO2n_;
    Teuchos::RCP<Epetra_Vector> scatraCO2np_;

    Teuchos::RCP<Epetra_Vector> junctionVolumeInMix_;
    Teuchos::RCP<Epetra_Vector> jVDofRowMix_;
    Teuchos::RCP<Epetra_Vector> junVolMix_Corrector_;
    Teuchos::RCP<Epetra_Vector> diffusionArea_;
    Teuchos::RCP<Epetra_Vector> cfls_;

    bool solveScatra_;
    bool compAwAcInter_;

    //@}

    //! @name state saving vectors

    // saving vector for pressure
    Teuchos::RCP<Epetra_Vector> saved_pnm_;
    Teuchos::RCP<Epetra_Vector> saved_pn_;
    Teuchos::RCP<Epetra_Vector> saved_pnp_;

    // saving vector for inflow rate
    Teuchos::RCP<Epetra_Vector> saved_qin_nm_;
    Teuchos::RCP<Epetra_Vector> saved_qin_n_;
    Teuchos::RCP<Epetra_Vector> saved_qin_np_;

    // saving vector for trajectory
    Teuchos::RCP<Epetra_Vector> saved_x_n_;
    Teuchos::RCP<Epetra_Vector> saved_x_np_;

    // saving vector for outflow rates
    Teuchos::RCP<Epetra_Vector> saved_qout_nm_;
    Teuchos::RCP<Epetra_Vector> saved_qout_n_;
    Teuchos::RCP<Epetra_Vector> saved_qout_np_;

    // saving vector for acinar volume
    Teuchos::RCP<Epetra_Vector> saved_acini_e_volumenm_;
    Teuchos::RCP<Epetra_Vector> saved_acini_e_volumen_;
    Teuchos::RCP<Epetra_Vector> saved_acini_e_volumenp_;

    // saving vector for element volume
    Teuchos::RCP<Epetra_Vector> saved_elemVolumenm_;
    Teuchos::RCP<Epetra_Vector> saved_elemVolumen_;
    Teuchos::RCP<Epetra_Vector> saved_elemVolumenp_;

    // saving vector for nodal O2 concentration
    Teuchos::RCP<Epetra_Vector> saved_scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> saved_scatraO2n_;
    Teuchos::RCP<Epetra_Vector> saved_scatraO2np_;

    // saving vector for element inlet O2 concentration
    Teuchos::RCP<Epetra_Vector> saved_e1scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> saved_e1scatraO2n_;
    Teuchos::RCP<Epetra_Vector> saved_e1scatraO2np_;

    // saving vector for element outlet O2 concentration
    Teuchos::RCP<Epetra_Vector> saved_e2scatraO2nm_;
    Teuchos::RCP<Epetra_Vector> saved_e2scatraO2n_;
    Teuchos::RCP<Epetra_Vector> saved_e2scatraO2np_;
    //@}

    //! Error vector that shows the convergenece of the nonlinear problem
    Teuchos::RCP<Epetra_Vector> residual_;
    Teuchos::RCP<Epetra_Vector> bc_residual_;
    //@}

    //! connection between master and slave nodes on this proc
    Teuchos::RCP<std::map<int, std::vector<int>>> pbcmapmastertoslave_;
    //@}

    //!
    bool coupledTo3D_;
    //@}

    //! @name Nonlinear solution parameters
    int maxiter_;
    double non_lin_tol_;
    //@}


    //! @name Additional stuff related to coupling with 3D tissue models

    /// map between coupling ID and conditions on structure
    std::map<int, Core::Conditions::Condition*> coupcond_;

    /// map of coupling IDs
    Teuchos::RCP<Epetra_Map> coupmap_;

    std::map<int, double> pres_;

    //@}

  };  // class RedAirwayImplicitTimeInt

}  // namespace Airway


FOUR_C_NAMESPACE_CLOSE

#endif
