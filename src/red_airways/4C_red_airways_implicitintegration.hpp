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

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>
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
    RedAirwayImplicitTimeInt(Teuchos::RCP<Core::FE::Discretization> dis,
        std::unique_ptr<Core::LinAlg::Solver> solver, Teuchos::ParameterList& params,
        Core::IO::DiscretizationWriter& output);


    /*!
    \brief Destructor

    */
    virtual ~RedAirwayImplicitTimeInt() = default;

    /*!
    \brief start time loop for startingalgo, normal problems and restarts

    */
    void integrate();

    /*!
    \brief start time loop for startingalgo, normal problems and restarts

    */
    void integrate(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

    /*!
    \brief Do time integration (time loop)

    */
    void time_loop(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
    \brief Do one step time integration (time loop)

    */
    void time_step(bool CoupledTo3D = false,
        Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams = Teuchos::null);

    /*!
    \brief Integrate one step

    */
    void integrate_step(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams = Teuchos::null);


    /// setup the variables to do a new time step
    void prepare_time_step();


    /*!
    \brief solve of airways

    */
    void solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
    \brief nonlinear solver of airways

    */
    void non_lin_solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    /*!
    \brief solve the scalar transport

    */
    void solve_scatra(Teuchos::ParameterList& CouplingTo3DParams);

    /*!
      \brief build linear system matrix and rhs

      \param vel new guess at velocity, cross-sectional area, and pressure
    */
    void evaluate(const Core::LinAlg::Vector<double>& vel){};

    /*!
    \brief Update the solution after convergence of the linear
           iteration. Current solution becomes old solution of next
           timestep.
    */
    void time_update();

    /*!
    \brief initialize state saving vector

    */
    void init_save_state();

    /*! Save state vectors into saving vectors
    \brief update configuration and output to file/screen

    */
    void save_state();

    /*! load saved vectors into state vectors
    \brief update configuration and output to file/screen

    */
    void load_state();

    /*!
    \brief update configuration and output to file/screen

    */
    void output(bool CoupledTo3D = false,
        Teuchos::RCP<Teuchos::ParameterList> CouplingParams = Teuchos::null);

    /*!
    \brief Output for UQ problemtype

    */
    void output_uq(Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

    /*!
    \brief Adjust acini_volume with prestress

    */
    void compute_vol0_for_pre_stress();

    /*!
    \brief Assembling AIRWAY_ACINUS_DEP Vector for getting pext of nearest acinus

    */
    void compute_nearest_acinus(const Core::FE::Discretization& search_discret,
        std::set<int>* elecolset, std::set<int>* nodecolset,
        Teuchos::RCP<Core::LinAlg::Vector<double>> airway_acinus_dep);

    /*!
    \brief Assembling of the RHS Vector and the LHS Matrix

    */
    void assemble_mat_and_rhs();

    /*!
    \brief Evaluate the error residual

    */
    void eval_residual(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);


    /*!
    \brief Summing all the values in an Element ColMap
           We use this to find the total lung volume by summing up
           the acinar volumes
    */
    bool sum_all_col_elem_val(
        Core::LinAlg::Vector<double>& vec, Core::LinAlg::Vector<double>& sumCond, double& sum);


    /*!
    \brief read restart data

    */
    void read_restart(int step, bool coupledTo3D = false);

    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test();


    //! @name access methods for composite algorithms
    /// Return nodal values
    Teuchos::RCP<Core::LinAlg::Vector<double>> pnp() { return pnp_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> pn() { return pn_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> pnm() { return pnm_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_np() { return qin_np_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_np() { return qout_np_; }

    /// provide access to the Dirichlet map
    Teuchos::RCP<const Core::LinAlg::MapExtractor> dirich_maps() { return dbcmaps_; }

    /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
    ///
    /// This method provides backward compatability only. Formerly, the Dirichlet conditions
    /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
    /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> dirichlet();

    /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
    ///
    /// This method provides backward compatability only. Formerly, the Dirichlet conditions
    /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
    /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> inv_dirichlet();

    Teuchos::RCP<Core::LinAlg::SparseMatrix> mass_matrix()
    {
      return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(massmat_);
    }
    Teuchos::RCP<Core::FE::Discretization> discretization() { return discret_; }

    double dt() const { return dta_; }
    double time() const { return time_; }
    int step() const { return step_; }

    int itemax() const { return params_.get<int>("max nonlin iter steps"); }
    void set_itemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

    Teuchos::RCP<Core::LinAlg::Vector<double>> pext_np() { return p_extnp_; }

    /// Return elemental acini volume
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_volume() { return acini_e_volumenp_; }

    /// Return elemental airway volume
    Teuchos::RCP<Core::LinAlg::Vector<double>> airway_volume() { return elemVolumenp_; }

    /// Return elemental airway status regarding opening
    Teuchos::RCP<Core::LinAlg::Vector<double>> open() { return open_; }

    /// Return elemental airway trajectory describing opening tendency
    Teuchos::RCP<Core::LinAlg::Vector<double>> opening_trajectory() { return x_np_; }

    //@}


    //! @name methods related to coupling with 3D tissue models

    void setup_for_coupling();
    void set_airway_flux_from_tissue(Core::LinAlg::Vector<double>& coupflux);
    void extract_pressure(Core::LinAlg::Vector<double>& couppres);

    /// Hand over outputwriter to redairway_tissue
    Core::IO::DiscretizationWriter& get_output_writer() { return output_; }

    /// Hand over restartreader to redairway_tissue
    Teuchos::RCP<Core::IO::DiscretizationReader> get_restart_reader(int step)
    {
      return Teuchos::make_rcp<Core::IO::DiscretizationReader>(
          discret_, Global::Problem::instance()->input_control_file(), step);
    }

    //@}

   protected:
    //! @name general algorithm parameters
    //! reduced dimensional airway network discretization
    Teuchos::RCP<Core::FE::Discretization> discret_;
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
    Teuchos::RCP<Core::LinAlg::Vector<double>> rhs_;


    //! @name pressures at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> pnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> pn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> pnm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> p_nonlin_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> n_intr_ac_ln_;
    //@}

    //! @name inlet volumetric flow rates at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_nm_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> qi_nl_np_;
    //@}

    //! @trajectory vector x at time n+1 and n
    Teuchos::RCP<Core::LinAlg::Vector<double>> x_np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> x_n_;

    //! @state of airway open or closed
    Teuchos::RCP<Core::LinAlg::Vector<double>> open_;

    //! @neighbouring acinus of airway
    Teuchos::RCP<Core::LinAlg::Vector<double>> airway_acinus_dep_;

    //! @external pressure of airway
    Teuchos::RCP<Core::LinAlg::Vector<double>> p_extnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> p_extn_;
    // p vector in colmap format for computing p_extnp_ and p_extn_
    Teuchos::RCP<Core::LinAlg::Vector<double>> pnp_colmap_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> pn_colmap_;

    //! @name outlet volumetric flow rates at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_nm_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> qo_nl_np_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> qexp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qexp2_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> pexp_;
    //@}

    //! @name element volume at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolumenp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolumen_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolumenm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolume0_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemArea0_;

    //! @name element radius at time n+1
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemRadiusnp_;
    //@}

    //@}
    //! @name node Id vector
    Teuchos::RCP<Core::LinAlg::Vector<double>> nodeIds_;
    //@}

    //! @name radii vector
    Teuchos::RCP<Core::LinAlg::Vector<double>> radii_;
    //@}

    //! @name generations vector
    Teuchos::RCP<Core::LinAlg::Vector<double>> generations_;
    //@}

    //! @name Dirichlet boundary condition vectors
    Teuchos::RCP<Core::LinAlg::Vector<double>> bcval_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> dbctog_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_bc_;

    //! @name acinar elementel values
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_e_volume0_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_e_volumenm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_e_volumen_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_e_volumenp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_e_volume_strain_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_max_strain_location_;

    bool calcV0PreStress_;
    double transpulmpress_;

    //@}

    //! @name scalar transport variables inside airways
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e1scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e2scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e1scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e2scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e1scatraO2np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e2scatraO2np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> dscatraO2_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> dVolumeO2_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acinarDO2_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2np_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> junctionVolumeInMix_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> jVDofRowMix_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> junVolMix_Corrector_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> diffusionArea_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> cfls_;

    bool solveScatra_;
    bool compAwAcInter_;

    //@}

    //! @name state saving vectors

    // saving vector for pressure
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_pnm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_pn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_pnp_;

    // saving vector for inflow rate
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qin_nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qin_n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qin_np_;

    // saving vector for trajectory
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_x_n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_x_np_;

    // saving vector for outflow rates
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qout_nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qout_n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qout_np_;

    // saving vector for acinar volume
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_acini_e_volumenm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_acini_e_volumen_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_acini_e_volumenp_;

    // saving vector for element volume
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_elemVolumenm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_elemVolumen_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_elemVolumenp_;

    // saving vector for nodal O2 concentration
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_scatraO2np_;

    // saving vector for element inlet O2 concentration
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_e1scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_e1scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_e1scatraO2np_;

    // saving vector for element outlet O2 concentration
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_e2scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_e2scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_e2scatraO2np_;
    //@}

    //! Error vector that shows the convergenece of the nonlinear problem
    Teuchos::RCP<Core::LinAlg::Vector<double>> residual_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> bc_residual_;
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
