
/*! \file
\brief Associated with control routine for artery solvers,

     including instationary solvers based on

     o two-step Taylor-Galerkin

\level 2




*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_EXPLICITINTEGRATION_HPP
#define FOUR_C_ART_NET_EXPLICITINTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_art_net_art_junction.hpp"
#include "4C_art_net_art_write_gnuplot.hpp"
#include "4C_art_net_timint.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

namespace Arteries
{
  /*!
  \brief time integration for arterial network problems

  */
  class ArtNetExplicitTimeInt : public TimInt
  {
    // friend class ArtNetResultTest;

   public:
    /*!
    \brief Standard Constructor

    */
    ArtNetExplicitTimeInt(Teuchos::RCP<Core::FE::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        Core::IO::DiscretizationWriter& output);



    /*!
    \brief Initialization

    */
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) override;

    // create field test
    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() override;


    /*!
    \brief solve linearised artery and bifurcation

    */
    void solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams) override;

    void solve_scatra() override;

    /*!
      \brief build linear system matrix and rhs


      \param vel new guess at velocity, cross-sectional area, and pressure
    */
    void evaluate(Teuchos::RCP<const Core::LinAlg::Vector<double>> vel){};

    /*!
    \brief Update the solution after convergence of the linear
           iteration. Current solution becomes old solution of next
           timestep.
    */
    void time_update() override;

    /*!
    \brief Initialize the saving state vectors
    */
    void init_save_state() override;

    /*!
    \brief Save the current vectors into the saving state vectors
    */
    void save_state() override;

    /*!
    \brief Load the currently saved state vectors into the currently used vectors
    */
    void load_state() override;

    /*!
    \brief update configuration and output to file/screen

    */
    void output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) override;

    /*!
    \brief Test results

    */
    void test_results() override;

    /*!
    \brief calculate values that could be used for postprocessing
           such as pressure and flowrate.
    */
    void calc_postprocessing_values();


    void calc_scatra_from_scatra_fw(Teuchos::RCP<Core::LinAlg::Vector<double>> scatra,
        Teuchos::RCP<Core::LinAlg::Vector<double>> scatra_fb);

    /*!
    \brief read restart data

    */
    void read_restart(int step, bool CoupledTo3D = false) override;

    //! @name access methods for composite algorithms

    //  Teuchos::RCP<Core::LinAlg::Vector<double>> Residual() { return residual_; } //This variable
    //  might be needed in future!
    Teuchos::RCP<Core::LinAlg::Vector<double>> qnp() { return qnp_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> q_anp() { return qanp_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> areanp() { return areanp_; }
    // Teuchos::RCP<Core::LinAlg::Vector<double>> Presnp() { return presnp_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> qn() { return qn_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> q_an() { return qan_; }
    Teuchos::RCP<Core::LinAlg::Vector<double>> arean() { return arean_; }
    // Teuchos::RCP<Core::LinAlg::Vector<double>> Presn()  { return presn_; }

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

    //@}


   protected:
    /// (standard) mass matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> massmat_;

    /// maps for scatra Dirichlet and free DOF sets
    Teuchos::RCP<Core::LinAlg::Vector<double>> nodeIds_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatra_bcval_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatra_dbctog_;

    //! @name Volumetric Flow rate and Cross-Sectional area at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> qanp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qan_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qanm_;
    //@}

    //! @name Volumetric Flow rate and Cross-Sectional area at time n before solving Fluid 3D
    Teuchos::RCP<Core::LinAlg::Vector<double>> qan_3D_;
    //@}

    //! @name Volumetric Flow rate at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> qnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qnm_;
    //@}

    //! @name Pressure at time n
    Teuchos::RCP<Core::LinAlg::Vector<double>> pn_;
    //@}

    //! @name Area at time n
    Teuchos::RCP<Core::LinAlg::Vector<double>> an_;
    //@}

    //! @name Forward and backwar characteristic wave speeds at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wfo_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wbo_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wfnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wfn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wfnm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wbnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wbn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> Wbnm_;
    //@}

    //! @name scalar transport vectors at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2nm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2wfn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2wfnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2wbn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraO2wbnp_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2wfn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2wfnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2wbn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatraCO2wbnp_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> export_scatra_;
    //@}

    //! @name saving state vectors
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qanp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qan_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_qanm_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_Wfnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_Wfn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_Wfnm_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_Wbnp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_Wbn_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_Wbnm_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_scatraO2np_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_scatraO2n_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> saved_scatraO2nm_;
    //@}

    //! @name cross-sectional area at time n+1, n and n-1
    Teuchos::RCP<Core::LinAlg::Vector<double>> arean_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> areanp_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> areanm_;
    //@}

    //! @name Dirichlet boundary condition vectors
    Teuchos::RCP<Core::LinAlg::Vector<double>> bcval_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> dbctog_;
    //@}

    //! @name Junction boundary condition
    Teuchos::RCP<UTILS::ArtJunctionWrapper> artjun_;
    //@}

    //! @name 1D artery values at the junctions
    Teuchos::RCP<std::map<const int, Teuchos::RCP<Arteries::UTILS::JunctionNodeParams>>>
        junc_nodal_vals_;
    //@}

    //! @name A condition to export 1D arteries as a gnuplot format
    Teuchos::RCP<UTILS::ArtWriteGnuplotWrapper> artgnu_;
    //@}

    //@}

  };  // class ArtNetExplicitTimeInt

}  // namespace Arteries


FOUR_C_NAMESPACE_CLOSE

#endif
