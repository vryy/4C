// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>

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
    ArtNetExplicitTimeInt(std::shared_ptr<Core::FE::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        Core::IO::DiscretizationWriter& output);



    /*!
    \brief Initialization

    */
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) override;

    // create field test
    std::shared_ptr<Core::Utils::ResultTest> create_field_test() override;


    /*!
    \brief solve linearised artery and bifurcation

    */
    void solve(std::shared_ptr<Teuchos::ParameterList> CouplingTo3DParams) override;

    void solve_scatra() override;

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
    void output(bool CoupledTo3D, std::shared_ptr<Teuchos::ParameterList> CouplingParams) override;

    /*!
    \brief Test results

    */
    void test_results() override;

    /*!
    \brief calculate values that could be used for postprocessing
           such as pressure and flowrate.
    */
    void calc_postprocessing_values();


    void calc_scatra_from_scatra_fw(std::shared_ptr<Core::LinAlg::Vector<double>> scatra,
        std::shared_ptr<Core::LinAlg::Vector<double>> scatra_fb);

    /*!
    \brief read restart data

    */
    void read_restart(int step, bool CoupledTo3D = false) override;

    //! @name access methods for composite algorithms

    //  std::shared_ptr<Core::LinAlg::Vector<double>> Residual() { return residual_; } //This
    //  variable might be needed in future!
    std::shared_ptr<Core::LinAlg::Vector<double>> qnp() { return qnp_; }
    std::shared_ptr<Core::LinAlg::Vector<double>> q_anp() { return qanp_; }
    std::shared_ptr<Core::LinAlg::Vector<double>> areanp() { return areanp_; }
    // std::shared_ptr<Core::LinAlg::Vector<double>> Presnp() { return presnp_; }
    std::shared_ptr<Core::LinAlg::Vector<double>> qn() { return qn_; }
    std::shared_ptr<Core::LinAlg::Vector<double>> q_an() { return qan_; }
    std::shared_ptr<Core::LinAlg::Vector<double>> arean() { return arean_; }
    // std::shared_ptr<Core::LinAlg::Vector<double>> Presn()  { return presn_; }

    /// provide access to the Dirichlet map
    std::shared_ptr<const Core::LinAlg::MapExtractor> dirich_maps() { return dbcmaps_; }

    /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
    ///
    /// This method provides backward compatability only. Formerly, the Dirichlet conditions
    /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
    /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
    const std::shared_ptr<const Core::LinAlg::Vector<double>> dirichlet();

    /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
    ///
    /// This method provides backward compatability only. Formerly, the Dirichlet conditions
    /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
    /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
    const std::shared_ptr<const Core::LinAlg::Vector<double>> inv_dirichlet();

    std::shared_ptr<Core::LinAlg::SparseMatrix> mass_matrix()
    {
      return std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(massmat_);
    }

    //@}


   protected:
    /// (standard) mass matrix
    std::shared_ptr<Core::LinAlg::SparseOperator> massmat_;

    /// maps for scatra Dirichlet and free DOF sets
    std::shared_ptr<Core::LinAlg::Vector<double>> nodeIds_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatra_bcval_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatra_dbctog_;

    //! @name Volumetric Flow rate and Cross-Sectional area at time n+1, n and n-1
    std::shared_ptr<Core::LinAlg::Vector<double>> qanp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> qan_;
    std::shared_ptr<Core::LinAlg::Vector<double>> qanm_;
    //@}

    //! @name Volumetric Flow rate and Cross-Sectional area at time n before solving Fluid 3D
    std::shared_ptr<Core::LinAlg::Vector<double>> qan_3D_;
    //@}

    //! @name Volumetric Flow rate at time n+1, n and n-1
    std::shared_ptr<Core::LinAlg::Vector<double>> qnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> qn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> qnm_;
    //@}

    //! @name Pressure at time n
    std::shared_ptr<Core::LinAlg::Vector<double>> pn_;
    //@}

    //! @name Area at time n
    std::shared_ptr<Core::LinAlg::Vector<double>> an_;
    //@}

    //! @name Forward and backwar characteristic wave speeds at time n+1, n and n-1
    std::shared_ptr<Core::LinAlg::Vector<double>> Wfo_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wbo_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wfnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wfn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wfnm_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wbnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wbn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> Wbnm_;
    //@}

    //! @name scalar transport vectors at time n+1, n and n-1
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2nm_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2n_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2np_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2wfn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2wfnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2wbn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraO2wbnp_;

    std::shared_ptr<Core::LinAlg::Vector<double>> scatraCO2n_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraCO2np_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraCO2wfn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraCO2wfnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraCO2wbn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scatraCO2wbnp_;

    std::shared_ptr<Core::LinAlg::Vector<double>> export_scatra_;
    //@}

    //! @name saving state vectors
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_qanp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_qan_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_qanm_;

    std::shared_ptr<Core::LinAlg::Vector<double>> saved_Wfnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_Wfn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_Wfnm_;

    std::shared_ptr<Core::LinAlg::Vector<double>> saved_Wbnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_Wbn_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_Wbnm_;

    std::shared_ptr<Core::LinAlg::Vector<double>> saved_scatraO2np_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_scatraO2n_;
    std::shared_ptr<Core::LinAlg::Vector<double>> saved_scatraO2nm_;
    //@}

    //! @name cross-sectional area at time n+1, n and n-1
    std::shared_ptr<Core::LinAlg::Vector<double>> arean_;
    std::shared_ptr<Core::LinAlg::Vector<double>> areanp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> areanm_;
    //@}

    //! @name Dirichlet boundary condition vectors
    std::shared_ptr<Core::LinAlg::Vector<double>> bcval_;
    std::shared_ptr<Core::LinAlg::Vector<double>> dbctog_;
    //@}

    //! @name Junction boundary condition
    std::shared_ptr<Utils::ArtJunctionWrapper> artjun_;
    //@}

    //! @name 1D artery values at the junctions
    std::shared_ptr<std::map<const int, std::shared_ptr<Arteries::Utils::JunctionNodeParams>>>
        junc_nodal_vals_;
    //@}

    //! @name A condition to export 1D arteries as a gnuplot format
    std::shared_ptr<Utils::ArtWriteGnuplotWrapper> artgnu_;
    //@}

    //@}

  };  // class ArtNetExplicitTimeInt

}  // namespace Arteries


FOUR_C_NAMESPACE_CLOSE

#endif
