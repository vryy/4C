// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_XFEM_XFACOUPLING_MANAGER_HPP
#define FOUR_C_FSI_XFEM_XFACOUPLING_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_fsi_xfem_coupling_comm_manager.hpp"
#include "4C_fsi_xfem_coupling_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class XFluid;
}

namespace Adapter
{
  class AleFpsiWrapper;
  class Structure;
}  // namespace Adapter

namespace XFEM
{
  class XfaCouplingManager : public CouplingManager, public CouplingCommManager
  {
   public:
    /// constructor
    // in idx ... idx[0] structureal discretization index , idx[1] fluid discretization index in the
    // blockmatrix
    explicit XfaCouplingManager(std::shared_ptr<FLD::XFluid> xfluid,
        std::shared_ptr<Adapter::AleFpsiWrapper> ale, std::vector<int> idx,
        std::shared_ptr<Adapter::Structure> structure = nullptr);
    //! @name Destruction
    //@{

    //! predict states in the coupling object
    void predict_coupling_states() override;

    //! Set required displacement & velocity states in the coupling object
    void set_coupling_states() override;

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n) -- not yet done here
    void init_coupling_states() override { return; }

    //! Add the coupling matrixes to the global systemmatrix
    // in ... scaling between xfluid evaluated coupling matrixes and coupled systemmatrix
    void add_coupling_matrix(
        Core::LinAlg::BlockSparseMatrixBase& systemmatrix, double scaling) override;

    //! Add the coupling rhs

    // in scaling ... scaling between xfluid evaluated coupling rhs and coupled rhs
    // in me ... global map extractor of coupled problem (same index used as for idx)
    void add_coupling_rhs(std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
        const Core::LinAlg::MultiMapExtractor& me, double scaling) override;

    //! Update (Perform after Each Timestep) -- nothing to do here
    void update(double scaling) override { return; }

    //! Write Output -- nothing to do here
    void output(Core::IO::DiscretizationWriter& writer) override { return; }

    //! Read Restart -- nothing to do here
    void read_restart(Core::IO::DiscretizationReader& reader) override { return; }

   private:
    //! Ale Object
    std::shared_ptr<Adapter::AleFpsiWrapper> ale_;
    //! eXtendedFluid
    std::shared_ptr<FLD::XFluid> xfluid_;

    // Global Index in the blockmatrix of the coupled system [0] = fluid-, [1] = ale- block,
    // [2]struct- block
    std::vector<int> idx_;

    //! Structural Object (just set if Ale is coupled to a structure)
    std::shared_ptr<Adapter::Structure> structure_;

    //! ALE-Structure coupling object on the matching interface
    std::shared_ptr<XFEM::CouplingCommManager> ale_struct_coupling_;
  };
}  // namespace XFEM
FOUR_C_NAMESPACE_CLOSE

#endif
