/*----------------------------------------------------------------------*/
/*! \file
\brief Coupling Manager for eXtended Fluid Ale Coupling

\level 3


*----------------------------------------------------------------------*/

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

namespace ADAPTER
{
  class AleFpsiWrapper;
  class Structure;
}  // namespace ADAPTER

namespace XFEM
{
  class XfaCouplingManager : public CouplingManager, public CouplingCommManager
  {
   public:
    /// constructor
    // in idx ... idx[0] structureal discretization index , idx[1] fluid discretization index in the
    // blockmatrix
    explicit XfaCouplingManager(Teuchos::RCP<FLD::XFluid> xfluid,
        Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale, std::vector<int> idx,
        Teuchos::RCP<ADAPTER::Structure> structure = Teuchos::null);
    //! @name Destruction
    //@{

    //! predict states in the coupling object
    void predict_coupling_states() override;

    //! Set required displacement & velocity states in the coupling object
    void SetCouplingStates() override;

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n) -- not yet done here
    void InitCouplingStates() override { return; }

    //! Add the coupling matrixes to the global systemmatrix
    // in ... scaling between xfluid evaluated coupling matrixes and coupled systemmatrix
    void AddCouplingMatrix(
        CORE::LINALG::BlockSparseMatrixBase& systemmatrix, double scaling) override;

    //! Add the coupling rhs

    // in scaling ... scaling between xfluid evaluated coupling rhs and coupled rhs
    // in me ... global map extractor of coupled problem (same index used as for idx)
    void AddCouplingRHS(Teuchos::RCP<Epetra_Vector> rhs, const CORE::LINALG::MultiMapExtractor& me,
        double scaling) override;

    //! Update (Perform after Each Timestep) -- nothing to do here
    void Update(double scaling) override { return; }

    //! Write Output -- nothing to do here
    void Output(CORE::IO::DiscretizationWriter& writer) override { return; }

    //! Read Restart -- nothing to do here
    void read_restart(CORE::IO::DiscretizationReader& reader) override { return; }

   private:
    //! Ale Object
    Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale_;
    //! eXtendedFluid
    Teuchos::RCP<FLD::XFluid> xfluid_;

    // Global Index in the blockmatrix of the coupled sytem [0] = fluid-, [1] = ale- block,
    // [2]struct- block
    std::vector<int> idx_;

    //! Structural Object (just set if Ale is coupled to a structure)
    Teuchos::RCP<ADAPTER::Structure> structure_;

    //! ALE-Structure coupling object on the matching interface
    Teuchos::RCP<XFEM::CouplingCommManager> ale_struct_coupling_;
  };
}  // namespace XFEM
FOUR_C_NAMESPACE_CLOSE

#endif
