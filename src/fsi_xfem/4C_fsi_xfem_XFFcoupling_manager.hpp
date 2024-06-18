/*----------------------------------------------------------------------*/
/*! \file
\brief Coupling Manager for eXtended Fluid Fluid Coupling

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_XFEM_XFFCOUPLING_MANAGER_HPP
#define FOUR_C_FSI_XFEM_XFFCOUPLING_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_fsi_xfem_coupling_comm_manager.hpp"
#include "4C_fsi_xfem_coupling_manager.hpp"

FOUR_C_NAMESPACE_OPEN

// #define FREESURFFLOW

namespace FLD
{
  class XFluid;
}

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFluidFluid;

  class XffCouplingManager : public CouplingManager, public CouplingCommManager
  {
   public:
    /// constructor
    // in idx ... idx[0] fluid discretization index , idx[1] fluid discretization index in the
    // blockmatrix
    explicit XffCouplingManager(Teuchos::RCP<ConditionManager> condmanager,
        Teuchos::RCP<FLD::XFluid> xfluid, Teuchos::RCP<FLD::XFluid> fluid, std::vector<int> idx);
    //! @name Destruction
    //@{

    //! predict states in the coupling object
    void predict_coupling_states() override {}

    //! Set required displacement & velocity states in the coupling object
    void SetCouplingStates() override;

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n) -- not yet done here
    void InitCouplingStates() override;

    //! Add the coupling matrixes to the global systemmatrix
    // in ... scaling between xfluid evaluated coupling matrixes and coupled systemmatrix
    void AddCouplingMatrix(
        Core::LinAlg::BlockSparseMatrixBase& systemmatrix, double scaling) override;

    //! Add the coupling rhs

    // in scaling ... scaling between xfluid evaluated coupling rhs and coupled rhs
    // in me ... global map extractor of coupled problem (same index used as for idx)
    void AddCouplingRHS(Teuchos::RCP<Epetra_Vector> rhs, const Core::LinAlg::MultiMapExtractor& me,
        double scaling) override;

    //! we need to think if inserting the ale matrixes are modifications (might conflict with other
    //! modifications)
    virtual bool ModifySysmatandRHS(Core::LinAlg::BlockSparseMatrixBase& systemmatrix,
        Teuchos::RCP<Epetra_Vector> rhs, const Core::LinAlg::MultiMapExtractor& me)
    {
      return false;
    }

    //! nothing to do
    virtual void PostLinearSolve(
        Teuchos::RCP<Epetra_Vector> inc, const Core::LinAlg::MultiMapExtractor& me)
    {
      return;
    }

    //! Update (Perform after Each Timestep) -- nothing to do here
    void Update(double scaling) override { return; }

    //! Write Output -- nothing to do here
    void output(Core::IO::DiscretizationWriter& writer) override { return; }

    //! Read Restart -- nothing to do here
    void read_restart(Core::IO::DiscretizationReader& reader) override { return; }

   private:
    //! FFI Mesh Coupling Object
    Teuchos::RCP<MeshCouplingFluidFluid> mcffi_;

    //! embeddedFluid Object
    Teuchos::RCP<FLD::XFluid> fluid_;
    //! eXtendedFluid
    Teuchos::RCP<FLD::XFluid> xfluid_;

    //"XFEMSurfFluidFluid"
    const std::string cond_name_;

    // Global Index in the blockmatrix of the coupled sytem [0] = fluid-, [1] = ale- block,
    // [2]struct- block
    std::vector<int> idx_;
  };
}  // namespace XFEM
FOUR_C_NAMESPACE_CLOSE

#endif
