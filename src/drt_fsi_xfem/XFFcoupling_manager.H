/*!----------------------------------------------------------------------
\brief Coupling Manager for eXtended Fluid Fluid Coupling

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

*----------------------------------------------------------------------*/

#ifndef XFFCOUPLING_MANAGER_H
#define XFFCOUPLING_MANAGER_H

#include "coupling_manager.H"
#include "coupling_comm_manager.H"

//#define FREESURFFLOW

namespace FLD
{
  class XFluid;
}

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFluidFluid;

  class XFFCoupling_Manager : public Coupling_Manager, public Coupling_Comm_Manager
  {
   public:
    /// constructor
    // in idx ... idx[0] fluid discretization index , idx[1] fluid discretization index in the
    // blockmatrix
    explicit XFFCoupling_Manager(Teuchos::RCP<ConditionManager> condmanager,
        Teuchos::RCP<FLD::XFluid> xfluid, Teuchos::RCP<FLD::XFluid> fluid, std::vector<int> idx);
    //! @name Destruction
    //@{

    /// virtual to get polymorph destruction
    virtual ~XFFCoupling_Manager() {}

    //! predict states in the coupling object
    virtual void PredictCouplingStates() {}

    //! Set required displacement & velocity states in the coupling object
    virtual void SetCouplingStates();

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n) -- not yet done here
    virtual void InitCouplingStates();

    //! Add the coupling matrixes to the global systemmatrix
    // in ... scaling between xfluid evaluated coupling matrixes and coupled systemmatrix
    virtual void AddCouplingMatrix(LINALG::BlockSparseMatrixBase& systemmatrix, double scaling);

    //! Add the coupling rhs

    // in scaling ... scaling between xfluid evaluated coupling rhs and coupled rhs
    // in me ... global map extractor of coupled problem (same index used as for idx)
    virtual void AddCouplingRHS(
        Teuchos::RCP<Epetra_Vector> rhs, const LINALG::MultiMapExtractor& me, double scaling);

    //! we need to think if inserting the ale matrixes are modifications (might conflict with other
    //! modifications)
    virtual bool ModifySysmatandRHS(LINALG::BlockSparseMatrixBase& systemmatrix,
        Teuchos::RCP<Epetra_Vector> rhs, const LINALG::MultiMapExtractor& me)
    {
      return false;
    }

    //! nothing to do
    virtual void PostLinearSolve(
        Teuchos::RCP<Epetra_Vector> inc, const LINALG::MultiMapExtractor& me)
    {
      return;
    }

    //! Update (Perform after Each Timestep) -- nothing to do here
    virtual void Update(double scaling) { return; }

    //! Write Output -- nothing to do here
    virtual void Output(IO::DiscretizationWriter& writer) { return; }

    //! Read Restart -- nothing to do here
    virtual void ReadRestart(IO::DiscretizationReader& reader) { return; }

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
#endif
