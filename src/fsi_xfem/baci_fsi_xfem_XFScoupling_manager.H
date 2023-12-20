/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Manager for eXtended Fluid Structural Coupling

\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_FSI_XFEM_XFSCOUPLING_MANAGER_H
#define BACI_FSI_XFEM_XFSCOUPLING_MANAGER_H

#include "baci_config.H"

#include "baci_fsi_xfem_coupling_comm_manager.H"
#include "baci_fsi_xfem_coupling_manager.H"

BACI_NAMESPACE_OPEN

namespace FLD
{
  class XFluid;
}

namespace ADAPTER
{
  class Structure;
}

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFSI;

  class XFSCoupling_Manager : public Coupling_Manager, public Coupling_Comm_Manager
  {
   public:
    /// constructor
    // in idx ... idx[0] structureal discretization index , idx[1] fluid discretization index in the
    // blockmatrix
    explicit XFSCoupling_Manager(Teuchos::RCP<ConditionManager> condmanager,
        Teuchos::RCP<ADAPTER::Structure> structure, Teuchos::RCP<FLD::XFluid> xfluid,
        std::vector<int> idx);

    //! @name Destruction
    //@{

    //! predict states in the coupling object
    void PredictCouplingStates() override{};

    //! init...
    void InitCouplingStates() override;

    //! Set required displacement & velocity states in the coupling object
    void SetCouplingStates() override;

    //! Add the coupling matrixes to the global systemmatrix
    // in ... scaling between xfluid evaluated coupling matrixes and coupled systemmatrix
    void AddCouplingMatrix(
        CORE::LINALG::BlockSparseMatrixBase& systemmatrix, double scaling) override;

    //! Add the coupling rhs

    // in scaling ... scaling between xfluid evaluated coupling rhs and coupled rhs
    // in me ... global map extractor of coupled problem (same index used as for idx)
    void AddCouplingRHS(Teuchos::RCP<Epetra_Vector> rhs, const CORE::LINALG::MultiMapExtractor& me,
        double scaling) override;

    //! Update (Perform after Each Timestep)
    void Update(double scaling) override;

    //! Write Output (For restart or write results on the interface)
    void Output(IO::DiscretizationWriter& writer) override;

    //! Read Restart (For lambda_)
    void ReadRestart(IO::DiscretizationReader& reader) override;

   private:
    //! Get Timeface on the interface (for OST this is 1/(theta dt))
    double GetInterfaceTimefac();

    //! FSI Mesh Coupling Object
    Teuchos::RCP<MeshCouplingFSI> mcfsi_;

    //! Structural Object
    Teuchos::RCP<ADAPTER::Structure> struct_;
    //! eXtendedFluid
    Teuchos::RCP<FLD::XFluid> xfluid_;

    //"XFEMSurfFSIMono"
    const std::string cond_name_;

    // Global Index in the blockmatrix of the coupled sytem [0] = structure-, [1] = fluid- block
    std::vector<int> idx_;

    bool interface_second_order_;

    //--------------------------------------------------------------------------//
    //! @name Store the Coupling RHS of the Old Timestep in lambda

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie forces onto the structure,
    //! Robin-type forces consisting of fluid forces and the Nitsche penalty term contribution)
    //! evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;
  };
}  // namespace XFEM
BACI_NAMESPACE_CLOSE

#endif
