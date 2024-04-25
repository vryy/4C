/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Manager for eXtended Fluid Poro Coupling

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_XFEM_XFPCOUPLING_MANAGER_HPP
#define FOUR_C_FSI_XFEM_XFPCOUPLING_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_fsi_xfem_coupling_comm_manager.hpp"
#include "4C_fsi_xfem_coupling_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class XFluid;
}

namespace POROELAST
{
  class PoroBase;
}

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFPI;

  class XfpCouplingManager : public CouplingManager, public CouplingCommManager
  {
   public:
    /// constructor
    explicit XfpCouplingManager(Teuchos::RCP<ConditionManager> condmanager,
        Teuchos::RCP<POROELAST::PoroBase> poro, Teuchos::RCP<FLD::XFluid> xfluid,
        std::vector<int> idx);

    //! @name Destruction
    //@{

    //! predict states in the coupling object
    void PredictCouplingStates() override {}

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n)
    void InitCouplingStates() override;

    void SetCouplingStates() override;

    void AddCouplingMatrix(
        CORE::LINALG::BlockSparseMatrixBase& systemmatrix, double scaling) override;

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



    Teuchos::RCP<MeshCouplingFPI> mcfpi_ps_ps_;
    Teuchos::RCP<MeshCouplingFPI> mcfpi_ps_pf_;
    Teuchos::RCP<MeshCouplingFPI> mcfpi_pf_ps_;
    Teuchos::RCP<MeshCouplingFPI> mcfpi_pf_pf_;

    Teuchos::RCP<POROELAST::PoroBase> poro_;
    Teuchos::RCP<FLD::XFluid> xfluid_;

    std::string cond_name_ps_ps_;
    std::string cond_name_ps_pf_;
    std::string cond_name_pf_ps_;
    std::string cond_name_pf_pf_;

    // Global Index in the blockmatrix of the coupled sytem [0] = structure-, [1] = fluid- block,
    // [2] = porofluid-block
    std::vector<int> idx_;

    bool interface_second_order_;

    //--------------------------------------------------------------------------//
    //! @name Store the Coupling RHS of the Old Timestep in lambda

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie forces onto the structure,
    //! Robin-type forces consisting of fluid forces and the Nitsche penalty term contribution)
    //! evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_ps_;
    Teuchos::RCP<Epetra_Vector> lambda_pf_;
  };
}  // namespace XFEM
FOUR_C_NAMESPACE_CLOSE

#endif
