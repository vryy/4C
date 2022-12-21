/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Manager for eXtended Fluid Poro Coupling

\level 3


*----------------------------------------------------------------------*/

#ifndef XFPCOUPLING_MANAGER_H
#define XFPCOUPLING_MANAGER_H

#include "coupling_manager.H"
#include "coupling_comm_manager.H"

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

  class XFPCoupling_Manager : public Coupling_Manager, public Coupling_Comm_Manager
  {
   public:
    /// constructor
    explicit XFPCoupling_Manager(Teuchos::RCP<ConditionManager> condmanager,
        Teuchos::RCP<POROELAST::PoroBase> poro, Teuchos::RCP<FLD::XFluid> xfluid,
        std::vector<int> idx);

    //! @name Destruction
    //@{

    /// virtual to get polymorph destruction
    virtual ~XFPCoupling_Manager() {}

    //! predict states in the coupling object
    virtual void PredictCouplingStates() {}

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n)
    virtual void InitCouplingStates();

    virtual void SetCouplingStates();

    virtual void AddCouplingMatrix(LINALG::BlockSparseMatrixBase& systemmatrix, double scaling);

    virtual void AddCouplingRHS(
        Teuchos::RCP<Epetra_Vector> rhs, const LINALG::MultiMapExtractor& me, double scaling);

    //! Update (Perform after Each Timestep)
    virtual void Update(double scaling);

    //! Write Output (For restart or write results on the interface)
    virtual void Output(IO::DiscretizationWriter& writer);

    //! Read Restart (For lambda_)
    virtual void ReadRestart(IO::DiscretizationReader& reader);

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
#endif
