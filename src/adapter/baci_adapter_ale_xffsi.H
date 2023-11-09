/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 2


*/
/*----------------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_ALE_XFFSI_H
#define BACI_ADAPTER_ALE_XFFSI_H


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_adapter_ale_fsi.H"
#include "baci_ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  /*! \brief ALE Wrapper for XFFSI Problems
   *
   *  Provide XFFSI specific ALE functionalities here by overloading the
   *  respective routines from ADAPTER::AleFsiWrapper and ADAPTER::AleWrapper
   *
   *  \sa ADAPTER::Ale, ADAPTER::AleWrapper, ADAPTER::AleFsiWrapper
   *
   *  \author mayr.mt \date 10/2014
   */
  class AleXFFsiWrapper : public AleFsiWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleXFFsiWrapper(Teuchos::RCP<Ale> ale);

    //@}

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override;

    /// evaluate and assemble residual and jacobian matrix
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc  ///< step increment such that \f$ x_{n+1}^{k+1} =
                                                   ///< x_{n}^{converged}+ stepinc \f$
        ) override;

    /*! \brief Relax ALE after some deformation has occured
     *
     *  At the beginning of the fluid-fluid-fsi step the xfem-dofs are
     *  dirichlet values so that they can not change in the next
     *  iterations. After the fsi step we put the ALE FSI-dofs to
     *  dirichlet and we solve the ALE again to find the real ALE
     *  displacement.
     */
    int Solve() override;

   private:
    /// map extractor for the fluid-fluid interface dof
    Teuchos::RCP<ALE::UTILS::XFluidFluidMapExtractor> xff_interface_;

  };  // class AleXFFsiWrapper
}  // namespace ADAPTER

#endif  // ADAPTER_ALE_XFFSI_H
