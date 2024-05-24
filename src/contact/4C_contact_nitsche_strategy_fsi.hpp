/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FSI

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_FSI_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_FSI_HPP

#include "4C_config.hpp"

#include "4C_contact_nitsche_strategy.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class Element;

  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   */
  class NitscheStrategyFsi : public NitscheStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategyFsi(const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::Interface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(
              dof_row_map, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          pen_n_(params.get<double>("PENALTYPARAM")),
          weighting_(CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING"))
    {
      if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION") !=
          INPAR::CONTACT::friction_none)
        FOUR_C_THROW("NitscheStrategyFsi: No frictional contact implemented for Nitsche FSCI!");
    }

    //! Shared data constructor
    NitscheStrategyFsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(data_ptr, dof_row_map, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof),
          pen_n_(params.get<double>("PENALTYPARAM")),
          weighting_(CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING"))
    {
      if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION") !=
          INPAR::CONTACT::friction_none)
        FOUR_C_THROW("NitscheStrategyFsi: No frictional contact implemented for Nitsche FSCI!");
    }
    //! Evaluate and apply RHS and Stiffness Matrix for Contact
    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor) override;

    //! Set Contact State and update search tree and normals
    void set_state(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    //! The the contact state at local coord of Element cele and compare to the fsi_traction,
    //! return true if contact is evaluated, reture false if FSI is evaluated
    bool check_nitsche_contact_state(CONTACT::Element* cele,
        const CORE::LINALG::Matrix<2, 1>& xsi,  ///< local coord on the ele element
        const double& full_fsi_traction,        ///< stressfluid + penalty
        double& gap                             ///< gap
    );

   protected:
    //! Update search tree and normals
    void DoContactSearch();

   private:
    //! Nitsche normal penalty parameter
    double pen_n_;
    //! Nitsche weighting strategy
    INPAR::CONTACT::NitscheWeighting weighting_;
  };

  namespace UTILS
  {
    //! The the contact state at local coord of Element cele and compare to the fsi_traction,
    //! return true if contact is evaluated, reture false if FSI is evaluated
    bool check_nitsche_contact_state(CONTACT::Interface& contactinterface,  ///< Contact interface
        const double& pen_n,                         ///< Nitsche normal penalty parameter
        INPAR::CONTACT::NitscheWeighting weighting,  ///< Nitsche interface stress weighting
        CONTACT::Element* cele,                      ///< the contact element
        const CORE::LINALG::Matrix<2, 1>& xsi,       ///< local coord on the ele element
        const double& full_fsi_traction,             ///< stressfluid + penalty
        double& gap                                  ///< gap
    );
  }  // namespace UTILS
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
