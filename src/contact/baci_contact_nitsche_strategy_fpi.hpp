/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FPI

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_FPI_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_FPI_HPP

#include "baci_config.hpp"

#include "baci_contact_nitsche_strategy_poro.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

#include <utility>

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  class Element;

  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   */
  class NitscheStrategyFpi : public NitscheStrategyPoro
  {
   public:
    //! Standard constructor
    NitscheStrategyFpi(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::Interface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategyPoro(
              DofRowMap, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          pen_n_(params.get<double>("PENALTYPARAM")),
          weighting_(CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING"))
    {
      if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION") !=
          INPAR::CONTACT::friction_none)
        dserror("NitscheStrategyFpi: No frictional contact implemented for Nitsche FPSCI!");
    }

    //! Shared data constructor
    NitscheStrategyFpi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategyPoro(data_ptr, DofRowMap, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof),
          pen_n_(params.get<double>("PENALTYPARAM")),
          weighting_(CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING"))
    {
      if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION") !=
          INPAR::CONTACT::friction_none)
        dserror("NitscheStrategyFpi: No frictional contact implemented for Nitsche FPSCI!");
    }
    //! Set Contact State and update search tree and normals
    void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    //! The the contact state at local coord of Element cele and compare to the fsi_traction,
    //! return true if contact is evaluated, reture false if FSI is evaluated
    bool CheckNitscheContactState(CONTACT::Element* cele,
        const CORE::LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
        const double& full_fsi_traction,        // stressfluid + penalty
        double& gap                             // gap
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
}  // namespace CONTACT
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_STRATEGY_FPI_H
