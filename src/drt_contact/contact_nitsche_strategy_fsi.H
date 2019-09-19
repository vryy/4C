/*---------------------------------------------------------------------*/
/*! \file

\brief Nitsche contact solving strategy for problems with FSI

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_FSI_H
#define CONTACT_NITSCHE_STRATEGY_FSI_H

#include "contact_nitsche_strategy.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace CONTACT
{
  class CoElement;

  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
   For a more general documentation of the involved functions refer to CoAbstract Strategy.

   */
  class CoNitscheStrategyFsi : public CoNitscheStrategy
  {
   public:
    //! Standard constructor
    CoNitscheStrategyFsi(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : CoNitscheStrategy(DofRowMap, NodeRowMap, params, interface, dim, comm, alphaf, maxdof),
          pen_n_(params.get<double>("PENALTYPARAM")),
          weighting_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING"))
    {
      if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION") !=
          INPAR::CONTACT::friction_none)
        dserror("CoNitscheStrategyFsi: No frictional contact implemented for Nitsche FSCI!");
    }

    //! Shared data constructor
    CoNitscheStrategyFsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : CoNitscheStrategy(
              data_ptr, DofRowMap, NodeRowMap, params, interface, dim, comm, alphaf, maxdof),
          pen_n_(params.get<double>("PENALTYPARAM")),
          weighting_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING"))
    {
      if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION") !=
          INPAR::CONTACT::friction_none)
        dserror("CoNitscheStrategyFsi: No frictional contact implemented for Nitsche FSCI!");
    }

    /*!
     \brief Destructor

     */
    virtual ~CoNitscheStrategyFsi() {}

    //! Evaluate and apply RHS and Stiffness Matrix for Contact
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor);

    //! Set Contact State and update search tree and normals
    void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec);

    //! The the contact state at local coord of CoElement cele and compare to the fsi_traction,
    //! return true if contact is evaluated, reture false if FSI is evaluated
    bool CheckNitscheContactState(CONTACT::CoElement* cele,
        const LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
        const double& full_fsi_traction,  // stressfluid + penalty
        double& gap                       // gap
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
#endif  //  CONTACT_NITSCHE_STRATEGY_FSI_H
