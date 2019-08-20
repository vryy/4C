/*----------------------------------------------------------------------------*/
/** \file
\brief xcontact level-set one step theta algorithm

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------------*/
#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_TIMINT_OST_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_TIMINT_OST_H_

#include "xcontact_levelset_algorithm.H"
#include "../drt_scatra/scatra_timint_ost.H"

namespace XCONTACT
{
  namespace LEVELSET
  {
    namespace TIMINT
    {
      class OneStepTheta : public virtual XCONTACT::LEVELSET::Algorithm,
                           public virtual SCATRA::TimIntOneStepTheta
      {
       public:
        /// constructor
        OneStepTheta(const Teuchos::RCP<DRT::Discretization>& actdis,
            const Teuchos::RCP<LINALG::Solver>& solver,
            const Teuchos::RCP<Teuchos::ParameterList>& params,
            const Teuchos::RCP<Teuchos::ParameterList>& sctratimintparams,
            const Teuchos::RCP<Teuchos::ParameterList>& extraparams,
            const Teuchos::RCP<IO::DiscretizationWriter>& output);

        /// destructor
        virtual ~OneStepTheta(){};

        /// call the initialize routines of the base classes
        virtual void Init();

        /// call the setup routines of the base classes
        virtual void Setup();

        /// read restart data
        virtual void ReadRestart(
            const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

        /// interpolate phi to intermediate time n+theta with 0<theta<1
        Teuchos::RCP<Epetra_Vector> Phinptheta(const double theta_inter) const;

        /// interpolate phidt to intermediate time n+theta with 0<theta<1
        Teuchos::RCP<Epetra_Vector> Phidtnptheta(const double theta_inter) const;

       protected:
        /// Print information about current time step to screen [derived]
        virtual void PrintTimeStepInfo();

        /// calculate consistent initial scalar time derivatives in compliance with initial scalar
        /// field
        virtual void CalcInitialTimeDerivative();

        /// additional predictor not intended for level-set methods
        virtual inline void ExplicitPredictor() { return; };

        /// Set the part of the righthandside belonging to the last timestep.
        virtual void SetOldPartOfRighthandside();

        /** \brief update state vectors
         *
         *  current solution becomes old solution of next time step */
        virtual void UpdateState();

        /** \brief update the solution after Solve()
         *
         *  extended version for coupled level-set problems including
         *  reinitialization if required */
        virtual void Update(const int num = 0);

        /// update phi within the reinitialization loop
        virtual void UpdateReinit();

        /// set the solution state in the current active time integration
        virtual void SetState(const Epetra_Vector& sol_state);

      };  // class OneStepTheta
    }     // namespace TIMINT
  }       // namespace LEVELSET
}  // namespace XCONTACT


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_TIMINT_OST_H_ */
