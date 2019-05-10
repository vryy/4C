/*!----------------------------------------------------------------------

\brief Structural time integration with one-step-theta for cell migration

\level 1

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#ifndef STRTIMINT_OST_IMMERSED_H
#define STRTIMINT_OST_IMMERSED_H

#include "../drt_structure/strtimint_ost.H"

namespace STR
{
  /*!
   * \brief Structural time integration with one-step-theta for cell migration
   *
   * This time integrator allows damping and nonlinear mass together without
   * touching the core implementation of the (old) structural time integration.
   * This is necessary because the elements in artificial solid phase of the ECM
   * do not assemble to the pseudo-boundary between asrtifical domain and physical
   * domain. So we have to allow to determine the mass matrix in case we perform
   * selective assembly for all other rhs contributions.
   *
   * Note that a poroelast background discretization needs both DAMPING=Material
   * AND MASSLIN=Standard in --/STRUCTURAL DYNAMIC section.
   *
   * \author rauch
   * \date 05/17
   */
  class TimIntOneStepThetaImmersed : public TimIntOneStepTheta
  {
   public:
    //! Constructor
    TimIntOneStepThetaImmersed(const Teuchos::ParameterList& timeparams,  //!< time parameters
        const Teuchos::ParameterList& ioparams,                           //!< ioflags
        const Teuchos::ParameterList& sdynparams,                         //!< input parameters
        const Teuchos::ParameterList& xparams,                            //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,      //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,           //!< the solver
        Teuchos::RCP<LINALG::Solver> contactsolver,    //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Init(const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver);

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

   protected:
    //! create fields, based on dofrowmap, whose previous time step values are unimportant
    virtual void createFields();

    //! Clear mass matrix and evaluate mass matrix again
    virtual void DetermineMass();

    //! check whether the initial conditions are fulfilled */
    virtual void NonlinearMassSanityCheck(
        Teuchos::RCP<const Epetra_Vector> fext,          ///< external forces
        Teuchos::RCP<const Epetra_Vector> dis,           ///< displacements
        Teuchos::RCP<const Epetra_Vector> vel,           ///< velocities
        Teuchos::RCP<const Epetra_Vector> acc,           ///< accelerations
        const Teuchos::ParameterList* sdynparams = NULL  ///< structural dynamics parameter list
        ) const;

  };  // class TimIntOneStepThetaCellmigration
}  // namespace STR



#endif  // #ifndef STRTIMINT_OST_IMMERSED_H
