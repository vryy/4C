/*----------------------------------------------------------------------*/
/*!
\brief Static Prestress  analysis
\level 2

\maintainer Fabian Braeu
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_PRESTRESS_H
#define STRTIMINT_PRESTRESS_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_statics.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */

namespace STR
{
  /*====================================================================*/
  /*!
   * \brief Static Prestress analysis
   *
   * This is the prestress version of the static analysis inside the structural dynamics section.
   *
   * Regarding this matter, please direct any complaints to Michael Gee.
   *
   * \authors amaier/kehl
   * \date 03/12
   */

  class TimIntPrestress : public TimIntStatics
  {
   public:
    //! @name Construction
    //@{
    //! Constructor
    TimIntPrestress(const Teuchos::ParameterList& timeparams,  //!< ioflags
        const Teuchos::ParameterList& ioparams,                //!< ioflags
        const Teuchos::ParameterList& sdynparams,              //!< input parameters
        const Teuchos::ParameterList& xparams,                 //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,              //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,                   //!< the solver
        Teuchos::RCP<LINALG::Solver> contactsolver,            //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output          //!< the output
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
    \return bool
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

    //! Update element
    void UpdateStepElement();

    //! Write external and internal forces for dynamic restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);
    //@}

  };  // class TimIntPrestress
}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_STATICS_H
