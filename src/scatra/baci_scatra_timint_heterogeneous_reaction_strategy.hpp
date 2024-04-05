/*----------------------------------------------------------------------*/
/*! \file

 \brief Solution strategy for heterogeneous reactions

   \level 3


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_HPP
#define FOUR_C_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_HPP

#include "baci_config.hpp"

#include "baci_scatra_timint_meshtying_strategy_std.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace SCATRA
{
  /*!
  \brief Solution strategy for heterogeneous reactions

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the solution strategy for heterogeneous reactions.

  */

  class HeterogeneousReactionStrategy : public MeshtyingStrategyStd
  {
   public:
    //! constructor
    explicit HeterogeneousReactionStrategy(SCATRA::ScaTraTimIntImpl* scatratimint);


    //! compute residual terms and their linearizations
    void EvaluateMeshtying() override;

    //! setup meshtying objects
    void SetupMeshtying() override;

    //! initialize meshtying objects
    void InitMeshtying() override;

    /*!
    \brief Evaluate a given condition

     Evaluate terms of your weak formulation on elements marked with a given condition.

    \note The implementation of EvaluateCondition in this class, calls
          \ref DRT::Discretization::EvaluateCondition on the auxiliary
          discretization. Since, this discretization has the dofs of both,
          the volume-bound scalars and the surface bound scalars, every
          term involving one or both kinds of scalars can easily be evaluated.

          This discretization has all scatra (surface+volume) in one dofset.
          This is required by the advanced reaction framework. This is also
          the reason for the fact, that in your input file first the
          surface species have to be listed and then the volume species, since
          here we fix the order of dofs to be first surface, and then volume
          dofs.

    \return void
    \date 08/16
    \author rauch
    */
    void EvaluateCondition(Teuchos::ParameterList& params,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring,
        const int condid) override;

    /*!
    \brief Set state on an discretization hidden in any MeshtyingStrategy.

     This method should encapsulate a standard SetState(...) of the discretizations
     you want to consider. See \ref HeterogeneousReactionStrategy for an example.

    \param nds (in): number of dofset
    \param name (in): Name of data
    \param state (in): vector of some data

    \return void
    \date 12/16
    \author rauch
    */
    void SetState(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state) override;

   private:
    //! the discretization for the reaction
    Teuchos::RCP<DRT::Discretization> discret_;

    //! private copy constructor
    HeterogeneousReactionStrategy(const HeterogeneousReactionStrategy& old);

    //! check assumptions and conventions for heterogeneous reaction simulations
    void HeterogeneousReactionSanityCheck();

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! check if \ref Init() was called
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class MeshtyingStrategyStd
}  // namespace SCATRA



BACI_NAMESPACE_CLOSE

#endif  // SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_H
