/*----------------------------------------------------------------------*/
/*! \file

 \brief Solution strategy for heterogeneous reactions

   \level 3


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_HPP
#define FOUR_C_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_meshtying_strategy_std.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace ScaTra
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
    explicit HeterogeneousReactionStrategy(ScaTra::ScaTraTimIntImpl* scatratimint);


    //! compute residual terms and their linearizations
    void evaluate_meshtying() override;

    //! setup meshtying objects
    void setup_meshtying() override;

    //! initialize meshtying objects
    void init_meshtying() override;

    /*!
    \brief Evaluate a given condition

     Evaluate terms of your weak formulation on elements marked with a given condition.

    \note The implementation of evaluate_condition in this class, calls
          \ref Core::FE::Discretization::evaluate_condition on the auxiliary
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
    void evaluate_condition(Teuchos::ParameterList& params,
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring,
        const int condid) override;

    /*!
    \brief Set state on an discretization hidden in any MeshtyingStrategy.

     This method should encapsulate a standard set_state(...) of the discretizations
     you want to consider. See \ref HeterogeneousReactionStrategy for an example.

    \param nds (in): number of dofset
    \param name (in): Name of data
    \param state (in): vector of some data

    \return void
    \date 12/16
    \author rauch
    */
    void set_state(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state) override;

   private:
    //! the discretization for the reaction
    Teuchos::RCP<Core::FE::Discretization> discret_;

    //! private copy constructor
    HeterogeneousReactionStrategy(const HeterogeneousReactionStrategy& old);

    //! check assumptions and conventions for heterogeneous reaction simulations
    void heterogeneous_reaction_sanity_check();

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if setup() was called and is still valid
    bool is_setup() { return issetup_; };

    //! returns true if init(..) was called and is still valid
    bool is_init() { return isinit_; };

    //! check if \ref setup() was called
    void check_is_setup()
    {
      if (not is_setup()) FOUR_C_THROW("setup() was not called.");
    };

    //! check if \ref init() was called
    void check_is_init()
    {
      if (not is_init()) FOUR_C_THROW("init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void set_is_setup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void set_is_init(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class MeshtyingStrategyStd
}  // namespace ScaTra



FOUR_C_NAMESPACE_CLOSE

#endif
