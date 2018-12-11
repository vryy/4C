/*!----------------------------------------------------------------------
\file MueLu_ContactSPAggregationFactory_decl.hpp

\brief MueLu contact aggregation factory class for saddle point formulations
\level 2
\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_CONTACTSPAGGREGATIONFACTORY_DECL_HPP_

#ifdef HAVE_MueLu

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"
#include <MueLu_AmalgamationFactory_fwd.hpp>

namespace MueLu
{
  /*!
    @class ContactSPAggregationFactory class.
    @brief Special factory for aggregation of contact problems in saddle-point formulation

    Starting from existing aggregates for the displacement DOFs, additional aggregates for the
    Lagrange multiplier DOFs are formed.

    Input matrix is assumed to be a block matrix of type
    \f[
    A = \left[\begin{array}{cc}
        A_{00} & A_{01}\\
        A_{10} & A_{11}
        \end{array}\right]
    \f]
    where the first and second row represent the displacement DOFs and the constraint equations,
    respectively.
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactSPAggregationFactory : public SingleLevelFactoryBase
  {
#undef MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

   public:
    //! @name Constructors/Destructors.
    //!@{

    /*! \brief Constructor
     *
     * @param[in] aggregatesFact
     * @param[in] amalgFact
     */
    ContactSPAggregationFactory(Teuchos::RCP<const FactoryBase> aggregatesFact = Teuchos::null,
        Teuchos::RCP<const FactoryBase> amalgFact = Teuchos::null);

    //! Destructor.
    virtual ~ContactSPAggregationFactory();
    //!@}

    //! Input
    //!@{

    /*! \brief Validate input parameters
     *
     * \param[in] paramList Parameter list to be validated
     *
     * \note The input \paramList is not used right now. Hence, we just create a new
     * Teuchos::ParameterListm, fill it with valid parameters, and then return it.
     *
     * \return Validated parameter list where missing parameters are set to default values
     */
    Teuchos::RCP<const Teuchos::ParameterList> GetValidParameterList(
        const Teuchos::ParameterList& paramList = Teuchos::ParameterList()) const;

    /*! \brief Specify data that this class needs and the factories that generate that data
     *
     * To create aggregates for a contact problem in saddle-point formulation, we need:
     * - the matrix A (a saddle-point type 2x2 block matrix)
     * - the slave interface DofRowMap containing all displacement DOFs on the slave interface
     *
     * @param[in/out] currentLevel Level to be dealt with
     */
    void DeclareInput(Level& currentLevel) const;

    //!@}

    //! @name Build methods
    //!@{

    /*! \brief Build an object with this factory
     *
     * \param[in/out] currentLevel Level to be processed by this factory
     */
    void Build(Level& currentLevel) const;

    //!@}

   private:
    //! Factory that creates aggregates
    Teuchos::RCP<const FactoryBase> aggregatesFact_;

    //! Factory that creates (Un)Amalgamation info from A
    Teuchos::RCP<const FactoryBase> amalgFact_;

    //! Define which matrix A is used in this factory
    Teuchos::RCP<const FactoryBase> AFact_;

  };  // class ContactSPAggregationFactory

}  // namespace MueLu

// Re-define as it has been undefined at the top
#define MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT

#endif  // HAVE_MueLu

#endif /* MUELU_CONTACTSPAGGREGATIONFACTORY_DECL_HPP_ */
