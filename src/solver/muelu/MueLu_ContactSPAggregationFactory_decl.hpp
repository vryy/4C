/*!----------------------------------------------------------------------
\file MueLu_ContactSPAggregationFactory_decl.hpp

\brief MueLu contact aggregation factory class
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_CONTACTSPAGGREGATIONFACTORY_DECL_HPP_

#ifdef HAVE_MueLu

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"
#include <MueLu_AmalgamationFactory_fwd.hpp>


//#include <Xpetra_MapExtractorFactory.hpp> // why no forward declarations in Xpetra?

namespace MueLu {

  /*!
    @class ContactSPAggregationFactory class.
    @brief special factory

  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactSPAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

    //typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass; // TODO move me to ShortNames...

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactSPAggregationFactory(Teuchos::RCP<const FactoryBase> aggregatesFact = Teuchos::null, Teuchos::RCP<const FactoryBase> amalgFact = Teuchos::null);

    //! Destructor.
    virtual ~ContactSPAggregationFactory();
    //@}

    //! Input
    //@{

    Teuchos::RCP<const Teuchos::ParameterList> GetValidParameterList(const Teuchos::ParameterList& paramList = Teuchos::ParameterList()) const;

    void DeclareInput(Level &currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}

  private:

    Teuchos::RCP<const FactoryBase> aggregatesFact_; //! Factory that creates aggregates
    Teuchos::RCP<const FactoryBase> amalgFact_;      //! Factory that (Un)Amalgamation info from A
    Teuchos::RCP<const FactoryBase> AFact_;          //! Define which matrix A is used in this factory

    //bool IsGlobalId(Teuchos::RCP<const Map> & map, GlobalOrdinal gid) const;

    //Teuchos::RCP<const MapExtractorClass> mapextractor_;   ///< user given map extractor (for finest level only)


  }; // class ContactSPAggregationFactory

} // namespace MueLu

#define MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
#endif // HAVE_MueLu

#endif /* MUELU_CONTACTSPAGGREGATIONFACTORY_DECL_HPP_ */
