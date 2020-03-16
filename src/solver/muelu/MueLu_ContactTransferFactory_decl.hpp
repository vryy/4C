/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu transfer factory class for contact
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTTRANSFERFACTORY_DECL_HPP_
#define MUELU_CONTACTTRANSFERFACTORY_DECL_HPP_


#ifdef HAVE_MueLuContact


#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
//#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

#include <Xpetra_MapExtractorFactory.hpp>  // why no forward declarations in Xpetra?

namespace MueLu
{
  /*!
    @class ContactAFilterFactory class.
    @brief special factory for segregation master/slave Dofs in matrix A for contact/meshtying
    problems

  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactTransferFactory : public TwoLevelFactoryBase
  {
#undef MUELU_CONTACTTRANSFERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        MapExtractorClass;  // TODO move me to ShortNames...

   public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactTransferFactory(Teuchos::RCP<FactoryBase> PtentFact = Teuchos::null);

    //! Destructor.
    virtual ~ContactTransferFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level &fineLevel, Level &coarseLevel) const;

    //@}

   private:
    std::string varName_;         ///< name of input and output variable
    const FactoryBase *factory_;  ///< generating factory of input variable

    Teuchos::RCP<FactoryBase> PtentFact_;  ///< tentative P Factory
    Teuchos::RCP<FactoryBase> AFact_;      ///< A factory (needed for maps)

    Teuchos::RCP<const MapExtractorClass>
        mapextractor_;  ///< user given map extractor (for finest level only)


  };  // class ContactTransferFactory

}  // namespace MueLu

#define MUELU_CONTACTTRANSFERFACTORY_SHORT
#endif  // HAVE_MueLuContact

#endif /* MUELU_CONTACTTRANSFERFACTORY_DECL_HPP_ */
