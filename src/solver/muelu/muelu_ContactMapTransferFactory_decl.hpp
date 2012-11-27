/*
 * MueLu_ContactMapTransferFactory_decl.hpp
 *
 *  Created on: Aug 2, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTMAPTRANSFERFACTORY_DECL_HPP_
#define MUELU_CONTACTMAPTRANSFERFACTORY_DECL_HPP_

#ifdef HAVE_MueLu


#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

#include <Xpetra_MapExtractorFactory.hpp> // why no forward declarations in Xpetra?

namespace MueLu {

  /*!
    @class ContactMapFilterFactory class.
    @brief special factory

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ContactMapTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_CONTACTMAPTRANSFERFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass; // TODO move me to ShortNames...

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactMapTransferFactory(std::string mapName, Teuchos::RCP<const FactoryBase> PtentFact = Teuchos::null, Teuchos::RCP<const FactoryBase> mapFact = Teuchos::null);

    //! Destructor.
    virtual ~ContactMapTransferFactory();
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

    std::string              mapName_;   ///< name of input and output variable
    RCP<const FactoryBase>   PtentFact_; ///< tentative P Factory
    RCP<const FactoryBase>   mapFact_;   ///< generating factory of input variable

  }; // class ContactMapTransferFactory

} // namespace MueLu

#define MUELU_CONTACTMAPTRANSFERFACTORY_SHORT
#endif // HAVE_MueLu

#endif /* MUELU_CONTACTMAPTRANSFERFACTORY_DECL_HPP_ */
