/*
 * muelu_ContactASlaveDofFilter_decl.hpp
 *
 *  Created on: Aug 2, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTASLAVEDOFFILTERFACTORY_DECL_HPP_
#define MUELU_CONTACTASLAVEDOFFILTERFACTORY_DECL_HPP_


#ifdef HAVE_MueLu

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

#include <Xpetra_MapExtractorFactory.hpp> // why no forward declarations in Xpetra?

namespace MueLu {

  /*!
    @class ContactASlaveDofFilterFactory class.
    @brief special factory

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ContactASlaveDofFilterFactory : public SingleLevelFactoryBase {
#undef MUELU_CONTACTASLAVEDOFFILTERFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass; // TODO move me to ShortNames...

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactASlaveDofFilterFactory(RCP<const FactoryBase> AFact = Teuchos::null);

    //! Destructor.
    virtual ~ContactASlaveDofFilterFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}

  private:

    RCP<const FactoryBase> AFact_;
    //bool IsGlobalId(Teuchos::RCP<const Map> & map, GlobalOrdinal gid) const;

    //RCP<const MapExtractorClass> mapextractor_;   ///< user given map extractor (for finest level only)


  }; // class ContactASlaveDofFilterFactory

} // namespace MueLu

#define MUELU_CONTACTASLAVEDOFFILTERFACTORY_SHORT
#endif // HAVE_MueLu

#endif /* MUELU_CONTACTASLAVEDOFFILTERFACTORY_DECL_HPP_ */
