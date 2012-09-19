/*
 * MueLu_ContactAFilterFactory_decl.hpp
 *
 *  Created on: Jan 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTAFILTERFACTORY_DECL_HPP_
#define MUELU_CONTACTAFILTERFACTORY_DECL_HPP_

#ifdef HAVE_MueLu

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

#include <Xpetra_MapExtractorFactory.hpp> // why no forward declarations in Xpetra?

namespace MueLu {


#if 1
/*!
  @class ContactAFilterFactory class.
  @brief special factory for segregation master/slave Dofs in matrix A for contact/meshtying problems

*/

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
class ContactAFilterFactory : public SingleLevelFactoryBase {
#undef MUELU_CONTACTAFILTERFACTORY_SHORT
  #include "MueLu_UseShortNames.hpp"

  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass; // TODO move me to ShortNames...

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  ContactAFilterFactory(const std::string& ename, const FactoryBase* fac);

  //! Destructor.
  virtual ~ContactAFilterFactory();
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

  //bool IsGlobalId(Teuchos::RCP<const Map> & map, GlobalOrdinal gid) const;

  std::string        varName_;   ///< name of input and output variable
  const FactoryBase* factory_;   ///< generating factory of input variable
  //const Scalar       threshold_; ///< threshold parameter

  //RCP<const MapExtractorClass> mapextractor_;   ///< user given map extractor (for finest level only)


}; // class ContactAFilterFactory

#else
  /*!
    @class ContactAFilterFactory class.
    @brief special factory for segregation master/slave Dofs in matrix A for contact/meshtying problems

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ContactAFilterFactory : public SingleLevelFactoryBase {
#undef MUELU_CONTACTAFILTERFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass; // TODO move me to ShortNames...

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactAFilterFactory(const std::string& ename, const FactoryBase* fac, Teuchos::RCP<const MapExtractorClass>& rangeMaps);

    //! Destructor.
    virtual ~ContactAFilterFactory();
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

    //bool IsGlobalId(Teuchos::RCP<const Map> & map, GlobalOrdinal gid) const;

    std::string        varName_;   ///< name of input and output variable
    const FactoryBase* factory_;   ///< generating factory of input variable
    const Scalar       threshold_; ///< threshold parameter

    RCP<const MapExtractorClass> mapextractor_;   ///< user given map extractor (for finest level only)


  }; // class ContactAFilterFactory
#endif

} // namespace MueLu


#define MUELU_CONTACTAFILTERFACTORY_SHORT
#endif // HAVE_MueLu
#endif /* MUELU_CONTACTAFILTERFACTORY_DECL_HPP_ */
