/*!----------------------------------------------------------------------

\brief MueLu contact info factory class
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTINFOFACTORY_DECL_HPP_
#define MUELU_CONTACTINFOFACTORY_DECL_HPP_

#ifdef HAVE_MueLu


#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

#include <Xpetra_MapExtractorFactory.hpp>  // why no forward declarations in Xpetra?

namespace MueLu
{
  /*!
    @class ContactAFilterFactory class.
    @brief special factory for exporting nullspace inforation in vtk format. Can be used together
    with the MueLuAggregationExportFactory. Extend aggregation information by using e.g. "cat
    output0.vtk agg_info_0.vtk >> aggregation0.vtk"


  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactInfoFactory : public TwoLevelFactoryBase
  {
#undef MUELU_CONTACTINFOFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        MapExtractorClass;  // TODO move me to ShortNames...

   public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactInfoFactory(std::string filename_prototype,
        Teuchos::RCP<FactoryBase> AFact = Teuchos::null,
        Teuchos::RCP<FactoryBase> nspFact = Teuchos::null);

    //! Destructor.
    virtual ~ContactInfoFactory();
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
    std::string replaceAll(std::string result, const std::string &replaceWhat,
        const std::string &replaceWithWhat) const;

    std::string filename_prototype_;  ///< prototype string for output filename

    Teuchos::RCP<FactoryBase> AFact_;    ///< A factory (needed for maps)
    Teuchos::RCP<FactoryBase> nspFact_;  ///< Nullspace factory

    Teuchos::RCP<const MapExtractorClass>
        mapextractor_;  ///< user given map extractor (for finest level only)


  };  // class ContactInfoFactory

}  // namespace MueLu

#define MUELU_CONTACTINFOFACTORY_SHORT
#endif  // HAVE_MueLu


#endif /* MUELU_CONTACTINFOFACTORY_DECL_HPP_ */
