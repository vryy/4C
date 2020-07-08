/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu contact filter factory class
\level 2

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTASLAVEDOFFILTERFACTORY_DECL_HPP_
#define MUELU_CONTACTASLAVEDOFFILTERFACTORY_DECL_HPP_


#ifdef TRILINOS_Q1_2015

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

//#include <Xpetra_MapExtractorFactory.hpp> // why no forward declarations in Xpetra?

namespace MueLu
{
  /*!
    @class ContactASlaveDofFilterFactory class.
    @brief special factory

  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactASlaveDofFilterFactory : public SingleLevelFactoryBase
  {
#undef MUELU_CONTACTASLAVEDOFFILTERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    // typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass; //
    // TODO move me to ShortNames...

   public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactASlaveDofFilterFactory(/*Teuchos::RCP<const FactoryBase> AFact = Teuchos::null*/);

    //! Destructor.
    virtual ~ContactASlaveDofFilterFactory();

    RCP<const ParameterList> GetValidParameterList() const;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level &currentLevel) const;

    //@}

   private:
  };  // class ContactASlaveDofFilterFactory

}  // namespace MueLu

#define MUELU_CONTACTASLAVEDOFFILTERFACTORY_SHORT
#endif  // TRILINOS_Q1_2015

#endif /* MUELU_CONTACTASLAVEDOFFILTERFACTORY_DECL_HPP_ */
