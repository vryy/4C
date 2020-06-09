/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu factory class for BACI
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTAFILTERFACTORY_DECL_HPP_
#define MUELU_CONTACTAFILTERFACTORY_DECL_HPP_

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

#include <Xpetra_MapExtractorFactory.hpp>  // why no forward declarations in Xpetra?

namespace MueLu
{
  /*!
    @class ContactAFilterFactory class.
    @brief special factory for segregation master/slave Dofs in matrix A for contact/meshtying
    problems

  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactAFilterFactory : public SingleLevelFactoryBase
  {
#undef MUELU_CONTACTAFILTERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

   public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ContactAFilterFactory(/*const std::string& ename, const FactoryBase* fac*/);

    //! Destructor.
    virtual ~ContactAFilterFactory();

    //! define valid factory parameters
    Teuchos::RCP<const Teuchos::ParameterList> GetValidParameterList(
        const Teuchos::ParameterList& paramList = Teuchos::ParameterList()) const;

    //@}

    //! Input
    //@{

    void DeclareInput(Level& currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level& currentLevel) const;

    //@}

   private:
  };  // class ContactAFilterFactory

}  // namespace MueLu

#define MUELU_CONTACTAFILTERFACTORY_SHORT

#endif /* MUELU_CONTACTAFILTERFACTORY_DECL_HPP_ */
