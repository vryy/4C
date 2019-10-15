/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu amalgamation factory for meshtying
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DECL_HPP_
#define MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DECL_HPP_

#ifdef HAVE_MueLu

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_StridedMap_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"
#include <MueLu_AmalgamationInfo_fwd.hpp>
#include <MueLu_AmalgamationFactory_fwd.hpp>


namespace MueLu
{
  /*!
    @class MeshtyingSPAmalgamationFactory class.
    @brief special factory

  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class MeshtyingSPAmalgamationFactory : public SingleLevelFactoryBase
  {
#undef MUELU_MESHTYINGSPAMALGAMATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"


   public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    MeshtyingSPAmalgamationFactory();

    //! Destructor.
    virtual ~MeshtyingSPAmalgamationFactory();
    //@}

    //! Input
    //@{

    Teuchos::RCP<const ParameterList> GetValidParameterList(
        const ParameterList& paramList = ParameterList()) const;

    void DeclareInput(Level& currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level& currentLevel) const;

    //@}

   private:
    // amalgamation information
    mutable Teuchos::RCP<std::map<GlobalOrdinal, std::vector<GlobalOrdinal>>> nodegid2dofgids_;


  };  // class MeshtyingSPAmalgamationFactory

}  // namespace MueLu

#define MUELU_MESHTYINGSPAMALGAMATIONFACTORY_SHORT
#endif  // HAVE_MueLu


#endif /* MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DECL_HPP_ */
