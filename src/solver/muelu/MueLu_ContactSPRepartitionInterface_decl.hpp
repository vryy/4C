/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu repartition algorithm for contact
\level 2

*----------------------------------------------------------------------*/
#ifndef MUELU_CONTACTSPREPARTITIONINTERFACE_DECL_HPP_
#define MUELU_CONTACTSPREPARTITIONINTERFACE_DECL_HPP_

#ifdef TRILINOS_Q1_2015

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include "MueLu_NodeDefinition.hpp"

namespace MueLu
{
  /*!
    @class ContactSPRepartitionInterface
    @brief Helper class which transforms an "AmalgamatedPartition" array to an unamalgamated
    "Partition".

    This class is meant to be used with IsorropiaInterface which in general provides the amalgamated
    partition information and an AmalgamationFactory which defines the amalgamation/unamalgamation
    process. The output is a "Partition" (unamalgamated) which can be used by the RepartitionFactory
    class.

    Input: matrix A, unamalgamation information (that corresponds to matrix A)
  */

  // FIXME: this class should not be templated
  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal,
      class Node = KokkosSerialNode,
      class LocalMatOps =
          typename KokkosClassic::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class ContactSPRepartitionInterface : public SingleLevelFactoryBase
  {
#undef MUELU_CONTACTSPREPARTITIONINTERFACE_SHORT
#include "MueLu_UseShortNames.hpp"

   public:
    //! @name Constructors/Destructors
    //@{

    //! Constructor
    ContactSPRepartitionInterface() {}

    //! Destructor
    virtual ~ContactSPRepartitionInterface() {}
    //@}

    Teuchos::RCP<const ParameterList> GetValidParameterList(
        const ParameterList& paramList = ParameterList()) const;

    //! @name Input
    //@{
    void DeclareInput(Level& level) const;
    //@}

    //! @name Build methods.
    //@{
    void Build(Level& level) const;

    //@}



   private:
  };  // class ContactSPRepartitionInterface

}  // namespace MueLu

#define MUELU_CONTACTSPREPARTITIONINTERFACE_SHORT

#endif  // TRILINOS_Q1_2015
#endif  /* MUELU_CONTACTSPREPARTITIONINTERFACE_DECL_HPP_ */
