/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu repartition algorithm for contact
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_
#define MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

#include "MueLu_ContactSPRepartitionInterface_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu
{
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<const ParameterList> ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node,
      LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const
  {
    Teuchos::RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set<Teuchos::RCP<const FactoryBase>>(
        "A", Teuchos::null, "Factory of the matrix A");
    validParamList->set<Teuchos::RCP<const FactoryBase>>("AmalgamatedPartition", Teuchos::null,
        "(advanced) Factory generating the AmalgamatedPartition (e.g. an IsorropiaInterface)");

    return validParamList;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(
      Level& currentLevel) const
  {
    Input(currentLevel, "A");
    Input(currentLevel, "AmalgamatedPartition");
  }  // DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(
      Level& level) const
  {
    FactoryMonitor m(*this, "Build", level);
    level.print(GetOStream(Statistics0, 0));
    // extract blocked operator A from current level
    Teuchos::RCP<Matrix> A = Get<Teuchos::RCP<Matrix>>(level, "A");
    Teuchos::RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
    const int myRank = comm->getRank();

    Teuchos::RCP<Xpetra::Vector<GO, LO, GO, NO>> decomposition =
        Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

    // fill decomposition vector
    for (LO i = 0; i < Teuchos::as<LO>(decomposition->getMap()->getNodeNumElements()); i++)
    {
      decompEntries[i] = myRank;  // no repartitioning for Lagr. multipliers
    }

    // Teuchos::RCP<Teuchos::FancyOStream> fos =
    // Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout)); decomposition->describe(*fos,
    // Teuchos::VERB_EXTREME);

    Set(level, "Partition", decomposition);

  }  // Build()



}  // namespace MueLu

#endif  // HAVE_MueLuContact

#endif /* MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_ */
