/*
 * MueLu_ContactSPRepartitionInterface_def.hpp
 *
 *  Created on: 11 Sep 2013
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_
#define MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q3_2013

#include "MueLu_ContactSPRepartitionInterface_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu {

 template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",                    Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("AmalgamatedPartition", Teuchos::null, "(advanced) Factory generating the AmalgamatedPartition (e.g. an IsorropiaInterface)");

    return validParamList;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "AmalgamatedPartition");
  } //DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &level) const {
    FactoryMonitor m(*this, "Build", level);
    level.print(GetOStream(Statistics0,0));
    // extract blocked operator A from current level
    RCP<Matrix> A = Get< RCP<Matrix> >     (level, "A");
    RCP<const Teuchos::Comm< int > > comm = A->getRowMap()->getComm();
    const int myRank = comm->getRank();

    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

    // fill decomposition vector
    for(LO i = 0; i<Teuchos::as<LO>(decomposition->getMap()->getNodeNumElements()); i++) {
      decompEntries[i] = myRank; // no repartitioning for Lagr. multipliers
    }

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    //decomposition->describe(*fos, Teuchos::VERB_EXTREME);

    Set(level, "Partition", decomposition);

  } //Build()



} //namespace MueLu

#endif
#endif

#endif /* MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_ */
