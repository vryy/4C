/*
 * MueLu_IterationAFactory_def.hpp
 *
 *  Created on: Jan 10, 2013
 *      Author: tobias
 */

#ifndef MUELU_ITERATIONAFACTORY_DEF_HPP_
#define MUELU_ITERATIONAFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "MueLu_IterationAFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::IterationAFactory(const std::string mapName,const Teuchos::RCP<const FactoryBase> & mapFact)
      : mapName_(mapName), mapFact_(mapFact)
  {

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~IterationAFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");

    currentLevel.DeclareInput(mapName_, mapFact_.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {

    Monitor m(*this, "IterationAFactory");

    // extract dof map from Level class
    // the dof map contains all global DOF ids that shall be handled as Dirichlet boundaries
    // i.e. zero out the corresponding lines in matrix A and put a 1.0 on the diagonal
    Teuchos::RCP<const Map> fixeddofmap = currentLevel.Get<Teuchos::RCP<const Map> >(mapName_, mapFact_.get());

    // extract the original matrix A
    Teuchos::RCP<Matrix> Ain = Get< Teuchos::RCP<Matrix> >(currentLevel, "A");  // corresponding Get call to Input

    // create new empty Operator
    Teuchos::RCP<CrsMatrixWrap> Aout = Teuchos::rcp(new CrsMatrixWrap(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile));

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++) {
        // get global row id
        GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id

        if(fixeddofmap->isNodeGlobalElement(grid) == false) {

          // extract information from current row
          //size_t nnz = Ain->getNumEntriesInLocalRow(row);
          Teuchos::ArrayView<const LocalOrdinal> indices;
          Teuchos::ArrayView<const Scalar> vals;
          Ain->getLocalRowView(row, indices, vals);

          //TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::IterationAFactory::Build: number of nonzeros not equal to number of indices? Error.");

          // just copy all values in output
          Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
          Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

          for(size_t i=0; i<(size_t)indices.size(); i++) {
              GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
              indout [i] = gcid;
              valout[i]  = vals[i];
          }

          Aout->insertGlobalValues(grid, indout.view(0,indout.size()), valout.view(0,valout.size()));
        } else {
          Aout->insertGlobalValues(grid, Teuchos::tuple<GlobalOrdinal>(grid),Teuchos::tuple<Scalar>(1.0));
        }
    }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    // copy block size information
    Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

    GetOStream(Statistics0, 0) << "Nonzeros in A  (input): " << Ain->getGlobalNumEntries() << ", Nonzeros introducing artificial DC bdries " << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(Aout), this);
  }

} // namespace MueLu

#endif // HAVE_MueLu

#endif /* MUELU_ITERATIONAFACTORY_DEF_HPP_ */
