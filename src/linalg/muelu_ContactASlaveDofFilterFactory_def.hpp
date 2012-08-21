/*
 * muelu_ContactASlaveDofFilger_def.hpp
 *
 *  Created on: Aug 2, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTASLAVEDOFFILTERFACTORY_DEF_HPP_
#define MUELU_CONTACTASLAVEDOFFILTERFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "muelu_ContactASlaveDofFilterFactory_decl.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ContactASlaveDofFilterFactory(RCP<const FactoryBase> AFact)
    : AFact_(AFact)
  {

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ContactASlaveDofFilterFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A",AFact_.get(),this);
    currentLevel.DeclareInput("SlaveDofMap", MueLu::NoFactory::get(),this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO
    typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactoryClass;
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;

    Monitor m(*this, "ContactASlaveDofFilter factory");

    // extract slave dof map from Level class
    // the slave dof map contains all global DOF ids that shall be handled as Dirichlet boundaries
    // i.e. zero out the corresponding lines in matrix A and put a 1.0 on the diagonal
    RCP<const MapClass> slavedofmap = currentLevel.Get<RCP<const MapClass> >("SlaveDofMap", MueLu::NoFactory::get());

    // extract the original matrix A
    RCP<OOperator> Ain = currentLevel.Get< RCP<OOperator> >("A", AFact_.get());

    /*RCP<Vector> blockVectorRowMap = VectorFactoryClass::Build(Ain->getRowMap());
    blockVectorRowMap->putScalar(-1.0);         // -1.0 denotes that this Dof is not a master or slave dof*/

    // use master map as source map (since all GIDs are uniquely owned by its corresponding proc
    // use column map of current matrix Ain as target map
    // define Xpetra::Import object
/*
    RCP<Vector> blockVectorMaster = VectorFactoryClass::Build(mastermap); blockVectorMaster->putScalar(0);
    RCP<Vector> blockVectorSlave  = VectorFactoryClass::Build(slavemap);  blockVectorSlave->putScalar(1);
    mapextractor->InsertVector(blockVectorMaster, 0, blockVectorRowMap);
    mapextractor->InsertVector(blockVectorSlave,  1, blockVectorRowMap);*/


    /*RCP<const Import> importer = ImportFactory::Build(Ain->getRowMap(), Ain->getColMap());
    RCP<Vector> blockVectorColMap = VectorFactoryClass::Build(Ain->getColMap());
    blockVectorColMap->putScalar(-1.0);         // -1.0 denotes that this Dof is not a master or slave dof
    blockVectorColMap->doImport(*blockVectorRowMap,*importer,Xpetra::INSERT);*/

    // create new empty Operator
    RCP<CrsOOperator> Aout = rcp(new CrsOOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile)); //FIXME

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++) {
        // get global row id
        GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id

        if(slavedofmap->isNodeGlobalElement(grid) == false) {

          size_t nnz = Ain->getNumEntriesInLocalRow(row);
          Teuchos::ArrayView<const LocalOrdinal> indices;
          Teuchos::ArrayView<const Scalar> vals;
          Ain->getLocalRowView(row, indices, vals);

          TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");

          // just copy all values in output
          Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
          Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

          for(size_t i=0; i<(size_t)indices.size(); i++) {
              GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
              indout [i] = gcid;
              valout [i] = vals[i];
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

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<OOperator>(Aout), this);
  }

} // namespace MueLu

#endif // HAVE_MueLu

#endif /* MUELU_CONTACTASLAVEDOFFILTERFACTORY_DEF_HPP_ */
