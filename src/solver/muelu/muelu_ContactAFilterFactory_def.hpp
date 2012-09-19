/*
 * MueLu_ContactAFilterFactory_def.hpp
 *
 *  Created on: Jan 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTAFILTERFACTORY_DEF_HPP_
#define MUELU_CONTACTAFILTERFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "muelu_ContactAFilterFactory_decl.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#if 1
namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ContactAFilterFactory(const std::string& ename, const FactoryBase* fac)
    : varName_(ename), factory_(fac)//, threshold_(0.0)
  {

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ContactAFilterFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput(varName_,factory_,this);
    //currentLevel.DeclareInput("SegAMapExtractor", MueLu::NoFactory::get(),this);
    currentLevel.DeclareInput("SlaveDofMap", MueLu::NoFactory::get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO
    typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactoryClass;

    Monitor m(*this, "A filtering (contact)");

    // fetch map with slave Dofs from Level
    RCP<const Map> slaveDofMap = currentLevel.Get< RCP<const Map> >("SlaveDofMap",MueLu::NoFactory::get());

    RCP<OOperator> Ain = currentLevel.Get< RCP<OOperator> >(varName_, factory_);

    RCP<Vector> blockVectorRowMap = VectorFactoryClass::Build(Ain->getRowMap());
    blockVectorRowMap->putScalar(-1.0);         // -1.0 denotes that this Dof is not slave DOF

    // use master map as source map (since all GIDs are uniquely owned by its corresponding proc
    // use column map of current matrix Ain as target map
    // define Xpetra::Import object
    RCP<Vector> blockVectorSlave  = VectorFactoryClass::Build(slaveDofMap);  blockVectorSlave->putScalar(1);

    RCP<const Import> importer = ImportFactory::Build(/*Ain->getRowMap()*/slaveDofMap, Ain->getColMap());
    RCP<Vector> blockVectorColMapData = VectorFactoryClass::Build(Ain->getColMap());
    blockVectorColMapData->putScalar(-1.0);         // -1.0 denotes that this Dof is not slave DOF
    blockVectorColMapData->doImport(*blockVectorSlave,*importer,Xpetra::INSERT);

    // create new empty Operator
    RCP<CrsOOperator> Aout = rcp(new CrsOOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile)); //FIXME

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++) {
        // get global row id
        GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id

        // check in which submap of mapextractor grid belongs to
        bool isSlaveDofRow = slaveDofMap->isNodeGlobalElement(grid);

        size_t nnz = Ain->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ContactAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");

        // just copy all values in output
        Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());
        size_t nNonzeros = 0;

        for(size_t i=0; i<(size_t)indices.size(); i++) {
            GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)

            //bool isSlaveDofCol = slaveDofMap->isNodeGlobalElement(gcid);

            Teuchos::ArrayRCP< const Scalar > colBlockData = blockVectorColMapData->getData(0);

            LocalOrdinal colBlockId = Teuchos::as<LocalOrdinal>(colBlockData[indices[i]]); // LID -> colBlockID

            // colBlockId can be
            // -1:  indices[i] is not a slavel dof
            //  1:  indices[i] is a slave dof
            //  else: error
            if( ( colBlockId == -1 && isSlaveDofRow == false) ||
                ( colBlockId ==  1 && isSlaveDofRow == true )) {
              indout [nNonzeros] = gcid;
              valout [nNonzeros] = vals[i];
              nNonzeros++;
            }
        }
        indout.resize(nNonzeros);
        valout.resize(nNonzeros);

        Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
    }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    // copy block size information
    Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

    GetOStream(Statistics0, 0) << "Nonzeros in " << varName_ << "(input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering " << varName_ << ": " << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set(varName_, Teuchos::rcp_dynamic_cast<OOperator>(Aout), this);
  }
} // namespace MueLu

#else  // OLD version with map extractor

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ContactAFilterFactory(const std::string& ename, const FactoryBase* fac, Teuchos::RCP<const MapExtractorClass>& rangeMaps)
    : varName_(ename), factory_(fac), threshold_(0.0), mapextractor_(rangeMaps)
  {

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ContactAFilterFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput(varName_,factory_,this);
    currentLevel.DeclareInput("SegAMapExtractor", MueLu::NoFactory::get(),this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO
    typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactoryClass;

    Monitor m(*this, "A filtering (contact)");

    if (currentLevel.IsAvailable("SegAMapExtractor", MueLu::NoFactory::get())==false) {
      GetOStream(Runtime0, 0) << "ContactAFilterFactory::Build(): Use user provided map extractor with " << mapextractor_->NumMaps() << " submaps for segregation filter for " << varName_ << std::endl;
      currentLevel.Set("SegAMapExtractor", mapextractor_, MueLu::NoFactory::get());
    }

    // fetch map extractor from level
    RCP<const MapExtractorClass> mapextractor = currentLevel.Get< RCP<const MapExtractorClass> >("SegAMapExtractor",MueLu::NoFactory::get());


#if 1

    // extract maps from mapextractor:                        // TODO: write a more common class for n distinct submaps in mapextractor (not just 2)
    RCP<const Map> mastermap = mapextractor->getMap(0);
    RCP<const Map> slavemap  = mapextractor->getMap(1);

    RCP<OOperator> Ain = currentLevel.Get< RCP<OOperator> >(varName_, factory_);

    RCP<Vector> blockVectorRowMap = VectorFactoryClass::Build(Ain->getRowMap());
    blockVectorRowMap->putScalar(-1.0);         // -1.0 denotes that this Dof is not a master or slave dof

    // use master map as source map (since all GIDs are uniquely owned by its corresponding proc
    // use column map of current matrix Ain as target map
    // define Xpetra::Import object
    RCP<Vector> blockVectorMaster = VectorFactoryClass::Build(mastermap); blockVectorMaster->putScalar(0);
    RCP<Vector> blockVectorSlave  = VectorFactoryClass::Build(slavemap);  blockVectorSlave->putScalar(1);
    mapextractor->InsertVector(blockVectorMaster, 0, blockVectorRowMap);
    mapextractor->InsertVector(blockVectorSlave,  1, blockVectorRowMap);

    RCP<const Import> importer = ImportFactory::Build(Ain->getRowMap(), Ain->getColMap());
    RCP<Vector> blockVectorColMap = VectorFactoryClass::Build(Ain->getColMap());
    blockVectorColMap->putScalar(-1.0);         // -1.0 denotes that this Dof is not a master or slave dof
    blockVectorColMap->doImport(*blockVectorRowMap,*importer,Xpetra::INSERT);

    // create new empty Operator
    RCP<CrsOOperator> Aout = rcp(new CrsOOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile)); //FIXME

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++) {
        // get global row id
        GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id

        // check in which submap of mapextractor grid belongs to
        LocalOrdinal rowBlockId = -Teuchos::ScalarTraits<LocalOrdinal>::one();
        for (size_t bb=0; bb<mapextractor->NumMaps(); bb++) {
            const RCP<const Map> cmap = mapextractor->getMap(bb);
            if (cmap->isNodeGlobalElement(grid)) {
                rowBlockId = bb;
                break;
            }
        }

        // rowBlockId can be
        // -1:  grid is not in one of the submaps of mapextractor
        // >=0: grid is in one of the submaps of mapextractor


        size_t nnz = Ain->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");


        // just copy all values in output
        Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());
        size_t nNonzeros = 0;

        for(size_t i=0; i<(size_t)indices.size(); i++) {
            GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)

            Teuchos::ArrayRCP< const Scalar > colBlockData = blockVectorColMap->getData(0);

            LocalOrdinal colBlockId = Teuchos::as<LocalOrdinal>(colBlockData[indices[i]]); // LID -> colBlockID

            // colBlockId can be
            // -1:  indices[i] is not in one of the submaps of mapextractor
            // >=0: indices[i] is in one of the submaps of mapextractor

            // check if mapextractor has grid
            if (rowBlockId >= 0) {
                // check if we have to filter columns
                if((colBlockId > -1 && colBlockId == rowBlockId) ||  // skip all entries with rowBlockId != colBlockId if colBlockId >= 0
                    colBlockId == -1) {
                    indout [nNonzeros] = gcid;
                    valout [nNonzeros] = vals[i];
                    nNonzeros++;
                }
            } else {
                indout [nNonzeros] = gcid;
                valout [nNonzeros] = vals[i];
                nNonzeros++;
            }
        }
        indout.resize(nNonzeros);
        valout.resize(nNonzeros);
	
        Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
    }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    // copy block size information
    Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());
    
    GetOStream(Statistics0, 0) << "Nonzeros in " << varName_ << "(input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering " << varName_ << " (parameter: " << threshold_ << "): " << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set(varName_, Teuchos::rcp_dynamic_cast<OOperator>(Aout), this);


#else
    // extract maps from mapextractor:
    RCP<const Map> mastermap = mapextractor->getMap(0);
    RCP<const Map> slavemap  = mapextractor->getMap(1);
    RCP<const Map> innermap  = mapextractor->getMap(2);

    RCP<OOperator> Ain = currentLevel.Get< RCP<OOperator> >(varName_, factory_);

    Monitor m(*this, "A filtering (contact)");


    // use master map as source map (since all GIDs are uniquely owned by its corresponding proc
    // use column map of current matrix Ain as target map
    // define Xpetra::Import object
    RCP<Vector> blockVectorRowMap = VectorFactoryClass::Build(Ain->getRowMap());
    RCP<Vector> blockVectorMaster = VectorFactoryClass::Build(mastermap); blockVectorMaster->putScalar(0);
    RCP<Vector> blockVectorSlave  = VectorFactoryClass::Build(slavemap);  blockVectorSlave->putScalar(1);
    RCP<Vector> blockVectorInner  = VectorFactoryClass::Build(innermap);  blockVectorInner->putScalar(2);
    mapextractor->InsertVector(blockVectorMaster, 0, blockVectorRowMap);
    mapextractor->InsertVector(blockVectorSlave,  1, blockVectorRowMap);
    mapextractor->InsertVector(blockVectorInner,  2, blockVectorRowMap);

    RCP<const Import> importer = ImportFactory::Build(Ain->getRowMap(), Ain->getColMap());
    RCP<Vector> blockVectorColMap = VectorFactoryClass::Build(Ain->getColMap());
    blockVectorColMap->doImport(*blockVectorRowMap,*importer,Xpetra::INSERT);


    // create new empty Operator
    RCP<CrsOOperator> Aout = rcp(new CrsOOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile)); //FIXME

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++)
      {

        // get global row id
        GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id

        // check in which submap of mapextractor grid belongs to
        LocalOrdinal rowBlockId = -Teuchos::ScalarTraits<LocalOrdinal>::one();
        for (size_t bb=0; bb<mapextractor->NumMaps(); bb++) {
          const RCP<const Map> cmap = mapextractor->getMap(bb);
          if (cmap->isNodeGlobalElement(grid)) {
            rowBlockId = bb;
            break;
          }
        }

        // rowBlockId is the map-number the current row belongs to
        TEUCHOS_TEST_FOR_EXCEPTION(rowBlockId == -1, Exceptions::RuntimeError, "MueLu::ContactAFilterFactory::Build: cannot determine in which submap the current row is (either master, slave or innermap). Error.");

        size_t nnz = Ain->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");


        // just copy all values in output
        Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());
        size_t nNonzeros = 0;

        for(size_t i=0; i<(size_t)indices.size(); i++) {

            GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)

            Teuchos::ArrayRCP< const Scalar > colBlockData = blockVectorColMap->getData(0);

            LocalOrdinal colBlockId = colBlockData[indices[i]]; // LID -> colBlockID

            bool bIsMasterDof = false;
            bool bIsSlaveDof  = false;
            bool bIsInnerDof  = false;
            if(colBlockId == 0) bIsMasterDof = true;
            if(colBlockId == 1) bIsSlaveDof  = true;
            if(colBlockId == 2) bIsInnerDof  = true;

            if (! ( (bIsMasterDof == true && bIsSlaveDof == false && bIsInnerDof == false) ||
                    (bIsMasterDof == false && bIsSlaveDof == true  && bIsInnerDof == false) ||
                    (bIsMasterDof == false && bIsSlaveDof == false && bIsInnerDof == true) ))
              TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "MueLu::ContactAFilterFactory::Build: Dof listed in more than one map? Cannot be. Error.");

            if(rowBlockId == 0 && (bIsMasterDof || bIsInnerDof)) {  // master block row (MM and MI blocks)
                indout [nNonzeros] = gcid;
                valout [nNonzeros] = vals[i];
                nNonzeros++;
            } else if (rowBlockId == 1 && (bIsSlaveDof || bIsInnerDof)) { // slave block row (SS and SI blocks)
                indout [nNonzeros] = gcid;
                valout [nNonzeros] = vals[i];
                nNonzeros++;
            } else if (rowBlockId == 2) {
                indout [nNonzeros] = gcid;
                valout [nNonzeros] = vals[i];
                nNonzeros++;
            }

        }

        indout.resize(nNonzeros);
        valout.resize(nNonzeros);

        Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
      }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    GetOStream(Statistics0, 0) << "Nonzeros in " << varName_ << "(input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering " << varName_ << " (parameter: " << threshold_ << "): " << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set(varName_, Teuchos::rcp_dynamic_cast<OOperator>(Aout), this);

#endif
  }

  /*template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::IsGlobalId(Teuchos::RCP<const Map> & map, GlobalOrdinal gid) const {
    std::vector<GlobalOrdinal> global_id(1);
    global_id[0] = gid;
    const Teuchos::ArrayView<const LocalOrdinal> global_id_view(&global_id[0],global_id.size());

    std::vector<int> node_id(global_id.size(),-1);
    const Teuchos::ArrayView<int> node_id_view(&node_id[0],node_id.size());

    map->getRemoteIndexList(global_id_view, node_id_view);
    if (node_id[0] != -1) {
        return true;
    }
    return false;
  }*/

} // namespace MueLu
#endif


#endif // HAVE_MueLu

#endif /* MUELU_CONTACTAFILTERFACTORY_DEF_HPP_ */
