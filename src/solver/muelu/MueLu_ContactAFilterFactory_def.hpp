/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu factory class for BACI
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/
#ifndef MUELU_CONTACTAFILTERFACTORY_DEF_HPP_
#define MUELU_CONTACTAFILTERFACTORY_DEF_HPP_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

#include "MueLu_ContactAFilterFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ContactAFilterFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~ContactAFilterFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::ParameterList>
  ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList(
      const Teuchos::ParameterList& paramList) const
  {
    Teuchos::RCP<Teuchos::ParameterList> validParamList =
        Teuchos::rcp(new Teuchos::ParameterList());

    validParamList->set<std::string>(
        "Input matrix name", "A", "Name of input matrix. (default='A')");
    validParamList->set<Teuchos::RCP<const FactoryBase>>(
        "Input matrix factory", Teuchos::null, "Generating factory of the input matrix.");

    validParamList->set<std::string>(
        "Map block 1 name", "SlaveDofMap", "Name of part 1 of map to be splitted.");
    validParamList->set<std::string>(
        "Map block 2 name", "MasterDofMap", "Name of part 2 of map to be splitted.");
    validParamList->set<std::string>(
        "Map block 1 factory", "null", "Name of generating factory for 'Map block 1 name'");
    validParamList->set<std::string>(
        "Map block 2 factory", "null", "Name of generating factory for 'Map block 2 name'");

    // This does not work
    // validParamList->set< Teuchos::RCP<const FactoryBase> >("Map block 1 factory",
    // MueLu::NoFactory::getRCP(), "Generating factory of part 1 of map to be segregated.");
    // validParamList->set< Teuchos::RCP<const FactoryBase> >("Map block 2 factory",
    // MueLu::NoFactory::getRCP(), "Generating factory of part 2 of map to be segregated.");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(
      Level& currentLevel) const
  {
    const ParameterList& pL = GetParameterList();
    std::string inputName = pL.get<std::string>("Input matrix name");
    Teuchos::RCP<const FactoryBase> inputFactory = GetFactory("Input matrix factory");

    currentLevel.DeclareInput(inputName, inputFactory.get(), this);

    std::string blockName1 = pL.get<std::string>("Map block 1 name");
    std::string blockFactName1 = pL.get<std::string>("Map block 1 factory");
    if (blockFactName1 == "NoFactory")
    {
      currentLevel.DeclareInput(blockName1, MueLu::NoFactory::getRCP().get(), this);
    }
    else if (blockFactName1 != "null")
    {
      Teuchos::RCP<const FactoryBase> blockFact1 =
          currentLevel.GetFactoryManager()->GetFactory(blockFactName1);
      currentLevel.DeclareInput(blockName1, blockFact1.get(), this);
    }
    else
    {
      currentLevel.DeclareInput(blockName1, NULL, this);
    }

    // Teuchos::RCP<const FactoryBase> blockFactory1= GetFactory          ("Map block 1 factory");
    // currentLevel.DeclareInput(blockName1,blockFactory1.get(),this);

    std::string blockName2 = pL.get<std::string>("Map block 2 name");
    std::string blockFactName2 = pL.get<std::string>("Map block 2 factory");
    if (blockFactName2 == "NoFactory")
    {
      currentLevel.DeclareInput(blockName2, MueLu::NoFactory::getRCP().get(), this);
    }
    else if (blockFactName2 != "null")
    {
      Teuchos::RCP<const FactoryBase> blockFact2 =
          currentLevel.GetFactoryManager()->GetFactory(blockFactName2);
      currentLevel.DeclareInput(blockName2, blockFact2.get(), this);
    }
    else
    {
      currentLevel.DeclareInput(blockName2, NULL, this);
    }

    // Teuchos::RCP<const FactoryBase> blockFactory2= GetFactory          ("Map block 2 factory");
    // currentLevel.DeclareInput(blockName2,blockFactory2.get(),this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
      Level& currentLevel) const
  {
    Monitor m(*this, "A filtering (contact)");

    const Teuchos::ParameterList& pL = GetParameterList();
    std::string inputName = pL.get<std::string>("Input matrix name");
    std::string blockName1 = pL.get<std::string>("Map block 1 name");
    std::string blockName2 = pL.get<std::string>("Map block 2 name");
    Teuchos::RCP<const FactoryBase> inputFactory = GetFactory("Input matrix factory");
    // Teuchos::RCP<const FactoryBase> blockFactory1 = GetFactory ("Map block 1 factory");
    // Teuchos::RCP<const FactoryBase> blockFactory2 = GetFactory ("Map block 2 factory");

    Teuchos::RCP<const FactoryBase> blockFactory1 = Teuchos::null;
    Teuchos::RCP<const FactoryBase> blockFactory2 = Teuchos::null;
    std::string blockFactName1 = pL.get<std::string>("Map block 1 factory");
    if (blockFactName1 == "NoFactory")
    {
      blockFactory1 = MueLu::NoFactory::getRCP();
    }
    else if (blockFactName1 != "null")
    {
      blockFactory1 = currentLevel.GetFactoryManager()->GetFactory(blockFactName1);
    }
    std::string blockFactName2 = pL.get<std::string>("Map block 2 factory");
    if (blockFactName2 == "NoFactory")
    {
      blockFactory2 = MueLu::NoFactory::getRCP();
    }
    else if (blockFactName2 != "null")
    {
      blockFactory2 = currentLevel.GetFactoryManager()->GetFactory(blockFactName2);
    }

    // fetch map with slave Dofs from Level
    Teuchos::RCP<const Map> DofMap1 =
        currentLevel.Get<Teuchos::RCP<const Map>>(blockName1, blockFactory1.get());
    Teuchos::RCP<const Map> DofMap2 =
        currentLevel.Get<Teuchos::RCP<const Map>>(blockName2, blockFactory2.get());

    Teuchos::RCP<Matrix> Ain =
        currentLevel.Get<Teuchos::RCP<Matrix>>(inputName, inputFactory.get());

    Teuchos::RCP<Vector> blockVectorRowMap = VectorFactory::Build(Ain->getRowMap());
    blockVectorRowMap->putScalar(-1.0);  // -1.0 denotes that this Dof is not slave DOF

    // define (sub) block vectors
    Teuchos::RCP<Vector> blockVector1 = VectorFactory::Build(DofMap1);
    blockVector1->putScalar(1);
    Teuchos::RCP<Vector> blockVector2 = VectorFactory::Build(DofMap2);
    blockVector2->putScalar(2);

    Teuchos::RCP<const Import> importer1 = ImportFactory::Build(DofMap1, Ain->getColMap());
    Teuchos::RCP<const Import> importer2 = ImportFactory::Build(DofMap2, Ain->getColMap());
    Teuchos::RCP<Vector> blockVectorColMapData = VectorFactory::Build(Ain->getColMap());
    blockVectorColMapData->putScalar(-1.0);  // -1.0 denotes that this Dof is not slave DOF
    blockVectorColMapData->doImport(*blockVector1, *importer1, Xpetra::INSERT);
    blockVectorColMapData->doImport(*blockVector2, *importer2, Xpetra::INSERT);

    // create new empty Operator
    Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(
        Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries(), Xpetra::StaticProfile);

    // loop over local rows
    Teuchos::ArrayRCP<const Scalar> colBlockData = blockVectorColMapData->getData(0);
    size_t numLocalRows = Ain->getNodeNumRows();
    bool ba = false;
    bool bb = false;
    bool bc = false;
    bool bd = false;
    bool be = false;
    bool bf = false;
    bool bg = false;
    bool bh = false;
    bool bi = false;
    bool bj = false;

    // declare helper variables
    bool isBlock1, isBlock2;
    size_t nNonzeros = 0;
    LocalOrdinal colBlockId = -1;

    GetOStream(Statistics0, 0) << "Filtering: ";
    for (size_t row = 0; row < numLocalRows; row++)
    {
      // get global row id
      GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row);  // global row id

      isBlock1 = DofMap1->isNodeGlobalElement(grid);
      isBlock2 = DofMap2->isNodeGlobalElement(grid);

      // this can happen due to the stupid permutation strategy which mixes up slave and master dofs
      // or interface and inner dofs
      // TEUCHOS_TEST_FOR_EXCEPTION(isBlock1 && isBlock2 == true, Exceptions::RuntimeError,
      // "MueLu::ContactAFilterFactory::Build: row is in subblock 1 and subblock 2? Error.");

      // size_t nnz = Ain->getNumEntriesInLocalRow(row);
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      Ain->getLocalRowView(row, indices, vals);

      // TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz,
      // Exceptions::RuntimeError, "MueLu::ContactAFilterFactory::Build: number of nonzeros not
      // equal to number of indices? Error.");

      // just copy all values in output
      Teuchos::ArrayRCP<GlobalOrdinal> indout(
          indices.size(), Teuchos::ScalarTraits<GlobalOrdinal>::zero());
      Teuchos::ArrayRCP<Scalar> valout(indices.size(), Teuchos::ScalarTraits<Scalar>::zero());

      nNonzeros = 0;
      colBlockId = -1;
      for (size_t i = 0; i < (size_t)indices.size(); i++)
      {
        colBlockId = Teuchos::as<LocalOrdinal>(colBlockData[indices[i]]);  // LID -> colBlockID

        // colBlockId can be
        //  -1: indices[i] is neither in DofMap1 nor in DofMap2
        //   1: indices[i] is in DofMap1
        //   2: indices[i] is in DofMap2
        bool bCopy = false;
        if (isBlock1 == true && isBlock2 == true)
          isBlock1 = false;  // if a row is in both submaps put it to the master side here
        if (isBlock1 == false && isBlock2 == false)
          bCopy = true;  // row is neither in block 1 or block 2 -> copy
        if (isBlock1 == true && colBlockId == 1)
          bCopy = true;  // row is block 1 and column is block 1 -> copy
        if (isBlock1 == true && colBlockId == -1)
          bCopy = true;  // row is block 1 and column is block -1-> copy
        if (isBlock2 == true && colBlockId == 2)
          bCopy = true;  // row is block 2 and column is block 2 -> copy
        if (isBlock2 == true && colBlockId == -1)
          bCopy = true;  // row is block 2 and column is block -1-> copy

        if (bCopy)
        {
          GlobalOrdinal gcid =
              Ain->getColMap()->getGlobalElement(indices[i]);  // LID -> GID (column)
          indout[nNonzeros] = gcid;
          valout[nNonzeros] = vals[i];
          nNonzeros++;
        }
      }
      indout.resize(nNonzeros);
      valout.resize(nNonzeros);

      Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row),
          indout.view(0, indout.size()), valout.view(0, valout.size()));

      // this is somewhat expensive, but do communication for debug output
      double pPerCent = 0.0;
      // double pPerCent = Teuchos::as<Scalar>(row) / Teuchos::as<Scalar>(numLocalRows);

      if (pPerCent > 0.9 && ba == false)
      {
        GetOStream(Statistics0, 0) << "+";
        ba = true;
      }
      if (pPerCent > 0.8 && bb == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bb = true;
      }
      if (pPerCent > 0.7 && bc == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bc = true;
      }
      if (pPerCent > 0.6 && bd == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bd = true;
      }
      if (pPerCent > 0.5 && be == false)
      {
        GetOStream(Statistics0, 0) << "+";
        be = true;
      }
      if (pPerCent > 0.4 && bf == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bf = true;
      }
      if (pPerCent > 0.3 && bg == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bg = true;
      }
      if (pPerCent > 0.2 && bh == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bh = true;
      }
      if (pPerCent > 0.1 && bi == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bi = true;
      }
      if (pPerCent > 0.0 && bj == false)
      {
        GetOStream(Statistics0, 0) << "+";
        bj = true;
      }

      GetOStream(Statistics0, 0).getOStream()->flush();
    }

    GetOStream(Statistics0, 0) << " complete ..." << std::endl;
    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    // copy block size information
    Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

    GetOStream(Statistics0, 0) << "Nonzeros in " << inputName
                               << "(input): " << Ain->getGlobalNumEntries()
                               << ", Nonzeros after filtering " << inputName << ": "
                               << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set(inputName, Teuchos::rcp_dynamic_cast<Matrix>(Aout), this);
  }
}  // namespace MueLu

#endif  // HAVE_MueLuContact

#endif /* MUELU_CONTACTAFILTERFACTORY_DEF_HPP_ */
