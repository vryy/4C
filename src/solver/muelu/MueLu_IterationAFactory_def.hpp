/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu iteration factory class
\level 2

*----------------------------------------------------------------------*/

#ifndef MUELU_ITERATIONAFACTORY_DEF_HPP_
#define MUELU_ITERATIONAFACTORY_DEF_HPP_

#ifdef TRILINOS_Q1_2015

#include "MueLu_IterationAFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IterationAFactory(
      /*const std::string mapName,const Teuchos::RCP<const FactoryBase> & mapFact*/)
  // : mapName_(mapName), mapFact_(mapFact)
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~IterationAFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList>
  IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const
  {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->setEntry("map: name", Teuchos::ParameterEntry(std::string("")));
    validParamList->setEntry("map: factory", Teuchos::ParameterEntry(std::string("null")));
    validParamList->set<RCP<const FactoryBase>>(
        std::string("A"), Teuchos::null, "Input matrix A the filter is applied to");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(
      Level &currentLevel) const
  {
    Input(currentLevel, "A");

    const ParameterList &pL = GetParameterList();
    std::string mapFactName = pL.get<std::string>("map: factory");
    std::string mapName = pL.get<std::string>("map: name");

    // check whether user has provided a specific name for the MapFactory
    Teuchos::RCP<const FactoryBase> mapFact = Teuchos::null;
    if (mapFactName == "NoFactory")
    {
      mapFact = MueLu::NoFactory::getRCP();
    }
    else if (mapFactName != "null")
    {
      mapFact = currentLevel.GetFactoryManager()->GetFactory(mapFactName);
    }

    currentLevel.DeclareInput(mapName, mapFact.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
      Level &currentLevel) const
  {
    Monitor m(*this, "IterationAFactory");

    const ParameterList &pL = GetParameterList();
    std::string mapFactName = pL.get<std::string>("map: factory");
    std::string mapName = pL.get<std::string>("map: name");

    // check whether user has provided a specific name for the MapFactory
    Teuchos::RCP<const FactoryBase> mapFact = Teuchos::null;
    if (mapFactName == "NoFactory")
    {
      mapFact = MueLu::NoFactory::getRCP();
    }
    else if (mapFactName != "null")
    {
      mapFact = currentLevel.GetFactoryManager()->GetFactory(mapFactName);
    }

    // extract dof map from Level class
    // the dof map contains all global DOF ids that shall be handled as Dirichlet boundaries
    // i.e. zero out the corresponding lines in matrix A and put a 1.0 on the diagonal
    Teuchos::RCP<const Map> fixeddofmap =
        currentLevel.Get<Teuchos::RCP<const Map>>(mapName, mapFact.get());

    // extract the original matrix A
    Teuchos::RCP<Matrix> Ain =
        Get<Teuchos::RCP<Matrix>>(currentLevel, "A");  // corresponding Get call to Input

    // create new empty Operator
    Teuchos::RCP<CrsMatrixWrap> Aout = Teuchos::rcp(new CrsMatrixWrap(
        Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries(), Xpetra::StaticProfile));

    // loop over local rows
    for (size_t row = 0; row < Ain->getNodeNumRows(); row++)
    {
      // get global row id
      GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row);  // global row id

      if (fixeddofmap->isNodeGlobalElement(grid) == false)
      {
        // extract information from current row
        // size_t nnz = Ain->getNumEntriesInLocalRow(row);
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        // TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz,
        // Exceptions::RuntimeError, "MueLu::IterationAFactory::Build: number of nonzeros not equal
        // to number of indices? Error.");

        // just copy all values in output
        Teuchos::ArrayRCP<GlobalOrdinal> indout(
            indices.size(), Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(), Teuchos::ScalarTraits<Scalar>::zero());

        for (size_t i = 0; i < (size_t)indices.size(); i++)
        {
          GlobalOrdinal gcid =
              Ain->getColMap()->getGlobalElement(indices[i]);  // LID -> GID (column)
          indout[i] = gcid;
          valout[i] = vals[i];
        }

        Aout->insertGlobalValues(
            grid, indout.view(0, indout.size()), valout.view(0, valout.size()));
      }
      else
      {
        Aout->insertGlobalValues(
            grid, Teuchos::tuple<GlobalOrdinal>(grid), Teuchos::tuple<Scalar>(1.0));
      }
    }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    // copy block size information
    Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

    GetOStream(Statistics0, 0) << "Nonzeros in A  (input): " << Ain->getGlobalNumEntries()
                               << ", Nonzeros introducing artificial DC bdries "
                               << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(Aout), this);
  }

}  // namespace MueLu

#endif  // TRILINOS_Q1_2015

#endif /* MUELU_ITERATIONAFACTORY_DEF_HPP_ */
