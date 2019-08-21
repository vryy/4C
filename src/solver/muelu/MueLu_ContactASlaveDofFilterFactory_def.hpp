/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu contact filter factory class
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTASLAVEDOFFILTERFACTORY_DEF_HPP_
#define MUELU_CONTACTASLAVEDOFFILTERFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "MueLu_ContactASlaveDofFilterFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal,
      Node>::ContactASlaveDofFilterFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal,
      Node>::~ContactASlaveDofFilterFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal,
      Node>::GetValidParameterList() const
  {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null,
        "Generating factory of the matrix A used for filtering slave-master coupling DOFs.");

    return validParamList;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(
      Level &currentLevel) const
  {
    RCP<const FactoryBase> AFact = GetFactory("A");
    currentLevel.DeclareInput("A", AFact.get(), this);
    currentLevel.DeclareInput("SlaveDofMap", MueLu::NoFactory::get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ContactASlaveDofFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
      Level &currentLevel) const
  {
    typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> OOperator;
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsOOperator;
    // typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactoryClass;
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;

    Monitor m(*this, "ContactASlaveDofFilter factory");

    // extract slave dof map from Level class
    // the slave dof map contains all global DOF ids that shall be handled as Dirichlet boundaries
    // i.e. zero out the corresponding lines in matrix A and put a 1.0 on the diagonal
    Teuchos::RCP<const MapClass> slavedofmap =
        currentLevel.Get<Teuchos::RCP<const MapClass>>("SlaveDofMap", MueLu::NoFactory::get());

    // extract the original matrix A
    RCP<const FactoryBase> AFact = GetFactory("A");
    Teuchos::RCP<OOperator> Ain = currentLevel.Get<Teuchos::RCP<OOperator>>("A", AFact.get());

    // create new empty Operator
    Teuchos::RCP<CrsOOperator> Aout = Teuchos::rcp(new CrsOOperator(
        Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries(), Xpetra::StaticProfile));  // FIXME

    // loop over local rows
    for (size_t row = 0; row < Ain->getNodeNumRows(); row++)
    {
      // get global row id
      GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row);  // global row id

      if (slavedofmap->isNodeGlobalElement(grid) == false)
      {
        size_t nnz = Ain->getNumEntriesInLocalRow(row);
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz,
            Exceptions::RuntimeError,
            "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of "
            "indices? Error.");

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

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<OOperator>(Aout), this);
  }

}  // namespace MueLu

#endif  // HAVE_MueLu

#endif /* MUELU_CONTACTASLAVEDOFFILTERFACTORY_DEF_HPP_ */
