/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu contact info factory class
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTINFOFACTORY_DEF_HPP_
#define MUELU_CONTACTINFOFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "MueLu_ContactInfoFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ContactInfoFactory(
      std::string filename_prototype, Teuchos::RCP<FactoryBase> AFact,
      Teuchos::RCP<FactoryBase> nspFact)
      : filename_prototype_(filename_prototype), AFact_(AFact), nspFact_(nspFact)
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~ContactInfoFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(
      Level& fineLevel, Level& coarseLevel) const
  {
    fineLevel.DeclareInput("A", AFact_.get(), this);
    fineLevel.DeclareInput("Nullspace", nspFact_.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
      Level& fineLevel, Level& coarseLevel) const
  {
    typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> OperatorClass;  // TODO
    typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorClass;
    typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
        MultiVectorFactoryClass;

    Monitor m(*this, "Contact info factory");

    Teuchos::RCP<OperatorClass> A = fineLevel.Get<Teuchos::RCP<OperatorClass>>("A", AFact_.get());
    Teuchos::RCP<MultiVectorClass> nsp =
        fineLevel.Get<Teuchos::RCP<MultiVectorClass>>("Nullspace", nspFact_.get());

    std::string outFile = filename_prototype_;
    std::stringstream streamLevel;
    streamLevel << fineLevel.GetLevelID();
    outFile = replaceAll(outFile, "%LEVELID", streamLevel.str());
    std::stringstream streamProc;
    streamProc << A->getRowMap()->getComm()->getRank();
    outFile = replaceAll(outFile, "%PROCID", streamProc.str());


    size_t nsdim = nsp->getNumVectors();
    int nPDE = 3;  // number of equations (usually 2 or 3)

    // matrix-vector multiplication
    Teuchos::RCP<MultiVectorClass> Ansp = MultiVectorFactoryClass::Build(A->getRowMap(), nsdim);

    A->apply(*nsp, *Ansp);



    std::ofstream os;
    os.open(outFile.c_str(), std::fstream::trunc);


    if (nPDE == 3)
    {
      // 3 dim problem (either fluid or structure), depending on number of null space vectors
      // (nsdim)

      if (nsdim != 3 && nsdim != 6)
        std::cout << "for 3d problems we expect the null space to be 3 or 6. It is " << nsdim
                  << std::endl;

      if (fineLevel.GetLevelID() == 0)
      {
        // loop over null space vectors
        for (int nsv = 0; nsv < Teuchos::as<int>(nsdim); nsv++)
        {
          os << "VECTORS nsp" << nsv << " float" << std::endl;
          Teuchos::ArrayRCP<Scalar> nullspacevec = nsp->getDataNonConst(nsv);
          // loop over all nodes
          int length = nsp->getLocalLength();
          for (int row = 0; row < Teuchos::as<int>(length / 3); row++)
          {
            for (int col = 0; col < 3; col++)
            {  // print components of vector (at max 3)
              os << std::setprecision(16);
              os << nullspacevec[3 * row + col];
              os << " ";
            }
            os << std::endl;
          }
        }
      }
      else
      {
        // loop over null space vectors
        for (int nsv = 0; nsv < Teuchos::as<int>(nsdim); nsv++)
        {
          // split information in vectors a 3 dofs (for visualization in vtk files)
          for (int k = 0; k < (int)(Teuchos::as<int>(nsdim) / 3); k++)
          {
            os << "VECTORS nsp" << nsv << "_" << k << " float" << std::endl;
            Teuchos::ArrayRCP<Scalar> nullspacevec = nsp->getDataNonConst(nsv);
            // loop over all nodes
            int length = nsp->getLocalLength();
            for (int row = 0; row < Teuchos::as<int>(length / nsdim); row++)
            {
              for (int col = 0; col < 3; col++)
              {  // print components of vector (at max 3)
                os << std::setprecision(16);
                os << nullspacevec[nsdim * row + 3 * k + col];
                os << " ";
              }
              os << std::endl;
            }
          }
        }
      }
    }
    else if (nPDE == 2)
    {
      std::cout << ("2d problems not supported yet.") << std::endl;
    }



    if (nPDE == 3)
    {
      // 3 dim problem (either fluid or structure), depending on number of null space vectors
      // (nsdim)

      if (nsdim != 3 && nsdim != 6)
        std::cout << "for 3d problems we expect the null space to be 3 or 6. It is " << nsdim
                  << std::endl;
      if (fineLevel.GetLevelID() == 0)
      {
        // loop over null space vectors
        for (int nsv = 0; nsv < Teuchos::as<int>(nsdim); nsv++)
        {
          os << "VECTORS Ansp" << nsv << " float" << std::endl;
          Teuchos::ArrayRCP<Scalar> nullspacevec = Ansp->getDataNonConst(nsv);
          // loop over all nodes
          int length = nsp->getLocalLength();
          for (int row = 0; row < Teuchos::as<int>(length / 3); row++)
          {
            for (int col = 0; col < 3; col++)
            {  // print components of vector (at max 3)
              os << std::setprecision(16);
              os << nullspacevec[3 * row + col];
              os << " ";
            }
            os << std::endl;
          }
        }
      }
      else
      {
        // loop over null space vectors
        for (int nsv = 0; nsv < Teuchos::as<int>(nsdim); nsv++)
        {
          // split information in vectors a 3 dofs (for visualization in vtk files)
          for (int k = 0; k < Teuchos::as<int>(nsdim / 3); k++)
          {
            os << "VECTORS Ansp" << nsv << "_" << k << " float" << std::endl;
            Teuchos::ArrayRCP<Scalar> Anullspacevec = Ansp->getDataNonConst(nsv);
            // loop over all nodes
            int length = nsp->getLocalLength();
            for (int row = 0; row < Teuchos::as<int>(length / nsdim); row++)
            {
              for (int col = 0; col < 3; col++)
              {  // print components of vector (at max 3)
                os << std::setprecision(16);
                os << Anullspacevec[nsdim * row + 3 * k + col];
                os << " ";
              }
              os << std::endl;
            }
          }
        }
      }
    }
    else if (nPDE == 2)
    {
      std::cout << ("2d problems not supported yet.") << std::endl;
    }

    os << std::flush;
    os.close();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceAll(
      std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const
  {
    while (1)
    {
      const int pos = result.find(replaceWhat);
      if (pos == -1) break;
      result.replace(pos, replaceWhat.size(), replaceWithWhat);
    }
    return result;
  }

}  // namespace MueLu

#endif  // HAVE_MueLu


#endif /* MUELU_CONTACTINFOFACTORY_DEF_HPP_ */
