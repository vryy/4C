/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu smoother interface
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_MYTRILINOSSMOOTHER_DEF_HPP_
#define MUELU_MYTRILINOSSMOOTHER_DEF_HPP_

#ifdef HAVE_MueLuContact

#include "MueLu_MyTrilinosSmoother_decl.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyTrilinosSmoother(
      std::string const &mapName, const Teuchos::RCP<const FactoryBase> &mapFact,
      std::string const &type, Teuchos::ParameterList const &paramList, LO const &overlap,
      Teuchos::RCP<const FactoryBase> AFact)
      : mapName_(mapName),
        mapFact_(mapFact),
        type_(type),
        paramList_(paramList),
        overlap_(overlap),
        AFact_(AFact)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(overlap_ < 0, Exceptions::RuntimeError, "overlap_ < 0");

    if (type_ == "ILU")
    {
      s_ = MueLu::GetIfpackSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
          type_, paramList_, overlap_);
      s_->SetFactory("A", AFact_);
    }
    else
    {
      s_ = Teuchos::rcp(new TrilinosSmoother(type_, paramList_, overlap_));
      s_->SetFactory("A", AFact_);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError,
        "MueLu::MyTrilinosSmoother: failed to create TrilinosSmoother.");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(
      Level &currentLevel) const
  {
    currentLevel.DeclareInput(mapName_, mapFact_.get());
    s_->DeclareInput(currentLevel);  // call explicitly DeclareInput
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel)
  {
    // FactoryMonitor m(*this, "Setup Smoother");
    Monitor m(*this, "Setup MyTrilinosSmoother");

    if (SmootherPrototype::IsSetup() == true)
      VerboseObject::GetOStream(Warnings0, 0)
          << "Warning: MueLu::MyTrilinosSmoother::Setup(): Setup() has already been called";
    // TEUCHOS_TEST_FOR_EXCEPTION(s_ != Teuchos::null, Exceptions::RuntimeError, "IsSetup() == false
    // but s_ != Teuchos::null. This does not make sense");

    if (currentLevel.IsAvailable(mapName_, mapFact_.get()))
    {
      map_ = currentLevel.Get<Teuchos::RCP<const Map>>(mapName_, mapFact_.get());
      TEUCHOS_TEST_FOR_EXCEPTION(map_ == Teuchos::null, Exceptions::RuntimeError,
          "MueLu::MyTrilinosSmoother::Setup: map is Teuchos::null.");
    }
    else
      map_ = Teuchos::null;  // no map with artificial Dirichlet boundaries available.

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(
      MultiVector &X, MultiVector const &B, bool InitialGuessIsZero) const
  {
    // void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X,
    // MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError,
        "MueLu::AmesosSmoother::Apply(): Setup() has not been called");
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError,
        "IsSetup() == true but s_ == Teuchos::null. This does not make sense");

    // create non-const copy of rhs vector B
    Teuchos::RCP<MultiVector> Btemp = MultiVectorFactory::Build(B.getMap(), 1, true);
    Btemp->update(1.0, B, 0.0);

    if (map_ != Teuchos::null)
    {
      size_t sz = map_->getNodeNumElements();
      for (size_t it = 0; it < sz; it++)
      {
        GlobalOrdinal gid = map_->getGlobalElement(Teuchos::as<LocalOrdinal>(it));
        if (X.getMap()->isNodeGlobalElement(gid))
        {
          LocalOrdinal xlid = X.getMap()->getLocalElement(gid);         //
          Teuchos::ArrayRCP<const Scalar> data = X.getDataNonConst(0);  // extract data
          Btemp->replaceGlobalValue(gid, 0, data[xlid]);
          X.replaceGlobalValue(gid, 0, data[xlid]);  // not necessary
        }
        else
        {
          // gid is not stored on current processor for (overlapping) map of X
          // this should not be the case
          std::cout << "GID not stored in variable X on current processor" << std::endl;
        }
      }
    }  // end if map_ != Teuchos::null
    s_->Apply(X, *Btemp, InitialGuessIsZero);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const
  {
    return Teuchos::rcp(new MyTrilinosSmoother(*this));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const
  {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(
      Teuchos::FancyOStream &out, const VerbLevel verbLevel) const
  {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
    {
      out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1)
    {
      out0 << "PrecType: " << type_ << std::endl;
      out0 << "Parameter list: " << std::endl;
      {
        Teuchos::OSTab tab2(out);
        out << paramList_;
      }
      out0 << "Overlap: " << overlap_ << std::endl;
    }

    if (verbLevel & Debug)
    {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

}  // namespace MueLu

#endif  // HAVE_MueLuContact

#endif /* MUELU_MYTRILINOSSMOOTHER_DEF_HPP_ */
