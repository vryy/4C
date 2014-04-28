// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_MyTrilinosSmoother_def.hpp
 *
 *  Created on: Aug 20, 2012
 *      Author: wiesner
 */

#ifndef MUELU_MYTRILINOSSMOOTHER_DEF_HPP_
#define MUELU_MYTRILINOSSMOOTHER_DEF_HPP_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q3_2013

#include "MueLu_MyTrilinosSmoother_decl.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyTrilinosSmoother(std::string const & mapName, const Teuchos::RCP<const FactoryBase> & mapFact, std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap, Teuchos::RCP<FactoryBase> AFact)
    : mapName_(mapName), mapFact_(mapFact), type_(type), paramList_(paramList), overlap_(overlap), AFact_(AFact)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(overlap_ < 0, Exceptions::RuntimeError, "overlap_ < 0");

    if(type_ == "ILU") {
      s_ = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(type_, paramList_,overlap_);
      s_->SetFactory("A", AFact_);
    } else {
      s_ = Teuchos::rcp(new TrilinosSmoother(type_, paramList_, overlap_));
      s_->SetFactory("A", AFact_);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::MyTrilinosSmoother: failed to create TrilinosSmoother.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput(mapName_, mapFact_.get());
    s_->DeclareInput(currentLevel); // call explicitely DeclareInput
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //FactoryMonitor m(*this, "Setup Smoother");
    Monitor m(*this, "Setup MyTrilinosSmoother");

    if (SmootherPrototype::IsSetup() == true) VerboseObject::GetOStream(Warnings0, 0) << "Warning: MueLu::MyTrilinosSmoother::Setup(): Setup() has already been called";
    //TEUCHOS_TEST_FOR_EXCEPTION(s_ != Teuchos::null, Exceptions::RuntimeError, "IsSetup() == false but s_ != Teuchos::null. This does not make sense");

    if(currentLevel.IsAvailable(mapName_,mapFact_.get())) {
      map_ = currentLevel.Get<Teuchos::RCP<const Map> >(mapName_,mapFact_.get());
      TEUCHOS_TEST_FOR_EXCEPTION(map_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::MyTrilinosSmoother::Setup: map is Teuchos::null.");
    }
    else map_ = Teuchos::null; // no map with artificial Dirichlet boundaries available.

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
#ifdef HAVE_Trilinos_Q3_2013
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool InitialGuessIsZero) const {
#else
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
#endif
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "IsSetup() == true but s_ == Teuchos::null. This does not make sense");

    // create non-const copy of rhs vector B
    Teuchos::RCP<MultiVector> Btemp = MultiVectorFactory::Build(B.getMap(),1,true);
    Btemp->update(1.0,B,0.0);

    if (map_ != Teuchos::null) {
      size_t sz = map_->getNodeNumElements();
      for(size_t it = 0; it < sz; it++) {
        GlobalOrdinal gid = map_->getGlobalElement(Teuchos::as<LocalOrdinal>(it));
        if(X.getMap()->isNodeGlobalElement(gid)) {
          LocalOrdinal xlid = X.getMap()->getLocalElement(gid);        //
          Teuchos::ArrayRCP<const Scalar> data = X.getDataNonConst(0); // extract data
          Btemp->replaceGlobalValue(gid,0,data[xlid]);
          X.replaceGlobalValue(gid,0,data[xlid]); // not necessary
        }
        else {
          // gid is not stored on current processor for (overlapping) map of X
          // this should not be the case
          std::cout << "GID not stored in variable X on current processor" << std::endl;
        }
      }
    } // end if map_ != Teuchos::null
    s_->Apply(X, *Btemp, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return Teuchos::rcp( new MyTrilinosSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {

    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1) {
      out0 << "PrecType: " << type_ << std::endl;
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      out0 << "Overlap: " << overlap_ << std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#endif // #ifdef HAVE_Trilinos_Q1_2013
#endif // HAVE_MueLu

#endif /* MUELU_MYTRILINOSSMOOTHER_DEF_HPP_ */
