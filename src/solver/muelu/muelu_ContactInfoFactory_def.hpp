/*
 * muelu_ContactInfoFactory_def.hpp
 *
 *  Created on: Feb 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTINFOFACTORY_DEF_HPP_
#define MUELU_CONTACTINFOFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "muelu_ContactInfoFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ContactInfoFactory(std::string filename_prototype, Teuchos::RCP<FactoryBase> AFact, Teuchos::RCP<FactoryBase> nspFact)
  : filename_prototype_(filename_prototype), AFact_(AFact), nspFact_(nspFact)
    {

    }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ContactInfoFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A", AFact_.get(), this);
    fineLevel.DeclareInput("Nullspace", nspFact_.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorClass; //TODO
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO
    typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactoryClass;
    typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorClass;
    typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorClass;
    typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactoryClass;
    typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactoryClass;

    Monitor m(*this, "Contact info factory");

    RCP<OperatorClass> A = fineLevel.Get<RCP<OperatorClass> >("A",AFact_.get());
    RCP<MultiVectorClass> nsp = fineLevel.Get<RCP<MultiVectorClass> >("Nullspace",nspFact_.get());

    std::string outFile = filename_prototype_;
    std::stringstream streamLevel; streamLevel << fineLevel.GetLevelID();
    outFile = replaceAll(outFile,"%LEVELID", streamLevel.str());
    std::stringstream streamProc; streamProc << A->getRowMap()->getComm()->getRank();
    outFile = replaceAll(outFile,"%PROCID", streamProc.str());

    
    size_t nsdim = nsp->getNumVectors();
    int    nPDE  = 3; // number of equations (usually 2 or 3)
    
    // matrix-vector multiplication
    RCP<MultiVectorClass> Ansp = MultiVectorFactoryClass::Build(A->getRowMap(),nsdim);
    
    A->apply(*nsp,*Ansp);


    
    std::ofstream os;
    os.open(outFile.c_str(),std::fstream::trunc);

    
    if(nPDE == 3) {
     // 3 dim problem (either fluid or structure), depending on number of null space vectors (nsdim)
     
     if(nsdim != 3 && nsdim != 6) std::cout << "for 3d problems we expect the null space to be 3 or 6. It is " << nsdim << std::endl;

     if(fineLevel.GetLevelID() == 0) {
       // loop over null space vectors
       for(int nsv = 0; nsv<Teuchos::as<int>(nsdim); nsv++) {
	    os << "VECTORS nsp" << nsv << " float" << std::endl;
	    Teuchos::ArrayRCP<Scalar> nullspacevec = nsp->getDataNonConst(nsv);
	    // loop over all nodes
	    int length = nsp->getLocalLength();
	    for (int row = 0; row < Teuchos::as<int>(length/3); row++) {
		for (int col = 0; col < 3; col++) { // print components of vector (at max 3)
		    os << std::setprecision(16);
		    os << nullspacevec[3*row+col];
		    os << " ";
		}
		os << std::endl;
	    }
	
       }

       
     }
     else {
      // loop over null space vectors
      for(int nsv = 0; nsv<Teuchos::as<int>(nsdim); nsv++) {
	// split information in vectors a 3 dofs (for visualization in vtk files)
	for (int k = 0; k<(int)(Teuchos::as<int>(nsdim)/3); k++) {
	    os << "VECTORS nsp" << nsv << "_" << k << " float" << std::endl;
	    Teuchos::ArrayRCP<Scalar> nullspacevec = nsp->getDataNonConst(nsv);
	    // loop over all nodes
	    int length = nsp->getLocalLength();
	    for (int row = 0; row < Teuchos::as<int>(length/nsdim); row++) {
		for (int col = 0; col < 3; col++) { // print components of vector (at max 3)
		    os << std::setprecision(16);
		    os << nullspacevec[nsdim*row+3*k+col];
		    os << " ";
		}
		os << std::endl;
	    }
	}
      }
     }
    }
    else if(nPDE == 2) {
	std::cout << ("2d problems not supported yet.") << std::endl;
    }



  
 
    if(nPDE == 3) {
     // 3 dim problem (either fluid or structure), depending on number of null space vectors (nsdim)
     
     if(nsdim != 3 && nsdim != 6) std::cout << "for 3d problems we expect the null space to be 3 or 6. It is " << nsdim << std::endl;
          if(fineLevel.GetLevelID() == 0) {
       // loop over null space vectors
       for(int nsv = 0; nsv<Teuchos::as<int>(nsdim); nsv++) {
	    os << "VECTORS Ansp" << nsv << " float" << std::endl;
	    Teuchos::ArrayRCP<Scalar> nullspacevec = Ansp->getDataNonConst(nsv);
	    // loop over all nodes
	    int length = nsp->getLocalLength();
	    for (int row = 0; row < Teuchos::as<int>(length/3); row++) {
		for (int col = 0; col < 3; col++) { // print components of vector (at max 3)
		    os << std::setprecision(16);
		    os << nullspacevec[3*row+col];
		    os << " ";
		}
		os << std::endl;
	    }
       } 
     }
     else {
     // loop over null space vectors
     for(int nsv = 0; nsv<Teuchos::as<int>(nsdim); nsv++) {
       // split information in vectors a 3 dofs (for visualization in vtk files)
       for (int k = 0; k<Teuchos::as<int>(nsdim/3); k++) {
	  os << "VECTORS Ansp" << nsv << "_" << k << " float" << std::endl;
	  Teuchos::ArrayRCP<Scalar> Anullspacevec = Ansp->getDataNonConst(nsv);
	  // loop over all nodes
	  int length = nsp->getLocalLength();
	  for (int row = 0; row < Teuchos::as<int>(length/nsdim); row++) {
	      for (int col = 0; col < 3; col++) { // print components of vector (at max 3)
		  os << std::setprecision(16);
		  os << Anullspacevec[nsdim*row+3*k+col];
		  os << " ";
	      }
	      os << std::endl;
	  }
       }
     }
     }
    }
    else if(nPDE == 2) {
	std::cout << ("2d problems not supported yet.") << std::endl;
    }

    os << flush;
    os.close();

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string ContactInfoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const {
    while(1)
    {
      const int pos = result.find(replaceWhat);
      if (pos==-1) break;
      result.replace(pos,replaceWhat.size(),replaceWithWhat);
    }
    return result;
  }

} // namespace MueLu

#endif // HAVE_MueLu


#endif /* MUELU_CONTACTINFOFACTORY_DEF_HPP_ */
