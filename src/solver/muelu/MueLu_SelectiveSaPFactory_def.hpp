/*
 * MueLu_SelectiveSaPFactory_def.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: wiesner
 */

#ifndef MUELU_SELECTIVESAPFACTORY_DEF_HPP_
#define MUELU_SELECTIVESAPFACTORY_DEF_HPP_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2013

#include <Epetra_RowMatrixTransposer.h>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_SelectiveSaPFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< Scalar >                         ("Damping factor",          4./3, "Smoothed-Aggregation damping factor");
    validParamList->set< const std::string >              ("Damping strategy", "Standard", "Damping strategy. Can be either \'Standard\' or \'User\'. (default: \'Standard\')");
    validParamList->set< Teuchos::RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    validParamList->set< Teuchos::RCP<const FactoryBase> >("P",              Teuchos::null, "Tentative prolongator factory");
    validParamList->set< const std::string >              ("NonSmoothRowMapName","", "Name of row map of Dofs which corresponding transfer operator basis functions are not to be smoothed.");
    validParamList->set< Teuchos::RCP<const FactoryBase> >("NonSmoothRowMapFactory",Teuchos::null, "Factory for row map of Dofs which corresponding transfer operator basis functions are not to be smoothed.");
    validParamList->set< const std::string >              ("NearZeroDiagMapName","", "Name of row map with nearly zero entries on diagonal of A");
    validParamList->set< Teuchos::RCP<const FactoryBase> >("NearZeroDiagMapFactory",Teuchos::null, "Factory for row map of Dofs generating \"NearZeroDiagMapName\"");

    //validParamList->set< bool >("bPerformSmoothing",true, "boolean flag to turn on/off smoothing of transfer operator basis functions. If set to false, the unsmoothed tentative transfer operator is used. This is to switch between smoothed and unsmoothed transfer operators on different multigrid levels and therefore it's more flexible than PA-AMG.");
    // validParamList->set                       ("Diagonal view",      "current", "Diagonal view used during the prolongator smoothing process");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P", initialPFact.get(), this); // --
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }

    // Level Get

    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    //Build final prolongator
    RCP<Matrix> finalP; // output

    // decide whether transfer operator basis functions shall be smoothed on current level
    bool bPerformSmoothing = true;
    const ParameterList & pL = GetParameterList();
    const std::string NearZeroDiagRowMapName = pL.get<std::string >("NearZeroDiagMapName");
    Teuchos::RCP<const FactoryBase> NearZeroDiagRowMapFact = GetFactory("NearZeroDiagMapFactory");
    if(NearZeroDiagRowMapName != "") {
      if(fineLevel.IsAvailable(NearZeroDiagRowMapName, NearZeroDiagRowMapFact.get())) {
        Teuchos::RCP<const Map> NearZeroDiagRowMap = fineLevel.Get< Teuchos::RCP<const Map> >(NearZeroDiagRowMapName, NearZeroDiagRowMapFact.get());
        // switch off transfer operator smoothing completely on current level
        // if there are some significant nearly zero entries on diagonal of matrix A
        // determined by solver object (e.g. PermutingAztecSolver)
        if(NearZeroDiagRowMap->getGlobalNumElements() != 0){
          bPerformSmoothing = false;
          std::cout << "               switch off transfer operator smoothing" << std::endl;
        }
      }
    }

    Scalar dampingFactor = pL.get<Scalar>("Damping factor");

    // perform basis function smoothing (only if requested/necessary)
    if ( bPerformSmoothing && dampingFactor != Teuchos::ScalarTraits<Scalar>::zero() ) {

      // extract matrix A from level (that shall be used for transfer operator smoothing)
      RCP<Matrix> A     = Get< RCP<Matrix> >(fineLevel, "A");
      if(restrictionMode_) {
        SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
        //A = Utils2::Transpose(A, true); // build transpose of A explicitely
        A = MyTranspose(A, true);
      }

      RCP<Matrix> AP;
      {
        SubFactoryMonitor m2(*this, "MxM: A x Ptentative", coarseLevel);
        //JJH -- If I switch doFillComplete to false, the resulting matrix seems weird when printed with describe.
        //JJH -- The final prolongator is wrong, to boot.  So right now, I fillComplete AP, but avoid fillComplete
        //JJH -- in the scaling.  Long story short, we're doing 2 fillCompletes, where ideally we'd do just one.
        bool doFillComplete=true;
        //FIXME Improved Epetra MM returns error code -1 optimizeStorage==true.  For now, do this
        bool optimizeStorage=true;  // false
        //FIXME but once fixed, reenable the next line.
        if (A->getRowMap()->lib() == Xpetra::UseTpetra) optimizeStorage=false;
#ifdef HAVE_Trilinos_Q3_2013
        RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        AP = Utils::Multiply(*A, false, *Ptent, false,*fos, doFillComplete, optimizeStorage);
#else
        AP = Utils::Multiply(*A, false, *Ptent, false, doFillComplete, optimizeStorage);
#endif
        //Utils::Write("AP.mat", *AP);
      }

      int bSkipDamping = 0;
      int globalbSkipDamping = 0;
      {
        SubFactoryMonitor m2(*this, "Scaling (A x Ptentative) by D^{-1}", coarseLevel);
        bool doFillComplete=true; //false;
        bool optimizeStorage=false;
        Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(*A);
#if 1 // TODO remove me. check is done outside (e.g. in MueLu_Contact preconditioner
        // needed if we do not provide the near zero map!!!
        if(NearZeroDiagRowMapName == "") {
          for(LocalOrdinal l=0; l<diag.size(); l++) {
            if(std::abs(diag[l])<1e-8) {  // was <1e-3 before
              std::cout << "               switch off transfer operator smoothing" << std::endl;
              bSkipDamping = 1;
              break;
            }
          }
        }

        // sum up all entries in multipleColRequests over all processors
        sumAll(A->getRowMap()->getComm(), (LocalOrdinal)bSkipDamping, globalbSkipDamping);
#endif

#ifdef HAVE_Trilinos_Q3_2013
        Utils::MyOldScaleMatrix(*AP, diag, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag
#else
        Utils::MyOldScaleMatrix(AP, diag, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag
#endif
        //Utils::Write("DinvAP.mat", *AP);
      }

      if (globalbSkipDamping == 0) {
        const std::string DampingStrategy = pL.get<std::string >("Damping strategy");
        Scalar lambdaMax;
        {
          if(DampingStrategy == "Standard") {
            SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
            Magnitude stopTol = 1e-4;
            lambdaMax = Utils::PowerMethod(*A, true, (LO) 10, stopTol);
            //Scalar lambdaMax = Utils::PowerMethod(*A, true, (LO) 50, (Scalar)1e-7, true);
            GetOStream(Statistics1, 0) << "Damping factor = " << dampingFactor/lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;
            if(lambdaMax < 0.0) {
              lambdaMax = 1.0;
              dampingFactor = 0.0; // no smoothing at all
              GetOStream(Statistics1, 0) << "Correct damping factor = " << dampingFactor/lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;
            }
          } else if(DampingStrategy == "User") {
            lambdaMax = 1.0;
            GetOStream(Statistics1, 0) << "Damping factor = " << dampingFactor/lambdaMax << " (user-given)" << std::endl;
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::SelectiveSaFactory::Build: wrong damping strategy. allowed values are Standard and User");
          }
        }


        {
          // fix product DinvAP
          // remove basis function smoothing for problematic rows
          Teuchos::RCP<Matrix> APtemp = FixAPproduct(fineLevel, coarseLevel, A, AP);
          if(APtemp != Teuchos::null) // TODO check this.
            AP = APtemp;
        }

        {
          SubFactoryMonitor m2(*this, "M+M: P = (Ptentative) + (D^{-1} x A x Ptentative)", coarseLevel);

          bool doTranspose=false;
          bool PtentHasFixedNnzPerRow=true;
#ifdef HAVE_Trilinos_Q3_2013
          RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
          Utils2::TwoMatrixAdd(*Ptent, doTranspose, Teuchos::ScalarTraits<Scalar>::one(), *AP, doTranspose, -dampingFactor/lambdaMax, finalP,*fos, PtentHasFixedNnzPerRow);
#else
          Utils2::TwoMatrixAdd(Ptent, doTranspose, Teuchos::ScalarTraits<Scalar>::one(), AP, doTranspose, -dampingFactor/lambdaMax, finalP, PtentHasFixedNnzPerRow);
#endif

        }

        {
          SubFactoryMonitor m2(*this, "FillComplete() of P", coarseLevel);
          finalP->fillComplete( Ptent->getDomainMap(), Ptent->getRangeMap() );
        }
      } else // bDoDamping == false // TODO fix this. move this to FixAP?
      {
        finalP = Ptent;
      }
    } else {
      finalP = Ptent;
    }

    // Level Set
    if(!restrictionMode_)
      {
        // prolongation factory is in prolongation mode
        Set(coarseLevel, "P", finalP);

        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) finalP->CreateView("stridedMaps", Ptent);
        ///////////////////////// EXPERIMENTAL
      }
    else
      {
        // prolongation factory is in restriction mode
#ifdef HAVE_Trilinos_Q3_2013
        RCP<Matrix> R = Utils2::Transpose(*finalP, true); // use Utils2 -> specialization for double
#else
        RCP<Matrix> R = Utils2::Transpose(finalP, true); // use Utils2 -> specialization for double
#endif
        //RCP<Matrix> R = MyTranspose(finalP, true);
        Set(coarseLevel, "R", R);

        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) R->CreateView("stridedMaps", Ptent, true);
        ///////////////////////// EXPERIMENTAL
      }

  } //Build()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::FixAPproduct(Level &fineLevel, Level &coarseLevel, Teuchos::RCP<Matrix> & A, Teuchos::RCP<Matrix> & DinvAP) const {
    // 1) extract NonSmoothRowMap from Level
    // TODO: move this code to Build function
    const ParameterList & pL = GetParameterList();

    const std::string NonSmoothedRowMapName = pL.get<std::string >("NonSmoothRowMapName");
    if(NonSmoothedRowMapName == "")
      return Teuchos::null; // the full transfer operator shall be smoothed

    Teuchos::RCP<const FactoryBase> NonSmoothedRowMapFact = GetFactory("NonSmoothRowMapFactory");

    Teuchos::RCP<const Map> NonSmoothedRowMap = fineLevel.Get< Teuchos::RCP<const Map> >(NonSmoothedRowMapName, NonSmoothedRowMapFact.get());

    // 2) find all columns in Ptent which correspond to rows in NonSmoothedRowMap

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }

    // Level Get
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    Teuchos::RCP<Vector> gPtentColVec = VectorFactory::Build(Ptent->getColMap());
    Teuchos::RCP<Vector> gPtentDomVec = VectorFactory::Build(Ptent->getDomainMap());
    gPtentColVec->putScalar(0.0);
    gPtentDomVec->putScalar(0.0);

    std::vector<GlobalOrdinal > coarseMapGids;

    // loop over all rows of Ptent
    for(size_t row=0; row<Ptent->getNodeNumRows(); row++) {
      // get global row id
      GlobalOrdinal grid = Ptent->getRowMap()->getGlobalElement(row); // global row id

      if(NonSmoothedRowMap->isNodeGlobalElement(grid)) {
        // extract data from current row (grid)
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ptent->getLocalRowView(row, indices, vals);

        for(size_t i=0; i<(size_t)indices.size(); i++) {
          if(vals[i] != 0.0) {  // TODO this is not necessary
            GlobalOrdinal gcid = Ptent->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
            gPtentColVec->sumIntoGlobalValue(gcid,1.0);
          }
        }
      }
    }

    // check which are necessary and which not
    Teuchos::RCP<Export> exporter = ExportFactory::Build(gPtentColVec->getMap(), gPtentDomVec->getMap());
    gPtentDomVec->doExport(*gPtentColVec,*exporter,Xpetra::ABSMAX);  // communicate blocked gcolids to all procs
    //gPtentColVec->doImport(*gPtentDomVec,*exporter,Xpetra::INSERT);

    Teuchos::RCP<Vector> gDinvAPColVec = VectorFactory::Build(DinvAP->getColMap());
    Teuchos::RCP<Import> importer = ImportFactory::Build(gPtentDomVec->getMap(), gDinvAPColVec->getMap());
    gDinvAPColVec->doImport(*gPtentDomVec,*importer,Xpetra::INSERT);

#if 0 // not sure about this...
    // 3) check diagonal of A for problematic entries
    Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(*A);
    for(LocalOrdinal row=0; row<diag.size(); row++) {
      if(std::abs(diag[row])<1e-3) {
        //std::cout << "fixed diagonal in row " << l << " from " << diag[l];
        //diag[l] = 0.001;
        //std::cout << " to be " << diag[l] << std::endl;

        // get global row id
        //GlobalOrdinal grid = DinvAP->getRowMap()->getGlobalElement(row); // global row id

        // extract data from current row (grid)
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        DinvAP->getLocalRowView(row, indices, vals);

        for(size_t i=0; i<(size_t)indices.size(); i++) {
          GlobalOrdinal gcid = DinvAP->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
          gDinvAPColVec->sumIntoGlobalValue(gcid,1.0);
        }
      }
    }

    Teuchos::RCP<Vector> gDinvAPDomVec = VectorFactory::Build(DinvAP->getDomainMap());
    gDinvAPDomVec->putScalar(0.0);
    Teuchos::RCP<Export> exporter2 = ExportFactory::Build(gDinvAPColVec->getMap(), gDinvAPDomVec->getMap());
    gDinvAPDomVec->doExport(*gDinvAPColVec,*exporter2,Xpetra::ABSMAX);  // communicate blocked gcolids to all procs
    gDinvAPColVec->doImport(*gDinvAPDomVec,*exporter2,Xpetra::INSERT);
#endif

    // 4) remove columns from DinvAP
#if 1

    // replace columns

    Teuchos::ArrayRCP<const Scalar > gDinvAPColVecData = gDinvAPColVec->getData(0);

    // arrays to extract data from AP matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;

    // empty arrays to store column indices to override
    size_t maxNNz = DinvAP->getNodeMaxNumRowEntries();
    Teuchos::Array<GlobalOrdinal> indout(maxNNz,0);
    Teuchos::Array<Scalar> valout(maxNNz,0.0);
    size_t curNNz = 0; // number nnz to be replaced in current row

    // loop over local rows
    for(size_t row=0; row<DinvAP->getNodeNumRows(); row++) {
        // extract data from current row (grid)
        DinvAP->getLocalRowView(row, indices, vals);

        // reset number of indices to override to zero
        curNNz = 0;

        for(size_t i=0; i<(size_t)indices.size(); i++) {
           if(gDinvAPColVecData[indices[i]] > 0.0) {
             // skip this column, explicitely set the value to zero
             indout[curNNz] = indices[i];
             valout[curNNz] = 0.0;
             curNNz++;
           }
         }
        DinvAP->replaceLocalValues(row, indout.view(0,curNNz), valout.view(0,curNNz));
    }

    return DinvAP;
#else
     // create new empty Operator
     RCP<CrsMatrixWrap> Aout = Teuchos::rcp(new CrsMatrixWrap(DinvAP->getRowMap(),DinvAP->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile));

     //GetOStream(Statistics1, 0) << "avoid transfer operator smoothing for " << coarseNonSmoothedMap->getGlobalNumElements() << "/" << DinvAP->getDomainMap()->getGlobalNumElements() << " transfer operator basis functions" << std::endl;

     Teuchos::ArrayRCP<const Scalar > gDinvAPColVecData = gDinvAPColVec->getData(0);

     // loop over local rows
     for(size_t row=0; row<DinvAP->getNodeNumRows(); row++) {
         // get global row id
         GlobalOrdinal grid = DinvAP->getRowMap()->getGlobalElement(row); // global row id

         // extract data from current row (grid)
         Teuchos::ArrayView<const LocalOrdinal> indices;
         Teuchos::ArrayView<const Scalar> vals;
         DinvAP->getLocalRowView(row, indices, vals);

         // just copy all column entries whose global id is not in cleardofgids
         Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
         Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

         for(size_t i=0; i<(size_t)indices.size(); i++) {
           GlobalOrdinal gcid = DinvAP->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)

           if(gDinvAPColVecData[indices[i]] > 0.0) {
             // skip this column, explicitely set the value to zero // TODO we can just skip this completely
             indout [i] = gcid;
             valout [i] = 0.0;
           }
           else {
             // keep the column entry
             indout [i] = gcid;
             valout [i] = vals[i];
           }
         }

         // fill Aout wiht new row/column entries
         Aout->insertGlobalValues(grid, indout.view(0,indout.size()), valout.view(0,valout.size()));

     }

     Aout->fillComplete(DinvAP->getDomainMap(), DinvAP->getRangeMap());

     ///////////////////////// EXPERIMENTAL
     // copy block size information
     if(DinvAP->IsView("stridedMaps")) Aout->CreateView("stridedMaps", DinvAP);
     ///////////////////////// EXPERIMENTAL

     return Aout;
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<Xpetra::Matrix<double, int, int> > SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyTranspose(Teuchos::RCP<Xpetra::Matrix<double, int, int> > const &Op, bool const & optimizeTranspose) const
  {
   Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("My Transpose"));

    Teuchos::RCP<Epetra_CrsMatrix> epetraOp;
    epetraOp = MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(Op);

    Epetra_RowMatrixTransposer et(&*epetraOp);
    Epetra_CrsMatrix *A;
    int rv = et.CreateTranspose(false,A);
    if(rv != 0) {
      std::ostringstream buf;
      buf << rv;
      std::string msg = "MueLu::Utils::Transpose: Epetra::RowMatrixTransposer returned value of " + buf.str();
      std::cout << msg << std::endl;
    }

    RCP<Epetra_CrsMatrix> rcpA(A);
    RCP<EpetraCrsMatrix> AA = rcp(new EpetraCrsMatrix(rcpA) );
    RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
    RCP<Xpetra::CrsMatrixWrap<SC> > AAAA = rcp( new Xpetra::CrsMatrixWrap<SC>(AAA) );
    AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
    return AAAA;

  } //Transpose

#if 0 // outdated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildUnsmoothedBasisFunctionMap(Level &fineLevel, Level &coarseLevel) const {

    const ParameterList & pL = GetParameterList();

    const std::string NonSmoothedRowMapName = pL.get<std::string >("NonSmoothRowMapName");
    if(NonSmoothedRowMapName == "")
      return Teuchos::null; // the full transfer operator shall be smoothed

    Teuchos::RCP<const FactoryBase> NonSmoothedRowMapFact = GetFactory("NonSmoothRowMapFactory");

    Teuchos::RCP<const Map> NonSmoothedRowMap = fineLevel.Get< Teuchos::RCP<const Map> >(NonSmoothedRowMapName, NonSmoothedRowMapFact.get());

    //Teuchos::RCP<const Map> NonSmoothedRowMap = pL.get<Teuchos::RCP<const Map> >("NonSmoothRowMap");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }

    // Level Get
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    std::vector<GlobalOrdinal > coarseMapGids;

    // loop over all rows of Ptent
    for(size_t row=0; row<Ptent->getNodeNumRows(); row++) {
        // get global row id
        GlobalOrdinal grid = Ptent->getRowMap()->getGlobalElement(row); // global row id

        if(NonSmoothedRowMap->isNodeGlobalElement(grid)) {
          // extract data from current row (grid)
          //size_t nnz = Ain->getNumEntriesInLocalRow(row);
          Teuchos::ArrayView<const LocalOrdinal> indices;
          Teuchos::ArrayView<const Scalar> vals;
          Ptent->getLocalRowView(row, indices, vals);

          for(size_t i=0; i<(size_t)indices.size(); i++) {
            if(vals[i] != 0.0) {  // TODO this is not necessary
              GlobalOrdinal gcid = Ptent->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
              coarseMapGids.push_back(gcid);
            }
          }
        }
    }

    Teuchos::ArrayView<GlobalOrdinal> coarseMapGidsView (&coarseMapGids[0],coarseMapGids.size());

    Teuchos::RCP<const Map> coarseNonSmoothedMap = MapFactory::Build(NonSmoothedRowMap->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        coarseMapGidsView,
        Ptent->getColMap()->getIndexBase(),
        Ptent->getColMap()->getComm());

    return coarseNonSmoothedMap;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RemoveColumnEntries(Teuchos::RCP<Matrix> & matrix, Teuchos::RCP<const Map> coarseNonSmoothedMap) const {
    // extract the original matrix A
    RCP<Matrix> Ain = matrix;

    // create new empty Operator
    RCP<CrsMatrixWrap> Aout = Teuchos::rcp(new CrsMatrixWrap(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile));

    GetOStream(Statistics1, 0) << "avoid transfer operator smoothing for " << coarseNonSmoothedMap->getGlobalNumElements() << "/" << Ain->getDomainMap()->getGlobalNumElements() << " transfer operator basis functions" << std::endl;

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++) {
        // get global row id
        GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id

        // extract data from current row (grid)
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        // just copy all column entries whose global id is not in cleardofgids
        Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

        Teuchos::ArrayRCP<GlobalOrdinal>  GIDList(1);
        Teuchos::ArrayRCP<int> nodeIDList(GIDList.size());

        for(size_t i=0; i<(size_t)indices.size(); i++) {
          GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)

          if(coarseNonSmoothedMap->isNodeGlobalElement(gcid)) {
            // skip this column, explicitely set the value to zero // TODO we can just skip this completely
            indout [i] = gcid;
            valout [i] = 0.0;
          }
          else {
            // in cleardofmap gcid is not stored on this proc
            // check whether gcid is in cleardofmap but on another proc
            /*GIDList[0] = gcid;
            Xpetra::LookupStatus checkOtherNodes = coarseNonSmoothedMap->getRemoteIndexList(GIDList(), nodeIDList());
            if(checkOtherNodes == Xpetra::AllIDsPresent) {
              // gcid is in cleardofmap, but on another processor
              // skip it, explicitely set the value to zero
              indout [i] = gcid;
              valout [i] = 0.0;
            } else {*/
              // keep the column entry
              indout [i] = gcid;
              valout [i] = vals[i];
            //}
          }
        }

        // fill Aout wiht new row/column entries
        Aout->insertGlobalValues(grid, indout.view(0,indout.size()), valout.view(0,valout.size()));

    }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    ///////////////////////// EXPERIMENTAL
    // copy block size information
    if(Ain->IsView("stridedMaps")) Aout->CreateView("stridedMaps", Ain);
    ///////////////////////// EXPERIMENTAL

    return Aout;
  }
#endif

  // deprecated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDampingFactor(Scalar dampingFactor) {
    SetParameter("Damping factor", ParameterEntry(dampingFactor)); // revalidate
  }

  // deprecated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDampingFactor() {
    const ParameterList & pL = GetParameterList();
    return pL.get<Scalar>("Damping factor");
  }

} //namespace MueLu

#endif //#ifdef HAVE_Trilinos_Q1_2013
#endif // HAVE_MueLu

#endif /* MUELU_SELECTIVESAPFACTORY_DEF_HPP_ */
