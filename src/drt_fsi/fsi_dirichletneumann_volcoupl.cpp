/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
       with volume coupling

\level 3

\maintainer Anh-Tu Vuong

*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | headers                                                 farah 03/16 |
 *---------------------------------------------------------------------*/
#include "fsi_dirichletneumann_volcoupl.H"
#include "fsi_debugwriter.H"

#include "../drt_inpar/inpar_fsi.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_adapter/ad_fld_fluid_xfem.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_fld_fluid_xfsi.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_ale_fluid.H"

#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include "../drt_mortar/mortar_calc_utils.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

// search
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannVolCoupl::DirichletNeumannVolCoupl(const Epetra_Comm& comm)
    : DirichletNeumannDisp(comm), coupsa_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::Setup()
{
  FSI::DirichletNeumann::Setup();

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  SetKinematicCoupling(
      DRT::INPUT::IntegralValue<int>(fsipart, "COUPVARIABLE") == INPAR::FSI::CoupVarPart::disp);

  if (!GetKinematicCoupling()) dserror("Currently only displacement coupling is supported!");

  SetupCouplingStructAle(fsidyn, Comm());

  SetupInterfaceCorrector(fsidyn, Comm());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::SetupCouplingStructAle(
    const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  const int ndim = DRT::Problem::Instance()->NDim();

  coupsa_ = Teuchos::rcp(new ADAPTER::MortarVolCoupl());

  // do a dynamic cast here
  Teuchos::RCP<ADAPTER::FluidAle> fluidale = Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(fluid_);

  // projection
  std::vector<int> coupleddof12 = std::vector<int>(ndim, 1);
  std::vector<int> coupleddof21 = std::vector<int>(ndim, 1);

  // define dof sets to be coupled for both projections
  std::pair<int, int> dofsets12(0, 0);
  std::pair<int, int> dofsets21(0, 0);

  // initialize coupling adapter
  coupsa_->Init(StructureField()->Discretization(),
      fluidale->AleField()->WriteAccessDiscretization(), &coupleddof12, &coupleddof21, &dofsets12,
      &dofsets21, Teuchos::null, false);

  // setup coupling adapter
  coupsa_->Setup();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::SetupInterfaceCorrector(
    const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  icorrector_ = Teuchos::rcp(new InterfaceCorrector());

  icorrector_->Setup(Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVolCoupl::FluidOp(
    Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag)
{
  FSI::Partitioned::FluidOp(idisp, fillFlag);

  // TODO cant this be done better?
  Teuchos::RCP<Epetra_Vector> vdisp = Teuchos::rcp(new Epetra_Vector(*StructureField()->Dispnp()));

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return FluidToStruct(MBFluidField()->RelaxationSolve(StructToFluid(idisp), Dt()));
  }
  else
  {
    // normal fluid solve
    // the displacement -> velocity conversion at the interface
    const Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idisp);

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) MBFluidField()->SetItemax(mfresitemax_ + 1);

    Teuchos::RCP<ADAPTER::FluidAle> fluidale =
        Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(MBFluidField());

    icorrector_->SetInterfaceDisplacements(idisp, StructureFluidCoupling());

    // important difference to dirichletneumann.cpp: vdisp is mapped from structure to ale here
    fluidale->NonlinearSolveVolCoupl(StructureToAle(vdisp), StructToFluid(ivel), icorrector_);

    MBFluidField()->SetItemax(itemax);

    return FluidToStruct(MBFluidField()->ExtractInterfaceForces());
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::ExtractPreviousInterfaceSolution()
{
  iveln_ = FluidToStruct(MBFluidField()->ExtractInterfaceVeln());
  idispn_ = StructureField()->ExtractInterfaceDispn();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVolCoupl::StructureToAle(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVolCoupl::AleToStructure(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::InterfaceCorrector::Setup(Teuchos::RCP<ADAPTER::FluidAle> fluidale)
{
  fluidale_ = fluidale;

  volcorrector_ = Teuchos::rcp(new VolCorrector);
  volcorrector_->Setup(DRT::Problem::Instance()->NDim(), fluidale);

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::InterfaceCorrector::SetInterfaceDisplacements(
    Teuchos::RCP<Epetra_Vector>& idisp_struct, ADAPTER::Coupling& icoupfs)
{
  idisp_ = idisp_struct;
  icoupfs_ = Teuchos::rcpFromRef(icoupfs);

  deltadisp_ = Teuchos::null;
  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::InterfaceCorrector::CorrectInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  if (icoupfs_ == Teuchos::null) dserror("Coupling adapter not set!");
  if (idisp_ == Teuchos::null) dserror("Interface displacements not set!");

  // std::cout<<*finterface->FullMap()<<std::endl;
  // std::cout<<*disp_fluid<<std::endl;
  deltadisp_ = LINALG::CreateVector(*finterface->FSICondMap(), true);

  LINALG::Export(*disp_fluid, *deltadisp_);
  // deltadisp_ = finterface->ExtractFSICondVector(disp_fluid);

  // dserror("stop");

  Teuchos::RCP<Epetra_Vector> idisp_fluid_corrected = icoupfs_->MasterToSlave(idisp_);

  deltadisp_->Update(1.0, *idisp_fluid_corrected, -1.0);

  LINALG::Export(*idisp_fluid_corrected, *disp_fluid);
  // finterface->InsertFSICondVector(idisp_fluid_corrected,disp_fluid);

  volcorrector_->CorrectVolDisplacements(fluidale_, deltadisp_, disp_fluid, finterface);

  // reset
  idisp_ = Teuchos::null;
  icoupfs_ = Teuchos::null;

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::CorrectVolDisplacements(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
    Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  if (fluidale->AleField()->Discretization()->Comm().MyPID() == 0)
    std::cout << "******************   FSI Volume Correction Step   **********************"
              << std::endl;

  // correction step in parameter space
  if (true) CorrectVolDisplacementsParaSpace(fluidale, deltadisp, disp_fluid, finterface);
  // correction step in physical space
  else
    CorrectVolDisplacementsPhysSpace(fluidale, deltadisp, disp_fluid, finterface);

  // output
  if (fluidale->AleField()->Discretization()->Comm().MyPID() == 0)
    std::cout << "******************FSI Volume Correction Step Done***********************"
              << std::endl;


  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::CorrectVolDisplacementsParaSpace(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
    Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  Teuchos::RCP<Epetra_Vector> correction = Teuchos::rcp(new Epetra_Vector(disp_fluid->Map(), true));
  Teuchos::RCP<Epetra_Vector> DofColMapDummy =
      Teuchos::rcp(new Epetra_Vector(*fluidale->FluidField()->Discretization()->DofColMap(), true));
  LINALG::Export(*deltadisp, *DofColMapDummy);

  const double tol = 1e-5;

  // loop over ale eles
  for (std::map<int, std::vector<int>>::iterator it = fluidalenodemap_.begin();
       it != fluidalenodemap_.end(); ++it)
  {
    DRT::Element* aleele = fluidale->AleField()->Discretization()->gElement(it->first);

    // loop over fluid volume nodes within one ale FSI element
    for (size_t i = 0; i < it->second.size(); ++i)
    {
      int gid = it->second[i];
      DRT::Node* fluidnode = fluidale->FluidField()->Discretization()->gNode(gid);

      if (fluidnode->Owner() != fluidale->AleField()->Discretization()->Comm().MyPID()) continue;

      double gpos[3] = {fluidnode->X()[0], fluidnode->X()[1], fluidnode->X()[2]};
      double lpos[3] = {0.0, 0.0, 0.0};
      if (aleele->Shape() == DRT::Element::quad4)
        MORTAR::UTILS::GlobalToLocal<DRT::Element::quad4>(*aleele, gpos, lpos);
      else if (aleele->Shape() == DRT::Element::hex8)
        MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*aleele, gpos, lpos);
      else
        dserror("ERROR: element type not implemented!");

      if (lpos[0] < -1.0 - tol || lpos[1] < -1.0 - tol || lpos[2] < -1.0 - tol ||
          lpos[0] > 1.0 + tol || lpos[1] > 1.0 + tol || lpos[2] > 1.0 + tol)
        continue;

      double dist = 1.0e12;
      int id = -1;
      for (size_t k = 0; k < fluidalenodeFSImap_[it->first].size(); ++k)
      {
        int gidfsi = fluidalenodeFSImap_[it->first][k];
        DRT::Node* fluidnodeFSI = fluidale->FluidField()->Discretization()->gNode(gidfsi);

        double gposFSI[3] = {fluidnodeFSI->X()[0], fluidnodeFSI->X()[1], fluidnodeFSI->X()[2]};
        double lposFSI[3] = {0.0, 0.0, 0.0};
        if (aleele->Shape() == DRT::Element::quad4)
          MORTAR::UTILS::GlobalToLocal<DRT::Element::quad4>(*aleele, gposFSI, lposFSI);
        else if (aleele->Shape() == DRT::Element::hex8)
          MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*aleele, gposFSI, lposFSI);
        else
          dserror("ERROR: element type not implemented!");

        if (lposFSI[0] < -1.0 - tol || lposFSI[1] < -1.0 - tol || lposFSI[2] < -1.0 - tol ||
            lposFSI[0] > 1.0 + tol || lposFSI[1] > 1.0 + tol || lposFSI[2] > 1.0 + tol)
          dserror("ERROR: wrong parameter space coordinates!");

        // valc distance to fsi node
        double vec0 = lposFSI[0] - lpos[0];
        double vec1 = lposFSI[1] - lpos[1];
        double vec2 = lposFSI[2] - lpos[2];
        double actdist = sqrt(vec0 * vec0 + vec1 * vec1 + vec2 * vec2);

        // check length
        if (actdist < dist)
        {
          id = fluidnodeFSI->Id();
          dist = actdist;
        }
      }  // end loop

      // safety
      if (id < 0) continue;

      double fac = 0.0;

      if (aleele->Shape() == DRT::Element::quad4 or aleele->Shape() == DRT::Element::hex8)
        fac = 1.0 - 0.5 * dist;
      else
        dserror("ERROR: element type not implemented!");

      // safety
      if (dist > 2.0) fac = 0.0;

      DRT::Node* fluidnodeFSI = fluidale->FluidField()->Discretization()->gNode(id);
      std::vector<int> temp = fluidale->FluidField()->Discretization()->Dof(fluidnodeFSI);
      std::vector<int> dofsFSI;
      for (int idof = 0; idof < dim_; idof++) dofsFSI.push_back(temp[idof]);

      // extract local values of the global vectors
      std::vector<double> FSIdisp(dofsFSI.size());
      DRT::UTILS::ExtractMyValues(*DofColMapDummy, FSIdisp, dofsFSI);

      std::vector<int> temp2 = fluidale->FluidField()->Discretization()->Dof(fluidnode);
      std::vector<int> dofs;
      for (int idof = 0; idof < dim_; idof++) dofs.push_back(temp2[idof]);

      Epetra_SerialDenseVector gnode(dim_);
      std::vector<int> lmowner(dim_);
      for (int idof = 0; idof < dim_; idof++)
      {
        gnode(idof) = fac * FSIdisp[idof];
        lmowner[idof] = fluidnode->Owner();
      }

      LINALG::Assemble(*correction, gnode, dofs, lmowner);
    }  // end fluid volume node loop
  }    // end ale fsi element loop

  // do correction
  disp_fluid->Update(1.0, *correction, 1.0);

  // calc norm
  double norm = 0.0;
  correction->Norm2(&norm);

  // output
  if (fluidale->AleField()->Discretization()->Comm().MyPID() == 0)
    std::cout << "Norm of correction (parameter space): " << norm << std::endl;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::CorrectVolDisplacementsPhysSpace(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
    Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  Teuchos::RCP<Epetra_Vector> correction = Teuchos::rcp(new Epetra_Vector(disp_fluid->Map(), true));
  Teuchos::RCP<Epetra_Vector> DofColMapDummy =
      Teuchos::rcp(new Epetra_Vector(*fluidale->FluidField()->Discretization()->DofColMap(), true));
  LINALG::Export(*deltadisp, *DofColMapDummy);

  std::map<int, LINALG::Matrix<9, 2>> CurrentDOPs =
      CalcBackgroundDops(fluidale->FluidField()->Discretization());

  Teuchos::RCP<std::set<int>> FSIaleeles =
      DRT::UTILS::ConditionedElementMap(*fluidale->AleField()->Discretization(), "FSICoupling");

  // evaluate search
  for (int i = 0; i < fluidale->AleField()->Discretization()->NumMyColElements(); ++i)
  {
    // 1 map node into bele
    int gid = fluidale->AleField()->Discretization()->ElementColMap()->GID(i);
    DRT::Element* aleele = fluidale->AleField()->Discretization()->gElement(gid);

    if (FSIaleeles->find(aleele->Id()) == FSIaleeles->end()) continue;
  }

  // do correction
  disp_fluid->Update(1.0, *correction, 1.0);

  // calc norm
  double norm = 0.0;
  correction->Norm2(&norm);

  // output
  if (fluidale->AleField()->Discretization()->Comm().MyPID() == 0)
    std::cout << "Norm of correction (physical space): " << norm << std::endl;

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::Setup(const int dim, Teuchos::RCP<ADAPTER::FluidAle> fluidale)
{
  if (fluidale->AleField()->Discretization()->Comm().MyPID() == 0)
    std::cout << "******************FSI Volume Correction Setup***********************"
              << std::endl;

  dim_ = dim;
  InitDopNormals();

  // init current positions
  std::map<int, LINALG::Matrix<3, 1>> currentpositions;

  for (int lid = 0; lid < fluidale->FluidField()->Discretization()->NumMyColElements(); ++lid)
  {
    DRT::Element* sele = fluidale->FluidField()->Discretization()->lColElement(lid);

    // calculate slabs for every node on every element
    for (int k = 0; k < sele->NumNode(); k++)
    {
      DRT::Node* node = sele->Nodes()[k];
      LINALG::Matrix<3, 1> currpos;

      currpos(0) = node->X()[0];
      currpos(1) = node->X()[1];
      currpos(2) = node->X()[2];

      currentpositions[node->Id()] = currpos;
    }
  }

  // init of 3D search tree
  searchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3, 2> rootBox =
      GEO::getXAABBofDis(*fluidale->FluidField()->Discretization(), currentpositions);
  searchTree_->initializeTree(
      rootBox, *fluidale->FluidField()->Discretization(), GEO::TreeType(GEO::OCTTREE));


  std::map<int, LINALG::Matrix<9, 2>> CurrentDOPs =
      CalcBackgroundDops(fluidale->FluidField()->Discretization());

  Teuchos::RCP<std::set<int>> FSIaleeles =
      DRT::UTILS::ConditionedElementMap(*fluidale->AleField()->Discretization(), "FSICoupling");

  // evaluate search
  for (int i = 0; i < fluidale->AleField()->Discretization()->NumMyColElements(); ++i)
  {
    // 1 map node into bele
    int gid = fluidale->AleField()->Discretization()->ElementColMap()->GID(i);
    DRT::Element* aleele = fluidale->AleField()->Discretization()->gElement(gid);

    if (FSIaleeles->find(aleele->Id()) == FSIaleeles->end()) continue;

    // get found elements from other discr.
    fluidaleelemap_[gid] = Search(*aleele, CurrentDOPs);
  }  // end node loop

  Teuchos::RCP<Epetra_Map> FSIfluidnodes =
      DRT::UTILS::ConditionNodeColMap(*fluidale->FluidField()->Discretization(), "FSICoupling");

  std::set<int> globalnodeids;
  // loop over ale eles
  for (std::map<int, std::vector<int>>::iterator it = fluidaleelemap_.begin();
       it != fluidaleelemap_.end(); ++it)
  {
    DRT::Element* aleele = fluidale->AleField()->Discretization()->gElement(it->first);

    std::vector<int> localnodeids;
    std::vector<int> localnodeidsFSI;

    // loop over fluid eles
    for (size_t i = 0; i < it->second.size(); ++i)
    {
      int gid = it->second[i];
      DRT::Element* fluidele = fluidale->FluidField()->Discretization()->gElement(gid);

      for (int j = 0; j < fluidele->NumNode(); ++j)
      {
        const int nodegid = fluidele->NodeIds()[j];
        DRT::Node* fluidnode = fluidele->Nodes()[j];

        double gpos[3] = {fluidnode->X()[0], fluidnode->X()[1], fluidnode->X()[2]};
        double lpos[3] = {0.0, 0.0, 0.0};
        if (aleele->Shape() == DRT::Element::quad4)
          MORTAR::UTILS::GlobalToLocal<DRT::Element::quad4>(*aleele, gpos, lpos);
        else if (aleele->Shape() == DRT::Element::hex8)
          MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*aleele, gpos, lpos);
        else
          dserror("ERROR: element type not implemented!");

        double tol = 1e-5;
        if (lpos[0] < -1.0 - tol || lpos[1] < -1.0 - tol || lpos[2] < -1.0 - tol ||
            lpos[0] > 1.0 + tol || lpos[1] > 1.0 + tol || lpos[2] > 1.0 + tol)
          continue;

        if (FSIfluidnodes->MyGID(nodegid))
        {
          localnodeidsFSI.push_back(nodegid);
          continue;
        }
        if (globalnodeids.find(nodegid) != globalnodeids.end()) continue;

        globalnodeids.insert(nodegid);
        localnodeids.push_back(nodegid);
      }
    }
    fluidalenodemap_[it->first] = localnodeids;
    fluidalenodeFSImap_[it->first] = localnodeidsFSI;
  }

  std::cout << "ALE elements found: " << fluidaleelemap_.size() << std::endl;

  if (fluidale->AleField()->Discretization()->Comm().MyPID() == 0)
    std::cout << "******************FSI Volume Correction Setup Done***********************"
              << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Init normals for Dop calculation                         farah 05/16|
 *----------------------------------------------------------------------*/
void FSI::VolCorrector::InitDopNormals()
{
  dopnormals_(0, 0) = 1.0;
  dopnormals_(0, 1) = 0.0;
  dopnormals_(0, 2) = 0.0;

  dopnormals_(1, 0) = 0.0;
  dopnormals_(1, 1) = 1.0;
  dopnormals_(1, 2) = 0.0;

  dopnormals_(2, 0) = 0.0;
  dopnormals_(2, 1) = 0.0;
  dopnormals_(2, 2) = 1.0;

  dopnormals_(3, 0) = 1.0;
  dopnormals_(3, 1) = 1.0;
  dopnormals_(3, 2) = 0.0;

  dopnormals_(4, 0) = 1.0;
  dopnormals_(4, 1) = 0.0;
  dopnormals_(4, 2) = 1.0;

  dopnormals_(5, 0) = 0.0;
  dopnormals_(5, 1) = 1.0;
  dopnormals_(5, 2) = 1.0;

  dopnormals_(6, 0) = 1.0;
  dopnormals_(6, 1) = 0.0;
  dopnormals_(6, 2) = -1.0;

  dopnormals_(7, 0) = 1.0;
  dopnormals_(7, 1) = -1.0;
  dopnormals_(7, 2) = 0.0;

  dopnormals_(8, 0) = 0.0;
  dopnormals_(8, 1) = 1.0;
  dopnormals_(8, 2) = -1.0;

  return;
}


/*----------------------------------------------------------------------*
 |  Calculate Dops for background mesh                       farah 05/16|
 *----------------------------------------------------------------------*/
std::map<int, LINALG::Matrix<9, 2>> FSI::VolCorrector::CalcBackgroundDops(
    Teuchos::RCP<DRT::Discretization> searchdis)
{
  std::map<int, LINALG::Matrix<9, 2>> currentKDOPs;

  for (int lid = 0; lid < searchdis->NumMyColElements(); ++lid)
  {
    DRT::Element* sele = searchdis->lColElement(lid);

    currentKDOPs[sele->Id()] = CalcDop(*sele);
  }

  return currentKDOPs;
}

/*----------------------------------------------------------------------*
 |  Calculate Dop for one Element                            farah 05/16|
 *----------------------------------------------------------------------*/
LINALG::Matrix<9, 2> FSI::VolCorrector::CalcDop(DRT::Element& ele)
{
  LINALG::Matrix<9, 2> dop;

  // calculate slabs
  for (int j = 0; j < 9; j++)
  {
    // initialize slabs
    dop(j, 0) = 1.0e12;
    dop(j, 1) = -1.0e12;
  }

  // calculate slabs for every node on every element
  for (int k = 0; k < ele.NumNode(); k++)
  {
    DRT::Node* node = ele.Nodes()[k];

    // get current node position
    double pos[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < dim_; ++j) pos[j] = node->X()[j];

    // calculate slabs
    for (int j = 0; j < 9; j++)
    {
      //= ax+by+cz=d/sqrt(aa+bb+cc)
      double num =
          dopnormals_(j, 0) * pos[0] + dopnormals_(j, 1) * pos[1] + dopnormals_(j, 2) * pos[2];
      double denom =
          sqrt((dopnormals_(j, 0) * dopnormals_(j, 0)) + (dopnormals_(j, 1) * dopnormals_(j, 1)) +
               (dopnormals_(j, 2) * dopnormals_(j, 2)));
      double dcurrent = num / denom;

      if (dcurrent > dop(j, 1)) dop(j, 1) = dcurrent;
      if (dcurrent < dop(j, 0)) dop(j, 0) = dcurrent;
    }
  }

  return dop;
}


/*----------------------------------------------------------------------*
 |  Perform searching procedure                              farah 05/16|
 *----------------------------------------------------------------------*/
std::vector<int> FSI::VolCorrector::Search(
    DRT::Element& ele, std::map<int, LINALG::Matrix<9, 2>>& currentKDOPs)
{
  // vector of global ids of found elements
  std::vector<int> gids;
  gids.clear();
  std::set<int> gid;
  gid.clear();

  LINALG::Matrix<9, 2> queryKDOP;

  // calc dop for considered element
  queryKDOP = CalcDop(ele);

  //**********************************************************
  // search for near elements to the background node's coord
  searchTree_->searchCollisions(currentKDOPs, queryKDOP, 0, gid);

  for (std::set<int>::iterator iter = gid.begin(); iter != gid.end(); ++iter) gids.push_back(*iter);

  return gids;
}
