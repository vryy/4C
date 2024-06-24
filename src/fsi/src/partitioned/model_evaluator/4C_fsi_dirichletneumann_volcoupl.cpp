/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
       with volume coupling

\level 3


*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
#include "4C_fsi_dirichletneumann_volcoupl.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_adapter_fld_fluid_xfem.hpp"
#include "4C_adapter_fld_fluid_xfsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_geometry_searchtree.hpp"
#include "4C_fem_geometry_searchtree_service.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_calc_utils.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannVolCoupl::DirichletNeumannVolCoupl(const Epetra_Comm& comm)
    : DirichletNeumannDisp(comm), coupsa_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::setup()
{
  FSI::DirichletNeumann::setup();

  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  set_kinematic_coupling(
      Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::disp);

  if (!get_kinematic_coupling()) FOUR_C_THROW("Currently only displacement coupling is supported!");

  setup_coupling_struct_ale(fsidyn, Comm());

  setup_interface_corrector(fsidyn, Comm());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::setup_coupling_struct_ale(
    const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  const int ndim = Global::Problem::Instance()->NDim();

  coupsa_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

  // do a dynamic cast here
  Teuchos::RCP<Adapter::FluidAle> fluidale = Teuchos::rcp_dynamic_cast<Adapter::FluidAle>(fluid_);

  // projection
  std::vector<int> coupleddof12 = std::vector<int>(ndim, 1);
  std::vector<int> coupleddof21 = std::vector<int>(ndim, 1);

  // define dof sets to be coupled for both projections
  std::pair<int, int> dofsets12(0, 0);
  std::pair<int, int> dofsets21(0, 0);

  // initialize coupling adapter
  coupsa_->init(ndim, structure_field()->discretization(),
      fluidale->ale_field()->write_access_discretization(), &coupleddof12, &coupleddof21,
      &dofsets12, &dofsets21, Teuchos::null, false);

  // setup coupling adapter
  coupsa_->setup(Global::Problem::Instance()->VolmortarParams());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::setup_interface_corrector(
    const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  icorrector_ = Teuchos::rcp(new InterfaceCorrector());

  icorrector_->setup(Teuchos::rcp_dynamic_cast<Adapter::FluidAle>(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVolCoupl::fluid_op(
    Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag)
{
  FSI::Partitioned::fluid_op(idisp, fillFlag);

  // TODO cant this be done better?
  Teuchos::RCP<Epetra_Vector> vdisp = Teuchos::rcp(new Epetra_Vector(*structure_field()->Dispnp()));

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return fluid_to_struct(MBFluidField()->RelaxationSolve(struct_to_fluid(idisp), Dt()));
  }
  else
  {
    // normal fluid solve
    // the displacement -> velocity conversion at the interface
    const Teuchos::RCP<Epetra_Vector> ivel = interface_velocity(idisp);

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) MBFluidField()->SetItemax(mfresitemax_ + 1);

    Teuchos::RCP<Adapter::FluidAle> fluidale =
        Teuchos::rcp_dynamic_cast<Adapter::FluidAle>(MBFluidField());

    icorrector_->set_interface_displacements(idisp, structure_fluid_coupling());

    // important difference to dirichletneumann.cpp: vdisp is mapped from structure to ale here
    fluidale->nonlinear_solve_vol_coupl(
        structure_to_ale(vdisp), struct_to_fluid(ivel), icorrector_);

    MBFluidField()->SetItemax(itemax);

    return fluid_to_struct(MBFluidField()->extract_interface_forces());
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::DirichletNeumannVolCoupl::extract_previous_interface_solution()
{
  iveln_ = fluid_to_struct(MBFluidField()->extract_interface_veln());
  idispn_ = structure_field()->extract_interface_dispn();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVolCoupl::structure_to_ale(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVolCoupl::ale_to_structure(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::InterfaceCorrector::setup(Teuchos::RCP<Adapter::FluidAle> fluidale)
{
  fluidale_ = fluidale;

  volcorrector_ = Teuchos::rcp(new VolCorrector);
  volcorrector_->setup(Global::Problem::Instance()->NDim(), fluidale);

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::InterfaceCorrector::set_interface_displacements(
    Teuchos::RCP<Epetra_Vector>& idisp_struct, Core::Adapter::Coupling& icoupfs)
{
  idisp_ = idisp_struct;
  icoupfs_ = Teuchos::rcpFromRef(icoupfs);

  deltadisp_ = Teuchos::null;
  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::InterfaceCorrector::correct_interface_displacements(
    Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  if (icoupfs_ == Teuchos::null) FOUR_C_THROW("Coupling adapter not set!");
  if (idisp_ == Teuchos::null) FOUR_C_THROW("Interface displacements not set!");

  // std::cout<<*finterface->FullMap()<<std::endl;
  // std::cout<<*disp_fluid<<std::endl;
  deltadisp_ = Core::LinAlg::CreateVector(*finterface->fsi_cond_map(), true);

  Core::LinAlg::Export(*disp_fluid, *deltadisp_);
  // deltadisp_ = finterface->extract_fsi_cond_vector(disp_fluid);

  // FOUR_C_THROW("stop");

  Teuchos::RCP<Epetra_Vector> idisp_fluid_corrected = icoupfs_->MasterToSlave(idisp_);

  deltadisp_->Update(1.0, *idisp_fluid_corrected, -1.0);

  Core::LinAlg::Export(*idisp_fluid_corrected, *disp_fluid);
  // finterface->insert_fsi_cond_vector(idisp_fluid_corrected,disp_fluid);

  volcorrector_->correct_vol_displacements(fluidale_, deltadisp_, disp_fluid, finterface);

  // reset
  idisp_ = Teuchos::null;
  icoupfs_ = Teuchos::null;

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::correct_vol_displacements(Teuchos::RCP<Adapter::FluidAle> fluidale,
    Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  if (fluidale->ale_field()->discretization()->Comm().MyPID() == 0)
    std::cout << "******************   FSI Volume Correction Step   **********************"
              << std::endl;

  // correction step in parameter space
  if (true) correct_vol_displacements_para_space(fluidale, deltadisp, disp_fluid, finterface);
  // correction step in physical space
  else
    correct_vol_displacements_phys_space(fluidale, deltadisp, disp_fluid, finterface);

  // output
  if (fluidale->ale_field()->discretization()->Comm().MyPID() == 0)
    std::cout << "******************FSI Volume Correction Step Done***********************"
              << std::endl;


  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::correct_vol_displacements_para_space(
    Teuchos::RCP<Adapter::FluidAle> fluidale, Teuchos::RCP<Epetra_Vector> deltadisp,
    Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  Teuchos::RCP<Epetra_Vector> correction = Teuchos::rcp(new Epetra_Vector(disp_fluid->Map(), true));
  Teuchos::RCP<Epetra_Vector> DofColMapDummy = Teuchos::rcp(
      new Epetra_Vector(*fluidale->fluid_field()->discretization()->DofColMap(), true));
  Core::LinAlg::Export(*deltadisp, *DofColMapDummy);

  const double tol = 1e-5;

  // loop over ale eles
  for (std::map<int, std::vector<int>>::iterator it = fluidalenodemap_.begin();
       it != fluidalenodemap_.end(); ++it)
  {
    Core::Elements::Element* aleele = fluidale->ale_field()->discretization()->gElement(it->first);

    // loop over fluid volume nodes within one ale FSI element
    for (size_t i = 0; i < it->second.size(); ++i)
    {
      int gid = it->second[i];
      Core::Nodes::Node* fluidnode = fluidale->fluid_field()->discretization()->gNode(gid);

      if (fluidnode->Owner() != fluidale->ale_field()->discretization()->Comm().MyPID()) continue;

      double gpos[3] = {fluidnode->X()[0], fluidnode->X()[1], fluidnode->X()[2]};
      double lpos[3] = {0.0, 0.0, 0.0};
      if (aleele->Shape() == Core::FE::CellType::quad4)
        Mortar::UTILS::GlobalToLocal<Core::FE::CellType::quad4>(*aleele, gpos, lpos);
      else if (aleele->Shape() == Core::FE::CellType::hex8)
        Mortar::UTILS::GlobalToLocal<Core::FE::CellType::hex8>(*aleele, gpos, lpos);
      else
        FOUR_C_THROW("ERROR: element type not implemented!");

      if (lpos[0] < -1.0 - tol || lpos[1] < -1.0 - tol || lpos[2] < -1.0 - tol ||
          lpos[0] > 1.0 + tol || lpos[1] > 1.0 + tol || lpos[2] > 1.0 + tol)
        continue;

      double dist = 1.0e12;
      int id = -1;
      for (size_t k = 0; k < fluidalenode_fs_imap_[it->first].size(); ++k)
      {
        int gidfsi = fluidalenode_fs_imap_[it->first][k];
        Core::Nodes::Node* fluidnodeFSI = fluidale->fluid_field()->discretization()->gNode(gidfsi);

        double gposFSI[3] = {fluidnodeFSI->X()[0], fluidnodeFSI->X()[1], fluidnodeFSI->X()[2]};
        double lposFSI[3] = {0.0, 0.0, 0.0};
        if (aleele->Shape() == Core::FE::CellType::quad4)
          Mortar::UTILS::GlobalToLocal<Core::FE::CellType::quad4>(*aleele, gposFSI, lposFSI);
        else if (aleele->Shape() == Core::FE::CellType::hex8)
          Mortar::UTILS::GlobalToLocal<Core::FE::CellType::hex8>(*aleele, gposFSI, lposFSI);
        else
          FOUR_C_THROW("ERROR: element type not implemented!");

        if (lposFSI[0] < -1.0 - tol || lposFSI[1] < -1.0 - tol || lposFSI[2] < -1.0 - tol ||
            lposFSI[0] > 1.0 + tol || lposFSI[1] > 1.0 + tol || lposFSI[2] > 1.0 + tol)
          FOUR_C_THROW("ERROR: wrong parameter space coordinates!");

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

      if (aleele->Shape() == Core::FE::CellType::quad4 or
          aleele->Shape() == Core::FE::CellType::hex8)
        fac = 1.0 - 0.5 * dist;
      else
        FOUR_C_THROW("ERROR: element type not implemented!");

      // safety
      if (dist > 2.0) fac = 0.0;

      Core::Nodes::Node* fluidnodeFSI = fluidale->fluid_field()->discretization()->gNode(id);
      std::vector<int> temp = fluidale->fluid_field()->discretization()->Dof(fluidnodeFSI);
      std::vector<int> dofsFSI;
      for (int idof = 0; idof < dim_; idof++) dofsFSI.push_back(temp[idof]);

      // extract local values of the global vectors
      std::vector<double> FSIdisp(dofsFSI.size());
      Core::FE::ExtractMyValues(*DofColMapDummy, FSIdisp, dofsFSI);

      std::vector<int> temp2 = fluidale->fluid_field()->discretization()->Dof(fluidnode);
      std::vector<int> dofs;
      for (int idof = 0; idof < dim_; idof++) dofs.push_back(temp2[idof]);

      Core::LinAlg::SerialDenseVector gnode(dim_);
      std::vector<int> lmowner(dim_);
      for (int idof = 0; idof < dim_; idof++)
      {
        gnode(idof) = fac * FSIdisp[idof];
        lmowner[idof] = fluidnode->Owner();
      }

      Core::LinAlg::Assemble(*correction, gnode, dofs, lmowner);
    }  // end fluid volume node loop
  }    // end ale fsi element loop

  // do correction
  disp_fluid->Update(1.0, *correction, 1.0);

  // calc norm
  double norm = 0.0;
  correction->Norm2(&norm);

  // output
  if (fluidale->ale_field()->discretization()->Comm().MyPID() == 0)
    std::cout << "Norm of correction (parameter space): " << norm << std::endl;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::correct_vol_displacements_phys_space(
    Teuchos::RCP<Adapter::FluidAle> fluidale, Teuchos::RCP<Epetra_Vector> deltadisp,
    Teuchos::RCP<Epetra_Vector> disp_fluid,
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface)
{
  Teuchos::RCP<Epetra_Vector> correction = Teuchos::rcp(new Epetra_Vector(disp_fluid->Map(), true));
  Teuchos::RCP<Epetra_Vector> DofColMapDummy = Teuchos::rcp(
      new Epetra_Vector(*fluidale->fluid_field()->discretization()->DofColMap(), true));
  Core::LinAlg::Export(*deltadisp, *DofColMapDummy);

  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPs =
      calc_background_dops(fluidale->fluid_field()->discretization());

  Teuchos::RCP<std::set<int>> FSIaleeles = Core::Conditions::conditioned_element_map(
      *fluidale->ale_field()->discretization(), "FSICoupling");

  // evaluate search
  for (int i = 0; i < fluidale->ale_field()->discretization()->NumMyColElements(); ++i)
  {
    // 1 map node into bele
    int gid = fluidale->ale_field()->discretization()->ElementColMap()->GID(i);
    Core::Elements::Element* aleele = fluidale->ale_field()->discretization()->gElement(gid);

    if (FSIaleeles->find(aleele->Id()) == FSIaleeles->end()) continue;
  }

  // do correction
  disp_fluid->Update(1.0, *correction, 1.0);

  // calc norm
  double norm = 0.0;
  correction->Norm2(&norm);

  // output
  if (fluidale->ale_field()->discretization()->Comm().MyPID() == 0)
    std::cout << "Norm of correction (physical space): " << norm << std::endl;

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::VolCorrector::setup(const int dim, Teuchos::RCP<Adapter::FluidAle> fluidale)
{
  if (fluidale->ale_field()->discretization()->Comm().MyPID() == 0)
    std::cout << "******************FSI Volume Correction Setup***********************"
              << std::endl;

  dim_ = dim;
  init_dop_normals();

  // init current positions
  std::map<int, Core::LinAlg::Matrix<3, 1>> currentpositions;

  for (int lid = 0; lid < fluidale->fluid_field()->discretization()->NumMyColElements(); ++lid)
  {
    Core::Elements::Element* sele = fluidale->fluid_field()->discretization()->lColElement(lid);

    // calculate slabs for every node on every element
    for (int k = 0; k < sele->num_node(); k++)
    {
      Core::Nodes::Node* node = sele->Nodes()[k];
      Core::LinAlg::Matrix<3, 1> currpos;

      currpos(0) = node->X()[0];
      currpos(1) = node->X()[1];
      currpos(2) = node->X()[2];

      currentpositions[node->Id()] = currpos;
    }
  }

  // init of 3D search tree
  search_tree_ = Teuchos::rcp(new Core::Geo::SearchTree(5));

  // find the bounding box of the elements and initialize the search tree
  const Core::LinAlg::Matrix<3, 2> rootBox =
      Core::Geo::getXAABBofDis(*fluidale->fluid_field()->discretization(), currentpositions);
  search_tree_->initializeTree(
      rootBox, *fluidale->fluid_field()->discretization(), Core::Geo::TreeType(Core::Geo::OCTTREE));


  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPs =
      calc_background_dops(fluidale->fluid_field()->discretization());

  Teuchos::RCP<std::set<int>> FSIaleeles = Core::Conditions::conditioned_element_map(
      *fluidale->ale_field()->discretization(), "FSICoupling");

  // evaluate search
  for (int i = 0; i < fluidale->ale_field()->discretization()->NumMyColElements(); ++i)
  {
    // 1 map node into bele
    int gid = fluidale->ale_field()->discretization()->ElementColMap()->GID(i);
    Core::Elements::Element* aleele = fluidale->ale_field()->discretization()->gElement(gid);

    if (FSIaleeles->find(aleele->Id()) == FSIaleeles->end()) continue;

    // get found elements from other discr.
    fluidaleelemap_[gid] = search(*aleele, CurrentDOPs);
  }  // end node loop

  Teuchos::RCP<Epetra_Map> FSIfluidnodes = Core::Conditions::ConditionNodeColMap(
      *fluidale->fluid_field()->discretization(), "FSICoupling");

  std::set<int> globalnodeids;
  // loop over ale eles
  for (std::map<int, std::vector<int>>::iterator it = fluidaleelemap_.begin();
       it != fluidaleelemap_.end(); ++it)
  {
    Core::Elements::Element* aleele = fluidale->ale_field()->discretization()->gElement(it->first);

    std::vector<int> localnodeids;
    std::vector<int> localnodeidsFSI;

    // loop over fluid eles
    for (size_t i = 0; i < it->second.size(); ++i)
    {
      int gid = it->second[i];
      Core::Elements::Element* fluidele = fluidale->fluid_field()->discretization()->gElement(gid);

      for (int j = 0; j < fluidele->num_node(); ++j)
      {
        const int nodegid = fluidele->NodeIds()[j];
        Core::Nodes::Node* fluidnode = fluidele->Nodes()[j];

        double gpos[3] = {fluidnode->X()[0], fluidnode->X()[1], fluidnode->X()[2]};
        double lpos[3] = {0.0, 0.0, 0.0};
        if (aleele->Shape() == Core::FE::CellType::quad4)
          Mortar::UTILS::GlobalToLocal<Core::FE::CellType::quad4>(*aleele, gpos, lpos);
        else if (aleele->Shape() == Core::FE::CellType::hex8)
          Mortar::UTILS::GlobalToLocal<Core::FE::CellType::hex8>(*aleele, gpos, lpos);
        else
          FOUR_C_THROW("ERROR: element type not implemented!");

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
    fluidalenode_fs_imap_[it->first] = localnodeidsFSI;
  }

  std::cout << "ALE elements found: " << fluidaleelemap_.size() << std::endl;

  if (fluidale->ale_field()->discretization()->Comm().MyPID() == 0)
    std::cout << "******************FSI Volume Correction Setup Done***********************"
              << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Init normals for Dop calculation                         farah 05/16|
 *----------------------------------------------------------------------*/
void FSI::VolCorrector::init_dop_normals()
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
std::map<int, Core::LinAlg::Matrix<9, 2>> FSI::VolCorrector::calc_background_dops(
    Teuchos::RCP<Core::FE::Discretization> searchdis)
{
  std::map<int, Core::LinAlg::Matrix<9, 2>> currentKDOPs;

  for (int lid = 0; lid < searchdis->NumMyColElements(); ++lid)
  {
    Core::Elements::Element* sele = searchdis->lColElement(lid);

    currentKDOPs[sele->Id()] = calc_dop(*sele);
  }

  return currentKDOPs;
}

/*----------------------------------------------------------------------*
 |  Calculate Dop for one Element                            farah 05/16|
 *----------------------------------------------------------------------*/
Core::LinAlg::Matrix<9, 2> FSI::VolCorrector::calc_dop(Core::Elements::Element& ele)
{
  Core::LinAlg::Matrix<9, 2> dop;

  // calculate slabs
  for (int j = 0; j < 9; j++)
  {
    // initialize slabs
    dop(j, 0) = 1.0e12;
    dop(j, 1) = -1.0e12;
  }

  // calculate slabs for every node on every element
  for (int k = 0; k < ele.num_node(); k++)
  {
    Core::Nodes::Node* node = ele.Nodes()[k];

    // get current node position
    std::array<double, 3> pos = {0.0, 0.0, 0.0};
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
std::vector<int> FSI::VolCorrector::search(
    Core::Elements::Element& ele, std::map<int, Core::LinAlg::Matrix<9, 2>>& currentKDOPs)
{
  // vector of global ids of found elements
  std::vector<int> gids;
  gids.clear();
  std::set<int> gid;
  gid.clear();

  Core::LinAlg::Matrix<9, 2> queryKDOP;

  // calc dop for considered element
  queryKDOP = calc_dop(ele);

  //**********************************************************
  // search for near elements to the background node's coord
  search_tree_->searchCollisions(currentKDOPs, queryKDOP, 0, gid);

  for (std::set<int>::iterator iter = gid.begin(); iter != gid.end(); ++iter) gids.push_back(*iter);

  return gids;
}

FOUR_C_NAMESPACE_CLOSE
