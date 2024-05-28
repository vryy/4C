/*----------------------------------------------------------------------------*/
/*! \file

\brief Utilities for FSI problems

\level 1

*/

/*----------------------------------------------------------------------------*/

#include "4C_fsi_utils.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_geometry_searchtree.hpp"
#include "4C_discretization_geometry_searchtree_service.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_surface.hpp"

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <NOX.H>
#include <NOX_Epetra.H>

#include <map>
#include <set>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DumpJacobian(::NOX::Epetra::Interface::Required& interface, double alpha,
    double beta, Teuchos::RCP<Epetra_Vector> soln, std::string filename)
{
  // that's really stupid again
  const Epetra_BlockMap& bmap = soln->Map();
  Epetra_Map map(
      bmap.NumGlobalElements(), bmap.NumMyElements(), bmap.MyGlobalElements(), 0, bmap.Comm());

  Teuchos::RCP<Epetra_CrsMatrix> jacobian =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, map.NumGlobalElements()));

  int nummyelements = map.NumMyElements();
  int mypos = CORE::LINALG::FindMyPos(nummyelements, map.Comm());
  double eta = 0.0;

  Epetra_Vector fo(*soln);
  Epetra_Vector fp(*soln);
  Epetra_Vector Jc(*soln);

  // Compute the RHS at the initial solution
  interface.computeF(*soln, fo, ::NOX::Epetra::Interface::Required::FD_Res);

  Epetra_Vector x_perturb = *soln;

  for (int i = 0; i < map.NumGlobalElements(); ++i)
  {
    if (map.Comm().MyPID() == 0) std::cout << "calculate column " << i << "\n";

    int proc = 0;
    int idx = 0;
    if (i >= mypos and i < mypos + nummyelements)
    {
      eta = alpha * (*soln)[i - mypos] + beta;
      x_perturb[i - mypos] += eta;
      idx = map.GID(i - mypos);
      proc = map.Comm().MyPID();
    }

    // Find what proc eta is on
    int broadcastProc = 0;
    map.Comm().SumAll(&proc, &broadcastProc, 1);

    // Send the perturbation variable, eta, to all processors
    map.Comm().Broadcast(&eta, 1, broadcastProc);

    map.Comm().Broadcast(&idx, 1, broadcastProc);

    // Compute the perturbed RHS
    interface.computeF(x_perturb, fp, ::NOX::Epetra::Interface::Required::FD_Res);

    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    Jc.Scale(1.0 / eta);

    // Insert nonzero column entries into the jacobian
    for (int j = 0; j < map.NumMyElements(); ++j)
    {
      int gid = map.GID(j);
      if (Jc[j] != 0.0)
      {
        int err = jacobian->SumIntoGlobalValues(gid, 1, &Jc[j], &idx);
        if (err > 0)
        {
          err = jacobian->InsertGlobalValues(gid, 1, &Jc[j], &idx);
        }
        if (err != 0) FOUR_C_THROW("Assembly failed");
      }
    }

    // Unperturb the solution vector
    x_perturb = *soln;
  }

  jacobian->FillComplete();

  EpetraExt::RowMatrixToMatlabFile(filename.c_str(), *jacobian);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::UTILS::FluidAleNodesDisjoint(
    Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> aledis)
{
  // flag indicating whether fluid and ALE node numbers have are non-overlapping or not
  bool isdisjoint = false;

  // try a simple check that should work for most cases
  if (fluiddis->NodeRowMap()->MaxAllGID() < aledis->NodeRowMap()->MinAllGID() or
      fluiddis->NodeRowMap()->MinAllGID() > aledis->NodeRowMap()->MaxAllGID())
  {
    // no overlap of node numbers
    isdisjoint = true;
  }
  else  // do a more sophisticated check
  {
    // get node row maps
    Teuchos::RCP<const Epetra_Map> fluidmap =
        Teuchos::rcp(new const Epetra_Map(*fluiddis->NodeRowMap()));
    Teuchos::RCP<const Epetra_Map> alemap =
        Teuchos::rcp(new const Epetra_Map(*aledis->NodeRowMap()));

    // Create intersection of fluid and ALE map
    std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
    intersectionmaps.push_back(fluidmap);
    intersectionmaps.push_back(alemap);
    Teuchos::RCP<Epetra_Map> intersectionmap =
        CORE::LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

    if (intersectionmap->NumGlobalElements() == 0) isdisjoint = true;
  }

  return isdisjoint;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// class SlideAleUtils
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::SlideAleUtils::SlideAleUtils(Teuchos::RCP<DRT::Discretization> structdis,
    Teuchos::RCP<DRT::Discretization> fluiddis, CORE::ADAPTER::CouplingMortar& coupsf,
    bool structcoupmaster, INPAR::FSI::SlideALEProj aleproj)
    : aletype_(aleproj)
{
  structcoupmaster_ = structcoupmaster;

  coupff_ = Teuchos::rcp(new CORE::ADAPTER::CouplingMortar(GLOBAL::Problem::Instance()->NDim(),
      GLOBAL::Problem::Instance()->mortar_coupling_params(),
      GLOBAL::Problem::Instance()->contact_dynamic_params(),
      GLOBAL::Problem::Instance()->spatial_approximation_type()));

  // declare struct objects in interface
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>> structelements;
  std::map<int, Teuchos::RCP<DRT::Element>> structmelements;
  std::map<int, Teuchos::RCP<DRT::Element>> structdelements;
  std::map<int, DRT::Node*> dummy1;                 // dummy map
  std::map<int, std::map<int, DRT::Node*>> dummy2;  // dummy map
  std::map<int, DRT::Node*> structmnodes;           // partial map of sticking structure nodes
  std::map<int, DRT::Node*> structdnodes;           // partial map of centerdisp structure nodes
  std::map<int, std::map<int, DRT::Node*>> structgnodes;  // complete map of strucutre nodes

  // initialize struct objects in interface
  CORE::Conditions::FindConditionObjects(
      *structdis, dummy2, structgnodes, structelements, "FSICoupling");
  CORE::Conditions::FindConditionObjects(
      *structdis, dummy1, structmnodes, structmelements, "FSICouplingNoSlide");
  CORE::Conditions::FindConditionObjects(
      *structdis, dummy1, structdnodes, structdelements, "FSICouplingCenterDisp");
  istructdispnodes_ = structdnodes;
  istructdispeles_ = structdelements;
  istructslideles_ = structelements;

  std::vector<int> slideeleidvector;

  std::map<int, Teuchos::RCP<DRT::Element>>::iterator eit;
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>::iterator meit;

  for (eit = structmelements.begin(); eit != structmelements.end(); eit++)
  {
    int err = 0;
    for (meit = istructslideles_.begin(); meit != istructslideles_.end(); meit++)
      err += meit->second.erase((*eit).first);
    if (!err) FOUR_C_THROW("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  int max_id = 0;
  // find max FSI condition ID
  for (meit = istructslideles_.begin(); meit != istructslideles_.end(); meit++)
  {
    //    for ( eit=meit->second.begin(); eit != meit->second.end(); eit++ )
    //    {
    //      //build slideeleidvector with unique distribution. Otherwise, AllreduceEMap() will
    //      complain in DEBUG if (structdis->Comm().MyPID()==(*eit).second->Owner())
    //        slideeleidvector.push_back((*eit).first);
    //    }
    //    const Epetra_Map slideelemap (-1, slideeleidvector.size(), slideeleidvector.data(), 0,
    //    structdis->Comm()); slideeleredmap_[meit->first] =
    //    CORE::LINALG::AllreduceEMap(slideelemap);
    if (meit->first > max_id) max_id = meit->first;
  }

  structdis->Comm().MaxAll(&max_id, &maxid_, 1);

  // declare fluid objects in interface
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>> fluidelements;
  std::map<int, Teuchos::RCP<DRT::Element>> fluidmelements;
  std::map<int, std::map<int, DRT::Node*>> fluidnodes;  // complete map of fluid nodes
  std::map<int, DRT::Node*> fluidmnodes;                // partial map of sticking fluid nodes

  // initialize struct objects in interface
  CORE::Conditions::FindConditionObjects(
      *fluiddis, fluidnodes, dummy2, fluidelements, "FSICoupling");
  CORE::Conditions::FindConditionObjects(
      *fluiddis, fluidmnodes, dummy1, fluidmelements, "FSICouplingNoSlide");
  ifluidconfnodes_ = fluidmnodes;
  ifluidslidnodes_ = fluidnodes;
  ifluidslideles_ = fluidelements;

  for (eit = fluidmelements.begin(); eit != fluidmelements.end(); eit++)
  {
    int err = 0;
    for (meit = ifluidslideles_.begin(); meit != ifluidslideles_.end(); meit++)
      err += meit->second.erase((*eit).first);
    if (!err) FOUR_C_THROW("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  std::map<int, DRT::Node*>::iterator nit;
  std::map<int, std::map<int, DRT::Node*>>::iterator mnit;
  for (nit = ifluidconfnodes_.begin(); nit != ifluidconfnodes_.end(); nit++)
  {
    int err = 0;
    for (mnit = ifluidslidnodes_.begin(); mnit != ifluidslidnodes_.end(); mnit++)
      err += mnit->second.erase((*nit).first);
    if (!err) FOUR_C_THROW("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  Teuchos::RCP<Epetra_Map> structdofrowmap;
  Teuchos::RCP<Epetra_Map> fluiddofrowmap;


  // useful displacement vectors
  if (structcoupmaster_)
  {
    structdofrowmap_ = coupsf.MasterDofMap();
    fluiddofrowmap_ = coupsf.SlaveDofMap();
  }
  else
  {
    structdofrowmap_ = coupsf.SlaveDofMap();
    fluiddofrowmap_ = coupsf.MasterDofMap();
  }

  Teuchos::RCP<Epetra_Map> dofrowmap =
      CORE::LINALG::MergeMap(*structdofrowmap_, *fluiddofrowmap_, true);
  idispms_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  iprojhist_ = Teuchos::rcp(new Epetra_Vector(*fluiddofrowmap_, true));


  centerdisptotal_.resize(GLOBAL::Problem::Instance()->NDim());

  redundant_elements(coupsf, structdis->Comm());

  maxmindist_ = 1.0e-1;

  // coupling condition at the fsi interface: displacements (=number spacial dimensions) are
  // coupled) e.g.: 3D: coupleddof = [1, 1, 1]
  std::vector<int> coupleddof(GLOBAL::Problem::Instance()->NDim(), 1);

  // this setup only initialize two sets of identical mortar elements (master and slave)
  // -> projection matrix is a unity matrix
  coupff_->Setup(
      fluiddis, fluiddis, Teuchos::null, coupleddof, "FSICoupling", fluiddis->Comm(), false, true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::Remeshing(ADAPTER::FSIStructureWrapper& structure,
    Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<Epetra_Vector> idispale,
    Teuchos::RCP<Epetra_Vector> iprojdispale, CORE::ADAPTER::CouplingMortar& coupsf,
    const Epetra_Comm& comm)
{
  Teuchos::RCP<Epetra_Vector> idisptotal = structure.extract_interface_dispnp();
  const int dim = GLOBAL::Problem::Instance()->NDim();

  // project sliding fluid nodes onto struct interface surface
  slide_projection(structure, fluiddis, idispale, iprojdispale, coupsf, comm);

  // For the NON sliding ALE Nodes, use standard ALE displacements

  std::map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = ifluidconfnodes_.begin(); nodeiter != ifluidconfnodes_.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    std::vector<int> lids(dim);
    for (int p = 0; p < dim; p++)
      // lids of gids of node
      lids[p] = fluiddofrowmap_->LID((fluiddis->Dof(node))[p]);

    // current coord of ale node = ref coord + ifluid_
    std::vector<double> finaldxyz(dim);

    for (int p = 0; p < dim; p++) finaldxyz[p] = (*idispale)[(lids[p])];

    int err = iprojdispale->ReplaceMyValues(dim, finaldxyz.data(), lids.data());
    if (err == 1) FOUR_C_THROW("error while replacing values");
  }

  // merge displacement values of interface nodes (struct+fluid) into idispms_ for mortar
  idispms_->PutScalar(0.0);

  Teuchos::RCP<Epetra_Map> dofrowmap =
      CORE::LINALG::MergeMap(*structdofrowmap_, *fluiddofrowmap_, true);
  Teuchos::RCP<Epetra_Import> msimpo =
      Teuchos::rcp(new Epetra_Import(*dofrowmap, *structdofrowmap_));
  Teuchos::RCP<Epetra_Import> slimpo =
      Teuchos::rcp(new Epetra_Import(*dofrowmap, *fluiddofrowmap_));

  idispms_->Import(*idisptotal, *msimpo, Add);
  idispms_->Import(*iprojdispale, *slimpo, Add);

  iprojhist_->Update(1.0, *iprojdispale, 0.0);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::EvaluateMortar(Teuchos::RCP<Epetra_Vector> idispstruct,
    Teuchos::RCP<Epetra_Vector> idispfluid, CORE::ADAPTER::CouplingMortar& coupsf)
{
  // merge displacement values of interface nodes (struct+fluid) into idispms_ for mortar
  idispms_->PutScalar(0.0);

  Teuchos::RCP<Epetra_Map> dofrowmap =
      CORE::LINALG::MergeMap(*structdofrowmap_, *fluiddofrowmap_, true);
  Teuchos::RCP<Epetra_Import> master_importer =
      Teuchos::rcp(new Epetra_Import(*dofrowmap, *structdofrowmap_));
  Teuchos::RCP<Epetra_Import> slave_importer =
      Teuchos::rcp(new Epetra_Import(*dofrowmap, *fluiddofrowmap_));

  if (idispms_->Import(*idispstruct, *master_importer, Add))
    FOUR_C_THROW("Import operation failed.");
  if (idispms_->Import(*idispfluid, *slave_importer, Add)) FOUR_C_THROW("Import operation failed.");

  // new D,M,Dinv out of disp of struct and fluid side
  coupsf.Evaluate(idispms_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::EvaluateFluidMortar(
    Teuchos::RCP<Epetra_Vector> ima, Teuchos::RCP<Epetra_Vector> isl)
{
  // new D,M,Dinv out of fluid disp before and after sliding
  coupff_->Evaluate(ima, isl);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::UTILS::SlideAleUtils::InterpolateFluid(
    Teuchos::RCP<const Epetra_Vector> uold)
{
  Teuchos::RCP<Epetra_Vector> unew = coupff_->MasterToSlave(uold);
  unew->ReplaceMap(uold->Map());

  return unew;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FSI::UTILS::SlideAleUtils::centerdisp(
    ADAPTER::FSIStructureWrapper& structure, const Epetra_Comm& comm)
{
  Teuchos::RCP<DRT::Discretization> structdis = structure.discretization();

  Teuchos::RCP<Epetra_Vector> idispn = structure.extract_interface_dispn();
  Teuchos::RCP<Epetra_Vector> idisptotal = structure.extract_interface_dispnp();
  Teuchos::RCP<Epetra_Vector> idispstep = structure.extract_interface_dispnp();

  int err = idispstep->Update(-1.0, *idispn, 1.0);
  if (err != 0) FOUR_C_THROW("ERROR");

  const int dim = GLOBAL::Problem::Instance()->NDim();
  // get structure and fluid discretizations  and set stated for element evaluation
  const Teuchos::RCP<Epetra_Vector> idisptotalcol =
      CORE::LINALG::CreateVector(*structdis->DofColMap(), true);
  CORE::LINALG::Export(*idisptotal, *idisptotalcol);
  const Teuchos::RCP<Epetra_Vector> idispstepcol =
      CORE::LINALG::CreateVector(*structdis->DofColMap(), true);
  CORE::LINALG::Export(*idispstep, *idispstepcol);

  structdis->set_state("displacementtotal", idisptotalcol);
  structdis->set_state("displacementincr", idispstepcol);

  // define stuff needed by the elements
  Teuchos::ParameterList params;
  CORE::LINALG::SerialDenseMatrix elematrix1;
  CORE::LINALG::SerialDenseMatrix elematrix2;
  CORE::LINALG::SerialDenseVector elevector1;
  CORE::LINALG::SerialDenseVector elevector2;
  CORE::LINALG::SerialDenseVector elevector3;

  // prepare variables for length (2D) or area (3D) of the interface
  std::vector<double> mycenterdisp(dim);
  std::vector<double> centerdisp(dim);
  double mylengthcirc = 0.0;
  double lengthcirc = 0.0;

  // calculating the center displacement by evaluating structure interface elements
  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator elemiter;
  for (elemiter = istructdispeles_.begin(); elemiter != istructdispeles_.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> iele = elemiter->second;
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    iele->LocationVector(*structdis, lm, lmowner, lmstride);
    elevector2.size(1);    // length of circ with gaussinteg
    elevector3.size(dim);  // centerdisp part of ele

    params.set<std::string>("action", "calc_struct_centerdisp");
    int err = iele->Evaluate(
        params, *structdis, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
    if (err) FOUR_C_THROW("error while evaluating elements");
    mylengthcirc += elevector2[0];

    // disp of the interface
    for (int i = 0; i < dim; i++)
    {
      mycenterdisp[i] += elevector3[i];
    }
  }  // end of ele loop
  structdis->ClearState();

  // Communicate to 'assemble' length and center displacements
  comm.SumAll(&mylengthcirc, &lengthcirc, 1);
  comm.SumAll(mycenterdisp.data(), centerdisp.data(), dim);

  if (lengthcirc <= 1.0E-6) FOUR_C_THROW("Zero interface length!");

  // calculating the final disp of the interface and summation over all time steps
  for (int i = 0; i < dim; i++)
  {
    centerdisp[i] = centerdisp[i] / lengthcirc;
    centerdisptotal_[i] += centerdisp[i];
  }

  return centerdisp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, CORE::LINALG::Matrix<3, 1>> FSI::UTILS::SlideAleUtils::current_struct_pos(
    Teuchos::RCP<Epetra_Vector> reddisp, DRT::Discretization& interfacedis,
    std::map<int, double>& maxcoord)
{
  std::map<int, CORE::LINALG::Matrix<3, 1>> currentpositions;
  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator eleiter;
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>::const_iterator meleiter;

  // map with fully reduced struct element distribution
  for (meleiter = structreduelements_.begin(); meleiter != structreduelements_.end(); meleiter++)
  {
    maxcoord[meleiter->first] = 0.0;
    for (eleiter = meleiter->second.begin(); eleiter != meleiter->second.end(); eleiter++)
    {
      Teuchos::RCP<DRT::Element> tmpele = eleiter->second;

      const int* n = tmpele->NodeIds();

      // fill currentpositions
      for (int j = 0; j < tmpele->num_node(); j++)
      {
        const int gid = n[j];
        const DRT::Node* node = interfacedis.gNode(gid);
        std::vector<int> lm;
        lm.reserve(3);
        // extract global dof ids
        interfacedis.Dof(node, lm);
        std::vector<double> mydisp(3);
        CORE::LINALG::Matrix<3, 1> currpos;

        CORE::FE::ExtractMyValues(*reddisp, mydisp, lm);

        for (int a = 0; a < 3; a++)
        {
          currpos(a, 0) = node->X()[a] + mydisp[a];
        }
        if (abs(currpos(2, 0)) > maxcoord[meleiter->first])
          maxcoord[meleiter->first] = abs(currpos(2, 0));
        currentpositions[node->Id()] = currpos;
      }
    }
  }

  return currentpositions;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::slide_projection(
    ADAPTER::FSIStructureWrapper& structure, Teuchos::RCP<DRT::Discretization> fluiddis,
    Teuchos::RCP<Epetra_Vector> idispale, Teuchos::RCP<Epetra_Vector> iprojdispale,
    CORE::ADAPTER::CouplingMortar& coupsf, const Epetra_Comm& comm

)
{
  const int dim = GLOBAL::Problem::Instance()->NDim();

  Teuchos::RCP<Epetra_Vector> idispnp = structure.extract_interface_dispnp();

  // Redistribute displacement of structnodes on the interface to all processors.
  Teuchos::RCP<Epetra_Import> interimpo =
      Teuchos::rcp(new Epetra_Import(*structfullnodemap_, *structdofrowmap_));
  Teuchos::RCP<Epetra_Vector> reddisp = CORE::LINALG::CreateVector(*structfullnodemap_, true);
  reddisp->Import(*idispnp, *interimpo, Add);

  DRT::Discretization& interfacedis = coupsf.Interface()->Discret();
  std::map<int, double> rotrat;
  // currentpositions of struct nodes for the search tree (always 3 coordinates)
  std::map<int, CORE::LINALG::Matrix<3, 1>> currentpositions =
      current_struct_pos(reddisp, interfacedis, rotrat);

  // calculate structural interface center of gravity
  std::vector<double> centerdisp_v = centerdisp(structure, comm);

  Teuchos::RCP<Epetra_Vector> frotfull = CORE::LINALG::CreateVector(*fluiddofrowmap_, true);
  if (aletype_ == INPAR::FSI::ALEprojection_rot_z ||
      aletype_ == INPAR::FSI::ALEprojection_rot_zsphere)
  {
    rotation(coupsf.Interface()->Discret(), idispale, comm, rotrat, frotfull);
  }


  std::map<int, std::map<int, DRT::Node*>>::iterator mnit;
  for (mnit = ifluidslidnodes_.begin(); mnit != ifluidslidnodes_.end(); ++mnit)
  {
    // translation + projection
    std::map<int, DRT::Node*>::const_iterator nodeiter;
    for (nodeiter = mnit->second.begin(); nodeiter != mnit->second.end(); ++nodeiter)
    {
      // Project fluid nodes onto the struct interface
      // init of search tree
      Teuchos::RCP<CORE::GEO::SearchTree> searchTree = Teuchos::rcp(new CORE::GEO::SearchTree(5));
      const CORE::LINALG::Matrix<3, 2> rootBox =
          CORE::GEO::getXAABBofEles(structreduelements_[mnit->first], currentpositions);

      if (dim == 2)
        searchTree->initialize_tree_slide_ale(
            rootBox, structreduelements_[mnit->first], CORE::GEO::TreeType(CORE::GEO::QUADTREE));
      else if (dim == 3)
        searchTree->initialize_tree_slide_ale(
            rootBox, structreduelements_[mnit->first], CORE::GEO::TreeType(CORE::GEO::OCTTREE));
      else
        FOUR_C_THROW("wrong dimension");


      DRT::Node* node = nodeiter->second;
      std::vector<int> lids(dim);
      for (int p = 0; p < dim; p++)
        // lids of gids of node
        lids[p] = (fluiddofrowmap_)->LID((fluiddis->Dof(node))[p]);

      // current coord of ale node.
      // Initialize as coordinates of current node, which is extremely important for 2D!
      CORE::LINALG::Matrix<3, 1> alenodecurr(node->X().data());

      // compute ALE position to project from
      if (aletype_ == INPAR::FSI::ALEprojection_curr)
      {
        // current coord of ale node = ref + centerdispincr + history
        for (int p = 0; p < dim; p++)
          alenodecurr(p, 0) = (node->X()[p]) + centerdisp_v[p] + 1.0 * (*iprojhist_)[(lids[p])];
      }
      else if (aletype_ == INPAR::FSI::ALEprojection_ref)
      {
        // current coord of ale node = ref + centerdisp
        for (int p = 0; p < dim; p++) alenodecurr(p, 0) = node->X()[p] + centerdisptotal_[p];
      }
      else if (aletype_ == INPAR::FSI::ALEprojection_rot_z ||
               aletype_ == INPAR::FSI::ALEprojection_rot_zsphere)
      {
        // current coord of ale node = ref + centerdisp
        for (int p = 0; p < dim; p++)
        {
          alenodecurr(p, 0) = node->X()[p] + (*idispale)[(lids[p])] -
                              1.0 * rotrat[mnit->first] * (*frotfull)[(lids[p])];
        }
      }
      else
        FOUR_C_THROW("you should not turn up here!");


      // final displacement of projection
      std::vector<double> finaldxyz(dim);

      // search for near elements next to the query point (ie within a radius of 2x maxmindist)
      std::map<int, std::set<int>> closeeles = searchTree->search_elements_in_radius(
          interfacedis, currentpositions, alenodecurr, maxmindist_, 0);
      // if no close elements could be found, try with a much larger radius and print a warning
      if (closeeles.empty())
      {
        const double enlarge_factor = 100;
        std::cout << "WARNING: no elements found in radius r=" << maxmindist_
                  << ". Will try once with a " << static_cast<int>(enlarge_factor)
                  << "-times bigger radius!" << std::endl;
        closeeles = searchTree->search_elements_in_radius(
            interfacedis, currentpositions, alenodecurr, enlarge_factor * maxmindist_, 0);
        maxmindist_ *= 10.0;

        // if still no element is found, complain about it!
        if (closeeles.empty()) FOUR_C_THROW("No elements in a large radius! Should not happen!");
      }
      // search for the nearest point to project on
      CORE::LINALG::Matrix<3, 1> minDistCoords;
      if (dim == 2)
      {
        CORE::GEO::nearest2DObjectInNode(Teuchos::rcp(&interfacedis, false),
            structreduelements_[mnit->first], currentpositions, closeeles, alenodecurr,
            minDistCoords);
        finaldxyz[0] = minDistCoords(0, 0) - node->X()[0];
        finaldxyz[1] = minDistCoords(1, 0) - node->X()[1];
      }
      else
      {
        CORE::GEO::nearest3DObjectInNode(Teuchos::rcp(&interfacedis, false),
            structreduelements_[mnit->first], currentpositions, closeeles, alenodecurr,
            minDistCoords);
        finaldxyz[0] = minDistCoords(0, 0) - node->X()[0];
        finaldxyz[1] = minDistCoords(1, 0) - node->X()[1];
        finaldxyz[2] = minDistCoords(2, 0) - node->X()[2];
      }

      // store displacement into parallel vector
      int err = iprojdispale->ReplaceMyValues(dim, finaldxyz.data(), lids.data());
      if (err == 1) FOUR_C_THROW("error while replacing values");
    }
  }
}

void FSI::UTILS::SlideAleUtils::redundant_elements(
    CORE::ADAPTER::CouplingMortar& coupsf, const Epetra_Comm& comm)
{
  // We need the structure elements (NOT THE MORTAR-ELEMENTS!) on every processor for the projection
  // of the fluid nodes. Furthermore we need the current position of the structnodes on every
  // processor. Elements provided by interface discretization, necessary maps provided by interface.

  int soffset = 0;
  int foffset = 0;
  if (structcoupmaster_)
  {
    structfullnodemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowDofs()));
    structfullelemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowElements()));
    fluidfullnodemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowDofs()));
    fluidfullelemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowElements()));
    soffset = 0;
    foffset = fluidfullelemap_->MinMyGID();
  }
  else
  {
    fluidfullnodemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowDofs()));
    fluidfullelemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowElements()));
    structfullnodemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowDofs()));
    structfullelemap_ = CORE::LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowElements()));
    soffset = structfullelemap_->MinMyGID();
    foffset = 0;
  }

  DRT::Discretization& interfacedis = coupsf.Interface()->Discret();

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>::iterator mapit;
  // build redundant version istructslideles_;
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator eit;
  int dim = GLOBAL::Problem::Instance()->NDim();

  for (int i = 0; i <= maxid_; ++i)
  {
    std::vector<int> vstruslideleids;  // vector for ele ids
    if (istructslideles_.find(i) != istructslideles_.end())
    {
      for (eit = istructslideles_[i].begin(); eit != istructslideles_[i].end(); eit++)
      {
        if (interfacedis.Comm().MyPID() == (*eit).second->Owner())
          vstruslideleids.push_back(eit->first);
      }
    }
    int globsum = 0;
    int partsum = (vstruslideleids.size());

    comm.SumAll(&partsum, &globsum, 1);
    // map with ele ids
    Epetra_Map mstruslideleids(globsum, vstruslideleids.size(), vstruslideleids.data(), 0, comm);
    // redundant version of it
    Epetra_Map redmstruslideleids(*CORE::LINALG::AllreduceEMap(mstruslideleids));

    for (int eleind = 0; eleind < redmstruslideleids.NumMyElements(); eleind++)
    {
      {
        DRT::Element* tmpele = interfacedis.gElement(redmstruslideleids.GID(eleind) + soffset);
        if (dim == 3)
        {
          structreduelements_[i][tmpele->Id()] =
              Teuchos::rcp(new DRT::ELEMENTS::StructuralSurface(tmpele->Id(), tmpele->Owner(),
                  tmpele->num_node(), tmpele->NodeIds(), tmpele->Nodes(), &(*tmpele), 0));
        }
        else if (dim == 2)
        {
          structreduelements_[i][tmpele->Id()] =
              Teuchos::rcp(new DRT::ELEMENTS::StructuralLine(tmpele->Id(), tmpele->Owner(),
                  tmpele->num_node(), tmpele->NodeIds(), tmpele->Nodes(), &(*tmpele), 0));
        }
      }
    }

    if (ifluidslideles_.find(i) != ifluidslideles_.end())
    {
      for (eit = ifluidslideles_[i].begin(); eit != ifluidslideles_[i].end(); eit++)
      {
        DRT::Element* tmpele = interfacedis.gElement(eit->first + foffset);
        if (dim == 3)
        {
          ifluidslidstructeles_[i][tmpele->Id()] =
              Teuchos::rcp(new DRT::ELEMENTS::StructuralSurface(tmpele->Id(), tmpele->Owner(),
                  tmpele->num_node(), tmpele->NodeIds(), tmpele->Nodes(), &(*tmpele), 0));
        }
        else if (dim == 2)
        {
          ifluidslidstructeles_[i][tmpele->Id()] =
              Teuchos::rcp(new DRT::ELEMENTS::StructuralLine(tmpele->Id(), tmpele->Owner(),
                  tmpele->num_node(), tmpele->NodeIds(), tmpele->Nodes(), &(*tmpele), 0));
        }
      }
    }
  }
}


void FSI::UTILS::SlideAleUtils::rotation(DRT::Discretization& mtrdis,  ///< fluid discretization
    Teuchos::RCP<Epetra_Vector> idispale,  ///< vector of ALE displacements
    const Epetra_Comm& comm,               ///< communicator
    std::map<int, double>& rotrat,         ///< rotation ratio of tangential displacements
    Teuchos::RCP<Epetra_Vector> rotfull  ///< vector of full displacements in tangential directions
)
{
  Teuchos::RCP<Epetra_Vector> idispstep = CORE::LINALG::CreateVector(*fluiddofrowmap_, false);
  idispstep->Update(1.0, *idispale, -1.0, *iprojhist_, 0.0);

  // get structure and fluid discretizations  and set state for element evaluation
  const Teuchos::RCP<Epetra_Vector> idispstepcol =
      CORE::LINALG::CreateVector(*mtrdis.DofColMap(), false);
  CORE::LINALG::Export(*idispstep, *idispstepcol);
  const Teuchos::RCP<Epetra_Vector> idispnpcol =
      CORE::LINALG::CreateVector(*mtrdis.DofColMap(), false);
  CORE::LINALG::Export(*idispale, *idispnpcol);

  mtrdis.set_state("displacementnp", idispnpcol);
  mtrdis.set_state("displacementincr", idispstepcol);

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>::iterator melit;
  for (int i = 0; i <= maxid_; ++i)
  {
    // prepare variables for length (2D) or area (3D) of the interface
    double myrotation = 0.0;
    double rotation = 0.0;
    double mylengthcirc = 0.0;
    double lengthcirc = 0.0;
    double maxcoord = 0.0;
    if (ifluidslidstructeles_.find(i) != ifluidslidstructeles_.end()) maxcoord = rotrat[i];

    std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator elemiter;
    for (elemiter = ifluidslidstructeles_[i].begin(); elemiter != ifluidslidstructeles_[i].end();
         elemiter++)
    {
      // define stuff needed by the elements
      CORE::LINALG::SerialDenseMatrix elematrix1;
      CORE::LINALG::SerialDenseMatrix elematrix2;
      CORE::LINALG::SerialDenseVector elevector1;
      CORE::LINALG::SerialDenseVector elevector2;
      CORE::LINALG::SerialDenseVector elevector3;
      Teuchos::ParameterList params;

      Teuchos::RCP<DRT::Element> iele = elemiter->second;
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      iele->LocationVector(mtrdis, lm, lmowner, lmstride);
      elevector2.size(1);  // circumference (2D) or surface area (3D) of the considered elements
      elevector3.size(1);  // normalized displacement in tangential direction ('rotation')

      params.set<std::string>("action", "calc_struct_rotation");
      params.set<double>("maxcoord", maxcoord);
      params.set<INPAR::FSI::SlideALEProj>("aletype", aletype_);
      int err = iele->Evaluate(
          params, mtrdis, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("error while evaluating elements");

      mylengthcirc += elevector2[0];
      // disp of the interface
      myrotation += elevector3[0];
    }  // end of ele loop

    // Communicate to 'assemble' length and center displacements
    comm.SumAll(&mylengthcirc, &lengthcirc, 1);
    comm.SumAll(&myrotation, &rotation, 1);

    if (lengthcirc >= 1.0E-6)
    {
      // calculating the final disp of the interface and summation over all time steps
      rotrat[i] = rotation / lengthcirc;
    }

    // second round!
    // compute correction displacement to account for rotation
    for (elemiter = ifluidslidstructeles_[i].begin(); elemiter != ifluidslidstructeles_[i].end();
         elemiter++)
    {
      // define stuff needed by the elements
      CORE::LINALG::SerialDenseMatrix elematrix1;
      CORE::LINALG::SerialDenseMatrix elematrix2;
      CORE::LINALG::SerialDenseVector elevector1;
      CORE::LINALG::SerialDenseVector elevector2;
      CORE::LINALG::SerialDenseVector elevector3;
      Teuchos::ParameterList params;

      Teuchos::RCP<DRT::Element> iele = elemiter->second;
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      iele->LocationVector(mtrdis, lm, lmowner, lmstride);
      elevector1.size(lm.size());

      params.set<std::string>("action", "calc_undo_struct_rotation");
      params.set<double>("maxcoord", maxcoord);
      params.set<INPAR::FSI::SlideALEProj>("aletype", aletype_);
      int err = iele->Evaluate(
          params, mtrdis, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      if (err) FOUR_C_THROW("error while evaluating elements");

      CORE::LINALG::Assemble(*rotfull, elevector1, lm, lmowner);
    }
  }
  mtrdis.ClearState();

  return;
}

void FSI::UTILS::SlideAleUtils::output_restart(IO::DiscretizationWriter& output)
{
  output.WriteVector("projhist", iprojhist_);

  return;
}

void FSI::UTILS::SlideAleUtils::read_restart(IO::DiscretizationReader& reader)
{
  reader.ReadVector(iprojhist_, "projhist");
}

FOUR_C_NAMESPACE_CLOSE
