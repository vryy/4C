/*----------------------------------------------------------------------*/
/*! \file
 \brief helper function/class for multiphase porous flow problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#include "porofluidmultiphase_utils.H"

#include "porofluidmultiphase_timint_ost.H"

#include "../drt_mat/material.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_adapter/ad_porofluidmultiphase.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_geometry/intersection_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::UTILS::SetupMaterial(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& fluid_disname)
{
  // get the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis(fluid_disname);

  // initialize material map
  std::map<int, int> matmap;
  {
    // get the cloning material map from the .dat file
    std::map<std::pair<std::string, std::string>, std::map<int, int>> clonefieldmatmap =
        DRT::Problem::Instance()->CloningMaterialMap();
    if (clonefieldmatmap.size() < 1)
      dserror("At least one material pairing required in --CLONING MATERIAL MAP.");

    // check if the current discretization is included in the material map
    std::pair<std::string, std::string> key(fluid_disname, struct_disname);
    matmap = clonefieldmatmap[key];
    if (matmap.size() < 1)
      dserror("Key pair '%s/%s' not defined in --CLONING MATERIAL MAP.", fluid_disname.c_str(),
          struct_disname.c_str());
  }


  // number of column elements within fluid discretization
  const int numelements = fluiddis->NumMyColElements();

  // loop over column elements
  for (int i = 0; i < numelements; ++i)
  {
    // get current element
    DRT::Element* ele = fluiddis->lColElement(i);

    // find the corresponding material in the matmap
    int src_matid = ele->Material()->Parameter()->Id();
    std::map<int, int>::iterator mat_iter = matmap.find(src_matid);
    if (mat_iter != matmap.end())
    {
      // get the ID of the secondary material
      const int tar_matid = mat_iter->second;
      // build the material usilng the factory
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(tar_matid);

      // add secondary material to poro fluid element
      if (ele->AddMaterial(mat) != 2) dserror("unexpected number of materials!");
    }
    else
    {
      // before we stop, print the material id map
      std::cout << "Material map on PROC " << comm.MyPID() << ":" << std::endl;
      for (mat_iter = matmap.begin(); mat_iter != matmap.end(); mat_iter++)
        std::cout << mat_iter->first << " -> " << mat_iter->second << std::endl;

      dserror("no matching material ID (%d) in map", src_matid);
    }

  }  // end loop over column elements

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> POROFLUIDMULTIPHASE::UTILS::ConvertDofVectorToNodeBasedMultiVector(
    const DRT::Discretization& dis, const Epetra_Vector& vector, const int nds,
    const int numdofpernode)
{
  // initialize multi vector
  Teuchos::RCP<Epetra_MultiVector> multi =
      Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), numdofpernode, true));

  // get maps
  const Epetra_BlockMap& vectormap = vector.Map();

  // loop over nodes of the discretization
  for (int inode = 0; inode < dis.NumMyRowNodes(); ++inode)
  {
    // get current node
    DRT::Node* node = dis.lRowNode(inode);
    // copy each dof value of node
    for (int idof = 0; idof < numdofpernode; ++idof)
      (*multi)[idof][inode] = vector[vectormap.LID(dis.Dof(nds, node, idof))];
  }

  return multi;
}

/*----------------------------------------------------------------------*
 | create algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::PoroFluidMultiphase> POROFLUIDMULTIPHASE::UTILS::CreateAlgorithm(
    INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme,
    Teuchos::RCP<DRT::Discretization> dis, const int linsolvernumber,
    const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams,
    FILE* errfile, Teuchos::RCP<IO::DiscretizationWriter> output)
{
  // Creation of Coupled Problem algortihm.
  Teuchos::RCP<ADAPTER::PoroFluidMultiphase> algo = Teuchos::null;

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------

  switch (timintscheme)
  {
    case INPAR::POROFLUIDMULTIPHASE::timeint_one_step_theta:
    {
      // create algorithm
      algo = Teuchos::rcp(new POROFLUIDMULTIPHASE::TimIntOneStepTheta(
          dis, linsolvernumber, probparams, poroparams, errfile, output));
      break;
    }
    default:
      dserror("Unknown time-integration scheme for multiphase poro fluid problem");
      break;
  }

  return algo;
}

/*--------------------------------------------------------------------------*
 | perform extended ghosting for artery dis                kremheller 03/19 |
 *--------------------------------------------------------------------------*/
std::map<int, std::set<int>> POROFLUIDMULTIPHASE::UTILS::ExtendedGhostingArteryDiscretization(
    Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis)
{
  // user output
  if (contdis->Comm().MyPID() == 0)
    std::cout
        << "\n<<<<<<<<<<<<<<< Starting extended ghosting of artery discretization >>>>>>>>>>>>>>>\n"
        << std::endl;

  artdis->FillComplete();
  if (!contdis->Filled()) contdis->FillComplete();

  // create the fully overlapping search discretization
  Teuchos::RCP<DRT::Discretization> artsearchdis = CreateSearchDiscretization(artdis);

  // to be filled with additional elements to be ghosted
  std::set<int> elecolset;
  const Epetra_Map* elecolmap = artdis->ElementColMap();
  for (int lid = 0; lid < elecolmap->NumMyElements(); ++lid)
  {
    int gid = elecolmap->GID(lid);
    elecolset.insert(gid);
  }

  // to be filled with additional nodes to be ghosted
  std::set<int> nodecolset;
  const Epetra_Map* nodecolmap = artdis->NodeColMap();
  for (int lid = 0; lid < nodecolmap->NumMyElements(); ++lid)
  {
    int gid = nodecolmap->GID(lid);
    nodecolset.insert(gid);
  }

  // search with the fully overlapping discretization
  std::map<int, std::set<int>> nearbyelepairs =
      OctTreeSearch(contdis, artdis, artsearchdis, &elecolset, &nodecolset);

  // extended ghosting for elements
  std::vector<int> coleles(elecolset.begin(), elecolset.end());
  Teuchos::RCP<const Epetra_Map> extendedelecolmap =
      Teuchos::rcp(new Epetra_Map(-1, coleles.size(), &coleles[0], 0, contdis->Comm()));

  artdis->ExportColumnElements(*extendedelecolmap);

  // extended ghosting for nodes
  std::vector<int> colnodes(nodecolset.begin(), nodecolset.end());
  Teuchos::RCP<const Epetra_Map> extendednodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, colnodes.size(), &colnodes[0], 0, contdis->Comm()));

  artdis->ExportColumnNodes(*extendednodecolmap);

  // fill and inform user
  artdis->FillComplete();
  DRT::UTILS::PrintParallelDistribution(*artdis);

  // user output
  if (contdis->Comm().MyPID() == 0)
    std::cout
        << "<<<<<<<<<<<<<<< Finished extended ghosting of artery discretization >>>>>>>>>>>>>>>\n"
        << std::endl;

  return nearbyelepairs;
}

/*--------------------------------------------------------------------------*
 | create the fully overlapping search discretization      kremheller 03/19 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> POROFLUIDMULTIPHASE::UTILS::CreateSearchDiscretization(
    Teuchos::RCP<DRT::Discretization> artdis)
{
  // we clone a search discretization of the artery discretization on which the search will be
  // performed in a brute force way fully overlapping
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase> discloner =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  Teuchos::RCP<DRT::Discretization> artsearchdis =
      discloner->CreateMatchingDiscretization(artdis, "artsearchdis", false, false, false, false);

  // ghost on all procs.
  DRT::UTILS::GhostDiscretizationOnAllProcs(artsearchdis);
  artsearchdis->FillComplete(false, false, false);

  return artsearchdis;
}

/*----------------------------------------------------------------------*
 | octtree search for nearby elements                  kremheller 10/19 |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> POROFLUIDMULTIPHASE::UTILS::OctTreeSearch(
    Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis,
    Teuchos::RCP<DRT::Discretization> artsearchdis, std::set<int>* elecolset,
    std::set<int>* nodecolset)
{
  // this map will be filled
  std::map<int, std::set<int>> nearbyelepairs;

  // search tree
  Teuchos::RCP<GEO::SearchTree> searchTree = Teuchos::rcp(new GEO::SearchTree(5));

  // nodal positions of 2D/3D-discretization
  std::map<int, LINALG::Matrix<3, 1>> my_positions_cont =
      GetNodalPositions(contdis, contdis->NodeColMap());
  // axis-aligned bounding boxes of all elements of 2D/3D discretization
  std::map<int, LINALG::Matrix<3, 2>> aabb_cont =
      GEO::getCurrentXAABBs(*contdis, my_positions_cont);

  // find the bounding box of the 2D/3D discretization
  const LINALG::Matrix<3, 2> sourceEleBox = GEO::getXAABBofDis(*contdis);
  searchTree->initializeTree(sourceEleBox, *contdis, GEO::TreeType(GEO::OCTTREE));

  // user info and timer
  if (contdis->Comm().MyPID() == 0)
    std::cout << "Starting with OctTree search for coupling ... " << std::endl;
  Epetra_Time timersearch(contdis->Comm());
  // reset timer
  timersearch.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timersearch.WallTime();
  // *********** time measurement ***********

  // nodal positions of artery-discretization (fully overlapping)
  std::map<int, LINALG::Matrix<3, 1>> positions_artery;
  // nodal positions of artery-discretization (row-map format)
  std::map<int, LINALG::Matrix<3, 1>> my_positions_artery =
      GetNodalPositions(artdis, artdis->NodeRowMap());

  // gather
  int procs[contdis->Comm().NumProc()];
  for (int i = 0; i < contdis->Comm().NumProc(); i++) procs[i] = i;
  LINALG::Gather<int, LINALG::Matrix<3, 1>>(
      my_positions_artery, positions_artery, contdis->Comm().NumProc(), &procs[0], contdis->Comm());

  // do the actual search on fully overlapping artery discretization
  for (int iart = 0; iart < artsearchdis->ElementColMap()->NumMyElements(); ++iart)
  {
    const int artelegid = artsearchdis->ElementColMap()->GID(iart);
    DRT::Element* artele = artsearchdis->gElement(artelegid);

    // axis-aligned bounding box of artery
    const LINALG::Matrix<3, 2> aabb_artery = GetAABB(artele, positions_artery);

    // get elements nearby
    std::set<int> closeeles;
    searchTree->searchCollisions(aabb_cont, aabb_artery, 0, closeeles);

    // nearby elements found
    if (closeeles.size() > 0)
    {
      nearbyelepairs[artelegid] = closeeles;

      // add elements and nodes for extended ghosting of artery discretization
      if (not artdis->HaveGlobalElement(artelegid))
      {
        elecolset->insert(artelegid);
        const int* nodeids = artele->NodeIds();
        for (int inode = 0; inode < artele->NumNode(); ++inode) nodecolset->insert(nodeids[inode]);
      }
    }

    // estimate of duration for search
    if (iart == (int)(0.05 * artsearchdis->ElementColMap()->NumMyElements()))
    {
      double mydtsearch = timersearch.WallTime() - dtcpu;
      double maxdtsearch = 0.0;
      contdis->Comm().MaxAll(&mydtsearch, &maxdtsearch, 1);
      if (contdis->Comm().MyPID() == 0)
        std::cout << "Estimated duration: " << 20.0 * (maxdtsearch) << "s" << std::endl;
    }
  }

  // *********** time measurement ***********
  double mydtsearch = timersearch.WallTime() - dtcpu;
  double maxdtsearch = 0.0;
  contdis->Comm().MaxAll(&mydtsearch, &maxdtsearch, 1);
  // *********** time measurement ***********
  if (contdis->Comm().MyPID() == 0) std::cout << "Completed in " << maxdtsearch << "s" << std::endl;

  return nearbyelepairs;
}

/*----------------------------------------------------------------------*
 | get axis-aligned bounding box of element            kremheller 03/19 |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> POROFLUIDMULTIPHASE::UTILS::GetAABB(
    DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& positions)
{
  const LINALG::SerialDenseMatrix xyze_element(GEO::getCurrentNodalPositions(ele, positions));
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  GEO::checkRoughGeoType(ele, xyze_element, eleGeoType);

  return GEO::computeFastXAABB(ele->Shape(), xyze_element, eleGeoType);
}

/*----------------------------------------------------------------------*
 | get nodal positions                                 kremheller 10/19 |
 *----------------------------------------------------------------------*/
std::map<int, LINALG::Matrix<3, 1>> POROFLUIDMULTIPHASE::UTILS::GetNodalPositions(
    Teuchos::RCP<DRT::Discretization> dis, const Epetra_Map* nodemap)
{
  std::map<int, LINALG::Matrix<3, 1>> positions;
  for (int lid = 0; lid < nodemap->NumMyElements(); ++lid)
  {
    const DRT::Node* node = dis->gNode(nodemap->GID(lid));
    LINALG::Matrix<3, 1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    positions[node->Id()] = currpos;
  }
  return positions;
}

/*----------------------------------------------------------------------*
 | get maximum nodal distance                          kremheller 05/18 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::UTILS::GetMaxNodalDistance(
    DRT::Element* ele, Teuchos::RCP<DRT::Discretization> dis)
{
  double maxdist = 0.0;

  for (int inode = 0; inode < ele->NumNode() - 1; inode++)
  {
    // get first node and its position
    int node0_gid = ele->NodeIds()[inode];
    DRT::Node* node0 = dis->gNode(node0_gid);

    static LINALG::Matrix<3, 1> pos0;
    pos0(0) = node0->X()[0];
    pos0(1) = node0->X()[1];
    pos0(2) = node0->X()[2];

    // loop over second node to numnode to compare distances with first node
    for (int jnode = inode + 1; jnode < ele->NumNode(); jnode++)
    {
      int node1_gid = ele->NodeIds()[jnode];
      DRT::Node* node1 = dis->gNode(node1_gid);

      static LINALG::Matrix<3, 1> pos1;
      pos1(0) = node1->X()[0];
      pos1(1) = node1->X()[1];
      pos1(2) = node1->X()[2];

      static LINALG::Matrix<3, 1> dist;
      dist.Update(1.0, pos0, -1.0, pos1, 0.0);

      maxdist = std::max(maxdist, dist.Norm2());
    }
  }

  return maxdist;
}

/*----------------------------------------------------------------------*
 | calculate vector norm                             kremheller 12/17   |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::UTILS::CalculateVectorNorm(
    const enum INPAR::POROFLUIDMULTIPHASE::VectorNorm norm,
    const Teuchos::RCP<const Epetra_Vector> vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::POROFLUIDMULTIPHASE::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_l1_scaled)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm / ((double)vect->GlobalLength());
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}  // CalculateVectorNorm()

/*----------------------------------------------------------------------*
 |                                                    kremheller 03/17  |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::PrintLogo()
{
  std::cout << "This is a Porous Media problem with multiphase flow" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "              +----------+" << std::endl;
  std::cout << "              |  Krebs-  |" << std::endl;
  std::cout << "              |  Modell  |" << std::endl;
  std::cout << "              +----------+" << std::endl;
  std::cout << "              |          |" << std::endl;
  std::cout << "              |          |" << std::endl;
  std::cout << " /\\           |          /\\" << std::endl;
  std::cout << "( /   @ @    (|)        ( /   @ @    ()" << std::endl;
  std::cout << " \\  __| |__  /           \\  __| |__  /" << std::endl;
  std::cout << "  \\/   \"   \\/             \\/   \"   \\/" << std::endl;
  std::cout << " /-|       |-\\           /-|       |-\\" << std::endl;
  std::cout << "/ /-\\     /-\\ \\         / /-\\     /-\\ \\" << std::endl;
  std::cout << " / /-`---'-\\ \\           / /-`---'-\\ \\" << std::endl;
  std::cout << "  /         \\             /         \\" << std::endl;

  return;
}
