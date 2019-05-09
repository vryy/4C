/*----------------------------------------------------------------------*/
/*!
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
      BruteForceSearch(contdis, artdis, artsearchdis, &elecolset, &nodecolset);

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
 | brute force search for near elements                kremheller 05/18 |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> POROFLUIDMULTIPHASE::UTILS::BruteForceSearch(
    Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis,
    Teuchos::RCP<DRT::Discretization> artsearchdis, std::set<int>* elecolset,
    std::set<int>* nodecolset)
{
  // user info and timer
  if (contdis->Comm().MyPID() == 0)
    std::cout << "Starting with brute force search for coupling ... " << std::endl;
  Epetra_Time timersearch(contdis->Comm());
  // reset timer
  timersearch.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timersearch.WallTime();
  // *********** time measurement ***********

  // this set will be filled
  std::map<int, std::set<int>> nearbyelepairs;

  // scheme is artery-dis: fully-overlapping elecolmap --> contdis: elecolmap
  for (int iart = 0; iart < artsearchdis->ElementColMap()->NumMyElements(); ++iart)
  {
    const int artelegid = artsearchdis->ElementColMap()->GID(iart);

    DRT::Element* artele = artsearchdis->gElement(artelegid);

    // AABB of artery
    std::vector<double> aabb_artery = GetAABB(artele);

    for (int icont = 0; icont < contdis->ElementColMap()->NumMyElements(); ++icont)
    {
      const int contelegid = contdis->ElementColMap()->GID(icont);
      DRT::Element* contele = contdis->gElement(contelegid);

      // AABB of 2D/3D-discretization
      std::vector<double> aabb_cont = GetAABB(contele);

      // check if bounding boxes intersect
      if (!(aabb_artery[0] >= aabb_cont[3] || aabb_artery[3] <= aabb_cont[0] ||
              aabb_artery[1] >= aabb_cont[4] || aabb_artery[4] <= aabb_cont[1] ||
              aabb_artery[2] >= aabb_cont[5] || aabb_artery[5] <= aabb_cont[2]))
      {
        nearbyelepairs[artelegid].insert(contelegid);
        if (not artdis->HaveGlobalElement(artelegid))
        {
          elecolset->insert(artelegid);
          const int* nodeids = artele->NodeIds();
          for (int inode = 0; inode < artele->NumNode(); ++inode)
            nodecolset->insert(nodeids[inode]);
        }
      }
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
std::vector<double> POROFLUIDMULTIPHASE::UTILS::GetAABB(DRT::Element* ele)
{
  // safety factor for brute force search
  const double bruteforcesafety = 1.1;
  // minimum length of AABB (in one coordinate direction) for brute force search
  const double bruteforceaabbminsize = 2.0e-10;

  // format is [xmin, ymin, zmin, xmax, ymax, zmax]
  std::vector<double> aabb(6);

  DRT::Node** nodes = ele->Nodes();

  // initialize with first node
  aabb[0] = nodes[0]->X()[0];
  aabb[1] = nodes[0]->X()[1];
  aabb[2] = nodes[0]->X()[2];
  aabb[3] = nodes[0]->X()[0];
  aabb[4] = nodes[0]->X()[1];
  aabb[5] = nodes[0]->X()[2];

  // loop over remaining nodes
  for (int inode = 1; inode < ele->NumNode(); inode++)
  {
    DRT::Node* node = nodes[inode];
    aabb[0] = std::min(aabb[0], node->X()[0]);
    aabb[1] = std::min(aabb[1], node->X()[1]);
    aabb[2] = std::min(aabb[2], node->X()[2]);

    aabb[3] = std::max(aabb[3], node->X()[0]);
    aabb[4] = std::max(aabb[4], node->X()[1]);
    aabb[5] = std::max(aabb[5], node->X()[2]);
  }

  // for safety: we increase the length of AABB by the factor bruteforcesafety (in x-, y- and
  // z-direction),
  // i.e. we add  bruteforcesafety/2.0*original length to the max values
  // and subtract bruteforcesafety/2.0*original length from the min values
  aabb[0] = (aabb[0] + aabb[3]) / 2.0 - (aabb[3] - aabb[0]) * bruteforcesafety / 2.0;
  aabb[1] = (aabb[1] + aabb[4]) / 2.0 - (aabb[4] - aabb[1]) * bruteforcesafety / 2.0;
  aabb[2] = (aabb[2] + aabb[5]) / 2.0 - (aabb[5] - aabb[2]) * bruteforcesafety / 2.0;
  aabb[3] = (aabb[0] + aabb[3]) / 2.0 + (aabb[3] - aabb[0]) * bruteforcesafety / 2.0;
  aabb[4] = (aabb[1] + aabb[4]) / 2.0 + (aabb[4] - aabb[1]) * bruteforcesafety / 2.0;
  aabb[5] = (aabb[2] + aabb[5]) / 2.0 + (aabb[5] - aabb[2]) * bruteforcesafety / 2.0;

  // take care of AABBs, that lie in surface orthogonal to one of the coordinate axes
  if (aabb[3] - aabb[0] < bruteforceaabbminsize)
  {
    aabb[3] += bruteforceaabbminsize / 2.0;
    aabb[0] -= bruteforceaabbminsize / 2.0;
  }
  if (aabb[4] - aabb[1] < bruteforceaabbminsize)
  {
    aabb[4] += bruteforceaabbminsize / 2.0;
    aabb[1] -= bruteforceaabbminsize / 2.0;
  }
  if (aabb[5] - aabb[2] < bruteforceaabbminsize)
  {
    aabb[5] += bruteforceaabbminsize / 2.0;
    aabb[2] -= bruteforceaabbminsize / 2.0;
  }

  return aabb;
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
