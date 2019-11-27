/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/


#include "drt_utils_parallel.H"

#include "../drt_binstrategy/binning_strategy.H"
#include "drt_utils_createdis.H"
#include "drt_node.H"
#include "drt_discret.H"
#include "drt_dserror.H"

#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"
#include "../drt_lib/drt_matchingoctree.H"
#include "../drt_lib/drt_condition.H"

#include <Epetra_IntVector.h>

/*----------------------------------------------------------------------*
 |  Redistribute using BinningStrategy                      rauch 08/16 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::RedistributeDiscretizationsByBinning(
    const std::vector<Teuchos::RCP<DRT::Discretization>>& vector_of_discretizations,
    bool revertextendedghosting)
{
  // safety check
  if (vector_of_discretizations.size() == 0) dserror("No discretizations provided for binning !");

  // get communicator
  const Epetra_Comm& comm = vector_of_discretizations[0]->Comm();

  // redistribute discr. with help of binning strategy
  if (comm.NumProc() > 1)
  {
    if (comm.MyPID() == 0)
    {
      IO::cout(IO::verbose) << "+---------------------------------------------------------------"
                            << IO::endl;
      IO::cout(IO::verbose) << "| Redistribute discretizations using Binning Strategy ...       "
                            << IO::endl;
      for (int dis_num = 0; dis_num < (int)(vector_of_discretizations.size()); dis_num++)
      {
        if (!vector_of_discretizations[dis_num]->Filled())
          dserror("FillComplete(false,false,false) was not called");
        IO::cout(IO::verbose) << "| Redistribute discretization " << std::setw(11)
                              << vector_of_discretizations[dis_num]->Name() << IO::endl;
      }
      IO::cout(IO::verbose) << "+---------------------------------------------------------------"
                            << IO::endl;
    }

    std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
    std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;

    // binning strategy is created and parallel redistribution is performed
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
        Teuchos::rcp(new BINSTRATEGY::BinningStrategy());

    binningstrategy->Init(vector_of_discretizations);

    binningstrategy->DoWeightedPartitioningOfBinsAndExtendGhostingOfDiscretToOneBinLayer(
        vector_of_discretizations, stdelecolmap, stdnodecolmap);

    // revert extended ghosting if requested
    if (revertextendedghosting)
      binningstrategy->RevertExtendedGhosting(
          vector_of_discretizations, stdelecolmap, stdnodecolmap);
  }  // if more than 1 proc
  else
  {
    for (int i = 0; i < static_cast<int>(vector_of_discretizations.size()); ++i)
      vector_of_discretizations[i]->FillComplete();
  }
  return;
}  // DRT::UTILS::RedistributeDiscretizationsByBinning

/*----------------------------------------------------------------------*
 |  Ghost input discr. redundantly on all procs             rauch 09/16 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::GhostDiscretizationOnAllProcs(
    const Teuchos::RCP<DRT::Discretization> distobeghosted)
{
  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(distobeghosted->Comm().Clone());

  if (com->MyPID() == 0)
  {
    IO::cout(IO::verbose)
        << "+-----------------------------------------------------------------------+" << IO::endl;
    IO::cout(IO::verbose) << "|   Ghost discretization " << std::setw(11) << distobeghosted->Name()
                          << " redundantly on all procs ...       |" << IO::endl;
    IO::cout(IO::verbose)
        << "+-----------------------------------------------------------------------+" << IO::endl;
  }

  std::vector<int> allproc(com->NumProc());
  for (int i = 0; i < com->NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = distobeghosted->NodeRowMap();
  std::vector<int> sdata;
  for (int lid = 0; lid < noderowmap->NumMyElements(); ++lid)
  {
    int gid = noderowmap->GID(lid);
    sdata.push_back(gid);
  }

  // gather all master row node gids redundantly in rdata
  std::vector<int> rdata;
  LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), &allproc[0], *com);

  // build new node column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, *com));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap = distobeghosted->ElementRowMap();
  sdata.resize(0);
  for (int i = 0; i < elerowmap->NumMyElements(); ++i)
  {
    int gid = elerowmap->GID(i);
    sdata.push_back(gid);
  }

  // gather all gids of elements redundantly
  rdata.resize(0);
  LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), &allproc[0], *com);

  // build new element column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, *com));
  sdata.clear();
  rdata.clear();
  allproc.clear();

  // redistribute the nodes and elements of the discr. according to the
  // new node / element column layout (i.e. master = full overlap)
  distobeghosted->ExportColumnNodes(*newnodecolmap);
  distobeghosted->ExportColumnElements(*newelecolmap);

  // Safety checks in DEBUG
#ifdef DEBUG
  int nummycolnodes = newnodecolmap->NumMyElements();
  int sizelist[com->NumProc()];
  com->GatherAll(&nummycolnodes, &sizelist[0], 1);
  com->Barrier();
  for (int k = 1; k < com->NumProc(); ++k)
  {
    if (sizelist[k - 1] != nummycolnodes)
      dserror(
          "Since whole dis. is ghosted every processor should have the same number of colnodes.\n"
          "This is not the case."
          "Fix this!");
  }
#endif
  return;
}  // DRT::UTILS::GhostDiscretizationOnAllProcs

/*---------------------------------------------------------------------*
|  Redistribute Nodes Matching Template Discretization     rauch 09/16 |
*----------------------------------------------------------------------*/
void DRT::UTILS::MatchNodalDistributionOfMatchingDiscretizations(
    DRT::Discretization& dis_template, DRT::Discretization& dis_to_redistribute)
{
  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(dis_template.Comm().Clone());
  if (com->NumProc() > 1)
  {
    // print to screen
    if (com->MyPID() == 0)
    {
      IO::cout(IO::verbose)
          << "+---------------------------------------------------------------------------+"
          << IO::endl;
      IO::cout(IO::verbose) << "|   Match nodal distribution of discr. " << std::setw(11)
                            << dis_to_redistribute.Name() << "to discr. " << std::setw(11)
                            << dis_template.Name() << " ... |" << IO::endl;
      IO::cout(IO::verbose)
          << "+---------------------------------------------------------------------------+"
          << IO::endl;
    }

    ////////////////////////////////////////
    // MATCH NODES
    ////////////////////////////////////////
    std::vector<int> redistribute_nodegid_vec(0);
    std::vector<int> redistribute_colnodegid_vec(0);

    // match nodes to be redistributed to template nodes and fill vectors
    // with desired row and col gids for redistribution.
    MatchNodalRowColDistribution(
        dis_template, dis_to_redistribute, redistribute_nodegid_vec, redistribute_colnodegid_vec);

    // construct redistributed node row map
    Teuchos::RCP<Epetra_Map> redistributed_noderowmap = Teuchos::rcp(
        new Epetra_Map(-1, redistribute_nodegid_vec.size(), &redistribute_nodegid_vec[0], 0, *com));

    // construct redistributed node col map
    Teuchos::RCP<Epetra_Map> redistributed_nodecolmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_colnodegid_vec.size(), &redistribute_colnodegid_vec[0], 0, *com));

    // we finally redistribute.
    // FillComplete(...) inside.
    dis_to_redistribute.Redistribute(*redistributed_noderowmap, *redistributed_nodecolmap,
        false,  // assigndegreesoffreedom
        false,  // initelements
        false,  // doboundaryconditions
        true,   // killdofs
        true    // killcond
    );

    // print to screen
    PrintParallelDistribution(dis_to_redistribute);
  }  // if more than one proc
  return;
}  // MatchDistributionOfMatchingDiscretizations


/*---------------------------------------------------------------------*
|  Redistribute Elements Matching Template Discretization  rauch 09/16 |
*----------------------------------------------------------------------*/
void DRT::UTILS::MatchElementDistributionOfMatchingDiscretizations(
    DRT::Discretization& dis_template, DRT::Discretization& dis_to_redistribute)
{
  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(dis_template.Comm().Clone());
  if (com->NumProc() > 1)
  {
    // print to screen
    if (com->MyPID() == 0)
    {
      IO::cout(IO::verbose)
          << "+-----------------------------------------------------------------------------+"
          << IO::endl;
      IO::cout(IO::verbose) << "|   Match element distribution of discr. " << std::setw(11)
                            << dis_to_redistribute.Name() << "to discr. " << std::setw(11)
                            << dis_template.Name() << " ... |" << IO::endl;
      IO::cout(IO::verbose)
          << "+-----------------------------------------------------------------------------+"
          << IO::endl;
    }

    ////////////////////////////////////////
    // MATCH ELEMENTS
    ////////////////////////////////////////
    std::vector<int> redistribute_rowelegid_vec(0);
    std::vector<int> redistribute_colelegid_vec(0);

    // match elements to be redistributed to template elements and fill vectors
    // with desired row and col gids for redistribution.
    MatchElementRowColDistribution(
        dis_template, dis_to_redistribute, redistribute_rowelegid_vec, redistribute_colelegid_vec);

    // construct redistributed element row map
    Teuchos::RCP<Epetra_Map> redistributed_elerowmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_rowelegid_vec.size(), &redistribute_rowelegid_vec[0], 0, *com));

    // construct redistributed element col map
    Teuchos::RCP<Epetra_Map> redistributed_elecolmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_colelegid_vec.size(), &redistribute_colelegid_vec[0], 0, *com));

    ////////////////////////////////////////
    // MATCH NODES
    ////////////////////////////////////////
    std::vector<int> redistribute_nodegid_vec(0);
    std::vector<int> redistribute_colnodegid_vec(0);

    // match nodes to be redistributed to template nodes and fill vectors
    // with desired row and col gids for redistribution.
    MatchNodalRowColDistribution(
        dis_template, dis_to_redistribute, redistribute_nodegid_vec, redistribute_colnodegid_vec);

    // construct redistributed node row map
    Teuchos::RCP<Epetra_Map> redistributed_noderowmap = Teuchos::rcp(
        new Epetra_Map(-1, redistribute_nodegid_vec.size(), &redistribute_nodegid_vec[0], 0, *com));

    // construct redistributed node col map
    Teuchos::RCP<Epetra_Map> redistributed_nodecolmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_colnodegid_vec.size(), &redistribute_colnodegid_vec[0], 0, *com));

    ////////////////////////////////////////
    // REDISTRIBUTE
    ////////////////////////////////////////
    // export the nodes
    dis_to_redistribute.ExportRowNodes(*redistributed_noderowmap, false, false);
    dis_to_redistribute.ExportColumnNodes(*redistributed_nodecolmap, false, false);
    // export the elements
    dis_to_redistribute.ExportRowElements(*redistributed_elerowmap, false, false);
    dis_to_redistribute.ExportColumnElements(*redistributed_elecolmap, false, false);

    ////////////////////////////////////////
    // FINISH
    ////////////////////////////////////////
    int err = dis_to_redistribute.FillComplete(false, false, false);

    if (err) dserror("FillComplete() returned err=%d", err);

    // print to screen
    PrintParallelDistribution(dis_to_redistribute);
  }  // if more than one proc
  return;
}  // DRT::UTILS::MatchElementDistributionOfMatchingDiscretizations


/*---------------------------------------------------------------------*
|  Redistribute Conditioned Elements Matching Template     rauch 10/16 |
*----------------------------------------------------------------------*/
void DRT::UTILS::MatchElementDistributionOfMatchingConditionedElements(
    DRT::Discretization& dis_template, DRT::Discretization& dis_to_redistribute,
    const std::string& condname_template, const std::string& condname_redistribute)
{
  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(dis_template.Comm().Clone());
  if (com->NumProc() > 1)
  {
    // print to screen
    if (com->MyPID() == 0)
    {
      IO::cout(IO::verbose)
          << "+-----------------------------------------------------------------------------+"
          << IO::endl;
      IO::cout(IO::verbose) << "|   Match element distribution of discr. " << std::setw(11)
                            << dis_to_redistribute.Name() << "                          |"
                            << IO::endl;
      IO::cout(IO::verbose) << "|   Condition : " << std::setw(35) << condname_redistribute
                            << "                           |" << IO::endl;
      IO::cout(IO::verbose) << "|   to template discr. " << std::setw(11) << dis_template.Name()
                            << "                                            |" << IO::endl;
      IO::cout(IO::verbose) << "|   Condition : " << std::setw(35) << condname_template
                            << "                           |" << IO::endl;
      IO::cout(IO::verbose)
          << "+-----------------------------------------------------------------------------+"
          << IO::endl;
    }

    // create vectors for element matching
    std::vector<int> my_template_colelegid_vec(0);
    std::vector<int> redistribute_rowelegid_vec(0);

    // create vectors for node matching
    std::vector<int> my_template_nodegid_vec(0);
    std::vector<int> redistribute_rownodegid_vec(0);
    std::vector<int> redistribute_colnodegid_vec(0);

    // geometry iterator
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator geom_it;

    // fill condition discretization by cloning scatra discretization
    Teuchos::RCP<DRT::Discretization> dis_from_template_condition;

    Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase> discreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
    std::vector<std::string> conditions_to_copy(0);
    dis_from_template_condition = discreator->CreateMatchingDiscretizationFromCondition(
        dis_template,        ///< discretization with condition
        condname_template,   ///< name of the condition, by which the derived discretization is
                             ///< identified
        "aux_dis",           ///< name of the new discretization
        "TRANSP",            ///< name/type of the elements to be created
        conditions_to_copy,  ///< list of conditions that will be copied to the new discretization
        -1  ///< coupling id, only elements conditioned with this coupling id are considered
    );

    // get element col map of conditioned template dis
    const Epetra_Map* template_cond_dis_elecolmap = dis_from_template_condition->ElementColMap();
    // get element row map of dis to be redistributed
    const Epetra_Map* redistribute_elerowmap = dis_to_redistribute.ElementRowMap();


    ////////////////////////////////////////
    // MATCH CONDITIONED ELEMENTS
    ////////////////////////////////////////
    // fill element gid vectors
    for (int lid = 0; lid < template_cond_dis_elecolmap->NumMyElements(); lid++)
      my_template_colelegid_vec.push_back(template_cond_dis_elecolmap->GID(lid));

    for (int lid = 0; lid < redistribute_elerowmap->NumMyElements(); lid++)
      redistribute_rowelegid_vec.push_back(redistribute_elerowmap->GID(lid));


    // initialize search tree for matching with template (source,master) elements
    DRT::UTILS::ElementMatchingOctree elementmatchingtree = DRT::UTILS::ElementMatchingOctree();
    elementmatchingtree.Init(*dis_from_template_condition, my_template_colelegid_vec, 150, 1e-06);
    elementmatchingtree.Setup();

    // map that will be filled with matched elements.
    // mapping: redistr. ele gid to (template ele gid, dist.).
    // note: 'FillSlaveToMasterGIDMapping' loops over all
    //        template eles and finds corresponding redistr. eles.
    std::map<int, std::vector<double>> matched_ele_map;
    // match target (slave) elements to source (master) elements using octtree
    elementmatchingtree.FillSlaveToMasterGIDMapping(
        dis_to_redistribute, redistribute_rowelegid_vec, matched_ele_map);

    // declare iterator
    std::map<int, std::vector<double>>::iterator it;

    // now we have a map matching the geometry ids of slave elements
    // to the geometry id of master elements (always starting from 0).
    // for redistribution we need to translate the geometry ids to the
    // actual element gids.
    // fill vectors with row and col gids for new distribution
    std::vector<int> redistribute_colelegid_vec;
    redistribute_rowelegid_vec.clear();
    for (it = matched_ele_map.begin(); it != matched_ele_map.end(); ++it)
    {
      // if this proc owns the template element we also want to own
      // the element of the redistributed discretization.
      // we also want to own all nodes of this element.
      if ((int)(it->second)[2] == 1) redistribute_rowelegid_vec.push_back(it->first);

      redistribute_colelegid_vec.push_back(it->first);
    }

    dis_to_redistribute.Comm().Barrier();
    for (it = matched_ele_map.begin(); it != matched_ele_map.end(); ++it)
    {
      std::cout << "ELEMENT : " << it->first << " ->  ( " << it->second[0] << ", " << it->second[1]
                << ", " << it->second[2] << " )"
                << " on PROC " << dis_to_redistribute.Comm().MyPID()
                << " map size = " << matched_ele_map.size() << std::endl;
    }


    ////////////////////////////////////////
    // ALSO APPEND UNCONDITIONED ELEMENTS
    ////////////////////////////////////////
    // add row elements
    for (int lid = 0; lid < dis_to_redistribute.ElementColMap()->NumMyElements(); lid++)
    {
      bool conditionedele = false;
      DRT::Element* ele =
          dis_to_redistribute.gElement(dis_to_redistribute.ElementColMap()->GID(lid));
      DRT::Node** nodes = ele->Nodes();
      for (int node = 0; node < ele->NumNode(); node++)
      {
        DRT::Condition* nodal_cond = nodes[node]->GetCondition(condname_redistribute);
        if (nodal_cond != NULL)
        {
          conditionedele = true;
          break;
        }
      }  // loop over nodes

      if (not conditionedele)
      {
        // append unconditioned ele id to col gid vec
        redistribute_colelegid_vec.push_back(ele->Id());

        // append unconditioned ele id to row gid vec
        if (ele->Owner() == com->MyPID()) redistribute_rowelegid_vec.push_back(ele->Id());
      }

    }  // loop over col elements


    // construct redistributed element row map
    Teuchos::RCP<Epetra_Map> redistributed_elerowmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_rowelegid_vec.size(), &redistribute_rowelegid_vec[0], 0, *com));

    // construct redistributed element col map
    Teuchos::RCP<Epetra_Map> redistributed_elecolmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_colelegid_vec.size(), &redistribute_colelegid_vec[0], 0, *com));


    ////////////////////////////////////////
    // MATCH CONDITIONED NODES
    ////////////////////////////////////////
    // fill vector with processor local conditioned node gids for template dis
    for (int lid = 0; lid < dis_template.NodeColMap()->NumMyElements(); ++lid)
    {
      if (dis_template.gNode(dis_template.NodeColMap()->GID(lid))
              ->GetCondition(condname_template) != NULL)
        my_template_nodegid_vec.push_back(dis_template.NodeColMap()->GID(lid));
    }

    // fill vec with processor local node gids of dis to be redistributed
    DRT::Condition* redistribute_cond = dis_to_redistribute.GetCondition(condname_redistribute);
    const std::vector<int>* redistribute_cond_nodes = redistribute_cond->Nodes();
    for (int lid = 0; lid < (int)redistribute_cond_nodes->size(); ++lid)
    {
      if (dis_to_redistribute.HaveGlobalNode(redistribute_cond_nodes->at(lid)))
        if (dis_to_redistribute.gNode(redistribute_cond_nodes->at(lid))->Owner() == com->MyPID())
          redistribute_rownodegid_vec.push_back(redistribute_cond_nodes->at(lid));
    }

    // initialize search tree for matching with template (source) nodes
    DRT::UTILS::NodeMatchingOctree nodematchingtree = DRT::UTILS::NodeMatchingOctree();
    nodematchingtree.Init(dis_template, my_template_nodegid_vec, 150, 1e-06);
    nodematchingtree.Setup();

    // map that will be filled with matched nodes.
    // mapping: redistr. node gid to (template node gid, dist.).
    // note: FindMatch loops over all template nodes
    //       and finds corresponding redistr. nodes.
    std::map<int, std::vector<double>> matched_node_map;
    // match target nodes to source nodes using octtree
    nodematchingtree.FillSlaveToMasterGIDMapping(
        dis_to_redistribute, redistribute_rownodegid_vec, matched_node_map);

    // fill vectors with row gids for new distribution
    redistribute_rownodegid_vec.clear();
    // std::vector<int> redistribute_colnodegid_vec;
    for (it = matched_node_map.begin(); it != matched_node_map.end(); ++it)
    {
      // if this proc owns the template node we also want to own
      // the node of the redistributed discretization
      if ((int)(it->second)[2] == 1) redistribute_rownodegid_vec.push_back(it->first);

      redistribute_colnodegid_vec.push_back(it->first);
    }

    dis_to_redistribute.Comm().Barrier();
    for (it = matched_node_map.begin(); it != matched_node_map.end(); ++it)
    {
      std::cout << "NODE : " << it->first << " ->  ( " << it->second[0] << ", " << it->second[1]
                << ", " << it->second[2] << " )"
                << " on PROC " << dis_to_redistribute.Comm().MyPID()
                << " map size = " << matched_node_map.size() << std::endl;
    }


    ////////////////////////////////////////
    // ALSO APPEND UNCONDITIONED  NODES
    ////////////////////////////////////////
    // add row nodes
    for (int lid = 0; lid < dis_to_redistribute.NodeRowMap()->NumMyElements(); lid++)
    {
      DRT::Condition* testcond =
          dis_to_redistribute.gNode(dis_to_redistribute.NodeRowMap()->GID(lid))
              ->GetCondition(condname_redistribute);
      if (testcond == NULL)
        redistribute_rownodegid_vec.push_back(
            dis_to_redistribute.gNode(dis_to_redistribute.NodeRowMap()->GID(lid))->Id());
    }
    // add col nodes
    for (int lid = 0; lid < dis_to_redistribute.NodeColMap()->NumMyElements(); lid++)
    {
      DRT::Condition* testcond =
          dis_to_redistribute.gNode(dis_to_redistribute.NodeColMap()->GID(lid))
              ->GetCondition(condname_redistribute);
      if (testcond == NULL)
        redistribute_colnodegid_vec.push_back(
            dis_to_redistribute.gNode(dis_to_redistribute.NodeColMap()->GID(lid))->Id());
    }

    // construct redistributed node row map
    Teuchos::RCP<Epetra_Map> redistributed_noderowmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_rownodegid_vec.size(), &redistribute_rownodegid_vec[0], 0, *com));

    // construct redistributed node col map
    Teuchos::RCP<Epetra_Map> redistributed_nodecolmap = Teuchos::rcp(new Epetra_Map(
        -1, redistribute_colnodegid_vec.size(), &redistribute_colnodegid_vec[0], 0, *com));


    ////////////////////////////////////////
    // REDISTRIBUTE
    ////////////////////////////////////////
    // export the nodes
    dis_to_redistribute.ExportRowNodes(*redistributed_noderowmap, false, false);
    dis_to_redistribute.ExportColumnNodes(*redistributed_nodecolmap, false, false);
    // export the elements
    dis_to_redistribute.ExportRowElements(*redistributed_elerowmap, false, false);
    dis_to_redistribute.ExportColumnElements(*redistributed_elecolmap, false, false);


    ////////////////////////////////////////
    // FINISH
    ////////////////////////////////////////
    int err = dis_to_redistribute.FillComplete(false, false, false);

    if (err) dserror("FillComplete() returned err=%d", err);

    // print to screen
    PrintParallelDistribution(dis_to_redistribute);

  }  // if more than one proc
  return;
}  // MatchElementDistributionOfMatchingConditionedElements


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PrintParallelDistribution(const DRT::Discretization& dis)
{
  const int numproc = dis.Comm().NumProc();

  if (numproc > 1)
  {
    const int myrank = dis.Comm().MyPID();

    std::vector<int> my_n_nodes(numproc, 0);
    std::vector<int> n_nodes(numproc, 0);
    std::vector<int> my_n_ghostnodes(numproc, 0);
    std::vector<int> n_ghostnodes(numproc, 0);
    std::vector<int> my_n_elements(numproc, 0);
    std::vector<int> n_elements(numproc, 0);
    std::vector<int> my_n_ghostele(numproc, 0);
    std::vector<int> n_ghostele(numproc, 0);

    my_n_nodes[myrank] = dis.NumMyRowNodes();
    my_n_ghostnodes[myrank] = dis.NumMyColNodes() - my_n_nodes[myrank];
    my_n_elements[myrank] = dis.NumMyRowElements();
    my_n_ghostele[myrank] = dis.NumMyColElements() - my_n_elements[myrank];

    dis.Comm().SumAll(&my_n_nodes[0], &n_nodes[0], numproc);
    dis.Comm().SumAll(&my_n_ghostnodes[0], &n_ghostnodes[0], numproc);
    dis.Comm().SumAll(&my_n_elements[0], &n_elements[0], numproc);
    dis.Comm().SumAll(&my_n_ghostele[0], &n_ghostele[0], numproc);

    if (myrank == 0)
    {
      std::cout << std::endl;
      std::cout << "   Discretization: " << dis.Name() << std::endl;
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      printf("   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |\n");
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      for (int npid = 0; npid < numproc; ++npid)
      {
        printf("   | %3d | %13d | %12d | %15d | %14d |\n", npid, n_nodes[npid], n_ghostnodes[npid],
            n_elements[npid], n_ghostele[npid]);
        printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      }
      std::cout << std::endl;
    }
  }
}  // PrintParallelDistribution

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> DRT::UTILS::GetColVersionOfRowVector(
    const Teuchos::RCP<const DRT::Discretization> dis,
    const Teuchos::RCP<const Epetra_Vector> state, const int nds)
{
  // note that this routine has the same functionality as SetState,
  // although here we do not store the new vector anywhere
  // maybe this routine can be used in SetState or become a member function of the discretization
  // class

  if (!dis->HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = dis->DofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  // if it's already in column map just set a reference
  // This is a rought test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap)) return state;
  // if it's not in column map export and allocate
  else
  {
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*colmap, false);
    LINALG::Export(*state, *tmp);
    return tmp;
  }
}  // GetColVersionOfRowVector

/*----------------------------------------------------------------------*
 |(private)                                                   tk 06/10  |
 |recompute nodecolmap of standard discretization to include all        |
 |nodes as of subdicretization                                          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::ComputeNodeColMap(
    const Teuchos::RCP<DRT::Discretization>
        sourcedis,  ///< standard discretization we want to redistribute
    const Teuchos::RCP<DRT::Discretization> subdis  ///< subdiscretization prescribing ghosting
)
{
  const Epetra_Map* oldcolnodemap = sourcedis->NodeColMap();

  std::vector<int> mycolnodes(oldcolnodemap->NumMyElements());
  oldcolnodemap->MyGlobalElements(&mycolnodes[0]);
  for (int inode = 0; inode != subdis->NumMyColNodes(); ++inode)
  {
    const DRT::Node* newnode = subdis->lColNode(inode);
    const int gid = newnode->Id();
    if (!(sourcedis->HaveGlobalNode(gid)))
    {
      mycolnodes.push_back(gid);
    }
  }

  // now reconstruct the extended colmap
  Teuchos::RCP<Epetra_Map> newcolnodemap =
      Teuchos::rcp(new Epetra_Map(-1, mycolnodes.size(), &mycolnodes[0], 0, sourcedis->Comm()));
  return newcolnodemap;
}  // DRT::UTILS::ComputeNodeColMap

/*----------------------------------------------------------------------*
 *                                                          rauch 10/16 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::MatchElementRowColDistribution(const DRT::Discretization& dis_template,
    const DRT::Discretization& dis_to_redistribute, std::vector<int>& row_id_vec_to_fill,
    std::vector<int>& col_id_vec_to_fill)
{
  // preliminary work
  const Epetra_Map* redistribute_elerowmap = dis_to_redistribute.ElementRowMap();
  const Epetra_Map* template_elecolmap = dis_template.ElementColMap();
  std::vector<int> my_template_elegid_vec(template_elecolmap->NumMyElements());
  std::vector<int> my_redistribute_elegid_vec(0);

  // fill vector with processor local ele gids for template dis
  for (int lid = 0; lid < template_elecolmap->NumMyElements(); ++lid)
    my_template_elegid_vec[lid] = template_elecolmap->GID(lid);

  // fill vec with processor local ele gids of dis to be redistributed
  for (int lid = 0; lid < redistribute_elerowmap->NumMyElements(); ++lid)
    my_redistribute_elegid_vec.push_back(redistribute_elerowmap->GID(lid));

  // initialize search tree for matching with template (source,master) elements
  DRT::UTILS::ElementMatchingOctree elementmatchingtree = DRT::UTILS::ElementMatchingOctree();
  elementmatchingtree.Init(dis_template, my_template_elegid_vec, 150, 1e-07);
  elementmatchingtree.Setup();

  // map that will be filled with matched elements.
  // mapping: redistr. ele gid to (template ele gid, dist.).
  // note: 'FillSlaveToMasterGIDMapping' loops over all
  //        template eles and finds corresponding redistr. eles.
  std::map<int, std::vector<double>> matched_ele_map;
  // match target (slave) nodes to source (master) nodes using octtree
  elementmatchingtree.FillSlaveToMasterGIDMapping(
      dis_to_redistribute, my_redistribute_elegid_vec, matched_ele_map);

  // declare iterator
  std::map<int, std::vector<double>>::iterator it;

  // fill vectors with row and col gids for new distribution
  for (it = matched_ele_map.begin(); it != matched_ele_map.end(); ++it)
  {
    // if this proc owns the template element we also want to own
    // the element of the redistributed discretization.
    // we also want to own all nodes of this element.
    if ((int)(it->second)[2] == 1) row_id_vec_to_fill.push_back(it->first);

    col_id_vec_to_fill.push_back(it->first);
  }
  return;
}  // DRT::UTILS::MatchElementRowColDistribution

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::MatchNodalRowColDistribution(const DRT::Discretization& dis_template,
    const DRT::Discretization& dis_to_redistribute, std::vector<int>& row_id_vec_to_fill,
    std::vector<int>& col_id_vec_to_fill)
{
  // temp sets
  std::set<int> temprowset;
  std::set<int> tempcolset;

  for (int i = 0; i < (int)row_id_vec_to_fill.size(); ++i)
  {
    temprowset.insert(row_id_vec_to_fill[i]);
  }

  for (int i = 0; i < (int)col_id_vec_to_fill.size(); ++i)
  {
    tempcolset.insert(col_id_vec_to_fill[i]);
  }

  // preliminary work
  const Epetra_Map* redistribute_noderowmap = dis_to_redistribute.NodeRowMap();
  const Epetra_Map* template_nodecolmap = dis_template.NodeColMap();
  std::vector<int> my_template_nodegid_vec(template_nodecolmap->NumMyElements());
  std::vector<int> my_redistribute_nodegid_vec(0);

  // fill vector with processor local node gids for template dis
  for (int lid = 0; lid < template_nodecolmap->NumMyElements(); ++lid)
    my_template_nodegid_vec[lid] = template_nodecolmap->GID(lid);

  // fill vec with processor local node gids of dis to be redistributed
  for (int lid = 0; lid < redistribute_noderowmap->NumMyElements(); ++lid)
    my_redistribute_nodegid_vec.push_back(redistribute_noderowmap->GID(lid));

  // initialize search tree for matching with template (source) nodes
  DRT::UTILS::NodeMatchingOctree nodematchingtree = DRT::UTILS::NodeMatchingOctree();
  nodematchingtree.Init(dis_template, my_template_nodegid_vec, 150, 1e-07);
  nodematchingtree.Setup();

  // map that will be filled with matched nodes.
  // mapping: redistr. node gid to (template node gid, dist.).
  // note: FindMatch loops over all template nodes
  //       and finds corresponding redistr. nodes.
  std::map<int, std::vector<double>> matched_node_map;
  // match target nodes to source nodes using octtree
  nodematchingtree.FillSlaveToMasterGIDMapping(
      dis_to_redistribute, my_redistribute_nodegid_vec, matched_node_map);

  // declare iterator
  std::map<int, std::vector<double>>::iterator it;

  // fill vectors with row gids for new distribution
  // std::vector<int> redistribute_colnodegid_vec;
  for (it = matched_node_map.begin(); it != matched_node_map.end(); ++it)
  {
    // if this proc owns the template node we also want to own
    // the node of the redistributed discretization
    if ((int)(it->second)[2] == 1) temprowset.insert(it->first);

    tempcolset.insert(it->first);
  }

  // assign temporary sets to vectors
  row_id_vec_to_fill.clear();
  row_id_vec_to_fill.reserve(temprowset.size());
  row_id_vec_to_fill.assign(temprowset.begin(), temprowset.end());

  col_id_vec_to_fill.clear();
  col_id_vec_to_fill.reserve(tempcolset.size());
  col_id_vec_to_fill.assign(tempcolset.begin(), tempcolset.end());
  return;
}  // DRT::UTILS::MatchNodalRowColDistribution

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::RedistributeInAccordanceWithReference(
    const Epetra_Map& ref_red_map, const Epetra_Map& unred_map)
{
  Teuchos::RCP<Epetra_Map> red_map = Teuchos::null;
  RedistributeInAccordanceWithReference(ref_red_map, unred_map, red_map);
  return red_map;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::RedistributeInAccordanceWithReference(
    const Epetra_Map& ref_red_map, const Epetra_Map& unred_map, Teuchos::RCP<Epetra_Map>& red_map)
{
  const Epetra_Comm& comm = unred_map.Comm();

  Epetra_IntVector source(unred_map);
  {
    const int num_myentries = unred_map.NumMyElements();
    const int* mygids = unred_map.MyGlobalElements();
    int* vals = source.Values();
    std::copy(mygids, mygids + num_myentries, vals);
  }

  Epetra_Import importer(ref_red_map, unred_map);

  Epetra_IntVector target(ref_red_map);
  target.PutValue(-1);

  const int err = target.Import(source, importer, Insert);
  if (err) dserror("Import failed with error %d!", err);

  std::vector<int> my_red_gids;
  {
    const int num_myentries = ref_red_map.NumMyElements();
    const int* my_received_gids = target.Values();
    my_red_gids.reserve(num_myentries);

    for (int i = 0; i < num_myentries; ++i)
    {
      if (my_received_gids[i] != -1) my_red_gids.push_back(my_received_gids[i]);
    }
  }

  // create new map
  red_map = Teuchos::rcp(new Epetra_Map(-1, my_red_gids.size(), my_red_gids.data(), 0, comm));

  //  red_map->Print( std::cout );
}
