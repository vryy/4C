/*!----------------------------------------------------------------------
\file drt_utils_parallel.cpp

\brief A collection of helper methods for namespace DRT

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/


#include "drt_utils_parallel.H"
#include "drt_node.H"
#include "drt_discret.H"
#include "drt_dserror.H"

#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"
#include "../drt_particle/binning_strategy.H"
#include "../drt_lib/drt_nodematchingoctree.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    )
{
  // redistribute nodes to column (ghost) map
  dis.ExportColumnNodes(nodecolmap);

  Teuchos::RCP< Epetra_Map > elerowmap;
  Teuchos::RCP< Epetra_Map > elecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  dis.BuildElementRowColumn(noderowmap, nodecolmap, elerowmap, elecolmap);

  // we can now export elements to resonable row element distribution
  dis.ExportRowElements(*elerowmap);

  // export to the column map / create ghosting of elements
  dis.ExportColumnElements(*elecolmap);

}

/*----------------------------------------------------------------------*
 |  Redistribute using BinningStrategy                      rauch 08/16 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::RedistributeDiscretizationsByBinning(
    const std::vector<Teuchos::RCP<DRT::Discretization> >& vector_of_discretizations,
    bool revertextendedghosting
    )
{
  // safety check
  if (vector_of_discretizations.size()==0)
    dserror("No discretizations provided for binning !");

  // get communicator
  const Epetra_Comm& comm = vector_of_discretizations[0]->Comm();

  // redistribute discr. with help of binning strategy
  if(comm.NumProc()>1)
  {
    if(comm.MyPID() == 0)
    {
      IO::cout(IO::verbose)<<"+---------------------------------------------------------------"<<IO::endl;
      IO::cout(IO::verbose)<<"| Redistribute discretizations using Binning Strategy ...       "<<IO::endl;
      for(int dis_num=0; dis_num < (int)(vector_of_discretizations.size()); dis_num++)
      {
        if (!vector_of_discretizations[dis_num]->Filled())
          dserror("FillComplete(false,false,false) was not called");
        IO::cout(IO::verbose)<<"| Redistribute discretization "<<std::setw(11)<<vector_of_discretizations[dis_num]->Name()<<IO::endl;
      }
      IO::cout(IO::verbose)<<"+---------------------------------------------------------------"<<IO::endl;
    }

    std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
    std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

    // binning strategy is created and parallel redistribution is performed
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
        Teuchos::rcp(new BINSTRATEGY::BinningStrategy(
            vector_of_discretizations,
            stdelecolmap,
            stdnodecolmap) );

    // revert extended ghosting if requested
    if(revertextendedghosting)
        binningstrategy->RevertExtendedGhosting(
            vector_of_discretizations,
            stdelecolmap,
            stdnodecolmap);

  }

  return;
}

/*----------------------------------------------------------------------*
 |  Ghost input discr. redundantly on all procs             rauch 09/16 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::GhostDiscretizationOnAllProcs(
    const Teuchos::RCP<DRT::Discretization> distobeghosted
    )
{
  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( distobeghosted->Comm().Clone());

  if(com->MyPID() == 0)
  {
    IO::cout(IO::verbose)<<"+-----------------------------------------------------------------------+"<<IO::endl;
    IO::cout(IO::verbose)<<"|   Ghost discretization "<<std::setw(11)<<distobeghosted->Name()<<
        " redundantly on all procs ...       |"<<IO::endl;
    IO::cout(IO::verbose)<<"+-----------------------------------------------------------------------+"<<IO::endl;
  }

  std::vector<int> allproc(com->NumProc());
  for (int i=0; i<com->NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = distobeghosted->NodeRowMap();
  std::vector<int> sdata;
  for (int lid=0; lid<noderowmap->NumMyElements(); ++lid)
  {
    int gid = noderowmap->GID(lid);
    sdata.push_back(gid);
#ifdef DEBUG
    DRT::Node* node = distobeghosted->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
#endif
  }

  // gather all master row node gids redundantly in rdata
  std::vector<int> rdata;
  LINALG::Gather<int>(sdata,rdata,(int)allproc.size(),&allproc[0],*com);

  // build new node column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap =
      Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,*com));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap  = distobeghosted->ElementRowMap();
  sdata.resize(0);
  for (int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    int gid = elerowmap->GID(i);
    sdata.push_back(gid);
#ifdef DEBUG
    DRT::Element* ele = distobeghosted->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %",gid);
#endif
  }

  // gather all gids of elements redundantly
  rdata.resize(0);
  LINALG::Gather<int>(sdata,rdata,(int)allproc.size(),&allproc[0],*com);

  // build new element column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap =
      Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,*com));
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
  com->GatherAll(&nummycolnodes,&sizelist[0],1);
  com->Barrier();
  for(int k=1;k<com->NumProc();++k)
  {
    if(sizelist[k-1]!=nummycolnodes)
      dserror("Since whole dis. is ghosted every processor should have the same number of colnodes.\n"
              "This is not the case."
              "Fix this!");
  }
#endif

  return;
} // end CreateGhosting

/*---------------------------------------------------------------------*
|  Redistribute Nodes Matching Template Discretization     rauch 09/16 |
*----------------------------------------------------------------------*/
void DRT::UTILS::MatchDistributionOfMatchingDiscretizations(
    DRT::Discretization&     dis_template,
    DRT::Discretization&     dis_to_redistribute
    )
{
  // clone communicator of target discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( dis_template.Comm().Clone());
  if(com->NumProc()>1)
  {
    // print to screen
    if(com->MyPID() == 0)
    {
      IO::cout(IO::verbose)<<"+-----------------------------------------------------------------------+"<<IO::endl;
      IO::cout(IO::verbose)<<"|   Match distribution of discr. "<<std::setw(11)<<
          dis_to_redistribute.Name()<<
          "to discr. "<<std::setw(11)<<dis_template.Name()<<" ...   |"<<IO::endl;
      IO::cout(IO::verbose)<<"+-----------------------------------------------------------------------+"<<IO::endl;
    }

    // get node map of template dis
    const Epetra_Map* template_nodecolmap =
        dis_template.NodeColMap();
    // get node map of dis to be redistributed
    const Epetra_Map* redistribute_noderowmap =
        dis_to_redistribute.NodeRowMap();

    // fill vector with processor local node gids for template dis
    std::vector<int> my_template_nodegid_vec(template_nodecolmap->NumMyElements());
    for(int lid=0; lid<template_nodecolmap->NumMyElements();++lid)
      my_template_nodegid_vec[lid] = template_nodecolmap->GID(lid);

    // fill vec with processor local node gids of dis to be redistributed
    std::vector<int> redistribute_nodegid_vec(redistribute_noderowmap->NumMyElements());
    for(int lid=0; lid<redistribute_noderowmap->NumMyElements();++lid)
      redistribute_nodegid_vec[lid] = redistribute_noderowmap->GID(lid);

    // initialize search tree for matching with template (source) nodes
    DRT::UTILS::NodeMatchingOctree tree(dis_template, my_template_nodegid_vec,150,1e-07);

    // map that will be filled with matched nodes.
    // mapping: redistr. node gid to (template node gid, dist.).
    // note: FindMatch loops over all template nodes
    //       and finds corresponding redistr. nodes.
    std::map<int,std::vector<double> > matched_node_map;
    // match target nodes to source nodes using octtree
    tree.FillSlaveToMasterGIDMapping(
        dis_to_redistribute,
        redistribute_nodegid_vec,
        matched_node_map );

    std::map<int,std::vector<double> >::iterator it;

//    // DEBUG
//    std::cout<<matched_node_map.size()<<std::endl;
//    std::cout <<" template gids" << " => " << "redistribute gids "<<'\n';
//    for (it=matched_node_map.begin(); it!=matched_node_map.end(); ++it)
//      std::cout <<"PROC "<<com->MyPID()<<" : " <<it->first << " => "
//      << (it->second)[0] << " , "<<(it->second)[1]<<" , "<<(it->second)[2]<< '\n';
//    // END DEBUG

    // fill vectors with row gids for new distribution
    redistribute_nodegid_vec.clear();
    std::vector<int> redistribute_colnodegid_vec;
    for (it=matched_node_map.begin(); it!=matched_node_map.end(); ++it)
    {
      // if this proc owns the template node we also want to own
      // the node of the redistributed discretization
      if((int)(it->second)[2]==1)
        redistribute_nodegid_vec.push_back(it->first);

      redistribute_colnodegid_vec.push_back(it->first);
    }

    // construct redistributed node row map
    Teuchos::RCP<Epetra_Map> redistributed_noderowmap =
        Teuchos::rcp(new Epetra_Map(-1, redistribute_nodegid_vec.size(), &redistribute_nodegid_vec[0], 0, *com));

    // construct redistributed node col map
    Teuchos::RCP<Epetra_Map> redistributed_nodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, redistribute_colnodegid_vec.size(), &redistribute_colnodegid_vec[0], 0, *com));

    // we finally redistribute.
    // FillComplete(...) inside.
    dis_to_redistribute.Redistribute(
        *redistributed_noderowmap,
        *redistributed_nodecolmap,
        false, // assigndegreesoffreedom
        false, // initelements
        false, // doboundaryconditions
        true,  // killdofs
        true   // killcond
    );

    // print to screen
    PrintParallelDistribution(dis_to_redistribute);
  } // if more than one proc
  return;
} // MatchDistributionOfMatchingDiscretizations


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PrintParallelDistribution(const DRT::Discretization& dis)
{
  const int numproc=dis.Comm().NumProc();

  if(numproc>1)
  {
    const int myrank=dis.Comm().MyPID();

    std::vector<int> my_n_nodes     (numproc,0);
    std::vector<int>    n_nodes     (numproc,0);
    std::vector<int> my_n_ghostnodes(numproc,0);
    std::vector<int>    n_ghostnodes(numproc,0);
    std::vector<int> my_n_elements  (numproc,0);
    std::vector<int>    n_elements  (numproc,0);
    std::vector<int> my_n_ghostele  (numproc,0);
    std::vector<int>    n_ghostele  (numproc,0);

    my_n_nodes     [myrank]=dis.NumMyRowNodes();
    my_n_ghostnodes[myrank]=dis.NumMyColNodes()-my_n_nodes[myrank];
    my_n_elements  [myrank]=dis.NumMyRowElements();
    my_n_ghostele  [myrank]=dis.NumMyColElements()-my_n_elements[myrank];

    dis.Comm().SumAll(&my_n_nodes     [0],&n_nodes     [0],numproc);
    dis.Comm().SumAll(&my_n_ghostnodes[0],&n_ghostnodes[0],numproc);
    dis.Comm().SumAll(&my_n_elements  [0],&n_elements  [0],numproc);
    dis.Comm().SumAll(&my_n_ghostele  [0],&n_ghostele  [0],numproc);

    if(myrank==0)
    {
      std::cout << std::endl;
      std::cout <<"   Discretization: " << dis.Name() << std::endl;
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      printf("   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |\n");
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      for(int npid=0;npid<numproc;++npid)
      {
        printf("   | %3d | %13d | %12d | %15d | %14d |\n",npid,n_nodes[npid],n_ghostnodes[npid],n_elements[npid],n_ghostele[npid]);
        printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      }
      std::cout << std::endl;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> DRT::UTILS::GetColVersionOfRowVector(
    const Teuchos::RCP<const DRT::Discretization> dis,
    const Teuchos::RCP<const Epetra_Vector> state,
    const int nds)
{
  // note that this routine has the same functionality as SetState,
  // although here we do not store the new vector anywhere
  // maybe this routine can be used in SetState or become a member function of the discretization class

  if (!dis->HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = dis->DofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  // if it's already in column map just set a reference
  // This is a rought test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
    return state;
  // if it's not in column map export and allocate
  else
  {
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state,*tmp);
    return tmp;
  }
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 06/10  |
 |recompute nodecolmap of standard discretization to include all        |
 |nodes as of subdicretization                                          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::ComputeNodeColMap(
         const Teuchos::RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
         const Teuchos::RCP<DRT::Discretization> subdis      ///< subdiscretization prescribing ghosting
         )
{
  const Epetra_Map* oldcolnodemap = sourcedis->NodeColMap();

  std::vector<int> mycolnodes(oldcolnodemap->NumMyElements());
  oldcolnodemap->MyGlobalElements (&mycolnodes[0]);
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
  Teuchos::RCP<Epetra_Map> newcolnodemap = Teuchos::rcp(new Epetra_Map(-1,
                                     mycolnodes.size(),
                                     &mycolnodes[0],
                                     0,
                                     sourcedis->Comm()));
  return newcolnodemap;
}
