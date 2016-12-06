/*----------------------------------------------------------------------*/
/*!
\file xfem_discretization_utils.cpp

\brief Basic discretization-related tools used in XFEM routines

\level 1

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include "xfem_discretization_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_lib/drt_dofset_fixed_size.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io_gmsh.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFEM::UTILS::PrintDiscretizationToStream(
  Teuchos::RCP<DRT::Discretization> dis,
  const std::string& disname,
  bool elements,
  bool elecol,
  bool nodes,
  bool nodecol,
  bool faces,
  bool facecol,
  std::ostream& s,
  std::map<int, LINALG::Matrix<3,1> >* curr_pos
)
{
  if(elements)
  {
    // draw bg elements with associated gid
    s << "View \" " << disname;
    if(elecol)
    {
      s << " col e->Id() \" {\n";
      for (int i=0; i<dis->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = dis->lColElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    else
    {
      s << " row e->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = dis->lRowElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    s << "};\n";
  }

  if(nodes)
  {
    s << "View \" " << disname;
    if(nodecol)
    {
      s << " col n->Id() \" {\n";
      for (int i=0; i<dis->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lColNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    else
    {
      s << " row n->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lRowNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    s << "};\n";
  }

  if(faces)
  {
    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(dis, true);
    if (xdis == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    s << "View \" " << disname;

    if( xdis->FilledExtension() == true )     // faces output
    {

      if(facecol)
      {
        s << " col f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyColFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lColFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      else
      {
        s << " row f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyRowFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lRowFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      s << "};\n";
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SetupXFEMDiscretization(
  const Teuchos::ParameterList&     xgen_params,
  Teuchos::RCP<DRT::Discretization> dis,
  int numdof)
{
  Teuchos::RCP<DRT::DiscretizationXFEM> xdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(dis, false);
  //
  if (xdis == Teuchos::null)
  {
    dserror("No XFEM discretization for XFEM problem available!");

    //REMARK: standard fluid could also step into this routine, as a special case! (remove dserror)
    if (!dis->Filled())
      dis->FillComplete();

    return;
  }

  if (!xdis->Filled())
    xdis->FillComplete();

  // now we can reserve dofs for xfem discretization
  int numglobalnodes = xdis->NumGlobalNodes();
  int maxNumMyReservedDofsperNode = (xgen_params.get<int>("MAX_NUM_DOFSETS"))*numdof;
    Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofsperNode,numglobalnodes));
  xdis->ReplaceDofSet(0, maxdofset, true ); //fluid dofset has nds = 0
  xdis->InitialFillComplete();

  // print all dofsets
  xdis->GetDofSetProxy()->PrintAllDofsets(xdis->Comm());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SetupXFEMDiscretization(
  const Teuchos::ParameterList&      xgen_params,
  Teuchos::RCP<DRT::Discretization>  dis,
  Teuchos::RCP<DRT::Discretization>  embedded_dis,
  int numdof)
{
  if (!embedded_dis->Filled())
    embedded_dis->FillComplete();

  Teuchos::RCP<DRT::DiscretizationXFEM> xdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(dis, true);
  if (!xdis->Filled())
    xdis->FillComplete();

  // get fluid mesh conditions: hereby we specify standalone fluid discretizations
  std::vector<DRT::Condition*> conditions;
  xdis->GetCondition("FluidMesh",conditions);

  std::vector<std::string> conditions_to_copy;
  xdis->GetConditionNames(conditions_to_copy);

  SplitDiscretizationByCondition(
      xdis,
      embedded_dis,
      conditions,
      "FLUID",
      conditions_to_copy);

  SetupXFEMDiscretization(xgen_params,xdis,numdof);

  DRT::UTILS::PrintParallelDistribution(*dis);
  DRT::UTILS::PrintParallelDistribution(*embedded_dis);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SplitDiscretizationByCondition(
  Teuchos::RCP<DRT::Discretization>  sourcedis,
  Teuchos::RCP<DRT::Discretization>  targetdis,
  std::vector<DRT::Condition*>&      conditions,
  const std::string&                 element_name,
  const std::vector<std::string>&    conditions_to_copy
)
{
  if (!sourcedis->Filled())
    dserror("sourcedis is not filled");
  const int myrank = targetdis->Comm().MyPID();

  // row node map (id -> pointer)
  std::map<int, DRT::Node*> sourcenodes;

  // column node map
  std::map<int, DRT::Node*> sourcegnodes;

  // element map
  std::map<int, Teuchos::RCP<DRT::Element> > sourceelements;

  const int numothernoderow = sourcedis->NumMyRowNodes();
  const int numothernodecol = sourcedis->NumMyColNodes();

  // find conditioned nodes (owned and ghosted) and elements
  DRT::UTILS::FindConditionObjects(*sourcedis, sourcenodes, sourcegnodes, sourceelements, conditions);

  // add the conditioned elements
  for (std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator sourceele_iter = sourceelements.begin();
      sourceele_iter != sourceelements.end();
       ++sourceele_iter)
  {
    if (sourceele_iter->second->Owner() == myrank)
    {
      targetdis->AddElement(sourceele_iter->second);
    }
  }

  // row/col vectors of conditioned node ids
  std::vector<int> condnoderowvec;
  condnoderowvec.reserve(sourcenodes.size());
  std::vector<int> condnodecolvec;
  condnodecolvec.reserve(sourcegnodes.size());

  // add conditioned nodes and fill the id vectors
  for (std::map<int, DRT::Node*>::const_iterator sourcegnode_iter = sourcegnodes.begin();
       sourcegnode_iter != sourcegnodes.end(); ++ sourcegnode_iter)
  {
    const int nid = sourcegnode_iter->first;
    if (sourcegnode_iter->second->Owner() == myrank)
    {
      Teuchos::RCP<DRT::Node> sourcegnode = Teuchos::rcp(new DRT::Node(nid, sourcegnode_iter->second->X(), myrank));
      targetdis->AddNode(sourcegnode);
      condnoderowvec.push_back(nid);
    }
    condnodecolvec.push_back(nid);
  }

  // row/col vectors of non-conditioned node ids
  std::vector<int> othernoderowvec;
  othernoderowvec.reserve(numothernoderow-condnoderowvec.size());
  std::vector<int> othernodecolvec;
  othernodecolvec.reserve(numothernodecol-condnodecolvec.size());

  // determine non-conditioned nodes
  for (int nlid = 0; nlid < sourcedis->NumMyColNodes(); ++nlid)
  {
    const int nid = sourcedis->lColNode(nlid)->Id();
    bool keep = true;
    for (size_t i= 0; i < condnodecolvec.size(); ++i)
    {
      if (condnodecolvec[i] == nid)
      {
        keep = false;
        break;
      }
    }

    if (keep)
      othernodecolvec.push_back(nid);

    if (sourcedis->NodeRowMap()->LID(nid) > -1 && keep)
    {
      othernoderowvec.push_back(nid);
    }
  }

  // delete conditioned nodes and elements from source discretizatio
  for (std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator sourceele_iter = sourceelements.begin();
      sourceele_iter != sourceelements.end();
       ++sourceele_iter)
    sourcedis->DeleteElement(sourceele_iter->first);

  for (size_t i= 0; i < condnodecolvec.size(); ++i)
    sourcedis->DeleteNode(condnodecolvec[i]);


  // copy selected conditions to the new discretization
  for (std::vector<std::string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end();
       ++conditername)
  {
    std::vector<DRT::Condition*> conds;
    sourcedis->GetCondition(*conditername, conds);
    for (unsigned i=0; i<conds.size(); ++i)
    {
      targetdis->SetCondition(*conditername, Teuchos::rcp(new DRT::Condition(*conds[i])));
    }
  }

  // re-partioning
  Redistribute(targetdis,condnoderowvec,condnodecolvec);
  Redistribute(sourcedis,othernoderowvec,othernodecolvec);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::Redistribute(
  Teuchos::RCP<DRT::Discretization> dis,
  std::vector<int>& noderowvec,
  std::vector<int>& nodecolvec)
{
  dis->CheckFilledGlobally();

  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(dis->Comm().Clone());

  Teuchos::RCP<Epetra_Map>  noderowmap = Teuchos::rcp(new Epetra_Map(-1,
                                                      noderowvec.size(),
                                                      &noderowvec[0],
                                                      0,
                                                      *comm));

  Teuchos::RCP<Epetra_Map>  nodecolmap = Teuchos::rcp(new Epetra_Map(-1,
                                                      nodecolvec.size(),
                                                      &nodecolvec[0],
                                                      0,
                                                      *comm));
  if (!dis->Filled())
    dis->Redistribute(*noderowmap,*nodecolmap);

  Teuchos::RCP<Epetra_Map> elerowmap = Teuchos::rcp( new Epetra_Map(*dis->ElementRowMap()));
  DRT::UTILS::PartUsingParMetis(dis,elerowmap,
      noderowmap,nodecolmap,comm,false);

  Teuchos::RCP<Epetra_Map> roweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> coleles  = Teuchos::null;
  dis->BuildElementRowColumn(*noderowmap,*nodecolmap,roweles,coleles);

  dis->ExportRowNodes(*noderowmap);
  dis->ExportRowElements(*roweles);

  dis->ExportColumnNodes(*nodecolmap);
  dis->ExportColumnElements(*coleles);

  dis->FillComplete();
}

