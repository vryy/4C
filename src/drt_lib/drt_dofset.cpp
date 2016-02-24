/*!----------------------------------------------------------------------
\file drt_dofset.cpp
\brief A set of degrees of freedom

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset.H"
#include "drt_dofset_proxy.H"
#include "drt_discret.H"
#include "drt_discret_hdg.H"
#include "drt_utils.H"

#include "../linalg/linalg_utils.H"

#include <Epetra_FECrsGraph.h>
#include <Ifpack_Graph_Epetra_CrsGraph.h>
#include <Ifpack_RCMReordering.h>
#include <Ifpack_METISReordering.h>
#include <Ifpack_AMDReordering.h>

// Bandwidth optimization is currently not working and to prevent it from
// generating compiler warnings, it is disabled by the following define.
#ifdef BW_OPT
  #undef BW_OPT
#endif

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet()
  : DRT::DofSetBase(), filled_(false), dspos_( 0 ), pccdofhandling_ (false)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::~DofSet()
{
  for (std::list<DofSetProxy*>::iterator i=proxies_.begin(); i!=proxies_.end(); ++i)
  {
    (*i)->Disconnect(this);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                               ukue 04/07|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::DofSet& dofset)
{
  dofset.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this  (public)                                      ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Print(std::ostream& os) const
{
  for (int proc=0; proc < numdfcolelements_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolelements_->Comm().MyPID())
    {
      if (numdfcolelements_->MyLength())
        os << "-------------------------- Proc " << proc << " :\n";
      for (int i=0; i<numdfcolelements_->MyLength(); ++i)
      {
        int numdf = (*numdfcolelements_)[i];
        int idx   = (*idxcolelements_)[i];
        os << i << ": ";
        for (int j=0; j<numdf; ++j)
          os << (idx+j) << " ";
        os << "\n";
      }
      os << std::endl;
    }
    numdfcolelements_->Comm().Barrier();
  }
  for (int proc=0; proc < numdfcolnodes_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolnodes_->Comm().MyPID())
    {
      if (numdfcolnodes_->MyLength())
        os << "-------------------------- Proc " << proc << " :\n";
      for (int i=0; i<numdfcolnodes_->MyLength(); ++i)
      {
        int numdf = (*numdfcolnodes_)[i];
        int idx   = (*idxcolnodes_)[i];
        os << i << ": ";
        for (int j=0; j<numdf; ++j)
          os << (idx+j) << " ";
        os << "\n";
      }
      os << std::endl;
    }
    numdfcolnodes_->Comm().Barrier();
  }
}


/*----------------------------------------------------------------------*
 |  reset everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Reset()
{
  dofrowmap_ = Teuchos::null;
  dofcolmap_ = Teuchos::null;
  numdfcolnodes_ = Teuchos::null;
  numdfcolelements_ = Teuchos::null;
  idxcolnodes_ = Teuchos::null;
  idxcolelements_ = Teuchos::null;

  filled_ = false;

  // tell all proxies
  NotifyReset();
}

/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int DRT::DofSet::AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
{
  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // A definite offset is currently not supported.
  // TODO (kronbichler) find a better solution for this
  //if (start!=0)
  //  dserror("right now user specified dof offsets are not supported");

  dspos_ = dspos;

  // Add DofSets in order of assignment to list. Once it is there it has its
  // place and will get its starting id from the previous DofSet.
  AddDofSettoList();

  // We assume that all dof sets before this one have been set up. Otherwise
  // we'd have to reorder the list.
  //
  // There is no test anymore to make sure that all prior dof sets have been
  // assigned. It seems people like to manipulate dof sets. People do create
  // dof sets that do not contain any dofs (on its first assignment), people
  // even shift dof set numbers to create overlapping dof sets. This is
  // perfectly fine.
  //
  // However if you rely on non-overlapping dof sets, you have to
  // FillComplete() your discretizations in the order of their creation. This
  // is guaranteed for all discretizations read from the input file since the
  // input reader calls FillComplete(). If you create your own discretizations
  // try to understand what you do.

  // Get highest GID used so far and add one
  int count = MaxGIDinList(dis.Comm()) + 1;

  // Check if we have a face discretization which supports degrees of freedom on faces
  Teuchos::RCP<const DiscretizationHDG> facedis =
    Teuchos::rcp_dynamic_cast<const DiscretizationHDG>(Teuchos::rcp(&dis,false));

  // Now this is tricky. We have to care for nodes, faces, and elements, both
  // row and column maps. In general both nodes, faces, and elements can have
  // dofs. In all cases these dofs might be shared with other nodes, faces,
  // or elements. (The very general case. For elements we'd probably
  // don't need that.)
  //
  // The point is that we have to make sure the dof numbering of a
  // mesh is independent of its parallel distribution. Otherwise we
  // could not redistribute a mesh. We would not be able to use old
  // distributed vectors afterwards.
  //
  // Each object (node or element) could have a different number of
  // dofs. The parallel distribution is arbitrary. So we fall back to
  // two redundant vectors here to gather the number of dofs per node
  // or element.

  // numdf for all nodes and elements
  numdfcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  numdfcolelements_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
  if (facedis != Teuchos::null && facedis->FaceColMap() != NULL)
    numdfcolfaces_ = Teuchos::rcp(new Epetra_IntVector(*facedis->FaceColMap()));

  // index of first dof for all nodes and elements
  idxcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  idxcolelements_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
  if (facedis != Teuchos::null && facedis->FaceColMap() != NULL)
    idxcolfaces_ = Teuchos::rcp(new Epetra_IntVector(*facedis->FaceColMap()));

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  int maxnodenumdf = 0;
  int maxelementnumdf = 0;
  std::map<int,std::vector<int> > nodedofset;
  std::map<int,std::vector<int> > nodeduplicatedofset;
  std::map<int,std::vector<int> > elementdofset;
  std::map<int,std::vector<int> > facedofset;

  // bandwidth optimization is implemented here and works. It though leads to
  // a totally different dof numbering (mixed between elements and nodes) and therefore
  // currently crashes the I/O of results.
  // Bandwidth opt seems to be also performed by most of our direct solvers
  // (umfpack, superlu) and does not matter that much in iterative methods.
  // So currently, we do not support bw, but might want to turn this on at some
  // later time.
  // Also, this method is also called from filters where the parameter list from input
  // does not exist. This is not resolved yet.
  // mwgee 01/13
#ifdef BW_OPT
  // TAW: Please do not include drt_problem.H!!!
  //bool bw = DRT::Problem::Instance()->BandWidthOpt();
  bool bw = false;
  if (bw)
  {
    int maxnodegid = dis.NodeRowMap()->MaxAllGID();
    // we have to invent some temporary element gids
    // that do not overlap with the node gids
    const Epetra_Map& elerowmap = *dis.ElementRowMap();
    Teuchos::RCP<Epetra_IntVector> elerowgids = Teuchos::rcp(new Epetra_IntVector(elerowmap));
    std::vector<int> ssizeelerowgids(elerowmap.Comm().NumProc(),0);
    std::vector<int> rsizeelerowgids(elerowmap.Comm().NumProc(),0);
    ssizeelerowgids[elerowmap.Comm().MyPID()] = elerowmap.NumMyElements();
    elerowmap.Comm().SumAll(&ssizeelerowgids[0],&rsizeelerowgids[0],(int)ssizeelerowgids.size());
    for (unsigned i=0; i<rsizeelerowgids.size(); ++i)
    {
      ssizeelerowgids[i] = maxnodegid+1;
      for (unsigned j=0; j<i; ++j)
        ssizeelerowgids[i] += rsizeelerowgids[j];
    }

    //if (elerowgids->Map().Comm().MyPID()==0)
    //  for (unsigned i=0; i<rsizeelerowgids.size(); ++i)
    //    printf("maxnodegid %d proc %d ssize %d rsize %d \n",maxnodegid,i,ssizeelerowgids[i],rsizeelerowgids[i]);

    for (int i=0; i<elerowgids->Map().NumMyElements(); ++i)
      (*elerowgids)[i] = ssizeelerowgids[elerowgids->Comm().MyPID()] + i;

    // export elerowgids to column map elecolgids
    Teuchos::RCP<Epetra_IntVector> elecolgids = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
    LINALG::Export(*elerowgids,*elecolgids);

    // Build a rowmap that contains node and element gids
    // Consider only elements that have dofs
    std::vector<int> mygids(0,0);
    for (int i=0; i<dis.NodeRowMap()->NumMyElements(); ++i)
      mygids.push_back(dis.NodeRowMap()->GID(i));
    for (int i=0; i<elerowgids->Map().NumMyElements(); ++i)
    {
      if (NumDofPerElement(*dis.lColElement(i)) != 0)
        mygids.push_back((*elerowgids)[i]);
    }

    Epetra_Map rowmapcombo(-1,(int)mygids.size(),&mygids[0],0,dis.Comm());

    // Loop Nodes and Elements and build a common nodal-elemental connectivity graph
    // Consider only elementsw rthat have dofs
    Teuchos::RCP<Epetra_FECrsGraph> graphcombo =
                        Teuchos::rcp( new Epetra_FECrsGraph(Copy,rowmapcombo,180,false));
    for (int i=0; i<dis.ElementColMap()->NumMyElements(); ++i)
    {
      const int  nnode   = dis.lColElement(i)->NumNode();
      int*       nodeids = const_cast<int*>(dis.lColElement(i)->NodeIds());
      for (int row=0; row<nnode; ++row)
      {
//        if (!dis.NodeRowMap()->MyGID(nodeids[row])) continue; // only add into my own nodal rows
        int err = graphcombo->InsertGlobalIndices(1,&nodeids[row],nnode,nodeids);
        if (err<0) dserror("graphcombo->InsertGlobalIndices returned err=%d",err);
      }
      // check whether this element is actually a row element
      //if (dis.ElementRowMap()->MyGID(dis.lColElement(i)->Id()) == false) continue;
      // check whether this element has elemental dofs
      if (NumDofPerElement(*dis.lColElement(i)) == 0) continue;
      // get the pseudo gid of this element
      int lid = dis.ElementColMap()->LID(dis.lColElement(i)->Id());
      if (lid==-1) dserror("Cannot find lid");
      int pseudogid = (*elecolgids)[lid];
      // check whether this pseudogid is actually in rowmapcombo on this proc
      //if (rowmapcombo.MyGID(pseudogid) != true) dserror("gid unknown on this proc");
      // the element dofs connect to itself
      int err = graphcombo->InsertGlobalIndices(1,&pseudogid,1,&pseudogid);
      if (err<0) dserror("graphcombo->InsertGlobalIndices returned err=%d",err);
      // the element dofs connect to the nodal dofs
      err = graphcombo->InsertGlobalIndices(1,&pseudogid,nnode,nodeids);
      if (err<0) dserror("graphcombo->InsertGlobalIndices returned err=%d",err);
      // the nodal dofs connect to the element dofs
      for (int row=0; row<nnode; ++row)
      {
        //if (!dis.NodeRowMap()->MyGID(nodeids[row])) continue; // only add into my own nodal rows
        int err = graphcombo->InsertGlobalIndices(1,&nodeids[row],1,&pseudogid);
        if (err<0) dserror("graphcombo->InsertGlobalIndices returned err=%d",err);
      }
    } // for (int i=0; i<dis.ElementColMap()->NumMyElements(); ++i)

    int err = graphcombo->GlobalAssemble(true);
    if (err) dserror("graphcombo->GlobalAssemble() returned err=%d",err);
    err = graphcombo->OptimizeStorage();
    if (err) dserror("graphcombo->OptimizeStorage() returned err=%d",err);
    //cout << *graphcombo; fflush(stdout);

    Ifpack_Graph_Epetra_CrsGraph ifgraph(graphcombo);
    Ifpack_RCMReordering reorderer;
//    Ifpack_AMDReordering reorderer;
    reorderer.Compute(ifgraph);
    std::vector<int> permute(rowmapcombo.NumMyElements(),-1);
    for (int i=0; i<rowmapcombo.NumMyElements(); ++i)
    {
      permute[i] = rowmapcombo.GID(reorderer.Reorder(i));
      if (permute[i]==-1) dserror("Cannot find gid for lid %d %d",i,reorderer.Reorder(i));
    }
    Epetra_IntVector permutation(::Copy,rowmapcombo,&permute[0]);
//    cout << rowmapcombo; fflush(stdout); rowmapcombo.Comm().Barrier();
//    fflush(stdout); rowmapcombo.Comm().Barrier(); cout << permutation; fflush(stdout);
//    rowmapcombo.Comm().Barrier();
//    exit(0);

#if 0
    for (int i=0; i<dis.Comm().NumProc(); ++i)
    {
      if (i==dis.Comm().MyPID())
      {
        for (int j=0; j<rowmapcombo.NumMyElements(); ++j)
          printf("proc %d oldgid %d newgid %d\n",i,rowmapcombo.GID(j),permute[j]);
      }
      fflush(stdout);
      dis.Comm().Barrier();
    }
#endif


    // Assign dofs to nodes using the permuted node gids
    // These two fields use the unpermuted gids and lids
    Epetra_IntVector numdfrownodes(*dis.NodeRowMap());
    Epetra_IntVector idxrownodes(*dis.NodeRowMap());

    int numrownodes = dis.NumMyRowNodes();
    for (int i=0; i<numrownodes; ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      numdfrownodes[i] = NumDofPerNode(*actnode);
    }

//    int minnodegid = dis.NodeRowMap()->MinAllGID();
    const int minnodegid = permutation.MinValue();
    maxnodenumdf = numdfrownodes.MaxValue();

    for (int i=0; i<numrownodes; ++i)
    {
      const int gid = dis.lRowNode(i)->Id();
      const int plid = permutation.Map().LID(gid);
      if (plid==-1) dserror("Cannot find lid for gid %d",gid);
      const int permutegid = permutation[plid];
      if (i!=plid) printf("Found an inconsistency %d vs. %d\n",i,plid);
      int numdf = numdfrownodes[i];
//      int dof = count + ( gid-minnodegid )*maxnodenumdf;
      int dof = count + ( permutegid-minnodegid )*maxnodenumdf;
      idxrownodes[i] = dof;
      std::vector<int> & dofs = nodedofset[gid]; // FIXME: permutegid? No!
      if (dofs.size()) dserror("Dofs for node-gid %d have been previously assigned",permutegid);
      dofs.reserve(numdf);
      if (!numdf) dserror("Node has numdf %d",numdf);
      for ( int j=0; j<numdf; ++j ) dofs.push_back(dof+j);
    }
    // These col vectors use the original gids
    Epetra_Import nodeimporter( numdfcolnodes_->Map(), numdfrownodes.Map() );
    err = numdfcolnodes_->Import( numdfrownodes, nodeimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );
    err = idxcolnodes_->Import( idxrownodes, nodeimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );

    //cout << idxrownodes;

    // do not increase the count, refer to the lowest gid in town which is nodal
    //count = idxrownodes.MaxValue() + maxnodenumdf + 1;

  //////////////////////////////////////////////////////////////////

  // Now do it again for the elements
  // These vectors use the original element gids
  Epetra_IntVector numdfrowelements(*dis.ElementRowMap());
  Epetra_IntVector idxrowelements(*dis.ElementRowMap());

  int numrowelements = dis.NumMyRowElements();
  for (int i=0; i<numrowelements; ++i)
  {
    DRT::Element* actele = dis.lRowElement(i);
    //const int gid = actele->Id();
    int numdf = NumDofPerElement(*actele);
    numdfrowelements[i] = numdf;
  }

  //int minelementgid = dis.ElementRowMap()->MinAllGID();
  // Since we have not increased count, we refer to the lowest gid in town
  const int minelementgid = minnodegid;
  maxelementnumdf = numdfrowelements.MaxValue();
  maxelementnumdf = std::max(maxelementnumdf,maxnodenumdf);

  for (int i=0; i<numrowelements; ++i)
  {
    const int gid = dis.lRowElement(i)->Id();
    const int uniquegid = (*elerowgids)[i];
    const int lid = permutation.Map().LID(uniquegid);
    const int permutegid = permutation[lid];
    int numdf = numdfrowelements[i];
    if (lid==-1 && numdf) dserror("Cannot find element lid in rowmapcombo");
    int dof = count + ( permutegid-minelementgid )*maxelementnumdf;
    // in case the element does not have any dofs the following is wrong,
    // but it won't matter, since numdf==0
    idxrowelements[i] = dof;
//    std::vector<int> & test = nodedofset[permutegid];
//    if (test.size()) dserror("Dofs for ele-gid %d have been previously assigned as nodal",permutegid);
    std::vector<int>& dofs = elementdofset[gid]; // FIXME permutegid? No!
    if (dofs.size()) dserror("Dofs for ele-gid %d have been previously assigned",permutegid);
    dofs.reserve(numdf);
    if (!numdf) dserror("Ele has numdf %d",numdf);
    for ( int j=0; j<numdf; ++j ) dofs.push_back(dof+j);
  }

  Epetra_Import elementimporter( numdfcolelements_->Map(), numdfrowelements.Map() );
  err = numdfcolelements_->Import( numdfrowelements, elementimporter, Insert );
  if (err) dserror( "Import using importer returned err=%d", err );
  err = idxcolelements_->Import( idxrowelements, elementimporter, Insert );
  if (err) dserror( "Import using importer returned err=%d", err );

  //cout << idxrowelements;
  } // if (bw)
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  else // don't do bandwidth optimization
#endif /* BW_OPT */
  {

    // get DoF coupling conditions
    std::vector<DRT::Condition*> couplingconditions(0);
    dis.GetCondition("PointCoupling", couplingconditions);
    if ((int)couplingconditions.size()>0)
      pccdofhandling_=true;

    // do the nodes first
    Epetra_IntVector numdfrownodes(*dis.NodeRowMap());
    Epetra_IntVector idxrownodes(*dis.NodeRowMap());

    int numrownodes = dis.NumMyRowNodes();
    for (int i=0; i<numrownodes; ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      numdfrownodes[i] = NumDofPerNode(*actnode);
    }

    int minnodegid = dis.NodeRowMap()->MinAllGID();
    maxnodenumdf = numdfrownodes.MaxValue();
    GetReservedMaxNumDofperNode(maxnodenumdf); //XFEM::XFEMDofSet set to const number!

    for (int i=0; i<numrownodes; ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const int gid = actnode->Id();

      // **********************************************************************
      // **********************************************************************
      // check for DoF coupling conditions                         popp 02/2016
      // **********************************************************************
      // **********************************************************************
      int relevantcondid = -1;
      for (int k=0; k<(int)couplingconditions.size();++k)
      {
        if (couplingconditions[k]->ContainsNode(gid))
        {
          if (relevantcondid != -1)
            dserror("ERROR: Two coupling conditions on one node");
          relevantcondid = k;
        }
      }

      // check for node coupling condition and slave/master status
      bool specialtreatment = false;
      if (relevantcondid>=0)
      {
        const std::vector<int>* nodeids = couplingconditions[relevantcondid]->Nodes();
        if (!nodeids) dserror("ERROR: Condition does not have Node Ids");

        // do nothing for first (master) node in coupling condition
        // do something for second, third, ... (slave) node
        if ((*nodeids)[0]!=gid)
        {
          // critical case
          specialtreatment = true;

          // check which dofs are to be coupled
          int numdofcond = couplingconditions[relevantcondid]->GetInt("numdof");
          if (numdofcond != numdfrownodes[i])
            dserror("ERROR: Number of DoFs in coupling condition does not match node");

          // check number of coupled DoFs
          int numdepdof = 0;
          const std::vector<int>* onoffcond = couplingconditions[relevantcondid]->Get<std::vector<int>>("onoff");
          for (int n=0;n<(int)(onoffcond->size());++n)
          {
            if ((*onoffcond)[n]==1)
              numdepdof++;
          }
          if (numdepdof==numdfrownodes[i])
            dserror("ERROR: All DoFs coupled at junction --> no junction needeed, just merge nodes!");

          // get master node of this condition
          int mgid = (*nodeids)[0];
          std::vector<int> & mdofs = nodedofset[mgid];
          if ((int)(mdofs.size())==0) dserror ("ERROR: Master node has not yet been initialized with DoFs");

          // special treatment
          int numdf = numdfrownodes[i];
          int dof = count + ( gid-minnodegid )*maxnodenumdf;
          idxrownodes[i] = dof;
          std::vector<int> & dofs = nodedofset[gid];
          std::vector<int> & duplicatedofs = nodeduplicatedofset[gid];
          dofs.reserve( numdf );
          duplicatedofs.reserve( numdf );
          for ( int j=0; j<numdf; ++j )
          {
            // push back master node DoF ID if coupled
            if ((*onoffcond)[j]==1)
            {
              dofs.push_back( mdofs[j] );
              duplicatedofs.push_back( 1 );
            }
            // push back new DoF ID if not coupled
            else
            {
              dofs.push_back( dof+j );
              duplicatedofs.push_back( 0 );
            }
          }
        }
      }

      // standard treatment for non-coupling nodes and master coupling nodes
      if (!specialtreatment)
      {
        int numdf = numdfrownodes[i];
        int dof = count + ( gid-minnodegid )*maxnodenumdf;
        idxrownodes[i] = dof;
        std::vector<int> & dofs = nodedofset[gid];
        std::vector<int> & duplicatedofs = nodeduplicatedofset[gid];
        dofs.reserve( numdf );
        duplicatedofs.reserve( numdf );
        for ( int j=0; j<numdf; ++j )
        {
          dofs.push_back( dof+j );
          duplicatedofs.push_back( 0 );
        }
      }
      // **********************************************************************
      // **********************************************************************
      // **********************************************************************
      // **********************************************************************
    }

    Epetra_Import nodeimporter( numdfcolnodes_->Map(), numdfrownodes.Map() );
    int err = numdfcolnodes_->Import( numdfrownodes, nodeimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );
    err = idxcolnodes_->Import( idxrownodes, nodeimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );

    count = maxnodenumdf>0 ? idxrownodes.MaxValue() + maxnodenumdf : 0;

    //////////////////////////////////////////////////////////////////

    // Now do it again for the faces
    if (facedis != Teuchos::null && facedis->FaceRowMap() != NULL)
    {
      Epetra_IntVector numdfrowfaces(*facedis->FaceRowMap());
      Epetra_IntVector idxrowfaces(*facedis->FaceRowMap());
      int numcolelements = dis.NumMyColElements();

      const int mypid = dis.Comm().MyPID();
      for (int i=0; i<numcolelements; ++i)
      {
        DRT::FaceElement** faces = dis.lColElement(i)->Faces();
        // If no faces are found, continue...
        if (faces == NULL)
          continue;
        for (int face=0; face<dis.lColElement(i)->NumFace(); ++face)
          if (faces[face]->Owner() == mypid) {
            const int mylid = facedis->FaceRowMap()->LID(faces[face]->Id());
            numdfrowfaces[mylid] = NumDofPerFace(*(dis.lColElement(i)),face);
          }
      }

      int minfacegid = facedis->FaceRowMap()->MinAllGID();
      int maxfacenumdf = numdfrowfaces.MaxValue();

      for (int i=0; i<numcolelements; ++i)
      {
        DRT::FaceElement** faces = dis.lColElement(i)->Faces();
        if (faces == NULL)
          continue;
        for (int face=0; face<dis.lColElement(i)->NumFace(); ++face)
          if (faces[face]->Owner() == mypid)
          {
            const int gid = faces[face]->Id();
            const int mylid = facedis->FaceRowMap()->LID(gid);
            int numdf = numdfrowfaces[mylid];
            int dof = count + ( gid-minfacegid )*maxfacenumdf;
            idxrowfaces[mylid] = dof;
            std::vector<int> & dofs = facedofset[gid];
            // do not visit the same face more than once
            if (dofs.empty())
            {
              dofs.reserve( numdf );
              for ( int j=0; j<numdf; ++j )
              {
                dofs.push_back( dof+j );
              }
            }
          }
      }

      Epetra_Import faceimporter( numdfcolfaces_->Map(), numdfrowfaces.Map() );
      err = numdfcolfaces_->Import( numdfrowfaces, faceimporter, Insert );
      if (err) dserror( "Import using importer returned err=%d", err );
      err = idxcolfaces_->Import( idxrowfaces, faceimporter, Insert );
      if (err) dserror( "Import using importer returned err=%d", err );

      count = idxrowfaces.MaxValue() + maxfacenumdf;
    }

    //////////////////////////////////////////////////////////////////

    // Now do it again for the elements
    Epetra_IntVector numdfrowelements(*dis.ElementRowMap());
    Epetra_IntVector idxrowelements(*dis.ElementRowMap());

    int numrowelements = dis.NumMyRowElements();
    for (int i=0; i<numrowelements; ++i)
    {
      DRT::Element* actele = dis.lRowElement(i);
      //const int gid = actele->Id();
      int numdf = NumDofPerElement(*actele);
      numdfrowelements[i] = numdf;
    }

    int minelementgid = dis.ElementRowMap()->MinAllGID();
    maxelementnumdf = numdfrowelements.MaxValue();

    for (int i=0; i<numrowelements; ++i)
    {
      DRT::Element* actelement = dis.lRowElement(i);
      const int gid = actelement->Id();
      int numdf = numdfrowelements[i];
      int dof = count + ( gid-minelementgid )*maxelementnumdf;
      idxrowelements[i] = dof;
      std::vector<int> & dofs = elementdofset[gid];
      dofs.reserve( numdf );
      for ( int j=0; j<numdf; ++j )
      {
        dofs.push_back( dof+j );
      }
    }

    Epetra_Import elementimporter( numdfcolelements_->Map(), numdfrowelements.Map() );
    err = numdfcolelements_->Import( numdfrowelements, elementimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );
    err = idxcolelements_->Import( idxrowelements, elementimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );

  } // end of else
  /////////////////////////////////////////////////////////////////

  // Now finally we have everything in place to build the maps.
  int numrownodes = dis.NumMyRowNodes();
  int numrowelements = dis.NumMyRowElements();

  std::vector<int> localrowdofs;
  std::vector<int> localcoldofs;
  localrowdofs.reserve( numrownodes*maxnodenumdf + numrowelements*maxelementnumdf );
  localcoldofs.reserve( numrownodes*maxnodenumdf + numrowelements*maxelementnumdf );

  std::vector<int> allnodelocalcoldofs;
  allnodelocalcoldofs.reserve( numrownodes*maxnodenumdf );

  for ( std::map<int,std::vector<int> >::iterator i=nodedofset.begin();
        i!=nodedofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::vector<int> & duplicatedofs = nodeduplicatedofset[i->first];
    std::vector<int> cleandofs;
    for (unsigned j=0; j<dofs.size(); ++j)
    {
      if (duplicatedofs[j]==0) cleandofs.push_back(dofs[j]);
    }
    std::copy( cleandofs.begin(), cleandofs.end(), std::back_inserter( localrowdofs ) );
    //printf("Proc %d nodal gid %d ndofs %d\n",proc,i->first,(int)dofs.size());
    //for (unsigned j=0; j<dofs.size(); ++j) printf(" %d ",dofs[j]);
    //printf("\n");
  }
  for ( std::map<int,std::vector<int> >::iterator i=facedofset.begin();
        i!=facedofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localrowdofs ) );
    //printf("Proc %d ele gid %d ndofs %d\n",dis.Comm().MyPID(),i->first,(int)dofs.size());
    //for (unsigned j=0; j<dofs.size(); ++j) printf(" %d ",dofs[j]);
    //printf("\n");
  }
  for ( std::map<int,std::vector<int> >::iterator i=elementdofset.begin();
        i!=elementdofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localrowdofs ) );
    //printf("Proc %d ele gid %d ndofs %d\n",dis.Comm().MyPID(),i->first,(int)dofs.size());
    //for (unsigned j=0; j<dofs.size(); ++j) printf(" %d ",dofs[j]);
    //printf("\n");
  }

  Exporter nodeexporter( *dis.NodeRowMap(), *dis.NodeColMap(), dis.Comm() );
  nodeexporter.Export( nodedofset );
  nodeexporter.Export( nodeduplicatedofset );

  Exporter elementexporter( *dis.ElementRowMap(), *dis.ElementColMap(), dis.Comm() );
  elementexporter.Export( elementdofset );

  if (facedis != Teuchos::null && facedis->FaceRowMap() != NULL)
  {
    Exporter faceexporter( *facedis->FaceRowMap(), *facedis->FaceColMap(), dis.Comm() );
    faceexporter.Export( facedofset );
  }

  for ( std::map<int,std::vector<int> >::iterator i=nodedofset.begin();
        i!=nodedofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::vector<int> & duplicatedofs = nodeduplicatedofset[i->first];
    std::vector<int> cleandofs;
    for (unsigned j=0; j<dofs.size(); ++j)
    {
      if (duplicatedofs[j]==0) cleandofs.push_back(dofs[j]);
    }
    std::copy( cleandofs.begin(), cleandofs.end(), std::back_inserter( localcoldofs ) );
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( allnodelocalcoldofs ) );
  }
  for ( std::map<int,std::vector<int> >::iterator i=facedofset.begin();
        i!=facedofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localcoldofs ) );
  }
  for ( std::map<int,std::vector<int> >::iterator i=elementdofset.begin();
        i!=elementdofset.end();
        ++i )
  {
    std::vector<int> & dofs = i->second;
    std::copy( dofs.begin(), dofs.end(), std::back_inserter( localcoldofs ) );
  }

  // in case of bandwidth optimization we have to sort localrowdofs and
  // localcoldofs in ascending order for the optimization to take effect.
  // the linear solver will operate in local ids, so it does not care for
  // gids to be renumbered if the lids do not change
#ifdef BW_OPT
  if (bw)
  {
    std::sort(localrowdofs.begin(),localrowdofs.end());
    std::sort(localcoldofs.begin(),localcoldofs.end());
  }
#endif

  dofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,localrowdofs.size(),&localrowdofs[0],0,dis.Comm()));
  if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");
  dofcolmap_ = Teuchos::rcp(new Epetra_Map(-1,localcoldofs.size(),&localcoldofs[0],0,dis.Comm()));

  // **********************************************************************
  // **********************************************************************
  // build map of all (non-unique) column DoFs
  dofscolnodes_ = Teuchos::rcp(new Epetra_Map(-1,allnodelocalcoldofs.size(),&allnodelocalcoldofs[0],0,dis.Comm()));

  // build shift vector
  shiftcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  int numcolnodes = dis.NumMyColNodes();
  for (int i=0; i<numcolnodes; ++i)
  {
    if (i==0)
    {
      (*shiftcolnodes_)[i] = 0;
    }
    else
    {
      DRT::Node* lastnode = dis.lColNode(i-1);
      (*shiftcolnodes_)[i] = (*shiftcolnodes_)[i-1] + NumDofPerNode(*lastnode);
    }
  }
  // **********************************************************************
  // **********************************************************************

  // degrees of freedom have now been assigned
  filled_ = true;

  // tell all proxies
  NotifyAssigned();

  // return maximum dof number of this dofset (+1)
  count = dofrowmap_->MaxAllGID() + 1;
  return count;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::RegisterProxy(DofSetProxy* proxy)
{
  proxies_.push_back(proxy);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::UnregisterProxy(DofSetProxy* proxy)
{
  proxies_.remove(proxy);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::NotifyAssigned()
{
  for (std::list<DofSetProxy*>::iterator i=proxies_.begin(); i!=proxies_.end(); ++i)
    (*i)->NotifyAssigned();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSet::NotifyReset()
{
  for (std::list<DofSetProxy*>::iterator i=proxies_.begin(); i!=proxies_.end(); ++i)
    (*i)->NotifyReset();
}


