/*!----------------------------------------------------------------------*
\file so_nstet_nodalstrain.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "so_nstet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_FECrsMatrix.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NStetType::Initialize(DRT::Discretization& dis)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::Initialize");

  const int myrank = dis.Comm().MyPID();
  const int numproc = dis.Comm().NumProc();

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // init elements, make maps of column elements and row nodes
  const int numele = dis.NumMyColElements();
  for (int i=0; i<numele; ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::NStet* actele = 
                    dynamic_cast<DRT::ELEMENTS::NStet*>(dis.lColElement(i));
    if (!actele) dserror("cast to NStet* failed");

    // init the element
    actele->InitElement();

    // register element in list of column nstet elements
    elecids_[actele->Id()] = actele;

    // compute a map of all row nodes adjacent to a NStet element
    for (int j=0; j<actele->NumNode(); ++j) 
    {
      DRT::Node* node = actele->Nodes()[j];
      if (myrank == node->Owner())
        noderids_[node->Id()] = node;
    }
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // compute adjacency for each row node
  // make patch of adjacent elements
  // make patch of adjacent nodes (including center node itself)
  // make lm for nodal patch
  // make lmowner for nodal patch
  std::map<int,DRT::Node*>::iterator node;
  for (node=noderids_.begin(); node != noderids_.end(); ++node)
  {
    DRT::Node* nodeL  = node->second;
    const int nodeidL = nodeL->Id();

    // list of adjacent elements
    vector<DRT::ELEMENTS::NStet*> adjele(0);
    for (int j=0; j<nodeL->NumElement(); ++j)
    {
      const int eleid = node->second->Elements()[j]->Id();
      std::map<int,DRT::ELEMENTS::NStet*>::iterator ele = elecids_.find(eleid);
      if (ele==elecids_.end()) continue;
      adjele.push_back(ele->second);
    }
    adjele_[nodeidL] = adjele;

    // patch of all nodes adjacent to adjacent elements
    map<int,DRT::Node*> nodepatch;
    for (int j=0; j<(int)adjele.size(); ++j) 
    {
      DRT::Node** nodes = adjele[j]->Nodes();
      for (int k=0; k<adjele[j]->NumNode(); ++k)
        nodepatch[nodes[k]->Id()] = nodes[k];
    }
    adjnode_[nodeidL] = nodepatch;

    // lm and lmowner arrays
    const int numnodepatch = (int)nodepatch.size();
    const int ndofperpatch = numnodepatch*3;

    // location and ownership vector of nodal patch
    vector<int> lm(ndofperpatch);
    std::map<int,DRT::Node*>::iterator pnode;
    int count=0;
    for (pnode=nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
    {
      const vector<int>& dofs = dis.Dof(pnode->second);
      for (int j=0; j<(int)dofs.size(); ++j)
        lm[count++]        = dofs[j];
    }
    adjlm_[nodeidL] = lm;
  } // for (node=noderids_.begin(); node != noderids_.end(); ++node)


  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // stuff related to MIS stabilization of pressure
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // build parallel maximum independent set of nodes (MIS-nodes)
  // this is done in pseudo-serial, as a true parallel MIS algorithm is not
  // easy and I'm too lazy to do one
  map<int,DRT::Node*> rnodes = noderids_; // working copy of row nodes
  vector<int> misnodes(0); // indicator which nodes are mis
  vector<int> smisnodes(0);// same but for communication
  for (int proc=0; proc<numproc; ++proc)
  {
    if (proc==myrank)
    {
      for (node=rnodes.begin(); node != rnodes.end(); ++node)
      {
        const int actnodeid = node->second->Id();
        // make this node a mis node 
        misnodes.push_back(actnodeid);
        
        cout << myrank << " MIS    " << *(node->second) << endl;
        
        // delete all neighboring nodes from its patch from map
        map<int,DRT::Node*>& nodepatch = adjnode_[actnodeid];
        std::map<int,DRT::Node*>::iterator neighbor;
        for (neighbor=nodepatch.begin(); neighbor != nodepatch.end(); ++neighbor)
        {
          if (neighbor->first == actnodeid) continue; // do not delete myself
          /*const int n =*/ rnodes.erase(neighbor->first);
          //if (n)
          //cout << myrank << " delon  " << *(neighbor->second) << endl;
        }
      }
    } // if (proc==mypid)

    // broadcast mis nodes of proc
    int size = misnodes.size();
    dis.Comm().Broadcast(&size,1,proc);
    if (proc==myrank) smisnodes = misnodes;
    else              smisnodes.resize(size);
    dis.Comm().Broadcast(&smisnodes[0],size,proc);

    // all other procs have to remove nodes adjacent to mis nodes from 
    // their potential list
    if (myrank>proc)
    {
      for (node=rnodes.begin(); node != rnodes.end();)
      {
        const int actnodeid = node->second->Id();
        map<int,DRT::Node*>& nodepatch = adjnode_[actnodeid];
        std::map<int,DRT::Node*>::iterator neighbor;
        bool foundit = false;
        for (neighbor=nodepatch.begin(); neighbor != nodepatch.end(); ++neighbor)
        {
          const int neighborid = neighbor->first;
          for (int i=0; i<size; ++i)
            if (neighborid == smisnodes[i])
            {
              foundit = true;
              break;
            }
          if (foundit) break;
        }
        if (foundit)
        {
          //cout << myrank << " deloff " << *(node->second) << endl;
          rnodes.erase(node++);
        }
        else node++;
      }
    }
    dis.Comm().Barrier();
    smisnodes.clear();
  } // for (proc=0; proc<numproc; ++proc)

  // convert the misnodes vector to a map because its easier to search
  map<int,int> misnodesmap;
  for (int i=0; i<(int)misnodes.size(); ++i)
    misnodesmap[misnodes[i]] = misnodes[i];
  misnodes.clear();

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // each MIS node is associated with a patch of column elements of which
  // it takes the full integration area
  // Additionally, there are leftover column elements (not adjacent to any MIS node)
  // they are taken by a greedy algorithm in a second phase

  //----------------------------------------------------------------------
  // make a working map of column elements
  map<int,DRT::ELEMENTS::NStet*> elecids = elecids_;
  dis.Comm().Barrier();
  
  //----------------------------------------------------------------------
  // assign mis nodes all surrounding row elements (greedy phase 1)
  for (int proc=0; proc<numproc; ++proc)
  {
    vector<int> sendeles(0);
    map<int,int>::iterator mis;
    if (proc==myrank)
    {
      for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        vector<DRT::ELEMENTS::NStet*> eles(0);
        DRT::Node* misnode = noderids_[mis->first];
        for (int i=0; i<misnode->NumElement(); ++i)
        {
          if (elecids.find(misnode->Elements()[i]->Id()) == elecids.end()) continue;
          eles.push_back((DRT::ELEMENTS::NStet*)misnode->Elements()[i]);
          elecids.erase(misnode->Elements()[i]->Id());
        } // i
        pstab_adjele_[mis->first] = eles;
      } 
      
      // Unassign mis nodes with patches that are too small
      for (mis=misnodesmap.begin(); mis != misnodesmap.end();)
      {
        vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele_[mis->first];
        const int patchsize = (int)adjele.size();
        cout << "Proc " << myrank << " MIS node " << mis->first << " patchsize " << patchsize << endl;
        if (patchsize >= MIS_MIN_PATCHSIZE) mis++;
        else
        {
          cout << "Proc " << myrank << " erased MIS node " << mis->first << " patchsize " << patchsize << endl;
          for (int i=0; i<patchsize; ++i) elecids[adjele[i]->Id()] = adjele[i];
          pstab_adjele_.erase(mis->first);
          misnodesmap.erase(mis++);
        }
      } 
      
      // make communication vector of all elements taken
      for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele_[mis->first];
        for (int i=0; i<(int)adjele.size(); ++i) sendeles.push_back(adjele[i]->Id());
      } 
      
    } // if (proc==myrank)
    
    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size,1,proc);
    if (proc != myrank) sendeles.resize(size);
    dis.Comm().Broadcast(&sendeles[0],size,proc);
    
    // all other procs remove the already taken elements from their list
    if (myrank != proc)
    {
      for (int i=0; i<size; ++i) elecids.erase(sendeles[i]);
    } // if (myrank != proc)
    
    sendeles.clear();
    dis.Comm().Barrier();
  } // for (int proc=0; proc<numproc; ++proc)

  //----------------------------------------------------------------------
  // assign all leftover elements to neighboring patches (greedy phase 2)
  // this is a distance 2 patch search
  for (int proc=0; proc<numproc; ++proc)
  {
    vector<int> sendeles(0);
    map<int,int>::iterator mis;
    if (proc==myrank)
    {
      for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        vector<DRT::ELEMENTS::NStet*> eles(0);
        //DRT::Node* misnode = noderids_[mis->first];
        vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele_[mis->first];
        for (int i=0; i<(int)adjele.size(); ++i)
        {
          DRT::ELEMENTS::NStet* actele = adjele[i];
          for (int j=0; j<actele->NumNode(); ++j)
          {
            DRT::Node* node = actele->Nodes()[j];
            if (node->Id() == mis->first) continue;
            for (int k=0; k<node->NumElement(); ++k)
            {
              DRT::Element* ele = node->Elements()[k];
              if (elecids.find(ele->Id()) == elecids.end()) continue;
              pstab_adjele_[mis->first].push_back((DRT::ELEMENTS::NStet*)ele);
              elecids.erase(ele->Id());
              sendeles.push_back(ele->Id());
              cout << "Proc " << myrank << " leftover NStet " << ele->Id() << 
                      " found MIS node " << mis->second << endl;
            } // k
          } // j
        } // i
      }
    } // if (proc==myrank)

    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size,1,proc);
    if (proc != myrank) sendeles.resize(size);
    dis.Comm().Broadcast(&sendeles[0],size,proc);

    // all other procs remove the already taken elements from their list
    if (myrank != proc)
    {
      for (int i=0; i<size; ++i) elecids.erase(sendeles[i]);
    } // if (myrank != proc)
    
    sendeles.clear();
    dis.Comm().Barrier();

  } // for (int proc=0; proc<numproc; ++proc)

  //----------------------------------------------------------------------
  // test whether all column elements on all procs have been assigned a patch
  if ((int)elecids.size() != 0)
  {
    map<int,DRT::ELEMENTS::NStet*>::iterator ele;
    for (ele=elecids.begin(); ele != elecids.end(); ++ele)
    {
      cout << "Proc " << myrank <<
              " leftover NStet " << ele->second->Id() << 
              " found NO MIS node (on any proc)" << endl << *ele->second << endl;
    }
    dserror("Proc %d has the above column elements leftover",myrank);
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // have to build adjnode and adjlm arrays for the patches
  map<int,vector<DRT::ELEMENTS::NStet*> >::iterator mis;
  for (mis=pstab_adjele_.begin(); mis != pstab_adjele_.end(); ++mis)
  {
    int id = mis->first;
    int mispatchsize = (int)mis->second.size();
    int patchsize = (int)adjele_[id].size();
    cout << "Proc " << myrank 
         << " MIS " << id
         << " mispatchsize " << mispatchsize
         << " patchsize " << patchsize << endl;
    if (mispatchsize==patchsize) 
    {
      pstab_ident_patch_[id] = true;
      mis->second.clear();
    }
    else // patch not identical
    {
      // adjnode
      pstab_ident_patch_[id] = false;
      map<int,DRT::Node*> nodepatch;
      vector<DRT::ELEMENTS::NStet*>& adjele = mis->second;
      for (int j=0; j<(int)adjele.size(); ++j) 
      {
        DRT::Node** nodes = adjele[j]->Nodes();
        for (int k=0; k<adjele[j]->NumNode(); ++k)
          nodepatch[nodes[k]->Id()] = nodes[k];
      }
      pstab_adjnode_[id] = nodepatch;
      
      // lm array
      const int numnodepatch = (int)nodepatch.size();
      const int ndofperpatch = numnodepatch*3;
      vector<int> lm(ndofperpatch);
      std::map<int,DRT::Node*>::iterator pnode;
      int count=0;
      for (pnode=nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
      {
        const vector<int>& dofs = dis.Dof(pnode->second);
        for (int j=0; j<(int)dofs.size(); ++j)
          lm[count++] = dofs[j];
      }
      pstab_adjlm_[id] = lm;
    } // else
  } // for (mis=pstab_adjele_.begin(); mis != pstab_adjele_.end(); ++mis)
  
  // the elements should know which MIS node they belong to and on what
  // processor that MIS node is on. This is neccessary to do stress output
  
  return 0;
}


/*----------------------------------------------------------------------*
 |  pre-evaluation of elements (public)                        gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::PreEvaluate(DRT::Discretization& dis,
                                          Teuchos::ParameterList& p,
                                          RCP<LINALG::SparseOperator> systemmatrix1,
                                          RCP<LINALG::SparseOperator> systemmatrix2,
                                          RCP<Epetra_Vector>          systemvector1,
                                          RCP<Epetra_Vector>          systemvector2,
                                          RCP<Epetra_Vector>          systemvector3)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::PreEvaluate");

  // nodal integration for nlnstiff and internal forces only
  // (this method does not compute stresses/strains/element updates/mass matrix)
  string& action = p.get<string>("action","none");
  if (action != "calc_struct_nlnstiffmass" &&
      action != "calc_struct_nlnstiff"     &&
      action != "calc_struct_stress") return;

  // These get filled in here, so remove old stuff
  nodestress_.clear();
  nodestrain_.clear();

  // see what we have for input
  bool assemblemat1 = systemmatrix1!=Teuchos::null;
  bool assemblevec1 = systemvector1!=Teuchos::null;
  bool assemblevec2 = systemvector2!=Teuchos::null;
  bool assemblevec3 = systemvector3!=Teuchos::null;
  if (assemblevec2 || assemblevec3) dserror("Wrong assembly expectations");

  // nodal stiffness and force (we don't do mass here)
  LINALG::SerialDenseMatrix stiff;
  LINALG::SerialDenseVector force;

  // nodal stiffness and force for pressure stabilization
  LINALG::SerialDenseMatrix pstab_stiff;
  LINALG::SerialDenseVector pstab_force;

  //-------------------------------------- construct F for each NStet
  // current displacement
  RCP<const Epetra_Vector> disp = dis.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  // loop elements
  std::map<int,DRT::ELEMENTS::NStet*>::iterator ele;
  for (ele=elecids_.begin(); ele != elecids_.end(); ++ele)
  {
    vector<int> lm;
    vector<int> lmowner;
    ele->second->LocationVector(dis,lm,lmowner);
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    ele->second->DeformationGradient(mydisp);
  }

  //-----------------------------------------------------------------
  // create a temporary matrix to assemble to in a baci-unusual way
  // (across-parallel-interface assembly)
  const Epetra_Map& rmap = *dis.DofRowMap();
  const Epetra_Map& dmap = rmap;

  RCP<Epetra_FECrsMatrix> stifftmp;
  RCP<LINALG::SparseMatrix> systemmatrix = rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix1);
  if (systemmatrix != null and systemmatrix->Filled())
    stifftmp = rcp(new Epetra_FECrsMatrix(Copy,systemmatrix->EpetraMatrix()->Graph()));
  else
    stifftmp = rcp(new Epetra_FECrsMatrix(Copy,rmap,256,false));

  // create temporary vector in column map to assemble to
  // FIXME: For MIS, this might NOT be sufficient!
  Epetra_Vector forcetmp1(*dis.DofColMap(),true);

  //================================================== do nodal stiffness
  std::map<int,DRT::Node*>::iterator node;
  for (node=noderids_.begin(); node != noderids_.end(); ++node)
  {
    DRT::Node* nodeL   = node->second;     // row node
    const int  nodeLid = nodeL->Id();

    // list of elements on standard patch
    vector<DRT::ELEMENTS::NStet*>& adjele = adjele_[nodeLid];

    // list of nodes on standard patch
    map<int,DRT::Node*>& adjnode = adjnode_[nodeLid];

    // location vector of nodal patch
    vector<int>& lm = adjlm_[nodeLid];

    // total number of degrees of freedom on patch
    const int ndofperpatch = (int)lm.size();

    // if node is MIS node it might have a second patch
    bool mis = (pstab_ident_patch_.find(nodeLid) != pstab_ident_patch_.end());
    bool samepatch = false;
    vector<DRT::ELEMENTS::NStet*>* ps_adjele = &adjele;
    map<int,DRT::Node*>*           ps_adjnode = &adjnode;
    vector<int>*                   ps_lm = &lm;
    int                            ps_ndofperpatch = ndofperpatch;
    LINALG::SerialDenseMatrix*     ps_stiff = &stiff;
    LINALG::SerialDenseVector*     ps_force = &force;
    if (mis)
    {
      samepatch = pstab_ident_patch_.find(nodeLid)->second;
      if (!samepatch)
      {
        ps_adjele = &(pstab_adjele_[nodeLid]);
        ps_adjnode = &(pstab_adjnode_[nodeLid]);
        ps_lm = &(pstab_adjlm_[nodeLid]);
        ps_ndofperpatch = (int)ps_lm->size();
        
        pstab_stiff.LightShape(ps_ndofperpatch,ps_ndofperpatch); pstab_stiff.Zero();
        pstab_force.LightSize(ps_ndofperpatch);                  pstab_force.Zero();
        
        ps_stiff = &pstab_stiff;
        ps_force = &pstab_force;
      }
    }

    if (action != "calc_struct_stress")
    {
      // do nodal integration of stiffness and internal force
      stiff.LightShape(ndofperpatch,ndofperpatch);
      force.LightSize(ndofperpatch);
      NodalIntegration(&stiff,&force,adjnode,adjele,
                       mis,samepatch,ps_stiff,ps_force,ps_adjnode,ps_adjele,
                       NULL,NULL,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    else
    {
      INPAR::STR::StressType iostress = p.get<INPAR::STR::StressType>("iostress",INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain = p.get<INPAR::STR::StrainType>("iostrain",INPAR::STR::strain_none);
      vector<double> nodalstress(6);
      vector<double> nodalstrain(6);
      NodalIntegration(NULL,NULL,adjnode,adjele,
                       mis,samepatch,NULL,NULL,ps_adjnode,ps_adjele,
                       &nodalstress,&nodalstrain,iostress,iostrain);
      nodestress_[nodeLid] = nodalstress;
      nodestrain_[nodeLid] = nodalstrain;
    }


    //---------------------- do assembly of stiffness and internal force
    // (note: this is non-standard-baci assembly and therefore a do it all yourself version!)
    if (assemblemat1)
    {
      for (int i=0; i<ndofperpatch; ++i)
      {
        const int rgid = lm[i];
        for (int j=0; j<ndofperpatch; ++j)
        {
          const int cgid = lm[j];
          int errone = stifftmp->SumIntoGlobalValues(1,&rgid,1,&cgid,&stiff(i,j));
          if (errone>0)
          {
            int errtwo = stifftmp->InsertGlobalValues(1,&rgid,1,&cgid,&stiff(i,j));
            if (errtwo<0) dserror("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d",errtwo);
          }
          else if (errone)
            dserror("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d",errone);
        }
      }
      if (mis && !samepatch)
      {
        for (int i=0; i<ps_ndofperpatch; ++i)
        {
          const int rgid = (*ps_lm)[i];
          for (int j=0; j<ps_ndofperpatch; ++j)
          {
            const int cgid = (*ps_lm)[j];
            int errone = stifftmp->SumIntoGlobalValues(1,&rgid,1,&cgid,&((*ps_stiff)(i,j)));
            if (errone>0)
            {
              int errtwo = stifftmp->InsertGlobalValues(1,&rgid,1,&cgid,&((*ps_stiff)(i,j)));
              if (errtwo<0) dserror("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d",errtwo);
            }
            else if (errone)
              dserror("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d",errone);
          }
        }
      }
    }
    
    if (assemblevec1)
    {
      for (int i=0; i<ndofperpatch; ++i)
      {
        const int rgid = lm[i];
        const int lid = forcetmp1.Map().LID(rgid);
        if (lid<0) dserror("global row %d does not exist in column map",rgid);
        forcetmp1[lid] += force[i];
      }
      if (mis && !samepatch)
      {
        for (int i=0; i<ps_ndofperpatch; ++i)
        {
          const int rgid = (*ps_lm)[i];
          const int lid = forcetmp1.Map().LID(rgid);
          if (lid<0) dserror("global row %d does not exist in column map",rgid);
          forcetmp1[lid] += (*ps_force)[i];
        }
      }
    }

  //=========================================================================
  } // for (node=noderids_.begin(); node != noderids_.end(); ++node)


  //-------------------------------------------------------------------------
  // need to export forcetmp to systemvector1 and insert stiffnesses from stifftmp
  // into systemmatrix1
  // Note that fillComplete is never called on stifftmp
  if (assemblevec1)
  {
    Epetra_Vector tmp(systemvector1->Map(),false);
    Epetra_Export exporter(forcetmp1.Map(),tmp.Map());
    int err = tmp.Export(forcetmp1,exporter,Add);
    if (err) dserror("Export using exporter returned err=%d",err);
    systemvector1->Update(1.0,tmp,1.0);
  }
  if (assemblemat1)
  {
    int err = stifftmp->GlobalAssemble(dmap,rmap,false);
    if (err) dserror("Epetra_FECrsMatrix::GlobalAssemble returned err=%d",err);
    const Epetra_Map& cmap = stifftmp->ColMap();
    for (int lrow=0; lrow<stifftmp->NumMyRows(); ++lrow)
    {
      int numentries;
      double* values;
      if (not stifftmp->Filled())
      {
        const int grow = stifftmp->RowMap().GID(lrow);
        int* gindices;
        int err = stifftmp->ExtractGlobalRowView(grow,numentries,values,gindices);
        if (err) dserror("Epetra_FECrsMatrix::ExtractGlobalRowView returned err=%d",err);
        for (int j=0; j<numentries; ++j)
          systemmatrix1->Assemble(values[j],grow,gindices[j]);
      }
      else
      {
        int* lindices;
        int err = stifftmp->ExtractMyRowView(lrow,numentries,values,lindices);
        if (err) dserror("Epetra_FECrsMatrix::ExtractMyRowView returned err=%d",err);
        if (systemmatrix != null and systemmatrix->Filled())
        {
          Epetra_CrsMatrix& matrix = *systemmatrix->EpetraMatrix();
          for (int j=0; j<numentries; ++j)
          {
            int err = matrix.SumIntoMyValues(lrow,1,&values[j],&lindices[j]);
            if (err!=0)
            {
              dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d",err);
            }
          }
        }
        else
        {
          const int grow = stifftmp->RowMap().GID(lrow);
          for (int j=0; j<numentries; ++j)
            systemmatrix1->Assemble(values[j],grow,cmap.GID(lindices[j]));
        }
      }
    }
  }

  if (action == "calc_struct_stress")
  {
    // we have to export the nodal stresses and strains to column map
    // so they can be written by the elements
    DRT::Exporter exporter(*dis.NodeRowMap(),*dis.NodeColMap(),dis.Comm());
    exporter.Export<double>(nodestress_);
    exporter.Export<double>(nodestrain_);
  }


  return;
}

/*----------------------------------------------------------------------*
 |  do nodal integration (public)                              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::NodalIntegration(Epetra_SerialDenseMatrix*       stiff,
                                                Epetra_SerialDenseVector*       force,
                                                map<int,DRT::Node*>&            adjnode,
                                                vector<DRT::ELEMENTS::NStet*>&  adjele,
                                                bool                            mis,
                                                bool                            samepatch,
                                                Epetra_SerialDenseMatrix*       ps_stiff,
                                                Epetra_SerialDenseVector*       ps_force,
                                                map<int,DRT::Node*>*            ps_adjnode,
                                                vector<DRT::ELEMENTS::NStet*>*  ps_adjele,
                                                vector<double>*                 nodalstress,
                                                vector<double>*                 nodalstrain,
                                                const INPAR::STR::StressType    iostress,
                                                const INPAR::STR::StrainType    iostrain)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::NodalIntegration");

  //-------------------------------------------------- standard quantities
  const int nnodeinpatch = (int)adjnode.size();
  const int ndofinpatch  = nnodeinpatch * 3;
  const int neleinpatch  = (int)adjele.size();

  // MIS quantities ------------------------------------------------------
  Epetra_SerialDenseMatrix ps_bop;
  LINALG::Matrix<3,3>      ps_FnodeL(true);
  LINALG::Matrix<3,3>      ps_cauchygreen;
  LINALG::Matrix<6,1>      ps_glstrain(false);
  LINALG::Matrix<6,6>      ps_cmat(true);
  LINALG::Matrix<6,1>      ps_stress(true);
  double                   ps_VnodeL = 0.0;
  int                      ps_nnodeinpatch = nnodeinpatch;
  int                      ps_neleinpatch = neleinpatch;
  if (mis && !samepatch) 
  {
    ps_nnodeinpatch = (*ps_adjnode).size() * 3;
    ps_neleinpatch  = (*ps_adjele).size();
  }
  int                      ps_ndofinpatch = ps_nnodeinpatch * 3;
  bool                     ps_matequal = true;

  //------------------------------ see whether materials in patch are equal
  bool matequal = true;
  {
    int mat = adjele[0]->material_;
    for (int i=1; i<neleinpatch; ++i)
      if (mat != adjele[i]->material_)
      {
        matequal = false;
        break;
      }
  }

  //-----------------------------------------------------------------------
  // MIS see whether materials are equal in patch
  if (mis && !samepatch)
  {
    int mat = (*ps_adjele)[0]->material_;
    for (int i=1; i<ps_neleinpatch; ++i)
      if (mat != (*ps_adjele)[i]->material_)
      {
        ps_matequal = false;
        break;
      }
  }
  else if (mis) ps_matequal = matequal;

  //-----------------------------------------------------------------------
  // build averaged deformation gradient and volume of node
  LINALG::Matrix<3,3> FnodeL(true);
  double VnodeL = 0.0;
  for (int i=0; i<neleinpatch; ++i)
  {
    const double V = adjele[i]->Volume()/4;
    VnodeL += V;
    FnodeL.Update(V,adjele[i]->F_,1.0);
  }
  FnodeL.Scale(1.0/VnodeL);

  //-----------------------------------------------------------------------
  // MIS build averaged deformation gradient and volume of node
  if (mis)
  {
    // MIS node takes entire element volume
    for (int i=0; i<ps_neleinpatch; ++i) 
    {
      const double V = (*ps_adjele)[i]->Volume();
      ps_VnodeL += V;
      if (!samepatch) ps_FnodeL.Update(V,(*ps_adjele)[i]->F_,1.0);
    }
    if (!samepatch) ps_FnodeL.Scale(1.0/ps_VnodeL);
    else            ps_FnodeL = FnodeL;
  }

  //-----------------------------------------------------------------------
  // do positioning map global nodes -> position in B-Operator
  map<int,int>  node_pos;
  {
    std::map<int,DRT::Node*>::iterator pnode;
    int count=0;
    for (pnode=adjnode.begin(); pnode != adjnode.end(); ++pnode)
    {
      node_pos[pnode->first] = count;
      count++;
    }
  }
  
  //-----------------------------------------------------------------------
  // MIS do positioning map global nodes -> position in B-Operator
  map<int,int>  ps_node_pos;
  if (mis && !samepatch)
  {
    std::map<int,DRT::Node*>::iterator pnode;
    int count=0;
    for (pnode=(*ps_adjnode).begin(); pnode != (*ps_adjnode).end(); ++pnode)
    {
      ps_node_pos[pnode->first] = count;
      count++;
    }
  }
  else if (mis) ps_node_pos = node_pos;

  //-----------------------------------------------------------------------
  // build B operator

  Epetra_SerialDenseMatrix bop(6,ndofinpatch);
  // loop elements in patch
  for (int ele=0; ele<neleinpatch; ++ele)
  {
    // current element
    DRT::ELEMENTS::NStet* actele = adjele[ele];

    // spatial deriv of that element
    LINALG::Matrix<4,3>& nxyz = actele->nxyz_;

    // volume of that element assigned to node L
    double V = actele->Volume()/4;

    // def-gradient of the element
    LINALG::Matrix<3,3>& F = actele->F_;
    
    // volume ratio of volume per node L of this element to
    // whole volume of node L
    V = V/VnodeL;

    // loop nodes of that element
    for (int i=0; i<actele->NumNode(); ++i)
    {
      DRT::Node* actnode = actele->Nodes()[i];
      const int  nodeid  = actnode->Id();

      // local  node index is i
      // global node index is nodeid
      // starting position in B-Operator is node_pos[nodeid]

      // find position in map of that node to determine place in bop
      int pos = node_pos[nodeid];

      bop(0,3*pos+0) += V * F(0,0)*nxyz(i,0);
      bop(0,3*pos+1) += V * F(1,0)*nxyz(i,0);
      bop(0,3*pos+2) += V * F(2,0)*nxyz(i,0);
      bop(1,3*pos+0) += V * F(0,1)*nxyz(i,1);
      bop(1,3*pos+1) += V * F(1,1)*nxyz(i,1);
      bop(1,3*pos+2) += V * F(2,1)*nxyz(i,1);
      bop(2,3*pos+0) += V * F(0,2)*nxyz(i,2);
      bop(2,3*pos+1) += V * F(1,2)*nxyz(i,2);
      bop(2,3*pos+2) += V * F(2,2)*nxyz(i,2);
      //
      bop(3,3*pos+0) += V * (F(0,0)*nxyz(i,1) + F(0,1)*nxyz(i,0) );
      bop(3,3*pos+1) += V * (F(1,0)*nxyz(i,1) + F(1,1)*nxyz(i,0) );
      bop(3,3*pos+2) += V * (F(2,0)*nxyz(i,1) + F(2,1)*nxyz(i,0) );
      bop(4,3*pos+0) += V * (F(0,1)*nxyz(i,2) + F(0,2)*nxyz(i,1) );
      bop(4,3*pos+1) += V * (F(1,1)*nxyz(i,2) + F(1,2)*nxyz(i,1) );
      bop(4,3*pos+2) += V * (F(2,1)*nxyz(i,2) + F(2,2)*nxyz(i,1) );
      bop(5,3*pos+0) += V * (F(0,2)*nxyz(i,0) + F(0,0)*nxyz(i,2) );
      bop(5,3*pos+1) += V * (F(1,2)*nxyz(i,0) + F(1,0)*nxyz(i,2) );
      bop(5,3*pos+2) += V * (F(2,2)*nxyz(i,0) + F(2,0)*nxyz(i,2) );

    } // for (int i=0; i<actele->NumNode(); ++i)

  } // for (int ele=0; ele<neleinpatch; ++ele)

  //-----------------------------------------------------------------------
  // MIS build B operator 

  if (mis && !samepatch)
  {
    ps_bop.Shape(6,ps_ndofinpatch);
    
    // loop elements in MIS patch
    for (int ele=0; ele<ps_neleinpatch; ++ele)
    {
      // current element
      DRT::ELEMENTS::NStet* actele = (*ps_adjele)[ele];

      // spatial deriv of that element
      LINALG::Matrix<4,3>& nxyz = actele->nxyz_;

      // entire volume of that element assigned to MIS node L
      double V = actele->Volume();

      // def-gradient of the element
      LINALG::Matrix<3,3>& F = actele->F_;

      // volume ratio of volume per node L of this element to
      // whole volume of node L
      double ratio = V/ps_VnodeL;

      // loop nodes of that element
      for (int i=0; i<actele->NumNode(); ++i)
      {
        DRT::Node* actnode = actele->Nodes()[i];
        const int  nodeid  = actnode->Id();

        // local  node index is i
        // global node index is nodeid
        // starting position in B-Operator is node_pos[nodeid]

        // find position in map of that node to determine place in bop
        int pos = ps_node_pos[nodeid];

        ps_bop(0,3*pos+0) += ratio * F(0,0)*nxyz(i,0);
        ps_bop(0,3*pos+1) += ratio * F(1,0)*nxyz(i,0);
        ps_bop(0,3*pos+2) += ratio * F(2,0)*nxyz(i,0);
        ps_bop(1,3*pos+0) += ratio * F(0,1)*nxyz(i,1);
        ps_bop(1,3*pos+1) += ratio * F(1,1)*nxyz(i,1);
        ps_bop(1,3*pos+2) += ratio * F(2,1)*nxyz(i,1);
        ps_bop(2,3*pos+0) += ratio * F(0,2)*nxyz(i,2);
        ps_bop(2,3*pos+1) += ratio * F(1,2)*nxyz(i,2);
        ps_bop(2,3*pos+2) += ratio * F(2,2)*nxyz(i,2);
        //
        ps_bop(3,3*pos+0) += ratio * (F(0,0)*nxyz(i,1) + F(0,1)*nxyz(i,0) );
        ps_bop(3,3*pos+1) += ratio * (F(1,0)*nxyz(i,1) + F(1,1)*nxyz(i,0) );
        ps_bop(3,3*pos+2) += ratio * (F(2,0)*nxyz(i,1) + F(2,1)*nxyz(i,0) );
        ps_bop(4,3*pos+0) += ratio * (F(0,1)*nxyz(i,2) + F(0,2)*nxyz(i,1) );
        ps_bop(4,3*pos+1) += ratio * (F(1,1)*nxyz(i,2) + F(1,2)*nxyz(i,1) );
        ps_bop(4,3*pos+2) += ratio * (F(2,1)*nxyz(i,2) + F(2,2)*nxyz(i,1) );
        ps_bop(5,3*pos+0) += ratio * (F(0,2)*nxyz(i,0) + F(0,0)*nxyz(i,2) );
        ps_bop(5,3*pos+1) += ratio * (F(1,2)*nxyz(i,0) + F(1,0)*nxyz(i,2) );
        ps_bop(5,3*pos+2) += ratio * (F(2,2)*nxyz(i,0) + F(2,0)*nxyz(i,2) );

      } // for (int i=0; i<actele->NumNode(); ++i)

    } // for (int ele=0; ele<neleinpatch; ++ele)
  }
  else if (mis)
  {
    ps_bop.Shape(6,ps_ndofinpatch);
    ps_bop = bop;
  }

  //----------------------------------------- averaged material and stresses
  LINALG::Matrix<6,6> cmat(true);
  LINALG::Matrix<6,1> stress(true);

  //----------------------------------------------------------------- strain
  // right cauchy green
  LINALG::Matrix<3,3> cauchygreen;
  cauchygreen.MultiplyTN(FnodeL,FnodeL);
  // Green-Lagrange ( 2x on offdiagonal!)
  LINALG::Matrix<6,1> glstrain(false);
  glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
  glstrain(3) = cauchygreen(0,1);
  glstrain(4) = cauchygreen(1,2);
  glstrain(5) = cauchygreen(2,0);

  //-----------------------------------------------------------------------
  // MIS strain
  if (mis && !samepatch)
  {
    ps_cauchygreen.MultiplyTN(ps_FnodeL,ps_FnodeL);
    ps_glstrain(0) = 0.5 * (ps_cauchygreen(0,0) - 1.0);
    ps_glstrain(1) = 0.5 * (ps_cauchygreen(1,1) - 1.0);
    ps_glstrain(2) = 0.5 * (ps_cauchygreen(2,2) - 1.0);
    ps_glstrain(3) = ps_cauchygreen(0,1);
    ps_glstrain(4) = ps_cauchygreen(1,2);
    ps_glstrain(5) = ps_cauchygreen(2,0);
  }
  else if (mis) 
  {
    ps_cauchygreen = cauchygreen;
    ps_glstrain    = glstrain;
  }

  //-------------------------------------------------------- output of strain
  // FIXME: MIS not considered yet
  switch (iostrain)
  {
  case INPAR::STR::strain_gl:
  {
    if (nodalstrain == NULL) dserror("no strain data available");
    for (int i = 0; i < 3; ++i) (*nodalstrain)[i] = glstrain(i);
    for (int i = 3; i < 6; ++i) (*nodalstrain)[i] = 0.5 * glstrain(i);
  }
  break;
  case INPAR::STR::strain_ea:
  {
    if (nodalstrain == NULL) dserror("no strain data available");

    // rewriting Green-Lagrange strains in matrix format
    LINALG::Matrix<3,3> gl;
    gl(0,0) = glstrain(0);
    gl(0,1) = 0.5*glstrain(3);
    gl(0,2) = 0.5*glstrain(5);
    gl(1,0) = gl(0,1);
    gl(1,1) = glstrain(1);
    gl(1,2) = 0.5*glstrain(4);
    gl(2,0) = gl(0,2);
    gl(2,1) = gl(1,2);
    gl(2,2) = glstrain(2);

    // inverse of deformation gradient
    LINALG::Matrix<3,3> invdefgrd;
    invdefgrd.Invert(FnodeL);

    LINALG::Matrix<3,3> temp;
    LINALG::Matrix<3,3> euler_almansi;
    temp.Multiply(gl,invdefgrd);
    euler_almansi.MultiplyTN(invdefgrd,temp);

    (*nodalstrain)[0] = euler_almansi(0,0);
    (*nodalstrain)[1] = euler_almansi(1,1);
    (*nodalstrain)[2] = euler_almansi(2,2);
    (*nodalstrain)[3] = euler_almansi(0,1);
    (*nodalstrain)[4] = euler_almansi(1,2);
    (*nodalstrain)[5] = euler_almansi(0,2);
  }
  break;
  case INPAR::STR::strain_none:
    break;
  default:
    dserror("requested strain type not available");
  }

  //-----------------------------------------------------------------------
  // material law and stresses
  if (matequal) // element patch has single material
  {
    double density; // just a dummy density
    RCP<MAT::Material> mat = adjele[0]->Material();
    SelectMaterial(mat,stress,cmat,density,glstrain,FnodeL,0);
  }
  else
  {
    double density; // just a dummy density
    LINALG::Matrix<6,6> cmatele;
    LINALG::Matrix<6,1> stressele;
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      cmatele = 0.0;
      stressele = 0.0;
      // current element
      DRT::ELEMENTS::NStet* actele = adjele[ele];
      // volume of that element assigned to node L
      const double V = actele->Volume()/4;
      // def-gradient of the element
      RCP<MAT::Material> mat = actele->Material();
      SelectMaterial(mat,stressele,cmatele,density,glstrain,FnodeL,0);
      cmat.Update(V,cmatele,1.0);
      stress.Update(V,stressele,1.0);
    } // for (int ele=0; ele<neleinpatch; ++ele)
    stress.Scale(1.0/VnodeL);
    cmat.Scale(1.0/VnodeL);
  }

  //-----------------------------------------------------------------------
  // MIS material law and stresses
  if (mis && !samepatch)
  {
    if (ps_matequal)
    {
      double density; // just a dummy density
      RCP<MAT::Material> mat = (*ps_adjele)[0]->Material();
      SelectMaterial(mat,ps_stress,ps_cmat,density,ps_glstrain,ps_FnodeL,0);
    }
    else
    {
      double density; // just a dummy density
      LINALG::Matrix<6,6> cmatele;
      LINALG::Matrix<6,1> stressele;
      for (int ele=0; ele<ps_neleinpatch; ++ele)
      {
        cmatele = 0.0;
        stressele = 0.0;
        // current element
        DRT::ELEMENTS::NStet* actele = (*ps_adjele)[ele];
        // volume of that element assigned to node L
        const double V = actele->Volume();
        // def-gradient of the element
        RCP<MAT::Material> mat = actele->Material();
        SelectMaterial(mat,stressele,cmatele,density,ps_glstrain,ps_FnodeL,0);
        ps_cmat.Update(V,cmatele,1.0);
        ps_stress.Update(V,stressele,1.0);
      } // for (int ele=0; ele<neleinpatch; ++ele)
      ps_stress.Scale(1.0/ps_VnodeL);
      ps_cmat.Scale(1.0/ps_VnodeL);
    }
  }
  else if (mis)
  {
    ps_stress = stress;
    ps_cmat = cmat;
  }

  //-----------------------------------------------------------------------
  // FIXME MIS not considered yet
  switch (iostress)
  {
  case INPAR::STR::stress_2pk:
  {
    if (nodalstress == NULL) dserror("no stress data available");
    for (int i = 0; i < 6; ++i) (*nodalstress)[i] = stress(i);
  }
  break;
  case INPAR::STR::stress_cauchy:
  {
    if (nodalstress == NULL) dserror("no stress data available");
    double detF = FnodeL.Determinant();

    LINALG::Matrix<3,3> pkstress;
    pkstress(0,0) = stress(0);
    pkstress(0,1) = stress(3);
    pkstress(0,2) = stress(5);
    pkstress(1,0) = pkstress(0,1);
    pkstress(1,1) = stress(1);
    pkstress(1,2) = stress(4);
    pkstress(2,0) = pkstress(0,2);
    pkstress(2,1) = pkstress(1,2);
    pkstress(2,2) = stress(2);

    LINALG::Matrix<3,3> temp;
    LINALG::Matrix<3,3> cauchystress;
    temp.Multiply(1.0/detF,FnodeL,pkstress);
    cauchystress.MultiplyNT(temp,FnodeL);

    (*nodalstress)[0] = cauchystress(0,0);
    (*nodalstress)[1] = cauchystress(1,1);
    (*nodalstress)[2] = cauchystress(2,2);
    (*nodalstress)[3] = cauchystress(0,1);
    (*nodalstress)[4] = cauchystress(1,2);
    (*nodalstress)[5] = cauchystress(0,2);
  }
  break;
  case INPAR::STR::stress_none:
    break;
  default:
    dserror("requested stress type not available");
  }

  //-----------------------------------------------------------------------
  // stress is plit as follows:
  // non-MIS node:
  //     stress = (1-beta) * vol_node + (1-alpha) * dev_node + alpha * dev_ele
  // MIS-node:
  //     stress = beta * vol_node(newpatch) + (1-beta) * vol_node + (1-alpha) * dev_node + alpha * dev_ele

  // define stuff we need in all cases
  LINALG::Matrix<6,6> cmatdev;
  LINALG::Matrix<6,6> cmatvol;
  LINALG::Matrix<6,1> stressdev;
  LINALG::Matrix<6,1> stressvol;
  
  // all nodes
  {
    // compute deviatoric stress and tangent from total stress and tangent
    DevStressTangent(stressdev,cmatdev,cmat,stress,cauchygreen);
    
    // compute volumetric stress and tangent
    stressvol.Update(-1.0,stressdev,1.0,stress,0.0);
    cmatvol.Update(-1.0,cmatdev,1.0,cmat,0.0);
    
    // compute nodal stress
    stress.Update(1-BETA_NSTET,stressvol,1-ALPHA_NSTET,stressdev,0.0);
    cmat.Update(1-BETA_NSTET,cmatvol,1-ALPHA_NSTET,cmatdev,0.0);
  }
  // MIS nodes
  if (mis)
  {
    // patch (and therefore stress/tangent) identical to standard patch
    if (samepatch)
    {
      ps_stress.Update(BETA_NSTET,stressvol,0.0);
      ps_cmat.Update(BETA_NSTET,cmatvol,0.0);
    }
    // patch (and therefore stress/tangent) different from standard patch
    else
    {
      // define stuff we need
      cmatdev = 0.0;
      stressdev = 0.0;

      // compute MIS deviatoric stress and tangent from MIS total stress and tangent
      DevStressTangent(stressdev,cmatdev,ps_cmat,ps_stress,ps_cauchygreen);
      
      // compute MIS volumetric stress and tangent
      stressvol.Update(-1.0,stressdev,1.0,ps_stress,0.0);
      cmatvol.Update(-1.0,cmatdev,1.0,ps_cmat,0.0);
      
      ps_stress.Update(BETA_NSTET,stressvol,0.0);
      ps_cmat.Update(BETA_NSTET,cmatvol,0.0);
    }
  }

  //----------------------------------------------------- internal forces
  if (force)
  {
    Epetra_SerialDenseVector stress_epetra(View,stress.A(),stress.Rows());
    force->Multiply('T','N',VnodeL,bop,stress_epetra,0.0);
    //------- MIS node part
    if (mis)
    {
      Epetra_SerialDenseVector ps_stress_epetra(View,ps_stress.A(),ps_stress.Rows());
      if (samepatch)
        ps_force->Multiply('T','N',ps_VnodeL,ps_bop,ps_stress_epetra,1.0);
      else
        ps_force->Multiply('T','N',ps_VnodeL,ps_bop,ps_stress_epetra,0.0);
    }
  }

  //--------------------------------------------------- elastic stiffness
  if (stiff)
  {
    Epetra_SerialDenseMatrix cmat_epetra(View,cmat.A(),cmat.Rows(),cmat.Rows(),cmat.Columns());
    LINALG::SerialDenseMatrix cb(6,ndofinpatch);
    cb.Multiply('N','N',1.0,cmat_epetra,bop,0.0);
    stiff->Multiply('T','N',VnodeL,bop,cb,0.0);
    //------- MIS node part
    if (mis)
    {
      Epetra_SerialDenseMatrix ps_cmat_epetra(View,ps_cmat.A(),ps_cmat.Rows(),ps_cmat.Rows(),ps_cmat.Columns());
      LINALG::SerialDenseMatrix ps_cb(6,ps_ndofinpatch);
      ps_cb.Multiply('N','N',1.0,ps_cmat_epetra,ps_bop,0.0);
      if (samepatch)
        ps_stiff->Multiply('T','N',ps_VnodeL,ps_bop,ps_cb,1.0);
      else
        ps_stiff->Multiply('T','N',ps_VnodeL,ps_bop,ps_cb,0.0);
    }
  }

  //----------------------------------------------------- geom. stiffness
  if (stiff)
  {
    // loop elements in patch
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      // current element
      DRT::ELEMENTS::NStet* actele = adjele[ele];
      // material deriv of that element
      LINALG::Matrix<4,3>& nxyz   = actele->nxyz_;
      // volume of actele assigned to node L
      double V = actele->Volume()/4;
      
      // loop nodes of that element
      double SmBL[3];
      DRT::Node** nodes = actele->Nodes();
      for (int i=0; i<4; ++i)
      {
        // row position of this node in matrix
        int ipos = node_pos[nodes[i]->Id()];
        SmBL[0] = V*(stress(0)*nxyz(i,0) + stress(3)*nxyz(i,1) + stress(5)*nxyz(i,2));
        SmBL[1] = V*(stress(3)*nxyz(i,0) + stress(1)*nxyz(i,1) + stress(4)*nxyz(i,2));
        SmBL[2] = V*(stress(5)*nxyz(i,0) + stress(4)*nxyz(i,1) + stress(2)*nxyz(i,2));
        for (int j=0; j<4; ++j)
        {
          // column position of this node in matrix
          int jpos = node_pos[nodes[j]->Id()];
          double bopstrbop = 0.0;
          for (int dim=0; dim<3; ++dim) bopstrbop += nxyz(j,dim) * SmBL[dim];
          (*stiff)(3*ipos+0,3*jpos+0) += bopstrbop;
          (*stiff)(3*ipos+1,3*jpos+1) += bopstrbop;
          (*stiff)(3*ipos+2,3*jpos+2) += bopstrbop;
        } // for (int j=0; j<4; ++j)
      } // for (int i=0; i<4; ++i)
    } // for (int ele=0; ele<neleinpatch; ++ele)
    
    //------- MIS node part
    if (mis && samepatch)
    {
      for (int ele=0; ele<ps_neleinpatch; ++ele)
      {
        // current element
        DRT::ELEMENTS::NStet* actele = (*ps_adjele)[ele];
        // material deriv of that element
        LINALG::Matrix<4,3>& nxyz   = actele->nxyz_;
        // all volume of actele assigned to MIS node L
        double V = actele->Volume();

        // loop nodes of that element
        double SmBL[3];
        DRT::Node** nodes = actele->Nodes();
        for (int i=0; i<4; ++i)
        {
          // row position of this node in matrix
          int ipos = ps_node_pos[nodes[i]->Id()];
          SmBL[0] = V*(ps_stress(0)*nxyz(i,0) + ps_stress(3)*nxyz(i,1) + ps_stress(5)*nxyz(i,2));
          SmBL[1] = V*(ps_stress(3)*nxyz(i,0) + ps_stress(1)*nxyz(i,1) + ps_stress(4)*nxyz(i,2));
          SmBL[2] = V*(ps_stress(5)*nxyz(i,0) + ps_stress(4)*nxyz(i,1) + ps_stress(2)*nxyz(i,2));
          for (int j=0; j<4; ++j)
          {
            // column position of this node in matrix
            int jpos = ps_node_pos[nodes[j]->Id()];
            double bopstrbop = 0.0;
            for (int dim=0; dim<3; ++dim) bopstrbop += nxyz(j,dim) * SmBL[dim];
            (*ps_stiff)(3*ipos+0,3*jpos+0) += bopstrbop;
            (*ps_stiff)(3*ipos+1,3*jpos+1) += bopstrbop;
            (*ps_stiff)(3*ipos+2,3*jpos+2) += bopstrbop;
          } // for (int j=0; j<4; ++j)
        } // for (int i=0; i<4; ++i)
      } // for (int ele=0; ele<ps_neleinpatch; ++ele)
    } // if (mis)
  } // if (stiff)

  return;
}

/*----------------------------------------------------------------------*
 | material laws for NStet (protected)                          gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::SelectMaterial(
                      RCP<MAT::Material> mat,
                      LINALG::Matrix<6,1>& stress,
                      LINALG::Matrix<6,6>& cmat,
                      double& density,
                      LINALG::Matrix<6,1>& glstrain,
                      LINALG::Matrix<3,3>& defgrd,
                      int gp)
{
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(glstrain,cmat,stress);
      density = stvk->Density();
    }
    break;
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain,cmat,stress);
      density = neo->Density();
    }
    break;
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast<MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(glstrain,cmat,stress);
      density = aaa->Density();
    }
    break;
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(&glstrain,&cmat,&stress);
      density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen = static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(&glstrain,&cmat,&stress);

      density = lungpen->Density();
      return;
      break;
    }
    default:
      dserror("Illegal type %d of material for element NStet tet4", mat->MaterialType());
    break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::NStet::SelectMaterial

/*----------------------------------------------------------------------*
 |  compute deviatoric tangent and stresses (private/static)   gee 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::DevStressTangent(
  LINALG::Matrix<6,1>& Sdev,
  LINALG::Matrix<6,6>& CCdev,
  LINALG::Matrix<6,6>& CC,
  const LINALG::Matrix<6,1>& S,
  const LINALG::Matrix<3,3>& C)
{

  //---------------------------------- things that we'll definitely need
  // inverse of C
  LINALG::Matrix<3,3> Cinv;
  const double detC = Cinv.Invert(C);

  // J = det(F) = sqrt(detC)
  const double J = sqrt(detC);

  // S as a 3x3 matrix
  LINALG::Matrix<3,3> Smat;
  Smat(0,0) = S(0);
  Smat(0,1) = S(3);
  Smat(0,2) = S(5);
  Smat(1,0) = Smat(0,1);
  Smat(1,1) = S(1);
  Smat(1,2) = S(4);
  Smat(2,0) = Smat(0,2);
  Smat(2,1) = Smat(1,2);
  Smat(2,2) = S(2);

  //--------------------------------------------- pressure p = -1/(3J) S:C
  double p = 0.0;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      p += Smat(i,j)*C(i,j);
  p *= (-1./(3.*J));

  //-------------------------------- compute volumetric PK2 Svol = -p J Cinv
  //-------------------------------------------------------- Sdev = S - Svol
  const double fac = -p*J;
  Sdev(0) = Smat(0,0) - fac*Cinv(0,0);
  Sdev(1) = Smat(1,1) - fac*Cinv(1,1);
  Sdev(2) = Smat(2,2) - fac*Cinv(2,2);
  Sdev(3) = Smat(0,1) - fac*Cinv(0,1);
  Sdev(4) = Smat(1,2) - fac*Cinv(1,2);
  Sdev(5) = Smat(0,2) - fac*Cinv(0,2);

  //======================================== volumetric tangent matrix CCvol
  LINALG::Matrix<6,6> CCvol(true); // fill with zeros

  //--------------------------------------- CCvol += 2pJ (Cinv boeppel Cinv)
  MAT::ElastSymTensor_o_Multiply(CCvol,-2.0*fac,Cinv,Cinv,0.0);

  //------------------------------------------ CCvol += 2/3 * Cinv dyad S
  MAT::ElastSymTensorMultiply(CCvol,2.0/3.0,Cinv,Smat,1.0);

  //-------------------------------------- CCvol += 1/3 Cinv dyad ( CC : C )
  {
    // C as Voigt vector
    LINALG::Matrix<6,1> Cvec;
    Cvec(0) = C(0,0);
    Cvec(1) = C(1,1);
    Cvec(2) = C(2,2);
    Cvec(3) = 2.0*C(0,1);
    Cvec(4) = 2.0*C(1,2);
    Cvec(5) = 2.0*C(0,2);

    LINALG::Matrix<6,1> CCcolonC;
    CCcolonC.Multiply(CC,Cvec);

    LINALG::Matrix<3,3> CCcC;
    CCcC(0,0) = CCcolonC(0);
    CCcC(0,1) = CCcolonC(3);
    CCcC(0,2) = CCcolonC(5);
    CCcC(1,0) = CCcC(0,1);
    CCcC(1,1) = CCcolonC(1);
    CCcC(1,2) = CCcolonC(4);
    CCcC(2,0) = CCcC(0,2);
    CCcC(2,1) = CCcC(1,2);
    CCcC(2,2) = CCcolonC(2);
    MAT::ElastSymTensorMultiply(CCvol,1./3.,Cinv,CCcC,1.0);
  }

  //----------------------------------------------------- CCdev = CC - CCvol
  CCdev.Update(1.0,CC,-1.0,CCvol);

  return;
}



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
