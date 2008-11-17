/*!----------------------------------------------------------------------*
\file so_ptet_nodalstrain.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_ptet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_FECrsMatrix.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/hyperpolyconvex.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PtetRegister::Initialize(DRT::Discretization& dis)
{
  const int myrank = dis.Comm().MyPID();

  const int numele = dis.NumMyColElements();
  for (int i=0; i<numele; ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_ptet) continue;

    DRT::ELEMENTS::Ptet* actele = dynamic_cast<DRT::ELEMENTS::Ptet*>(dis.lColElement(i));
    if (!actele) dserror("cast to Ptet* failed");

    // init the element (also set pointer to this register class)
    actele->InitElement(this);

    // register element in list of column Ptet elements
    elecids_[actele->Id()] = actele;

    // compute a map of all row nodes adjacent to a Ptet element
    for (int j=0; j<actele->NumNode(); ++j) {
      DRT::Node* node = actele->Nodes()[j];
      if (myrank == node->Owner())
        noderids_.insert(pair<int,DRT::Node*>(node->Id(),node));
    }
  }

#if 1 // efficiency (at the price of memory)
  // compute adjacency for each row node
  std::map<int,DRT::Node*>::iterator node;
  for (node=noderids_.begin(); node != noderids_.end(); ++node)
  {
    DRT::Node* nodeL  = node->second;
    const int nodeidL = nodeL->Id();

    // list of adjacent elements
    vector<DRT::ELEMENTS::Ptet*> adjele(0);
    for (int j=0; j<nodeL->NumElement(); ++j)
    {
      const int eleid = node->second->Elements()[j]->Id();
      std::map<int,DRT::ELEMENTS::Ptet*>::iterator ele = elecids_.find(eleid);
      if (ele==elecids_.end()) continue;
      adjele.push_back(ele->second);
    }
    adjele_[nodeidL] = adjele;

    // patch of all nodes adjacent to adjacent elements
    map<int,DRT::Node*> nodepatch;
    for (int j=0; j<(int)adjele.size(); ++j) {
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
      vector<int> dofs = dis.Dof(pnode->second);
      for (int j=0; j<(int)dofs.size(); ++j)
        lm[count++]        = dofs[j];
    }
    if (count != ndofperpatch) dserror("dimension mismatch");
    adjlm_[nodeidL] = lm;

  } // for (node=noderids_.begin(); node != noderids_.end(); ++node)
#endif

  return 0;
}


/*----------------------------------------------------------------------*
 |  pre-evaluation of elements (public)                        gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PtetRegister::PreEvaluate(DRT::Discretization& dis,
                                              Teuchos::ParameterList& p,
                                              RCP<LINALG::SparseOperator> systemmatrix1,
                                              RCP<LINALG::SparseOperator> systemmatrix2,
                                              RCP<Epetra_Vector>          systemvector1,
                                              RCP<Epetra_Vector>          systemvector2,
                                              RCP<Epetra_Vector>          systemvector3)
{
  // nodal integration for nlnstiff and internal forces only
  // (this method does not compute stresses/strains/element updates)
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
  LINALG::SerialDenseVector force1;

  //-------------------------------------- construct F for each Ptet
  // current displacement
  RCP<const Epetra_Vector> disp = dis.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  // loop elements
  std::map<int,DRT::ELEMENTS::Ptet*>::iterator ele;
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
  Epetra_FECrsMatrix stifftmp(Copy,rmap,256,false);
  // create temporary vector in column map to assemble to
  Epetra_Vector forcetmp1(*dis.DofColMap(),true);

  //-------------------------------------------------do nodal stiffness
  std::map<int,DRT::Node*>::iterator node;
  for (node=noderids_.begin(); node != noderids_.end(); ++node)
  {
    DRT::Node* nodeL   = node->second;     // row node
    const int  nodeLid = nodeL->Id();

    // list of adjacent elements
#if 0
    vector<DRT::ELEMENTS::Ptet*> adjele(0);
    for (int j=0; j<nodeL->NumElement(); ++j)
    {
      const int eleid = node->second->Elements()[j]->Id();
      std::map<int,DRT::ELEMENTS::Ptet*>::iterator ele = elecids_.find(eleid);
      if (ele==elecids_.end()) continue;
      adjele.push_back(ele->second);
    }
#else
    vector<DRT::ELEMENTS::Ptet*>& adjele = adjele_[nodeLid];
#endif


#if 0
    // patch of all nodes adjacent to adjacent elements
    map<int,DRT::Node*> nodepatch;
    nodepatch.clear();
    for (int j=0; j<(int)adjele.size(); ++j)
      for (int k=0; k<adjele[j]->NumNode(); ++k)
        nodepatch[adjele[j]->Nodes()[k]->Id()] = adjele[j]->Nodes()[k];
#else
    map<int,DRT::Node*>& nodepatch = adjnode_[nodeLid];
#endif

    // total number of nodes
    const int numnodepatch = (int)nodepatch.size();
    // total number of degrees of freedom on patch
    const int ndofperpatch = numnodepatch*3;

#if 0
    // location and ownership vector of nodal patch
    vector<int> lm(ndofperpatch);
    std::map<int,DRT::Node*>::iterator pnode;
    int count=0;
    for (pnode=nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
    {
      vector<int> dofs = dis.Dof(pnode->second);
      for (int j=0; j<(int)dofs.size(); ++j)
        lm[count++]        = dofs[j];
    }
    if (count != ndofperpatch) dserror("dimension mismatch");
#else
    vector<int>& lm      = adjlm_[nodeLid];
#endif

    if (action != "calc_struct_stress")
    {
      // do nodal integration of stiffness and internal force
      stiff.LightShape(ndofperpatch,ndofperpatch);
      force1.LightSize(ndofperpatch);
      NodalIntegration(&stiff,&force1,nodepatch,adjele,NULL,NULL,false,false);
    }
    else
    {
      bool cauchy     = p.get<bool>("cauchy",false);
      string iostrain = p.get<string>("iostrain","none");
      bool ea         = (iostrain=="euler_almansi");
      vector<double> nodalstress(6);
      vector<double> nodalstrain(6);
      NodalIntegration(NULL,NULL,nodepatch,adjele,&nodalstress,&nodalstrain,cauchy,ea);
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
          int errone = stifftmp.SumIntoGlobalValues(1,&rgid,1,&cgid,&stiff(i,j));
          if (errone>0)
          {
            int errtwo = stifftmp.InsertGlobalValues(1,&rgid,1,&cgid,&stiff(i,j));
            if (errtwo<0) dserror("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d",errtwo);
          }
          else if (errone)
            dserror("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d",errone);
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
        forcetmp1[lid] += force1[i];
      }
    }

  //---------------------------------------------------------------------
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
    int err = stifftmp.GlobalAssemble(dmap,rmap,false);
    if (err) dserror("Epetra_FECrsMatrix::GlobalAssemble returned err=%d",err);
    for (int lrow=0; lrow<stifftmp.NumMyRows(); ++lrow)
    {
      const int grow = stifftmp.RowMap().GID(lrow);
      int numentries;
      double* values;
      int* gindices;
      int err = stifftmp.ExtractGlobalRowView(grow,numentries,values,gindices);
      if (err) dserror("Epetra_FECrsMatrix::ExtractGlobalRowView returned err=%d",err);
      for (int j=0; j<numentries; ++j)
        systemmatrix1->Assemble(values[j],grow,gindices[j]);
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
void DRT::ELEMENTS::PtetRegister::NodalIntegration(Epetra_SerialDenseMatrix*     stiff,
                                                   Epetra_SerialDenseVector*     force,
                                                   map<int,DRT::Node*>&          nodepatch,
                                                   vector<DRT::ELEMENTS::Ptet*>& adjele,
                                                   vector<double>*               nodalstress,
                                                   vector<double>*               nodalstrain,
                                                   bool                          cauchy,
                                                   bool                          ea)
{
  const int nnodeinpatch = (int)nodepatch.size();
  const int ndofinpatch  = nnodeinpatch*3;
  const int neleinpatch  = (int)adjele.size();

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
  // build averaged deformation gradient and volume of node
  // this can make get rid of FnodeL and VnodeL
  LINALG::Matrix<3,3> FnodeL(true);
  double VnodeL = 0.0;
  for (int i=0; i<neleinpatch; ++i)
  {
    const double V = adjele[i]->Volume()/NUMNOD_PTET;
    VnodeL += V;
    FnodeL.Update(V,adjele[i]->F_,1.0);
  }
  FnodeL.Scale(1.0/VnodeL);

  //-----------------------------------------------------------------------
  // do positioning map global nodes -> position in B-Operator
  map<int,int>  node_pos;
  std::map<int,DRT::Node*>::iterator pnode;
  int count=0;
  for (pnode=nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
  {
    node_pos[pnode->first] = count;
    count++;
  }

  //------------------------------------------------------ build B operator
  Epetra_SerialDenseMatrix bop(NUMSTR_PTET,ndofinpatch);
  // loop elements in patch
  for (int ele=0; ele<neleinpatch; ++ele)
  {
    // current element
    DRT::ELEMENTS::Ptet*      actele = adjele[ele];
    // spatial deriv of that element
    LINALG::Matrix<NUMNOD_PTET,NUMDIM_PTET>& nxyz   = actele->nxyz_;
    // volume of that element assigned to node L
    double V = actele->Volume()/NUMNOD_PTET;
    // def-gradient of the element
    LINALG::Matrix<NUMDIM_PTET,NUMDIM_PTET>& F = actele->F_;

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

      bop(0,NODDOF_PTET*pos+0) += V * F(0,0)*nxyz(i,0);
      bop(0,NODDOF_PTET*pos+1) += V * F(1,0)*nxyz(i,0);
      bop(0,NODDOF_PTET*pos+2) += V * F(2,0)*nxyz(i,0);
      bop(1,NODDOF_PTET*pos+0) += V * F(0,1)*nxyz(i,1);
      bop(1,NODDOF_PTET*pos+1) += V * F(1,1)*nxyz(i,1);
      bop(1,NODDOF_PTET*pos+2) += V * F(2,1)*nxyz(i,1);
      bop(2,NODDOF_PTET*pos+0) += V * F(0,2)*nxyz(i,2);
      bop(2,NODDOF_PTET*pos+1) += V * F(1,2)*nxyz(i,2);
      bop(2,NODDOF_PTET*pos+2) += V * F(2,2)*nxyz(i,2);
      //
      bop(3,NODDOF_PTET*pos+0) += V * (F(0,0)*nxyz(i,1) + F(0,1)*nxyz(i,0) );
      bop(3,NODDOF_PTET*pos+1) += V * (F(1,0)*nxyz(i,1) + F(1,1)*nxyz(i,0) );
      bop(3,NODDOF_PTET*pos+2) += V * (F(2,0)*nxyz(i,1) + F(2,1)*nxyz(i,0) );
      bop(4,NODDOF_PTET*pos+0) += V * (F(0,1)*nxyz(i,2) + F(0,2)*nxyz(i,1) );
      bop(4,NODDOF_PTET*pos+1) += V * (F(1,1)*nxyz(i,2) + F(1,2)*nxyz(i,1) );
      bop(4,NODDOF_PTET*pos+2) += V * (F(2,1)*nxyz(i,2) + F(2,2)*nxyz(i,1) );
      bop(5,NODDOF_PTET*pos+0) += V * (F(0,2)*nxyz(i,0) + F(0,0)*nxyz(i,2) );
      bop(5,NODDOF_PTET*pos+1) += V * (F(1,2)*nxyz(i,0) + F(1,0)*nxyz(i,2) );
      bop(5,NODDOF_PTET*pos+2) += V * (F(2,2)*nxyz(i,0) + F(2,0)*nxyz(i,2) );

    } // for (int i=0; i<actele->NumNode(); ++i)

  } // for (int ele=0; ele<neleinpatch; ++ele)

  //----------------------------------------- averaged material and stresses
  LINALG::Matrix<6,6> cmat(true);
  LINALG::Matrix<6,1> stress(true);

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

  if (nodalstrain)
  {
    if (!ea)
    {
      for (int i = 0; i < 3; ++i) (*nodalstrain)[i] = glstrain(i);
      for (int i = 3; i < 6; ++i) (*nodalstrain)[i] = 0.5 * glstrain(i);
    }
    else
    {
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
  }

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
    LINALG::Matrix<6,6> cmatele(true);
    LINALG::Matrix<6,1> stressele(true);
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      // current element
      DRT::ELEMENTS::Ptet* actele = adjele[ele];
      // volume of that element assigned to node L
      const double V = actele->Volume()/NUMNOD_PTET;
      // def-gradient of the element
      RCP<MAT::Material> mat = actele->Material();
      SelectMaterial(mat,stressele,cmatele,density,glstrain,FnodeL,0);
      // here scaling of stresses/tangent with 1.0-alpha
      cmat.Update(V,cmatele,1.0);
      stress.Update(V,stressele,1.0);
    } // for (int ele=0; ele<neleinpatch; ++ele)
    stress.Scale(1.0/VnodeL);
    cmat.Scale(1.0/VnodeL);
  }

  if (nodalstress)
  {
    if (!cauchy)
      for (int i = 0; i < NUMSTR_PTET; ++i) (*nodalstress)[i] = stress(i);
    else
    {
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
  }

#if 1 // dev stab on cauchy stresses
  {
    // define stuff we need
    LINALG::Matrix<6,6> cmatdev;
    LINALG::Matrix<6,1> stressdev;

    // compute deviatoric stress and tangent
    DevStressTangent(stressdev,cmatdev,cmat,stress,cauchygreen);

    // reduce deviatoric stresses
    stress.Update(-ALPHA_PTET,stressdev,1.0);

    // reduce deviatoric tangent
    cmat.Update(-ALPHA_PTET,cmatdev,1.0);
  }
#endif

#if 0 // orig. puso tet
    stress.Scale(1.0-ALPHA_PTET);
    cmat.Scale(1.0-ALPHA_PTET);
#endif

  //----------------------------------------------------- internal forces
  if (force) 
  {
    Epetra_SerialDenseVector stress_epetra(View,stress.A(),stress.Rows());
    force->Multiply('T','N',VnodeL,bop,stress_epetra,0.0);
  }

  //--------------------------------------------------- elastic stiffness
  if (stiff)
  {
    Epetra_SerialDenseMatrix cmat_epetra(View,cmat.A(),cmat.Rows(),cmat.Rows(),cmat.Columns());
    LINALG::SerialDenseMatrix cb(NUMSTR_PTET,ndofinpatch);
    cb.Multiply('N','N',1.0,cmat_epetra,bop,0.0);
    stiff->Multiply('T','N',VnodeL,bop,cb,0.0);
  }

  //----------------------------------------------------- geom. stiffness
  if (stiff)
  {
    // loop elements in patch
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      // current element
      DRT::ELEMENTS::Ptet*      actele = adjele[ele];
      // spatial deriv of that element
      LINALG::Matrix<NUMNOD_PTET,NUMDIM_PTET>& nxyz   = actele->nxyz_;
      // volume of actele assigned to node L
      double V = actele->Volume()/NUMNOD_PTET;
      // loop nodes of that element
      double SmBL[3];
      DRT::Node** nodes = actele->Nodes();
      for (int i=0; i<NUMNOD_PTET; ++i)
      {
        // row position of this node in matrix
        int ipos = node_pos[nodes[i]->Id()];
        SmBL[0] = V*(stress(0)*nxyz(i,0) + stress(3)*nxyz(i,1) + stress(5)*nxyz(i,2));
        SmBL[1] = V*(stress(3)*nxyz(i,0) + stress(1)*nxyz(i,1) + stress(4)*nxyz(i,2));
        SmBL[2] = V*(stress(5)*nxyz(i,0) + stress(4)*nxyz(i,1) + stress(2)*nxyz(i,2));
        for (int j=0; j<NUMNOD_PTET; ++j)
        {
          // column position of this node in matrix
          int jpos = node_pos[nodes[j]->Id()];
          double bopstrbop = 0.0;
          for (int dim=0; dim<NUMDIM_PTET; ++dim)
            bopstrbop += nxyz(j,dim) * SmBL[dim];
          (*stiff)(NUMDIM_PTET*ipos+0,NUMDIM_PTET*jpos+0) += bopstrbop;
          (*stiff)(NUMDIM_PTET*ipos+1,NUMDIM_PTET*jpos+1) += bopstrbop;
          (*stiff)(NUMDIM_PTET*ipos+2,NUMDIM_PTET*jpos+2) += bopstrbop;
        }
      } // for (int i=0; i<actele->NumNode(); ++i)
    } // for (int ele=0; ele<neleinpatch; ++ele)
  }

  // there is no nodal mass matrix - this is done the conventional way in the elements
  return;
}

/*----------------------------------------------------------------------*
 | material laws for Ptet (protected)                          gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PtetRegister::SelectMaterial(
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
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(glstrain,cmat,stress);
      density = stvk->Density();
    }
    break;
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain,cmat,stress);
      density = neo->Density();
    }
    break;
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast<MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(glstrain,cmat,stress);
      density = aaa->Density();
    }
    break;
    default:
      dserror("Illegal type %d of material for element Ptet tet4", mat->MaterialType());
    break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::Ptet::SelectMaterial

/*----------------------------------------------------------------------*
 |  compute deviatoric tangent and stresses (private/static)   gee 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PtetRegister::DevStressTangent(
  LINALG::Matrix<NUMSTR_PTET,1>& Sdev,
  LINALG::Matrix<NUMSTR_PTET,NUMSTR_PTET>& CCdev,
  LINALG::Matrix<NUMSTR_PTET,NUMSTR_PTET>& CC,
  const LINALG::Matrix<NUMSTR_PTET,1>& S,
  const LINALG::Matrix<NUMDIM_PTET,NUMDIM_PTET>& C)
{

  //---------------------------------- things that we'll definitely need
  // inverse of C
  LINALG::Matrix<3,3> Cinv;
  const double detC = Cinv.Invert(C);

  // J = det(F) = sqrt(detC)
  const double J = sqrt(detC);

  // S as a 3x3 matrix
  LINALG::Matrix<NUMDIM_PTET,NUMDIM_PTET> Smat;
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
  LINALG::Matrix<NUMSTR_PTET,NUMSTR_PTET> CCvol(true); // fill with zeros

  //--------------------------------------- CCvol += 2pJ (Cinv boeppel Cinv)
  MAT::ElastSymTensor_o_Multiply(CCvol,-2.0*fac,Cinv,Cinv,0.0);

  //------------------------------------------ CCvol += 2/3 * Cinv dyad S
  MAT::ElastSymTensorMultiply(CCvol,2.0/3.0,Cinv,Smat,1.0);

  //-------------------------------------- CCvol += 1/3 Cinv dyad ( CC : C )
  {
    // C as Voigt vector
    LINALG::Matrix<NUMSTR_PTET,1> Cvec;
    Cvec(0) = C(0,0);
    Cvec(1) = C(1,1);
    Cvec(2) = C(2,2);
    Cvec(3) = 2.0*C(0,1);
    Cvec(4) = 2.0*C(1,2);
    Cvec(5) = 2.0*C(0,2);

    LINALG::Matrix<NUMSTR_PTET,1> CCcolonC;
    CCcolonC.Multiply(CC,Cvec);

    LINALG::Matrix<NUMDIM_PTET,NUMDIM_PTET> CCcC;
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
