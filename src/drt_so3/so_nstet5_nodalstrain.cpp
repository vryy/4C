/*!----------------------------------------------------------------------*
\file so_nstet5_nodalstrain.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_FECrsMatrix.h"
#include "Sacado.hpp"
#include "../linalg/linalg_serialdensevector.H"
//#include "../linalg/linalg_fixedsizematrix.H"


#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"
#include "../drt_mat/elasthyper.H"

#include "so_nstet5.H"

/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::ElementDeformationGradient(DRT::Discretization& dis)
{
  // current displacement
  RCP<const Epetra_Vector> disp = dis.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  // loop elements
  std::map<int,DRT::ELEMENTS::NStet5*>::iterator ele;
  for (ele=elecids_.begin(); ele != elecids_.end(); ++ele)
  {
    DRT::ELEMENTS::NStet5* e = ele->second;
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    e->LocationVector(dis,lm,lmowner,lmstride);
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

    //------------------------------------subelement F
    LINALG::Matrix<4,3> disp(false);
    LINALG::Matrix<5,3> subdisp(false);
    for (int j=0; j<3; ++j)
    {
      for (int i=0; i<4; ++i)
        subdisp(i,j) = mydisp[i*3+j];
      subdisp(4,j) = mydisp[4*3+j];
    }
    for (int k=0; k<4; ++k) 
    {
      for (int i=0; i<4; ++i)
        for (int j=0; j<3; ++j)
          disp(i,j) = subdisp(e->SubLM(k)[i],j);
      
      e->SubF(k) = e->BuildF(disp,e->SubNxyz(k));
      double J = e->SubF(k).Determinant();
      if (J<=0.0) dserror("det(F) of Element %d / Subelement %d %10.5e <= 0 !!\n",e->Id(),k,J);
    } // for (int k=0; k<4; ++k) 

  } // ele
  return;
}


/*----------------------------------------------------------------------*
 |  pre-evaluation of elements (public)                        gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::PreEvaluate(DRT::Discretization& dis,
                                          Teuchos::ParameterList& p,
                                          RCP<LINALG::SparseOperator> systemmatrix1,
                                          RCP<LINALG::SparseOperator> systemmatrix2,
                                          RCP<Epetra_Vector>          systemvector1,
                                          RCP<Epetra_Vector>          systemvector2,
                                          RCP<Epetra_Vector>          systemvector3)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::PreEvaluate");

  // nodal integration for nlnstiff and internal forces only
  // (this method does not compute stresses/strains/element updates/mass matrix)
  string& action = p.get<string>("action","none");
  if (action != "calc_struct_nlnstiffmass" &&
      action != "calc_struct_nlnstiff"     &&
      action != "calc_struct_stress") return;

  // These get filled in here, so remove old stuff
  if (action == "calc_struct_stress")
  {
    nstress_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(),6,false));
    nstrain_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(),6,false));
  }
  else
  {
    nstress_ = Teuchos::null;
    nstrain_ = Teuchos::null;
  }

  // see what we have for input
  bool assemblemat1 = systemmatrix1!=Teuchos::null;
  bool assemblevec1 = systemvector1!=Teuchos::null;
  bool assemblevec2 = systemvector2!=Teuchos::null;
  bool assemblevec3 = systemvector3!=Teuchos::null;
  if (assemblevec2 || assemblevec3) dserror("Wrong assembly expectations");

  //-----------------------------------------------------------------
  // nodal stiffness and force (we don't do mass here)
  LINALG::SerialDenseMatrix stiff;
  LINALG::SerialDenseVector force;

  //-------------------------------------- construct F for each NStet5
  ElementDeformationGradient(dis);

  //-----------------------------------------------------------------
  // create a temporary matrix to assemble to in a baci-unusual way
  // (across-parallel-interface assembly)
  const Epetra_Map* rmap = NULL;
  const Epetra_Map* dmap = NULL;

  RCP<Epetra_FECrsMatrix> stifftmp;
  RCP<LINALG::SparseMatrix> systemmatrix;
  if (systemmatrix1 != Teuchos::null)
  {
    rmap = &(systemmatrix1->OperatorRangeMap());
    dmap = rmap;
    systemmatrix = rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix1);
    if (systemmatrix != null && systemmatrix->Filled())
      stifftmp = rcp(new Epetra_FECrsMatrix(::Copy,systemmatrix->EpetraMatrix()->Graph()));
    else
      stifftmp = rcp(new Epetra_FECrsMatrix(::Copy,*rmap,256,false));
  }

  //-----------------------------------------------------------------
  // make some tests for fast assembly
  if (systemmatrix != null && systemmatrix->Filled())
  {
    Epetra_CrsMatrix& matrix = *(systemmatrix->EpetraMatrix());
    if (!matrix.StorageOptimized()) dserror("Matrix must be StorageOptimized() when Filled()");
  }

  //-----------------------------------------------------------------
  // create temporary vector in column map to assemble to
  Epetra_Vector forcetmp1(*dis.DofColMap(),true);

  //-----------------------------------------------------------------
  // current displacements
  RCP<const Epetra_Vector> disp = dis.GetState("displacement");

  //================================================== do nodal stiffness
  std::map<int,DRT::Node*>::iterator node;
  for (node=noderids_.begin(); node != noderids_.end(); ++node)
  {
    DRT::Node* nodeL   = node->second;     // row node
    const int  nodeLid = nodeL->Id();

    // standard quantities for all nodes
    vector<DRT::ELEMENTS::NStet5*>&  adjele    = adjele_[nodeLid];
    map<int,vector<int> >&           adjsubele = adjsubele_[nodeLid];
    map<int,DRT::Node*>&             adjnode   = adjnode_[nodeLid];
    vector<int>&                     lm        = adjlm_[nodeLid];
    vector<vector<vector<int> > >&   lmlm      = lmlm_[nodeLid];
    const int ndofperpatch = (int)lm.size();

    if (action != "calc_struct_stress")
    {
      // do nodal integration of stiffness and internal force
      stiff.LightShape(ndofperpatch,ndofperpatch); 
      force.LightSize(ndofperpatch);
      TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::NodalIntegration");
      NodalIntegration(&stiff,&force,adjnode,adjele,adjsubele,lm,lmlm,*disp,dis,
                       NULL,NULL,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    else
    {
      INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(p,"iostress",INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(p,"iostrain",INPAR::STR::strain_none);
      vector<double> nodalstress(6);
      vector<double> nodalstrain(6);
      NodalIntegration(NULL,NULL,adjnode,adjele,adjsubele,lm,lmlm,*disp,dis,
                       &nodalstress,&nodalstrain,iostress,iostrain);

      const int lid = dis.NodeRowMap()->LID(nodeLid);
      if (lid==-1) dserror("Cannot find local id for row node");
      for (int i=0; i<6; ++i)
      {
        (*(*nstress_)(i))[lid] = nodalstress[i];
        (*(*nstrain_)(i))[lid] = nodalstrain[i];
      }
    }


    //---------------------- do assembly of stiffness and internal force
    // (note: this is non-standard-baci assembly and therefore a do it all yourself version!)
    // there is no guarantee that systemmatrix exists
    // (e.g. if systemmatrix1 is actually a BlockSparseMatrix)
    bool fastassemble = false;
    if (systemmatrix != null) fastassemble = true;

    if (assemblemat1)
    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::PreEvaluate Assembly");
      vector<int> lrlm;
      vector<int> lclm;

      const Epetra_Map& dofrowmap = systemmatrix1->OperatorRangeMap();
      lrlm.resize(ndofperpatch);
      for (int i=0; i<ndofperpatch; ++i)
        lrlm[i] = dofrowmap.LID(lm[i]);
      if (fastassemble)
      {
        const Epetra_Map& dofcolmap = systemmatrix->ColMap();
        lclm.resize(ndofperpatch);
        for (int i=0; i<ndofperpatch; ++i)
          lclm[i] = dofcolmap.LID(lm[i]);
      }

      for (int i=0; i<ndofperpatch; ++i)
      {
        if (lrlm[i]==-1) // off-processor row
        {
          for (int j=0; j<ndofperpatch; ++j)
          {
            int errone = stifftmp->SumIntoGlobalValues(1,&lm[i],1,&lm[j],&stiff(i,j));
            if (errone>0)
            {
              int errtwo = stifftmp->InsertGlobalValues(1,&lm[i],1,&lm[j],&stiff(i,j));
              if (errtwo<0) dserror("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d",errtwo);
            }
            else if (errone) dserror("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d",errone);
          }
        }
        else // local row
        {
          if (systemmatrix != null && systemmatrix->Filled()) // matrix is SparseMatrix
          {
            Epetra_CrsMatrix& matrix = *(systemmatrix->EpetraMatrix());
            int length;
            double* values;
            int* indices;
            matrix.ExtractMyRowView(lrlm[i],length,values,indices);
            for (int j=0; j<ndofperpatch; ++j)
            {
              int* loc = std::lower_bound(indices,indices+length,lclm[j]);
//#ifdef DEBUG
              if (*loc != lclm[j]) dserror("Cannot find local column entry %d",lclm[j]);
//#endif
              int pos = loc-indices;

              // test physical continuity of nodal values inside the Epetra_CrsMatrix
              bool continuous = true;
              for (int k=1; k<3; ++k)
                if (indices[pos+k] == lclm[j+k]) continue;
                else
                {
                  continuous = false;
                  break;
                }

              if (continuous)
              {
                values[pos++] += stiff(i,j++);
                values[pos++] += stiff(i,j++);
                values[pos]   += stiff(i,j);
              }
              else
              {
                int err =  matrix.SumIntoMyValues(lrlm[i],1,&stiff(i,j),&lclm[j]); j++;
                    err += matrix.SumIntoMyValues(lrlm[i],1,&stiff(i,j),&lclm[j]); j++;
                    err += matrix.SumIntoMyValues(lrlm[i],1,&stiff(i,j),&lclm[j]);
                if (err) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d",err);
              }
            }
          }
          else // matrix not SparseMatrix (e.g. BlockMatrix) -> fall back to standard assembly
          {
            for (int j=0; j<ndofperpatch; ++j)
              systemmatrix1->Assemble(stiff(i,j),lm[i],lm[j]);
          }
        }
      }
    }

    //-----------------------------------------------------------------------------------
    if (assemblevec1)
    {
      for (int i=0; i<ndofperpatch; ++i)
      {
        const int rgid = lm[i];
        const int lid = forcetmp1.Map().LID(rgid);
        if (lid<0) dserror("global row %d does not exist in column map",rgid);
        forcetmp1[lid] += force[i];
      }
    }

  //=========================================================================
  } // for (node=noderids_.begin(); node != noderids_.end(); ++node)

  //-------------------------------------------------------------------------
  if (action == "calc_struct_stress")
  {
    // we have to export the nodal stresses and strains to column map
    // so they can be written by the elements
    RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeColMap(),6,false));
    LINALG::Export(*nstress_,*tmp);
    nstress_ = tmp;
    tmp = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeColMap(),6,false));
    LINALG::Export(*nstrain_,*tmp);
    nstrain_ = tmp;
  }


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
    int err = stifftmp->GlobalAssemble(*dmap,*rmap,false);
    if (err) dserror("Epetra_FECrsMatrix::GlobalAssemble returned err=%d",err);
    const Epetra_Map& cmap = stifftmp->ColMap();
    for (int lrow=0; lrow<stifftmp->NumMyRows(); ++lrow)
    {
      int numentries;
      double* values;
      if (!stifftmp->Filled())
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
        if (systemmatrix != null && systemmatrix->Filled())
        {
          Epetra_CrsMatrix& matrix = *systemmatrix->EpetraMatrix();
          for (int j=0; j<numentries; ++j)
          {
            int err = matrix.SumIntoMyValues(lrow,1,&values[j],&lindices[j]);
            if (err) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d",err);
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

  return;
}

/*----------------------------------------------------------------------*
 |  do nodal integration (public)                              gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::NodalIntegration(
                                                Epetra_SerialDenseMatrix*      stiff,
                                                Epetra_SerialDenseVector*       force,
                                                map<int,DRT::Node*>&            adjnode,
                                                vector<DRT::ELEMENTS::NStet5*>& adjele,
                                                map<int,vector<int> >&          adjsubele,
                                                vector<int>&                    lm,
                                                vector<vector<vector<int> > >&  lmlm,
                                                const Epetra_Vector&            disp,
                                                DRT::Discretization&            dis,
                                                vector<double>*                 nodalstress,
                                                vector<double>*                 nodalstrain,
                                                const INPAR::STR::StressType    iostress,
                                                const INPAR::STR::StrainType    iostrain)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::NodalIntegration");
  typedef Sacado::Fad::DFad<double> FAD; // for first derivs
  //typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD; // for second derivs

  //-------------------------------------------------- standard quantities
  const int ndofinpatch  = (int)lm.size();
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
  // get displacements of this patch
  vector<FAD> patchdisp(ndofinpatch);
  for (int i=0; i<ndofinpatch; ++i)
  {
    int lid = disp.Map().LID(lm[i]);
    if (lid==-1) dserror("Cannot find degree of freedom on this proc");
    patchdisp[i] = disp[disp.Map().LID(lm[i])];
    patchdisp[i].diff(i,ndofinpatch);
  }

  //-----------------------------------------------------------------------
  // build averaged F and volume of node (using sacado)
  double Vnode = 0.0;
  //FAD    TVnode = 0.0;
  LINALG::TMatrix<FAD,3,3> TFnode(true);
  for (int i=0; i<neleinpatch; ++i)
  {
    DRT::ELEMENTS::NStet5* ele = adjele[i];
    vector<int>& subele = adjsubele[ele->Id()];

    for (unsigned j=0; j<subele.size(); ++j)
    {
      const int   subeleid = subele[j];
      
      // copy subelement displacements to 4x3 format
     LINALG::TMatrix<FAD,4,3> eledispmat(false);
     for (int k=0; k<4; ++k)
       for (int l=0; l<3; ++l)
         eledispmat(k,l) = patchdisp[lmlm[i][j][k*3+l]];
    
     // add 1/3 of subelement volume to this node
     const double V = ele->SubV(subeleid)/3.0;
     Vnode += V;

     // build F from this subelement
     LINALG::TMatrix<FAD,3,3> F(true);
     F = TBuildF<FAD>(eledispmat,ele->SubNxyz(subeleid));

     //FAD detFe = Determinant3x3<FAD>(F);
     //TVnode += (detFe*V);
     
     // add to nodal deformation gradient
     F.Scale(V);
     TFnode += F;

    } // for (unsigned j=0; j<subele.size(); ++j)
  } // for (int i=0; i<neleinpatch; ++i)
  
  // do the actual averaging
  TFnode.Scale(1.0/Vnode);
  //FAD detFnode = Determinant3x3<FAD>(TFnode);
  
  // build corrected integration area of this node
  //TVnode /= detFnode;
  //printf("V %10.5e Vtilde %10.5e V/Vtilde %10.5e\n",Vnode,TVnode.val(),Vnode/TVnode.val());
  //Vnode = TVnode.val();
  
  // copy fad F to double F
  LINALG::Matrix<3,3> FnodeL(false);
  for (int j=0; j<3; ++j)
    for (int k=0; k<3; ++k)
      FnodeL(j,k) = TFnode(j,k).val();

  //-----------------------------------------------------------------------
  // build B operator
  Epetra_SerialDenseMatrix bop(6,ndofinpatch);
  // loop elements in patch
  for (int ele=0; ele<neleinpatch; ++ele)
  {
    // current element
    DRT::ELEMENTS::NStet5* actele = adjele[ele];
    vector<int>& subele = adjsubele[actele->Id()];
    
    // loop subelements in this element
    for (unsigned j=0; j<subele.size(); ++j)
    {
      const int   subeleid = subele[j];

      // volume of that subelement assigned to node L
      double V = actele->SubV(subeleid)/3;

      // volume ratio of volume per node of this element to
      // whole volume of node L
      const double ratio = V/Vnode;
      
      // get derivatives with respect to X
      const LINALG::Matrix<4,3>& nxyz = actele->SubNxyz(subeleid);
      
      // get defgrd
      LINALG::Matrix<3,3>& F = actele->SubF(subeleid);
      
      LINALG::Matrix<6,12> bele(false);
      
      for (int i=0; i<4; ++i)
      {
        bele(0,3*i+0) = F(0,0)*nxyz(i,0);
        bele(0,3*i+1) = F(1,0)*nxyz(i,0);
        bele(0,3*i+2) = F(2,0)*nxyz(i,0);
        bele(1,3*i+0) = F(0,1)*nxyz(i,1);
        bele(1,3*i+1) = F(1,1)*nxyz(i,1);
        bele(1,3*i+2) = F(2,1)*nxyz(i,1);
        bele(2,3*i+0) = F(0,2)*nxyz(i,2);
        bele(2,3*i+1) = F(1,2)*nxyz(i,2);
        bele(2,3*i+2) = F(2,2)*nxyz(i,2);

        bele(3,3*i+0) = F(0,0)*nxyz(i,1) + F(0,1)*nxyz(i,0);
        bele(3,3*i+1) = F(1,0)*nxyz(i,1) + F(1,1)*nxyz(i,0);
        bele(3,3*i+2) = F(2,0)*nxyz(i,1) + F(2,1)*nxyz(i,0);
        bele(4,3*i+0) = F(0,1)*nxyz(i,2) + F(0,2)*nxyz(i,1);
        bele(4,3*i+1) = F(1,1)*nxyz(i,2) + F(1,2)*nxyz(i,1);
        bele(4,3*i+2) = F(2,1)*nxyz(i,2) + F(2,2)*nxyz(i,1);
        bele(5,3*i+0) = F(0,2)*nxyz(i,0) + F(0,0)*nxyz(i,2);
        bele(5,3*i+1) = F(1,2)*nxyz(i,0) + F(1,0)*nxyz(i,2);
        bele(5,3*i+2) = F(2,2)*nxyz(i,0) + F(2,0)*nxyz(i,2);
      }

      for (int k=0; k<6; ++k)
        for (int i=0; i<12; ++i)
          bop(k,lmlm[ele][j][i]) += ratio * bele(k,i);

    } // for (unsigned j=0; j<subele.size(); ++j)
  } // for (int ele=0; ele<neleinpatch; ++ele)

  //-------------------------------------------------------------- averaged strain
  // right cauchy green
  LINALG::TMatrix<FAD,3,3> CG(false);
  CG.MultiplyTN(TFnode,TFnode);
  vector<FAD> Ebar(6);
  Ebar[0] = 0.5 * (CG(0,0) - 1.0);
  Ebar[1] = 0.5 * (CG(1,1) - 1.0);
  Ebar[2] = 0.5 * (CG(2,2) - 1.0);
  Ebar[3] =        CG(0,1);
  Ebar[4] =        CG(1,2);
  Ebar[5] =        CG(2,0);

  // for material law and output, copy to LINALG::Matrix object
  LINALG::Matrix<3,3> cauchygreen(false);
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      cauchygreen(i,j) = CG(i,j).val();
  LINALG::Matrix<6,1> glstrain(false);
  glstrain(0) = Ebar[0].val();
  glstrain(1) = Ebar[1].val();
  glstrain(2) = Ebar[2].val();
  glstrain(3) = Ebar[3].val();
  glstrain(4) = Ebar[4].val();
  glstrain(5) = Ebar[5].val();

  //-------------------------------------------------------- output of strain
  if (iostrain != INPAR::STR::strain_none)
  {
#ifndef PUSO_NSTET5
    StrainOutput(iostrain,*nodalstrain,FnodeL,FnodeL.Determinant(),1.0,1.0-ALPHA_NSTET5);
#else
    StrainOutput(iostrain,*nodalstrain,FnodeL,glstrain,1.0-ALPHA_NSTET5);
#endif
  }

  //-------------------------------------------------------------------------
  // build a second B-operator from the averaged strains that are based on
  // the averaged F
  Epetra_SerialDenseMatrix bopbar(6,ndofinpatch);
  for (int i=0; i<ndofinpatch; ++i)
    for (int k=0; k<6; ++k)
      bopbar(k,i) = Ebar[k].fastAccessDx(i);

  //----------------------------------------- averaged material and stresses
  LINALG::Matrix<6,6> cmat(true);
  LINALG::Matrix<6,1> stress(true);

  //-----------------------------------------------------------------------
  // material law
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
      DRT::ELEMENTS::NStet5* actele = adjele[ele];
      // volume of that element assigned to node L
      double V = 0.0;
      for (unsigned j=0; j<adjsubele[actele->Id()].size(); ++j)
        V += (actele->SubV(adjsubele[actele->Id()][j])/3.0);
      // material of the element
      RCP<MAT::Material> mat = actele->Material();
      SelectMaterial(mat,stressele,cmatele,density,glstrain,FnodeL,0);
      cmat.Update(V,cmatele,1.0);
      stress.Update(V,stressele,1.0);
    } // for (int ele=0; ele<neleinpatch; ++ele)
    stress.Scale(1.0/Vnode);
    cmat.Scale(1.0/Vnode);
  }

  //-----------------------------------------------------------------------
  // stress is split as follows:
  // stress = vol_node + (1-alpha) * dev_node + alpha * dev_ele
#ifndef PUSO_NSTET5
  {
    LINALG::Matrix<6,1> stressdev(true);
    LINALG::Matrix<6,6> cmatdev(true);
    LINALG::Matrix<6,1> stressvol(false);
    LINALG::Matrix<6,6> cmatvol(false);

    // compute deviatoric stress and tangent from total stress and tangent
    DevStressTangent(stressdev,cmatdev,cmat,stress,cauchygreen);

    // compute volumetric stress and tangent
    stressvol.Update(-1.0,stressdev,1.0,stress,0.0);
    cmatvol.Update(-1.0,cmatdev,1.0,cmat,0.0);

    // compute nodal stress
    stress.Update(1.0,stressvol,1-ALPHA_NSTET5,stressdev,0.0);
    cmat.Update(1.0,cmatvol,1-ALPHA_NSTET5,cmatdev,0.0);
  }
#else
  {
    stress.Scale(1.-ALPHA_NSTET5);
    cmat.Scale(1.-ALPHA_NSTET5);
  }
#endif
  //-----------------------------------------------------------------------
  // stress output
  if (iostress != INPAR::STR::stress_none)
  {
    StressOutput(iostress,*nodalstress,stress,FnodeL,FnodeL.Determinant());
  }

  //----------------------------------------------------- internal forces
  if (force)
  {
    Epetra_SerialDenseVector stress_epetra(::View,stress.A(),stress.Rows());
    force->Multiply('T','N',Vnode,bop,stress_epetra,0.0);
  }

  //--------------------------------------------------- elastic stiffness
  if (stiff)
  {
    Epetra_SerialDenseMatrix cmat_epetra(::View,cmat.A(),cmat.Rows(),cmat.Rows(),cmat.Columns());
    LINALG::SerialDenseMatrix cb(6,ndofinpatch);
    cb.Multiply('N','N',1.0,cmat_epetra,bopbar,0.0);
    stiff->Multiply('T','N',Vnode,bop,cb,0.0);// bop
  }

  //----------------------------------------------------- geom. stiffness
  // do not use sacado for second derivative of E as it is too expensive!
  // As long as the 2nd deriv is as easy as this, do it by hand
  if (stiff)
  {
    // loop elements in patch
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      // current element
      DRT::ELEMENTS::NStet5* actele = adjele[ele];
      vector<int>& subele = adjsubele[actele->Id()];
      
      // loop subelements in this element
      for (unsigned sub=0; sub<subele.size(); ++sub)
      {
        const int   subeleid = subele[sub];
      
        // material deriv of element
        const LINALG::Matrix<4,3>& nxyz = actele->SubNxyz(subeleid);
      
        // volume of element assigned to node L
        const double V = actele->SubV(subeleid)/3;

        // loop nodes of that element
        double SmBL[3];
        for (int i=0; i<4; ++i)
        {
          SmBL[0] = V*(stress(0)*nxyz(i,0) + stress(3)*nxyz(i,1) + stress(5)*nxyz(i,2));
          SmBL[1] = V*(stress(3)*nxyz(i,0) + stress(1)*nxyz(i,1) + stress(4)*nxyz(i,2));
          SmBL[2] = V*(stress(5)*nxyz(i,0) + stress(4)*nxyz(i,1) + stress(2)*nxyz(i,2));
          for (int j=0; j<4; ++j)
          {
            double bopstrbop = 0.0;
            for (int dim=0; dim<3; ++dim) bopstrbop += nxyz(j,dim) * SmBL[dim];
            (*stiff)(lmlm[ele][sub][i*3+0],lmlm[ele][sub][j*3+0]) += bopstrbop;
            (*stiff)(lmlm[ele][sub][i*3+1],lmlm[ele][sub][j*3+1]) += bopstrbop;
            (*stiff)(lmlm[ele][sub][i*3+2],lmlm[ele][sub][j*3+2]) += bopstrbop;
          } // for (int j=0; j<4; ++j)
        } // for (int i=0; i<4; ++i)
      } // for (unsigned sub=0; sub<subele.size(); ++sub)
    } // for (int ele=0; ele<neleinpatch; ++ele)

  } // if (stiff)

  return;
}









/*----------------------------------------------------------------------*
 | material laws for NStet5 (protected)                        gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::SelectMaterial(
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
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(glstrain,cmat,stress);
      density = hyper->Density();
      return;
      break;
    }
    default:
      dserror("Illegal type %d of material for element NStet5 tet4", mat->MaterialType());
    break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::NStet5::SelectMaterial

/*----------------------------------------------------------------------*
 |  compute deviatoric tangent and stresses (private/static)   gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::DevStressTangent(
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

/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::StrainOutput(
                    const INPAR::STR::StrainType iostrain,
                    vector<double>&              nodalstrain,
                    LINALG::Matrix<3,3>&         F,
                    const double&                detF,
                    const double                 volweight,
                    const double                 devweight)
{
  LINALG::Matrix<3,3> Fiso = F;
  Fiso.Scale(pow(detF,-1.0/3.0));

  LINALG::Matrix<3,3> Fvol(true);
  Fvol(0,0) = 1.0; Fvol(1,1) = 1.0; Fvol(2,2) = 1.0;
  Fvol.Scale(pow(detF,1.0/3.0));

  LINALG::Matrix<3,3> cauchygreeniso(false);
  cauchygreeniso.MultiplyTN(Fiso,Fiso);

  LINALG::Matrix<3,3> cauchygreenvol(false);
  cauchygreenvol.MultiplyTN(Fvol,Fvol);

  LINALG::Matrix<3,3> glstrainiso(false);
  glstrainiso(0,0) = 0.5 * (cauchygreeniso(0,0) - 1.0);
  glstrainiso(0,1) = 0.5 *  cauchygreeniso(0,1);
  glstrainiso(0,2) = 0.5 *  cauchygreeniso(0,2);
  glstrainiso(1,0) = glstrainiso(0,1);
  glstrainiso(1,1) = 0.5 * (cauchygreeniso(1,1) - 1.0);
  glstrainiso(1,2) = 0.5 *  cauchygreeniso(1,2);
  glstrainiso(2,0) = glstrainiso(0,2);
  glstrainiso(2,1) = glstrainiso(1,2);
  glstrainiso(2,2) = 0.5 * (cauchygreeniso(2,2) - 1.0);

  LINALG::Matrix<3,3> glstrainvol(false);
  glstrainvol(0,0) = 0.5 * (cauchygreenvol(0,0) - 1.0);
  glstrainvol(0,1) = 0.5 *  cauchygreenvol(0,1);
  glstrainvol(0,2) = 0.5 *  cauchygreenvol(0,2);
  glstrainvol(1,0) = glstrainvol(0,1);
  glstrainvol(1,1) = 0.5 * (cauchygreenvol(1,1) - 1.0);
  glstrainvol(1,2) = 0.5 *  cauchygreenvol(1,2);
  glstrainvol(2,0) = glstrainvol(0,2);
  glstrainvol(2,1) = glstrainvol(1,2);
  glstrainvol(2,2) = 0.5 * (cauchygreenvol(2,2) - 1.0);

  LINALG::Matrix<3,3> glstrainout = glstrainiso;
  glstrainout.Update(volweight,glstrainvol,devweight);

  switch (iostrain)
  {
  case INPAR::STR::strain_gl:
  {
    nodalstrain[0] = glstrainout(0,0);
    nodalstrain[1] = glstrainout(1,1);
    nodalstrain[2] = glstrainout(2,2);
    nodalstrain[3] = glstrainout(0,1);
    nodalstrain[4] = glstrainout(1,2);
    nodalstrain[5] = glstrainout(0,2);
  }
  break;
  case INPAR::STR::strain_ea:
  {
    // inverse of deformation gradient
    LINALG::Matrix<3,3> invdefgrd;
    invdefgrd.Invert(F);
    LINALG::Matrix<3,3> temp;
    LINALG::Matrix<3,3> euler_almansi;
    temp.Multiply(glstrainout,invdefgrd);
    euler_almansi.MultiplyTN(invdefgrd,temp);
    nodalstrain[0] = euler_almansi(0,0);
    nodalstrain[1] = euler_almansi(1,1);
    nodalstrain[2] = euler_almansi(2,2);
    nodalstrain[3] = euler_almansi(0,1);
    nodalstrain[4] = euler_almansi(1,2);
    nodalstrain[5] = euler_almansi(0,2);
  }
  break;
  case INPAR::STR::strain_none:
    break;
  default:
    dserror("requested strain type not available");
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::StrainOutput(
                    const INPAR::STR::StrainType iostrain,
                    vector<double>&              nodalstrain,
                    LINALG::Matrix<3,3>&         F,
                    LINALG::Matrix<6,1>&         glstrain,
                    const double                 weight)
{
  LINALG::Matrix<3,3> glstrainout;

  glstrainout(0,0) = weight * glstrain(0);
  glstrainout(1,1) = weight * glstrain(1);
  glstrainout(2,2) = weight * glstrain(2);
  glstrainout(0,1) = weight * glstrain(3);
  glstrainout(1,2) = weight * glstrain(4);
  glstrainout(0,2) = weight * glstrain(5);


  switch (iostrain)
  {
  case INPAR::STR::strain_gl:
  {
    nodalstrain[0] = glstrainout(0,0);
    nodalstrain[1] = glstrainout(1,1);
    nodalstrain[2] = glstrainout(2,2);
    nodalstrain[3] = glstrainout(0,1);
    nodalstrain[4] = glstrainout(1,2);
    nodalstrain[5] = glstrainout(0,2);
  }
  break;
  case INPAR::STR::strain_ea:
  {
    // inverse of deformation gradient
    LINALG::Matrix<3,3> invdefgrd;
    invdefgrd.Invert(F);
    LINALG::Matrix<3,3> temp;
    LINALG::Matrix<3,3> euler_almansi;
    temp.Multiply(glstrainout,invdefgrd);
    euler_almansi.MultiplyTN(invdefgrd,temp);
    nodalstrain[0] = euler_almansi(0,0);
    nodalstrain[1] = euler_almansi(1,1);
    nodalstrain[2] = euler_almansi(2,2);
    nodalstrain[3] = euler_almansi(0,1);
    nodalstrain[4] = euler_almansi(1,2);
    nodalstrain[5] = euler_almansi(0,2);
  }
  break;
  case INPAR::STR::strain_none:
    break;
  default:
    dserror("requested strain type not available");
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::StressOutput(
                    const INPAR::STR::StressType iostress,
                    vector<double>&              nodalstress,
                    LINALG::Matrix<6,1>&         stress,
                    LINALG::Matrix<3,3>&         F,
                    const double&                detF)
{
  switch (iostress)
  {
  case INPAR::STR::stress_2pk:
  {
    for (int i = 0; i < 6; ++i) nodalstress[i] = stress(i);
  }
  break;
  case INPAR::STR::stress_cauchy:
  {
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
    temp.Multiply(1.0/detF,F,pkstress);
    cauchystress.MultiplyNT(temp,F);
    nodalstress[0] = cauchystress(0,0);
    nodalstress[1] = cauchystress(1,1);
    nodalstress[2] = cauchystress(2,2);
    nodalstress[3] = cauchystress(0,1);
    nodalstress[4] = cauchystress(1,2);
    nodalstress[5] = cauchystress(0,2);
  }
  break;
  case INPAR::STR::stress_none:
    break;
  default:
    dserror("requested stress type not available");
  }
  return;
}




#endif  // #ifdef CCADISCRET
