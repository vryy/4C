/*!----------------------------------------------------------------------*

\brief Nstet element nodal strain implementation

\maintainer Christoph Meier

\level 3

*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_FECrsMatrix.h"
//#include "Sacado.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "../linalg/linalg_serialdensevector.H"
//#include "../linalg/linalg_fixedsizematrix.H"

#include "../drt_mat/material_service.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/elasthyper.H"

#include "so_nstet5.H"

/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::ElementDeformationGradient(DRT::Discretization& dis)
{
  // current displacement
  Teuchos::RCP<const Epetra_Vector> disp = dis.GetState("displacement");
  if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
  // loop elements
  std::map<int, DRT::ELEMENTS::NStet5*>::iterator ele;
  for (ele = elecids_.begin(); ele != elecids_.end(); ++ele)
  {
    DRT::ELEMENTS::NStet5* e = ele->second;
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    e->LocationVector(dis, lm, lmowner, lmstride);
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

    //------------------------------------subelement F
    LINALG::Matrix<5, 3> subdisp(false);
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 5; ++i) subdisp(i, j) = mydisp[i * 3 + j];

    LINALG::Matrix<4, 3> disp(false);
    for (int k = 0; k < 4; ++k)  // subelement k
    {
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j) disp(i, j) = subdisp(e->SubLM(k)[i], j);

      e->SubF(k) = e->BuildF(disp, e->SubNxyz(k));
      double J = e->SubF(k).Determinant();
      if (J <= 0.0) dserror("det(F) of Element %d / Subelement %d %10.5e <= 0 !!\n", e->Id(), k, J);
    }  // for (int k=0; k<4; ++k)

  }  // ele
  return;
}


/*----------------------------------------------------------------------*
 |  pre-evaluation of elements (public)                        gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::PreEvaluate");

  // nodal integration for nlnstiff and internal forces only
  // (this method does not compute stresses/strains/element updates/mass matrix)
  std::string& action = p.get<std::string>("action", "none");
  if (action != "calc_struct_nlnstiffmass" && action != "calc_struct_nlnstifflmass" &&
      action != "calc_struct_nlnstiff" && action != "calc_struct_stress" &&
      action != "calc_struct_internalforce")
    return;

  // These get filled in here, so remove old stuff
  if (action == "calc_struct_stress")
  {
    nstress_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), 6, false));
    nstrain_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), 6, false));
  }
  else
  {
    nstress_ = Teuchos::null;
    nstrain_ = Teuchos::null;
  }

  // see what we have for input
  bool assemblemat1 = systemmatrix1 != Teuchos::null;
  bool assemblevec1 = systemvector1 != Teuchos::null;
  bool assemblevec2 = systemvector2 != Teuchos::null;
  bool assemblevec3 = systemvector3 != Teuchos::null;
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

  Teuchos::RCP<Epetra_FECrsMatrix> stifftmp;
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix;
  if (systemmatrix1 != Teuchos::null)
  {
    rmap = &(systemmatrix1->OperatorRangeMap());
    dmap = rmap;
    systemmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix1);
    if (systemmatrix != Teuchos::null && systemmatrix->Filled())
      stifftmp =
          Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, systemmatrix->EpetraMatrix()->Graph()));
    else
      stifftmp = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, *rmap, 256, false));
  }

  //-----------------------------------------------------------------
  // make some tests for fast assembly
  if (systemmatrix != Teuchos::null && systemmatrix->Filled())
  {
    Epetra_CrsMatrix& matrix = *(systemmatrix->EpetraMatrix());
    if (!matrix.StorageOptimized()) dserror("Matrix must be StorageOptimized() when Filled()");
  }

  //-----------------------------------------------------------------
  // create temporary vector in column map to assemble to
  Epetra_Vector forcetmp1(*dis.DofColMap(), true);

  //-----------------------------------------------------------------
  // current displacements
  Teuchos::RCP<const Epetra_Vector> disp = dis.GetState("displacement");

  //================================================== do nodal stiffness
  std::map<int, DRT::Node*>::iterator node;
  for (node = noderids_.begin(); node != noderids_.end(); ++node)
  {
    DRT::Node* nodeL = node->second;  // row node
    const int nodeLid = nodeL->Id();

    // standard quantities for all nodes
    std::vector<DRT::ELEMENTS::NStet5*>& adjele = adjele_[nodeLid];
    std::map<int, std::vector<int>>& adjsubele = adjsubele_[nodeLid];
    std::map<int, DRT::Node*>& adjnode = adjnode_[nodeLid];
    std::vector<int>& lm = adjlm_[nodeLid];
    std::vector<std::vector<std::vector<int>>>& lmlm = lmlm_[nodeLid];
    const int ndofperpatch = (int)lm.size();

    if (action == "calc_struct_nlnstiffmass" || action == "calc_struct_nlnstifflmass" ||
        action == "calc_struct_nlnstiff" || action == "calc_struct_internalforce")
    {
      // do nodal integration of stiffness and internal force
      stiff.LightShape(ndofperpatch, ndofperpatch);
      force.LightSize(ndofperpatch);
      LINALG::SerialDenseMatrix* stiffptr = &stiff;
      if (action == "calc_struct_internalforce") stiffptr = NULL;
      TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::NodalIntegration");
      NodalIntegration(stiffptr, &force, adjnode, adjele, adjsubele, lm, lmlm, *disp, dis, NULL,
          NULL, INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    else if (action == "calc_struct_stress")
    {
      INPAR::STR::StressType iostress =
          DRT::INPUT::get<INPAR::STR::StressType>(p, "iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain =
          DRT::INPUT::get<INPAR::STR::StrainType>(p, "iostrain", INPAR::STR::strain_none);
      std::vector<double> nodalstress(6);
      std::vector<double> nodalstrain(6);
      NodalIntegration(NULL, NULL, adjnode, adjele, adjsubele, lm, lmlm, *disp, dis, &nodalstress,
          &nodalstrain, iostress, iostrain);

      const int lid = dis.NodeRowMap()->LID(nodeLid);
      if (lid == -1) dserror("Cannot find local id for row node");
      for (int i = 0; i < 6; ++i)
      {
        (*(*nstress_)(i))[lid] = nodalstress[i];
        (*(*nstrain_)(i))[lid] = nodalstrain[i];
      }
    }
    else
      dserror("Unknown action");


    //---------------------- do assembly of stiffness and internal force
    // (note: this is non-standard-baci assembly and therefore a do it all yourself version!)
    // there is no guarantee that systemmatrix exists
    // (e.g. if systemmatrix1 is actually a BlockSparseMatrix)
    bool fastassemble = false;
    if (systemmatrix != Teuchos::null) fastassemble = true;

    if (assemblemat1)
    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::PreEvaluate Assembly");
      std::vector<int> lrlm;
      std::vector<int> lclm;

      const Epetra_Map& dofrowmap = systemmatrix1->OperatorRangeMap();
      lrlm.resize(ndofperpatch);
      for (int i = 0; i < ndofperpatch; ++i) lrlm[i] = dofrowmap.LID(lm[i]);
      if (fastassemble)
      {
        const Epetra_Map& dofcolmap = systemmatrix->ColMap();
        lclm.resize(ndofperpatch);
        for (int i = 0; i < ndofperpatch; ++i) lclm[i] = dofcolmap.LID(lm[i]);
      }

      for (int i = 0; i < ndofperpatch; ++i)
      {
        if (lrlm[i] == -1)  // off-processor row
        {
          for (int j = 0; j < ndofperpatch; ++j)
          {
            int errone = stifftmp->SumIntoGlobalValues(1, &lm[i], 1, &lm[j], &stiff(i, j));
            if (errone > 0)
            {
              int errtwo = stifftmp->InsertGlobalValues(1, &lm[i], 1, &lm[j], &stiff(i, j));
              if (errtwo < 0)
                dserror("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d", errtwo);
            }
            else if (errone)
              dserror("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d", errone);
          }
        }
        else  // local row
        {
          if (systemmatrix != Teuchos::null && systemmatrix->Filled())  // matrix is SparseMatrix
          {
            Epetra_CrsMatrix& matrix = *(systemmatrix->EpetraMatrix());
            int length;
            double* values;
            int* indices;
            matrix.ExtractMyRowView(lrlm[i], length, values, indices);
            for (int j = 0; j < ndofperpatch; ++j)
            {
              int* loc = std::lower_bound(indices, indices + length, lclm[j]);
              //#ifdef DEBUG
              if (*loc != lclm[j]) dserror("Cannot find local column entry %d", lclm[j]);
              //#endif
              int pos = loc - indices;

              // test physical continuity of nodal values inside the Epetra_CrsMatrix
              bool continuous = true;
              for (int k = 1; k < 3; ++k)
                if (indices[pos + k] == lclm[j + k])
                  continue;
                else
                {
                  continuous = false;
                  break;
                }

              if (continuous)
              {
                values[pos++] += stiff(i, j++);
                values[pos++] += stiff(i, j++);
                values[pos] += stiff(i, j);
              }
              else
              {
                int err = matrix.SumIntoMyValues(lrlm[i], 1, &stiff(i, j), &lclm[j]);
                j++;
                err += matrix.SumIntoMyValues(lrlm[i], 1, &stiff(i, j), &lclm[j]);
                j++;
                err += matrix.SumIntoMyValues(lrlm[i], 1, &stiff(i, j), &lclm[j]);
                if (err) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d", err);
              }
            }
          }
          else  // matrix not SparseMatrix (e.g. BlockMatrix) -> fall back to standard assembly
          {
            for (int j = 0; j < ndofperpatch; ++j)
              systemmatrix1->Assemble(stiff(i, j), lm[i], lm[j]);
          }
        }
      }
    }

    //-----------------------------------------------------------------------------------
    if (assemblevec1)
    {
      for (int i = 0; i < ndofperpatch; ++i)
      {
        const int rgid = lm[i];
        const int lid = forcetmp1.Map().LID(rgid);
        if (lid < 0) dserror("global row %d does not exist in column map", rgid);
        forcetmp1[lid] += force[i];
      }
    }

    //=========================================================================
  }  // for (node=noderids_.begin(); node != noderids_.end(); ++node)

  //-------------------------------------------------------------------------
  if (action == "calc_struct_stress")
  {
    // we have to export the nodal stresses and strains to column map
    // so they can be written by the elements
    Teuchos::RCP<Epetra_MultiVector> tmp =
        Teuchos::rcp(new Epetra_MultiVector(*dis.NodeColMap(), 6, false));
    LINALG::Export(*nstress_, *tmp);
    nstress_ = tmp;
    tmp = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeColMap(), 6, false));
    LINALG::Export(*nstrain_, *tmp);
    nstrain_ = tmp;
  }


  //-------------------------------------------------------------------------
  // need to export forcetmp to systemvector1 and insert stiffnesses from stifftmp
  // into systemmatrix1
  // Note that fillComplete is never called on stifftmp
  if (assemblevec1)
  {
    Epetra_Vector tmp(systemvector1->Map(), false);
    Epetra_Export exporter(forcetmp1.Map(), tmp.Map());
    int err = tmp.Export(forcetmp1, exporter, Add);
    if (err) dserror("Export using exporter returned err=%d", err);
    systemvector1->Update(1.0, tmp, 1.0);
  }
  if (assemblemat1)
  {
    int err = stifftmp->GlobalAssemble(*dmap, *rmap, false);
    if (err) dserror("Epetra_FECrsMatrix::GlobalAssemble returned err=%d", err);
    const Epetra_Map& cmap = stifftmp->ColMap();
    for (int lrow = 0; lrow < stifftmp->NumMyRows(); ++lrow)
    {
      int numentries;
      double* values;
      if (!stifftmp->Filled())
      {
        const int grow = stifftmp->RowMap().GID(lrow);
        int* gindices;
        int err = stifftmp->ExtractGlobalRowView(grow, numentries, values, gindices);
        if (err) dserror("Epetra_FECrsMatrix::ExtractGlobalRowView returned err=%d", err);
        for (int j = 0; j < numentries; ++j) systemmatrix1->Assemble(values[j], grow, gindices[j]);
      }
      else
      {
        int* lindices;
        int err = stifftmp->ExtractMyRowView(lrow, numentries, values, lindices);
        if (err) dserror("Epetra_FECrsMatrix::ExtractMyRowView returned err=%d", err);
        if (systemmatrix != Teuchos::null && systemmatrix->Filled())
        {
          Epetra_CrsMatrix& matrix = *systemmatrix->EpetraMatrix();
          for (int j = 0; j < numentries; ++j)
          {
            int err = matrix.SumIntoMyValues(lrow, 1, &values[j], &lindices[j]);
            if (err) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d", err);
          }
        }
        else
        {
          const int grow = stifftmp->RowMap().GID(lrow);
          for (int j = 0; j < numentries; ++j)
            systemmatrix1->Assemble(values[j], grow, cmap.GID(lindices[j]));
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  do nodal integration (public)                              gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::NodalIntegration(Epetra_SerialDenseMatrix* stiff,
    Epetra_SerialDenseVector* force, std::map<int, DRT::Node*>& adjnode,
    std::vector<DRT::ELEMENTS::NStet5*>& adjele, std::map<int, std::vector<int>>& adjsubele,
    std::vector<int>& lm, std::vector<std::vector<std::vector<int>>>& lmlm,
    const Epetra_Vector& disp, DRT::Discretization& dis, std::vector<double>* nodalstress,
    std::vector<double>* nodalstrain, const INPAR::STR::StressType iostress,
    const INPAR::STR::StrainType iostrain)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::NodalIntegration");
  //  typedef Sacado::Fad::DFad<double> FAD; // for first derivs
  //  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD; // for second derivs

  //-------------------------------------------------- standard quantities
  const int ndofinpatch = (int)lm.size();
  const int neleinpatch = (int)adjele.size();

  //------------------------------ see whether materials in patch are equal
  bool matequal = true;
  {
    int mat = adjele[0]->material_;
    for (int i = 1; i < neleinpatch; ++i)
      if (mat != adjele[i]->material_)
      {
        matequal = false;
        break;
      }
  }

  //-----------------------------------------------------------------------
  // get displacements of this patch
  std::vector<double> patchdisp(ndofinpatch);
  //  std::vector<FADFAD> tpatchdisp(ndofinpatch);
  for (int i = 0; i < ndofinpatch; ++i)
  {
    int lid = disp.Map().LID(lm[i]);
    if (lid == -1) dserror("Cannot find degree of freedom on this proc");
    patchdisp[i] = disp[disp.Map().LID(lm[i])];
    //    tpatchdisp[i] = patchdisp[i];
    //    tpatchdisp[i].diff(i,ndofinpatch);
    //    tpatchdisp[i].val().diff(i,ndofinpatch);
  }

  //-----------------------------------------------------------------------
  // build averaged F and volume of node
  double Vnode = 0.0;
  LINALG::Matrix<3, 3> Fnode(true);
  //  LINALG::TMatrix<FADFAD,3,3> tFnode(true);

  for (int i = 0; i < neleinpatch; ++i)
  {
    DRT::ELEMENTS::NStet5* ele = adjele[i];
    std::vector<int>& subele = adjsubele[ele->Id()];

    for (unsigned j = 0; j < subele.size(); ++j)
    {
      const int subeleid = subele[j];

      // copy subelement displacements to 4x3 format
      //      LINALG::Matrix<4,3> eledispmat(false);
      //      LINALG::TMatrix<FADFAD,4,3> teledispmat(false);
      //      for (int k=0; k<4; ++k)
      //        for (int l=0; l<3; ++l)
      //        {
      //          eledispmat(k,l) = patchdisp[lmlm[i][j][k*3+l]];
      //          teledispmat(k,l) = tpatchdisp[lmlm[i][j][k*3+l]];
      //        }

      // build F from this subelement
      //     LINALG::Matrix<3,3> F(false);
      //     LINALG::TMatrix<FADFAD,3,3> tF(true);
      //     F = ele->BuildF(eledispmat,ele->SubNxyz(subeleid));
      LINALG::Matrix<3, 3> F = ele->SubF(subeleid);
      //     tF = TBuildF(teledispmat,ele->SubNxyz(subeleid));

      // add 1/3 of subelement material volume to this node
      const double V = ele->SubV(subeleid) / 3.0;
      Vnode += V;

      // add to nodal deformation gradient
      F.Scale(V);
      //     tF.Scale(V);
      Fnode += F;
      //     tFnode += tF;

    }  // for (unsigned j=0; j<subele.size(); ++j)
  }    // for (int i=0; i<neleinpatch; ++i)

  // do the actual averaging
  Fnode.Scale(1.0 / Vnode);
  //  tFnode.Scale(1.0/Vnode);


  //-----------------------------------------------------------------------
  // build B operator

  Epetra_SerialDenseMatrix bopbar(6, ndofinpatch);
  //  Teuchos::SerialDenseMatrix<int,FAD> tbopbar(6,ndofinpatch,true);
  Epetra_SerialDenseMatrix nxyzbar(9, ndofinpatch);

  // loop elements in patch
  for (int ele = 0; ele < neleinpatch; ++ele)
  {
    // current element
    DRT::ELEMENTS::NStet5* actele = adjele[ele];
    std::vector<int>& subele = adjsubele[actele->Id()];
    // loop subelements in this element
    for (unsigned j = 0; j < subele.size(); ++j)
    {
      const int subeleid = subele[j];
      double V = actele->SubV(subeleid) / 3;
      V = V / Vnode;

      // get derivatives with respect to X
      const LINALG::Matrix<4, 3>& nxyz = actele->SubNxyz(subeleid);
      for (int k = 0; k < 4; ++k)
      {
        for (int l = 0; l < 3; ++l)
        {
          int index = lmlm[ele][j][k * 3 + l];
          if (index % 3 == 0)
          {
            nxyzbar(l * 3 + 0, index + 0) += V * nxyz(k, l);
            nxyzbar(l * 3 + 1, index + 1) += V * nxyz(k, l);
            nxyzbar(l * 3 + 2, index + 2) += V * nxyz(k, l);
          }
          else if (index % 3 == 1)
          {
            nxyzbar(l * 3 + 0, index - 1) += V * nxyz(k, l);
            nxyzbar(l * 3 + 1, index + 0) += V * nxyz(k, l);
            nxyzbar(l * 3 + 2, index + 1) += V * nxyz(k, l);
          }
          else
          {
            nxyzbar(l * 3 + 0, index - 2) += V * nxyz(k, l);
            nxyzbar(l * 3 + 1, index - 1) += V * nxyz(k, l);
            nxyzbar(l * 3 + 2, index + 0) += V * nxyz(k, l);
          }
        }  // for (int l=0; l<3; ++l)
      }    // for (int k=0; k<4; ++k)
    }      // for (unsigned j=0; j<subele.size();++j)
  }        // for (int ele=0; ele<neleinpatch; ++ele)



  for (int j = 0; j < ndofinpatch; ++j)
  {
    bopbar(0, j) =
        Fnode(0, 0) * nxyzbar(0, j) + Fnode(1, 0) * nxyzbar(1, j) + Fnode(2, 0) * nxyzbar(2, j);
    bopbar(1, j) =
        Fnode(0, 1) * nxyzbar(3, j) + Fnode(1, 1) * nxyzbar(4, j) + Fnode(2, 1) * nxyzbar(5, j);
    bopbar(2, j) =
        Fnode(0, 2) * nxyzbar(6, j) + Fnode(1, 2) * nxyzbar(7, j) + Fnode(2, 2) * nxyzbar(8, j);
    bopbar(3, j) = Fnode(0, 1) * nxyzbar(0, j) + Fnode(0, 0) * nxyzbar(3, j) +
                   Fnode(1, 1) * nxyzbar(1, j) + Fnode(1, 0) * nxyzbar(4, j) +
                   Fnode(2, 1) * nxyzbar(2, j) + Fnode(2, 0) * nxyzbar(5, j);
    bopbar(4, j) = Fnode(0, 2) * nxyzbar(3, j) + Fnode(0, 1) * nxyzbar(6, j) +
                   Fnode(1, 2) * nxyzbar(4, j) + Fnode(1, 1) * nxyzbar(7, j) +
                   Fnode(2, 2) * nxyzbar(5, j) + Fnode(2, 1) * nxyzbar(8, j);
    bopbar(5, j) = Fnode(0, 2) * nxyzbar(0, j) + Fnode(0, 0) * nxyzbar(6, j) +
                   Fnode(1, 2) * nxyzbar(1, j) + Fnode(1, 0) * nxyzbar(7, j) +
                   Fnode(2, 2) * nxyzbar(2, j) + Fnode(2, 0) * nxyzbar(8, j);
  }



  //-------------------------------------------------------------- averaged strain
  // right cauchy green

  LINALG::Matrix<3, 3> cauchygreen;
  //  LINALG::TMatrix<FADFAD,3,3> tcauchygreen;
  cauchygreen.MultiplyTN(Fnode, Fnode);
  //  tcauchygreen.MultiplyTN(tFnode,tFnode);

  // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::Matrix<6, 1> glstrain(true);
  glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
  glstrain(3) = cauchygreen(0, 1);
  glstrain(4) = cauchygreen(1, 2);
  glstrain(5) = cauchygreen(2, 0);

#if 0
  LINALG::TMatrix<FADFAD,6,1> tglstrain(false);
  tglstrain(0) = 0.5 * (tcauchygreen(0,0) - 1.0);
  tglstrain(1) = 0.5 * (tcauchygreen(1,1) - 1.0);
  tglstrain(2) = 0.5 * (tcauchygreen(2,2) - 1.0);
  tglstrain(3) = tcauchygreen(0,1);
  tglstrain(4) = tcauchygreen(1,2);
  tglstrain(5) = tcauchygreen(2,0);

  for (int i=0; i<ndofinpatch; ++i)
    for (int k=0; k<6; ++k)
      tbopbar(k,i) = tglstrain(k).dx(i);
#endif

  //-------------------------------------------------------- output of strain
  if (iostrain != INPAR::STR::strain_none)
  {
#ifndef PUSO_NSTET5
    StrainOutput(iostrain, *nodalstrain, Fnode, Fnode.Determinant(), 1.0, 1.0 - ALPHA_NSTET5);
#else
    StrainOutput(iostrain, *nodalstrain, Fnode, glstrain, 1.0 - ALPHA_NSTET5);
#endif
  }

  //----------------------------------------- averaged material and stresses
  LINALG::Matrix<6, 6> cmat(true);
  LINALG::Matrix<6, 1> stress(true);

  //-----------------------------------------------------------------------
  // material law
  if (matequal)  // element patch has single material
  {
    double density;  // just a dummy density
    Teuchos::RCP<MAT::Material> mat = adjele[0]->Material();
    // EleGID is set to -1 errorcheck is performed in
    // MAT::Evaluate. I.e if we have elementwise mat params you will catch an error
    SelectMaterial(mat, stress, cmat, density, glstrain, Fnode, 0, -1);
  }
  else
  {
    double density;  // just a dummy density
    LINALG::Matrix<6, 6> cmatele;
    LINALG::Matrix<6, 1> stressele;
    for (int ele = 0; ele < neleinpatch; ++ele)
    {
      cmatele = 0.0;
      stressele = 0.0;
      // current element
      DRT::ELEMENTS::NStet5* actele = adjele[ele];
      // volume of that element assigned to node L
      double V = 0.0;
      for (unsigned j = 0; j < adjsubele[actele->Id()].size(); ++j)
        V += (actele->SubV(adjsubele[actele->Id()][j]) / 3.0);
      // material of the element
      Teuchos::RCP<MAT::Material> mat = actele->Material();
      // EleGID is set to -1 errorcheck is performed in
      // MAT::Evaluate. I.e if we have elementwise mat params you will catch an error
      SelectMaterial(mat, stressele, cmatele, density, glstrain, Fnode, 0, -1);
      cmat.Update(V, cmatele, 1.0);
      stress.Update(V, stressele, 1.0);
    }  // for (int ele=0; ele<neleinpatch; ++ele)
    stress.Scale(1.0 / Vnode);
    cmat.Scale(1.0 / Vnode);
  }

  //-----------------------------------------------------------------------
  // stress is split as follows:
  // stress = vol_node + (1-alpha) * dev_node + alpha * dev_ele
#ifndef PUSO_NSTET5
  {
    LINALG::Matrix<6, 1> stressdev(true);
    LINALG::Matrix<6, 6> cmatdev(true);
    LINALG::Matrix<6, 1> stressvol(true);
    LINALG::Matrix<6, 6> cmatvol(true);

    // compute deviatoric stress and tangent from total stress and tangent
    DevStressTangent(stressdev, cmatdev, cmat, stress, cauchygreen);

    // compute volumetric stress and tangent
    stressvol.Update(-1.0, stressdev, 1.0, stress, 0.0);
    cmatvol.Update(-1.0, cmatdev, 1.0, cmat, 0.0);

    // compute nodal stress
    stress.Update(1.0, stressvol, 1 - ALPHA_NSTET5, stressdev, 0.0);
    cmat.Update(1.0, cmatvol, 1 - ALPHA_NSTET5, cmatdev, 0.0);
  }
#else
  {
    stress.Scale(1. - ALPHA_NSTET5);
    cmat.Scale(1. - ALPHA_NSTET5);
  }
#endif
  //-----------------------------------------------------------------------
  // stress output
  if (iostress != INPAR::STR::stress_none)
  {
    StressOutput(iostress, *nodalstress, stress, Fnode, Fnode.Determinant());
  }

  //  Teuchos::SerialDenseMatrix<int,FAD> tforce(ndofinpatch,1,false);
  //----------------------------------------------------- internal forces
  if (force)
  {
    Epetra_SerialDenseVector stress_epetra(::View, stress.A(), stress.Rows());
    force->Multiply('T', 'N', Vnode, bopbar, stress_epetra, 0.0);

#if 0
    for (int i=0; i<ndofinpatch; ++i)
    {
      FAD sum = 0.0;
      for (int j=0; j<6; ++j)
        sum += tbopbar(j,i) * stress(j);
      tforce(i,0) = Vnode * sum;
    }
#endif
  }
  //--------------------------------------------------- elastic stiffness
  if (stiff)
  {
    Epetra_SerialDenseMatrix cmat_epetra(
        ::View, cmat.A(), cmat.Rows(), cmat.Rows(), cmat.Columns());
    LINALG::SerialDenseMatrix cb(6, ndofinpatch);
    cb.Multiply('N', 'N', 1.0, cmat_epetra, bopbar, 0.0);
    stiff->Multiply('T', 'N', Vnode, bopbar, cb, 0.0);
  }

  //----------------------------------------------------- geom. stiffness

  if (stiff)
  {
    if (!force) dserror("Cannot compute stiffness matrix without computing internal force");

    Teuchos::SerialDenseMatrix<int, double> kg(ndofinpatch, ndofinpatch, true);
    for (int m = 0; m < ndofinpatch; ++m)
    {
      for (int n = 0; n < ndofinpatch; ++n)
      {
        for (int i = 0; i < 3; ++i)
        {
          double sum = Vnode * (stress(0) * nxyzbar(i, n) * nxyzbar(i, m) +
                                   stress(1) * nxyzbar(i + 1 * 3, n) * nxyzbar(i + 1 * 3, m) +
                                   stress(2) * nxyzbar(i + 2 * 3, n) * nxyzbar(i + 2 * 3, m) +
                                   stress(3) * (nxyzbar(i, n) * nxyzbar(i + 1 * 3, m) +
                                                   nxyzbar(i + 1 * 3, n) * nxyzbar(i, m)) +
                                   stress(4) * (nxyzbar(i + 1 * 3, n) * nxyzbar(i + 2 * 3, m) +
                                                   nxyzbar(i + 2 * 3, n) * nxyzbar(i + 1 * 3, m)) +
                                   stress(5) * (nxyzbar(i, n) * nxyzbar(i + 2 * 3, m) +
                                                   nxyzbar(i + 2 * 3, n) * nxyzbar(i, m)));
          (*stiff)(m, n) += sum;
          kg(m, n) += sum;

        }  // for (int i=0; i<3; ++i)
      }    // for (int n=0; n<ndofinpatch; ++n)
    }      // for (int m=0; m<ndofinpatch; ++m)

#if 0
    Teuchos::SerialDenseMatrix<int,double> tkg(ndofinpatch,ndofinpatch,false);
    for (int i=0; i<ndofinpatch; ++i)
      for (int j=0; j<ndofinpatch; ++j)
        tkg(i,j) = tforce(i,0).dx(j);
#endif

  }  // if (stiff)

  return;
}



/*----------------------------------------------------------------------*
 | material laws for NStet5 (protected)                        gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::SelectMaterial(Teuchos::RCP<MAT::Material> mat,
    LINALG::Matrix<6, 1>& stress, LINALG::Matrix<6, 6>& cmat, double& density,
    LINALG::Matrix<6, 1>& glstrain, LINALG::Matrix<3, 3>& defgrd, int gp, const int eleGID)
{
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
      Teuchos::ParameterList params;
      LINALG::Matrix<3, 3> defgrd(true);
      stvk->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, eleGID);
      density = stvk->Density();
    }
    break;
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain, cmat, stress);
      density = neo->Density();
    }
    break;
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see
                                       Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast<MAT::AAAneohooke*>(mat.get());
      Teuchos::ParameterList params;
      aaa->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, eleGID);
      density = aaa->Density();
    }
    break;
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast<MAT::ElastHyper*>(mat.get());
      Teuchos::ParameterList params;
      hyper->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, eleGID);
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
void DRT::ELEMENTS::NStet5Type::DevStressTangent(LINALG::Matrix<6, 1>& Sdev,
    LINALG::Matrix<6, 6>& CCdev, LINALG::Matrix<6, 6>& CC, const LINALG::Matrix<6, 1>& S,
    const LINALG::Matrix<3, 3>& C)
{
  //---------------------------------- things that we'll definitely need
  // inverse of C
  LINALG::Matrix<3, 3> Cinv;
  const double detC = Cinv.Invert(C);
  // J = det(F) = sqrt(detC)
  const double J = sqrt(detC);

  // S as a 3x3 matrix
  LINALG::Matrix<3, 3> Smat;
  Smat(0, 0) = S(0);
  Smat(0, 1) = S(3);
  Smat(0, 2) = S(5);
  Smat(1, 0) = Smat(0, 1);
  Smat(1, 1) = S(1);
  Smat(1, 2) = S(4);
  Smat(2, 0) = Smat(0, 2);
  Smat(2, 1) = Smat(1, 2);
  Smat(2, 2) = S(2);

  //--------------------------------------------- pressure p = -1/(3J) S:C
  double p = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) p += Smat(i, j) * C(i, j);
  p *= (-1. / (3. * J));

  //-------------------------------- compute volumetric PK2 Svol = -p J Cinv
  //-------------------------------------------------------- Sdev = S - Svol
  const double fac = -p * J;
  Sdev(0) = Smat(0, 0) - fac * Cinv(0, 0);
  Sdev(1) = Smat(1, 1) - fac * Cinv(1, 1);
  Sdev(2) = Smat(2, 2) - fac * Cinv(2, 2);
  Sdev(3) = Smat(0, 1) - fac * Cinv(0, 1);
  Sdev(4) = Smat(1, 2) - fac * Cinv(1, 2);
  Sdev(5) = Smat(0, 2) - fac * Cinv(0, 2);

  //======================================== volumetric tangent matrix CCvol
  LINALG::Matrix<6, 6> CCvol(true);  // fill with zeros

  //--------------------------------------- CCvol += 2pJ (Cinv boeppel Cinv)
  MAT::ElastSymTensor_o_Multiply(CCvol, -2.0 * fac, Cinv, Cinv, 0.0);

  //------------------------------------------ CCvol += 2/3 * Cinv dyad S
  MAT::ElastSymTensorMultiply(CCvol, 2.0 / 3.0, Cinv, Smat, 1.0);

  //-------------------------------------- CCvol += 1/3 Cinv dyad ( CC : C )
  {
    // C as Voigt vector
    LINALG::Matrix<6, 1> Cvec;
    Cvec(0) = C(0, 0);
    Cvec(1) = C(1, 1);
    Cvec(2) = C(2, 2);
    Cvec(3) = 2.0 * C(0, 1);
    Cvec(4) = 2.0 * C(1, 2);
    Cvec(5) = 2.0 * C(0, 2);

    LINALG::Matrix<6, 1> CCcolonC;
    CCcolonC.Multiply(CC, Cvec);

    LINALG::Matrix<3, 3> CCcC;
    CCcC(0, 0) = CCcolonC(0);
    CCcC(0, 1) = CCcolonC(3);
    CCcC(0, 2) = CCcolonC(5);
    CCcC(1, 0) = CCcC(0, 1);
    CCcC(1, 1) = CCcolonC(1);
    CCcC(1, 2) = CCcolonC(4);
    CCcC(2, 0) = CCcC(0, 2);
    CCcC(2, 1) = CCcC(1, 2);
    CCcC(2, 2) = CCcolonC(2);
    MAT::ElastSymTensorMultiply(CCvol, 1. / 3., Cinv, CCcC, 1.0);
  }

  //----------------------------------------------------- CCdev = CC - CCvol
  CCdev.Update(1.0, CC, -1.0, CCvol);

  return;
}

/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::StrainOutput(const INPAR::STR::StrainType iostrain,
    std::vector<double>& nodalstrain, LINALG::Matrix<3, 3>& F, const double& detF,
    const double volweight, const double devweight)
{
  LINALG::Matrix<3, 3> Fiso = F;
  Fiso.Scale(pow(detF, -1.0 / 3.0));

  LINALG::Matrix<3, 3> Fvol(true);
  Fvol(0, 0) = 1.0;
  Fvol(1, 1) = 1.0;
  Fvol(2, 2) = 1.0;
  Fvol.Scale(pow(detF, 1.0 / 3.0));

  LINALG::Matrix<3, 3> cauchygreeniso(false);
  cauchygreeniso.MultiplyTN(Fiso, Fiso);

  LINALG::Matrix<3, 3> cauchygreenvol(false);
  cauchygreenvol.MultiplyTN(Fvol, Fvol);

  LINALG::Matrix<3, 3> glstrainiso(false);
  glstrainiso(0, 0) = 0.5 * (cauchygreeniso(0, 0) - 1.0);
  glstrainiso(0, 1) = 0.5 * cauchygreeniso(0, 1);
  glstrainiso(0, 2) = 0.5 * cauchygreeniso(0, 2);
  glstrainiso(1, 0) = glstrainiso(0, 1);
  glstrainiso(1, 1) = 0.5 * (cauchygreeniso(1, 1) - 1.0);
  glstrainiso(1, 2) = 0.5 * cauchygreeniso(1, 2);
  glstrainiso(2, 0) = glstrainiso(0, 2);
  glstrainiso(2, 1) = glstrainiso(1, 2);
  glstrainiso(2, 2) = 0.5 * (cauchygreeniso(2, 2) - 1.0);

  LINALG::Matrix<3, 3> glstrainvol(false);
  glstrainvol(0, 0) = 0.5 * (cauchygreenvol(0, 0) - 1.0);
  glstrainvol(0, 1) = 0.5 * cauchygreenvol(0, 1);
  glstrainvol(0, 2) = 0.5 * cauchygreenvol(0, 2);
  glstrainvol(1, 0) = glstrainvol(0, 1);
  glstrainvol(1, 1) = 0.5 * (cauchygreenvol(1, 1) - 1.0);
  glstrainvol(1, 2) = 0.5 * cauchygreenvol(1, 2);
  glstrainvol(2, 0) = glstrainvol(0, 2);
  glstrainvol(2, 1) = glstrainvol(1, 2);
  glstrainvol(2, 2) = 0.5 * (cauchygreenvol(2, 2) - 1.0);

  LINALG::Matrix<3, 3> glstrainout = glstrainiso;
  glstrainout.Update(volweight, glstrainvol, devweight);

  switch (iostrain)
  {
    case INPAR::STR::strain_gl:
    {
      nodalstrain[0] = glstrainout(0, 0);
      nodalstrain[1] = glstrainout(1, 1);
      nodalstrain[2] = glstrainout(2, 2);
      nodalstrain[3] = glstrainout(0, 1);
      nodalstrain[4] = glstrainout(1, 2);
      nodalstrain[5] = glstrainout(0, 2);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      // inverse of deformation gradient
      LINALG::Matrix<3, 3> invdefgrd;
      invdefgrd.Invert(F);
      LINALG::Matrix<3, 3> temp;
      LINALG::Matrix<3, 3> euler_almansi;
      temp.Multiply(glstrainout, invdefgrd);
      euler_almansi.MultiplyTN(invdefgrd, temp);
      nodalstrain[0] = euler_almansi(0, 0);
      nodalstrain[1] = euler_almansi(1, 1);
      nodalstrain[2] = euler_almansi(2, 2);
      nodalstrain[3] = euler_almansi(0, 1);
      nodalstrain[4] = euler_almansi(1, 2);
      nodalstrain[5] = euler_almansi(0, 2);
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
void DRT::ELEMENTS::NStet5Type::StrainOutput(const INPAR::STR::StrainType iostrain,
    std::vector<double>& nodalstrain, LINALG::Matrix<3, 3>& F, LINALG::Matrix<6, 1>& glstrain,
    const double weight)
{
  LINALG::Matrix<3, 3> glstrainout;

  glstrainout(0, 0) = weight * glstrain(0);
  glstrainout(1, 1) = weight * glstrain(1);
  glstrainout(2, 2) = weight * glstrain(2);
  glstrainout(0, 1) = weight * glstrain(3);
  glstrainout(1, 2) = weight * glstrain(4);
  glstrainout(0, 2) = weight * glstrain(5);


  switch (iostrain)
  {
    case INPAR::STR::strain_gl:
    {
      nodalstrain[0] = glstrainout(0, 0);
      nodalstrain[1] = glstrainout(1, 1);
      nodalstrain[2] = glstrainout(2, 2);
      nodalstrain[3] = glstrainout(0, 1);
      nodalstrain[4] = glstrainout(1, 2);
      nodalstrain[5] = glstrainout(0, 2);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      // inverse of deformation gradient
      LINALG::Matrix<3, 3> invdefgrd;
      invdefgrd.Invert(F);
      LINALG::Matrix<3, 3> temp;
      LINALG::Matrix<3, 3> euler_almansi;
      temp.Multiply(glstrainout, invdefgrd);
      euler_almansi.MultiplyTN(invdefgrd, temp);
      nodalstrain[0] = euler_almansi(0, 0);
      nodalstrain[1] = euler_almansi(1, 1);
      nodalstrain[2] = euler_almansi(2, 2);
      nodalstrain[3] = euler_almansi(0, 1);
      nodalstrain[4] = euler_almansi(1, 2);
      nodalstrain[5] = euler_almansi(0, 2);
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
void DRT::ELEMENTS::NStet5Type::StressOutput(const INPAR::STR::StressType iostress,
    std::vector<double>& nodalstress, LINALG::Matrix<6, 1>& stress, LINALG::Matrix<3, 3>& F,
    const double& detF)
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
      LINALG::Matrix<3, 3> pkstress;
      pkstress(0, 0) = stress(0);
      pkstress(0, 1) = stress(3);
      pkstress(0, 2) = stress(5);
      pkstress(1, 0) = pkstress(0, 1);
      pkstress(1, 1) = stress(1);
      pkstress(1, 2) = stress(4);
      pkstress(2, 0) = pkstress(0, 2);
      pkstress(2, 1) = pkstress(1, 2);
      pkstress(2, 2) = stress(2);
      LINALG::Matrix<3, 3> temp;
      LINALG::Matrix<3, 3> cauchystress;
      temp.Multiply(1.0 / detF, F, pkstress);
      cauchystress.MultiplyNT(temp, F);
      nodalstress[0] = cauchystress(0, 0);
      nodalstress[1] = cauchystress(1, 1);
      nodalstress[2] = cauchystress(2, 2);
      nodalstress[3] = cauchystress(0, 1);
      nodalstress[4] = cauchystress(1, 2);
      nodalstress[5] = cauchystress(0, 2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress type not available");
  }
  return;
}
