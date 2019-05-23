/*!----------------------------------------------------------------------*

\brief NStet element

\level 3

\maintainer Christoph Meier
*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_utils_densematrix_determinant.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_FECrsMatrix.h"
#include "Sacado.hpp"

#include "../drt_mat/material_service.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/elasthyper.H"

#include "so_nstet.H"

/*----------------------------------------------------------------------*
 |                                                             gee 10/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::ElementDeformationGradient(DRT::Discretization& dis)
{
  // current displacement
  Teuchos::RCP<const Epetra_Vector> disp = dis.GetState("displacement");
  if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
  // loop elements
  std::map<int, DRT::ELEMENTS::NStet*>::iterator ele;
  for (ele = elecids_.begin(); ele != elecids_.end(); ++ele)
  {
    DRT::ELEMENTS::NStet* e = ele->second;
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    e->LocationVector(dis, lm, lmowner, lmstride);
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

    // create dfad version of nxyz and mydisp
    LINALG::Matrix<4, 3> disp(false);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 3; ++j) disp(i, j) = mydisp[i * 3 + j];

    // create deformation gradient
    e->F() = e->BuildF(disp, e->Nxyz());

  }  // ele
  return;
}


void AutoDiffDemo();
void AutoDiffDemo(DRT::Discretization& dis);

/*----------------------------------------------------------------------*
 |  pre-evaluation of elements (public)                        gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::PreEvaluate");

  // nodal integration for nlnstiff and internal forces only
  // (this method does not compute stresses/strains/element updates/mass matrix)
  std::string& action = p.get<std::string>("action", "none");
  if (action != "calc_struct_nlnstiffmass" && action != "calc_struct_nlnstifflmass" &&
      action != "calc_struct_nlnstiff" && action != "calc_struct_internalforce" &&
      action != "calc_struct_stress")
    return;

  // These get filled in here, so remove old stuff
  if (action == "calc_struct_stress")
  {
    nstress_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), 6, false));
    nstrain_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), 6, false));
    pstab_nstress_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), 6, true));
    pstab_nstrain_ = Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(), 6, true));
  }
  else
  {
    nstress_ = Teuchos::null;
    nstrain_ = Teuchos::null;
    pstab_nstress_ = Teuchos::null;
    pstab_nstrain_ = Teuchos::null;
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
  LINALG::SerialDenseMatrix mis_stiff;
  LINALG::SerialDenseVector mis_force;

  //-------------------------------------- construct F for each NStet
  // AutoDiffDemo();
  // AutoDiffDemo(dis);
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
    bool mis = (pstab_adjele_.find(nodeLid) != pstab_adjele_.end());

    // standard quantities for all nodes
    std::vector<DRT::ELEMENTS::NStet*>& adjele = adjele_[nodeLid];
    std::map<int, DRT::Node*>& adjnode = adjnode_[nodeLid];
    std::vector<int>& lm = adjlm_[nodeLid];
    const int ndofperpatch = (int)lm.size();

    // quantities for mis nodes
    int mis_ndofperpatch = 0;
    std::vector<DRT::ELEMENTS::NStet*>* mis_adjele = NULL;
    std::map<int, DRT::Node*>* mis_adjnode = NULL;
    std::vector<double>* mis_weight = NULL;
    std::vector<int>* mis_lm = NULL;
    if (mis)
    {
      mis_adjele = &pstab_adjele_[nodeLid];
      mis_adjnode = &pstab_adjnode_[nodeLid];
      mis_weight = &pstab_adjele_weight_[nodeLid];
      mis_lm = &pstab_adjlm_[nodeLid];
      mis_ndofperpatch = (int)(*mis_lm).size();
    }

    if (action == "calc_struct_nlnstiffmass" || action == "calc_struct_nlnstifflmass" ||
        action == "calc_struct_nlnstiff" || action == "calc_struct_internalforce")
    {
      // do nodal integration of stiffness and internal force
      stiff.LightShape(ndofperpatch, ndofperpatch);
      force.LightSize(ndofperpatch);
      {
        TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::NodalIntegration");
        NodalIntegration(&stiff, &force, adjnode, adjele, lm, *disp, dis, NULL, NULL,
            INPAR::STR::stress_none, INPAR::STR::strain_none);
      }

#ifndef PUSOSOLBERG
      if (mis)
      {
        TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::MISNodalIntegration");
        mis_stiff.LightShape(
            mis_ndofperpatch, mis_ndofperpatch);  // put in Light once values are present
        mis_force.LightSize(mis_ndofperpatch);    // put in Light once values are present
        MISNodalIntegration(&mis_stiff, &mis_force, *mis_adjnode, *mis_adjele, *mis_weight, *mis_lm,
            *disp, dis, NULL, NULL, INPAR::STR::stress_none, INPAR::STR::strain_none);
      }  // mis
#endif
    }
    else if (action == "calc_struct_stress")
    {
      INPAR::STR::StressType iostress =
          DRT::INPUT::get<INPAR::STR::StressType>(p, "iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain =
          DRT::INPUT::get<INPAR::STR::StrainType>(p, "iostrain", INPAR::STR::strain_none);
      std::vector<double> nodalstress(6);
      std::vector<double> nodalstrain(6);
      NodalIntegration(NULL, NULL, adjnode, adjele, lm, *disp, dis, &nodalstress, &nodalstrain,
          iostress, iostrain);

      const int lid = dis.NodeRowMap()->LID(nodeLid);
      if (lid == -1) dserror("Cannot find local id for row node");
      for (int i = 0; i < 6; ++i)
      {
        (*(*nstress_)(i))[lid] = nodalstress[i];
        (*(*nstrain_)(i))[lid] = nodalstrain[i];
      }

#ifndef PUSOSOLBERG
      if (mis)
      {
        std::vector<double> mis_nodalstress(6);
        std::vector<double> mis_nodalstrain(6);
        MISNodalIntegration(NULL, NULL, *mis_adjnode, *mis_adjele, *mis_weight, *mis_lm, *disp, dis,
            &mis_nodalstress, &mis_nodalstrain, iostress, iostrain);

        for (int i = 0; i < 6; ++i)
        {
          (*(*pstab_nstress_)(i))[lid] = mis_nodalstress[i];
          (*(*pstab_nstrain_)(i))[lid] = mis_nodalstrain[i];
        }
      }  // mis
#endif
    }
    else
      dserror("Unknown type of action");


    //---------------------- do assembly of stiffness and internal force
    // (note: this is non-standard-baci assembly and therefore a do it all yourself version!)
    // there is no guarantee that systemmatrix exists
    // (e.g. if systemmatrix1 is actually a BlockSparseMatrix)
    bool fastassemble = false;
    if (systemmatrix != Teuchos::null) fastassemble = true;

    if (assemblemat1)
    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::PreEvaluate Assembly");
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
#ifdef DEBUG
              if (*loc != lclm[j]) dserror("Cannot find local column entry %d", lclm[j]);
#endif
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

#ifndef PUSOSOLBERG
      if (mis)
      {
        lrlm.resize(mis_ndofperpatch);
        for (int i = 0; i < mis_ndofperpatch; ++i) lrlm[i] = dofrowmap.LID((*mis_lm)[i]);
        if (fastassemble)
        {
          const Epetra_Map& dofcolmap = systemmatrix->ColMap();
          lclm.resize(mis_ndofperpatch);
          for (int i = 0; i < mis_ndofperpatch; ++i) lclm[i] = dofcolmap.LID((*mis_lm)[i]);
        }

        for (int i = 0; i < mis_ndofperpatch; ++i)
        {
          if (lrlm[i] == -1)  // off processor row
          {
            for (int j = 0; j < mis_ndofperpatch; ++j)
            {
              int errone = stifftmp->SumIntoGlobalValues(
                  1, &(*mis_lm)[i], 1, &(*mis_lm)[j], &mis_stiff(i, j));
              if (errone > 0)
              {
                int errtwo = stifftmp->InsertGlobalValues(
                    1, &(*mis_lm)[i], 1, &(*mis_lm)[j], &mis_stiff(i, j));
                if (errtwo < 0)
                  dserror("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d", errtwo);
              }
              else if (errone)
                dserror("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d", errone);
            }
          }
          else
          {
            if (systemmatrix != Teuchos::null && systemmatrix->Filled())  // matrix is SparseMatrix
            {
              Epetra_CrsMatrix& matrix = *(systemmatrix->EpetraMatrix());
              int length;
              double* values;
              int* indices;
              matrix.ExtractMyRowView(lrlm[i], length, values, indices);
              for (int j = 0; j < mis_ndofperpatch; ++j)
              {
                int* loc = std::lower_bound(indices, indices + length, lclm[j]);
#ifdef DEBUG
                if (*loc != lclm[j]) dserror("Cannot find local column entry %d", lclm[j]);
#endif
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
                  values[pos++] += mis_stiff(i, j++);
                  values[pos++] += mis_stiff(i, j++);
                  values[pos] += mis_stiff(i, j);
                }
                else
                {
                  int err = matrix.SumIntoMyValues(lrlm[i], 1, &mis_stiff(i, j), &lclm[j]);
                  j++;
                  err += matrix.SumIntoMyValues(lrlm[i], 1, &mis_stiff(i, j), &lclm[j]);
                  j++;
                  err += matrix.SumIntoMyValues(lrlm[i], 1, &mis_stiff(i, j), &lclm[j]);
                  if (err) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d", err);
                }
              }
            }
            else  // matrix not SparseMatrix (e.g. BlockMatrix) -> fall back to standard assembly
            {
              for (int j = 0; j < mis_ndofperpatch; ++j)
                systemmatrix1->Assemble(mis_stiff(i, j), (*mis_lm)[i], (*mis_lm)[j]);
            }
          }
        }
      }  // mis
#endif
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
#ifndef PUSOSOLBERG
      if (mis)
      {
        for (int i = 0; i < mis_ndofperpatch; ++i)
        {
          const int rgid = (*mis_lm)[i];
          const int lid = forcetmp1.Map().LID(rgid);
          if (lid < 0) dserror("global row %d does not exist in column map", rgid);
          forcetmp1[lid] += mis_force[i];
        }
      }  // mis
#endif
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

#ifndef PUSOSOLBERG
    // export mis stress and strains to mis overlapping map to allow for output
    tmp = Teuchos::rcp(new Epetra_MultiVector(*pstab_misstressout_, 6, true));
    LINALG::Export(*pstab_nstress_, *tmp);
    pstab_nstress_ = tmp;
    tmp = Teuchos::rcp(new Epetra_MultiVector(*pstab_misstressout_, 6, true));
    LINALG::Export(*pstab_nstrain_, *tmp);
    pstab_nstrain_ = tmp;
#endif
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
 |  do nodal integration (public)                              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::NodalIntegration(Epetra_SerialDenseMatrix* stiff,
    Epetra_SerialDenseVector* force, std::map<int, DRT::Node*>& adjnode,
    std::vector<DRT::ELEMENTS::NStet*>& adjele, std::vector<int>& lm, const Epetra_Vector& disp,
    DRT::Discretization& dis, std::vector<double>* nodalstress, std::vector<double>* nodalstrain,
    const INPAR::STR::StressType iostress, const INPAR::STR::StrainType iostrain)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::NodalIntegration");
  typedef Sacado::Fad::DFad<double> FAD;  // for first derivs
  // typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD; // for second derivs

  //-------------------------------------------------- standard quantities
  const int nnodeinpatch = (int)adjnode.size();
  const int ndofinpatch = nnodeinpatch * 3;
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
  std::vector<FAD> patchdisp(ndofinpatch);
  for (int i = 0; i < ndofinpatch; ++i)
  {
    int lid = disp.Map().LID(lm[i]);
    if (lid == -1) dserror("Cannot find degree of freedom on this proc");
    patchdisp[i] = disp[disp.Map().LID(lm[i])];
    patchdisp[i].diff(i, ndofinpatch);
  }

  //-----------------------------------------------------------------------
  // build averaged F and volume of node (using sacado)
  double VnodeL = 0.0;
  LINALG::TMatrix<FAD, 3, 3> fad_FnodeL(true);
  std::vector<std::vector<int>> lmlm(neleinpatch);
  for (int i = 0; i < neleinpatch; ++i)
  {
    const double V = adjele[i]->Vol() / 4;
    VnodeL += V;

    // get the element's displacements out of the patch' displacements
    std::vector<int> elelm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    adjele[i]->LocationVector(dis, elelm, lmowner, lmstride);

    // have to find position of elelm[i] in lm
    // lmlm[i][j] : element i, degree of freedom j, lmlm[i][j] position in patchdisp[0..ndofinpatch]
    lmlm[i].resize(12);
    for (int j = 0; j < 12; ++j)
    {
      std::vector<int>::iterator k = find(lm.begin(), lm.end(), elelm[j]);
      lmlm[i][j] = k - lm.begin();  // the position of elelm[j] in lm
    }

    // copy element disp to 4x3 format
    LINALG::TMatrix<FAD, 4, 3> eledispmat(false);
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 3; ++k) eledispmat(j, k) = patchdisp[lmlm[i][j * 3 + k]];

    // build F of this element
    LINALG::TMatrix<FAD, 3, 3> Fele(true);
    Fele = adjele[i]->BuildF<FAD>(eledispmat, adjele[i]->Nxyz());

    // add up to nodal deformation gradient
    Fele.Scale(V);
    fad_FnodeL += Fele;

  }  // for (int i=0; i<neleinpatch; ++i)

  // do averaging
  fad_FnodeL.Scale(1.0 / VnodeL);

  // copy values of fad to 'normal' values
  LINALG::Matrix<3, 3> FnodeL(false);
  for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k) FnodeL(j, k) = fad_FnodeL(j, k).val();

  //-----------------------------------------------------------------------
  // build B operator
  Epetra_SerialDenseMatrix bop(6, ndofinpatch);
  // loop elements in patch
  for (int ele = 0; ele < neleinpatch; ++ele)
  {
    // current element
    DRT::ELEMENTS::NStet* actele = adjele[ele];

    // volume of that element assigned to node L
    double V = actele->Vol() / 4;

    // volume ratio of volume per node of this element to
    // whole volume of node L
    const double ratio = V / VnodeL;

    // get derivatives with respect to X
    LINALG::Matrix<4, 3>& nxyz = actele->Nxyz();

    // get defgrd
    LINALG::Matrix<3, 3>& F = actele->F();

    LINALG::Matrix<6, 12> bele(false);
    for (int i = 0; i < 4; ++i)
    {
      bele(0, 3 * i + 0) = F(0, 0) * nxyz(i, 0);
      bele(0, 3 * i + 1) = F(1, 0) * nxyz(i, 0);
      bele(0, 3 * i + 2) = F(2, 0) * nxyz(i, 0);
      bele(1, 3 * i + 0) = F(0, 1) * nxyz(i, 1);
      bele(1, 3 * i + 1) = F(1, 1) * nxyz(i, 1);
      bele(1, 3 * i + 2) = F(2, 1) * nxyz(i, 1);
      bele(2, 3 * i + 0) = F(0, 2) * nxyz(i, 2);
      bele(2, 3 * i + 1) = F(1, 2) * nxyz(i, 2);
      bele(2, 3 * i + 2) = F(2, 2) * nxyz(i, 2);

      bele(3, 3 * i + 0) = F(0, 0) * nxyz(i, 1) + F(0, 1) * nxyz(i, 0);
      bele(3, 3 * i + 1) = F(1, 0) * nxyz(i, 1) + F(1, 1) * nxyz(i, 0);
      bele(3, 3 * i + 2) = F(2, 0) * nxyz(i, 1) + F(2, 1) * nxyz(i, 0);
      bele(4, 3 * i + 0) = F(0, 1) * nxyz(i, 2) + F(0, 2) * nxyz(i, 1);
      bele(4, 3 * i + 1) = F(1, 1) * nxyz(i, 2) + F(1, 2) * nxyz(i, 1);
      bele(4, 3 * i + 2) = F(2, 1) * nxyz(i, 2) + F(2, 2) * nxyz(i, 1);
      bele(5, 3 * i + 0) = F(0, 2) * nxyz(i, 0) + F(0, 0) * nxyz(i, 2);
      bele(5, 3 * i + 1) = F(1, 2) * nxyz(i, 0) + F(1, 0) * nxyz(i, 2);
      bele(5, 3 * i + 2) = F(2, 2) * nxyz(i, 0) + F(2, 0) * nxyz(i, 2);
    }

    for (int k = 0; k < 6; ++k)
      for (int j = 0; j < 12; ++j) bop(k, lmlm[ele][j]) += ratio * bele(k, j);
  }  // for (int ele=0; ele<neleinpatch; ++ele)

  //-------------------------------------------------------------- averaged strain
  // right cauchy green
  LINALG::TMatrix<FAD, 3, 3> CG(false);
  CG.MultiplyTN(fad_FnodeL, fad_FnodeL);
  std::vector<FAD> Ebar(6);
  Ebar[0] = 0.5 * (CG(0, 0) - 1.0);
  Ebar[1] = 0.5 * (CG(1, 1) - 1.0);
  Ebar[2] = 0.5 * (CG(2, 2) - 1.0);
  Ebar[3] = CG(0, 1);
  Ebar[4] = CG(1, 2);
  Ebar[5] = CG(2, 0);

  // for material law and output, copy to LINALG::Matrix object
  LINALG::Matrix<3, 3> cauchygreen(false);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) cauchygreen(i, j) = CG(i, j).val();
  LINALG::Matrix<6, 1> glstrain(false);
  glstrain(0) = Ebar[0].val();
  glstrain(1) = Ebar[1].val();
  glstrain(2) = Ebar[2].val();
  glstrain(3) = Ebar[3].val();
  glstrain(4) = Ebar[4].val();
  glstrain(5) = Ebar[5].val();

  //-------------------------------------------------------- output of strain
  if (iostrain != INPAR::STR::strain_none)
  {
#ifndef PUSOSOLBERG
    StrainOutput(
        iostrain, *nodalstrain, FnodeL, FnodeL.Determinant(), 1.0 - BETA_NSTET, 1.0 - ALPHA_NSTET);
#else
    StrainOutput(iostrain, *nodalstrain, FnodeL, glstrain, 1.0 - ALPHA_NSTET);
#endif
  }

  //-------------------------------------------------------------------------
  // build a second B-operator from the averaged strains that are based on
  // the averaged F
  Epetra_SerialDenseMatrix bopbar(6, ndofinpatch);
  for (int i = 0; i < ndofinpatch; ++i)
    for (int k = 0; k < 6; ++k) bopbar(k, i) = Ebar[k].fastAccessDx(i);

#if 0
  Epetra_SerialDenseMatrix bopbar2(6,ndofinpatch);
  {

    Epetra_SerialDenseMatrix Graddm(9,ndofinpatch);
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      LINALG::Matrix<4,3>& nxyz = adjele[ele]->Nxyz();
      const double V = (adjele[ele]->Vol()/4.);
      LINALG::Matrix<9,12> Gradde(true);
      for (int i=0; i<4; ++i) // FIXME: do += ratio * directly into Graddm for efficiency
      {
        Gradde(0,3*i+0) = nxyz(i,0);
        //Gradde(0,3*i+1) = 0.0;
        //Gradde(0,3*i+2) = 0.0;

        //Gradde(1,3*i+0) = 0.0;
        Gradde(1,3*i+1) = nxyz(i,1);
        //Gradde(1,3*i+2) = 0.0;

        //Gradde(2,3*i+0) = 0.0;
        //Gradde(2,3*i+1) = 0.0;
        Gradde(2,3*i+2) = nxyz(i,2);

        Gradde(3,3*i+0) = nxyz(i,1);
        //Gradde(3,3*i+1) = 0.0;
        //Gradde(3,3*i+2) = 0.0;
        Gradde(4,3*i+0) = nxyz(i,2);
        //Gradde(4,3*i+1) = 0.0;
        //Gradde(4,3*i+2) = 0.0;

        //Gradde(5,3*i+0) = 0.0;
        Gradde(5,3*i+1) = nxyz(i,0);
        //Gradde(5,3*i+2) = 0.0;
        //Gradde(6,3*i+0) = 0.0;
        Gradde(6,3*i+1) = nxyz(i,2);
        //Gradde(6,3*i+2) = 0.0;

        //Gradde(7,3*i+0) = 0.0;
        //Gradde(7,3*i+1) = 0.0;
        Gradde(7,3*i+2) = nxyz(i,0);
        //Gradde(8,3*i+0) = 0.0;
        //Gradde(8,3*i+1) = 0.0;
        Gradde(8,3*i+2) = nxyz(i,1);
      } // for (int i=0; i<4; ++i)
      for (int i=0; i<9; ++i)
        for (int j=0; j<12; ++j)
          Graddm(i,lmlm[ele][j]) += V * Gradde(i,j);
    } // for (int ele=0; ele<neleinpatch; ++ele)
    Graddm.Scale(1./VnodeL);
    for (int j=0; j<ndofinpatch; ++j)
    {
      bopbar2(0,j) = FnodeL(0,0)*Graddm(0,j) + FnodeL(1,0)*Graddm(5,j) + FnodeL(2,0)*Graddm(7,j);
      bopbar2(1,j) = FnodeL(0,1)*Graddm(3,j) + FnodeL(1,1)*Graddm(1,j) + FnodeL(2,1)*Graddm(8,j);
      bopbar2(2,j) = FnodeL(0,2)*Graddm(4,j) + FnodeL(1,2)*Graddm(6,j) + FnodeL(2,2)*Graddm(2,j);

      bopbar2(3,j) = FnodeL(0,1)*Graddm(0,j) + FnodeL(1,1)*Graddm(5,j) + FnodeL(2,1)*Graddm(7,j) +
                     FnodeL(0,0)*Graddm(3,j) + FnodeL(1,0)*Graddm(1,j) + FnodeL(2,0)*Graddm(8,j);

      bopbar2(4,j) = FnodeL(0,2)*Graddm(3,j) + FnodeL(1,2)*Graddm(1,j) + FnodeL(2,2)*Graddm(8,j) +
                     FnodeL(0,1)*Graddm(4,j) + FnodeL(1,1)*Graddm(6,j) + FnodeL(2,1)*Graddm(2,j);

      bopbar2(5,j) = FnodeL(0,0)*Graddm(4,j) + FnodeL(1,0)*Graddm(6,j) + FnodeL(2,0)*Graddm(2,j) +
                     FnodeL(0,2)*Graddm(0,j) + FnodeL(1,2)*Graddm(5,j) + FnodeL(2,2)*Graddm(7,j);
    }

  //bopbar = bopbar2;

  }
#endif

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
    SelectMaterial(mat, stress, cmat, density, glstrain, FnodeL, 0, -1);
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
      DRT::ELEMENTS::NStet* actele = adjele[ele];
      // volume of that element assigned to node L
      const double V = actele->Vol() / 4;
      // def-gradient of the element
      Teuchos::RCP<MAT::Material> mat = actele->Material();
      // EleGID is set to -1 errorcheck is performed in
      // MAT::Evaluate. I.e if we have elementwise mat params you will catch an error
      SelectMaterial(mat, stressele, cmatele, density, glstrain, FnodeL, 0, -1);
      cmat.Update(V, cmatele, 1.0);
      stress.Update(V, stressele, 1.0);
    }  // for (int ele=0; ele<neleinpatch; ++ele)
    stress.Scale(1.0 / VnodeL);
    cmat.Scale(1.0 / VnodeL);
  }

  //-----------------------------------------------------------------------
  // stress is split as follows:
  // stress = beta * vol_misnode + (1-beta) * vol_node + (1-alpha) * dev_node + alpha * dev_ele
#ifndef PUSOSOLBERG
  {
    LINALG::Matrix<6, 1> stressdev(true);
    LINALG::Matrix<6, 6> cmatdev(true);
    LINALG::Matrix<6, 1> stressvol(false);
    LINALG::Matrix<6, 6> cmatvol(false);

    // compute deviatoric stress and tangent from total stress and tangent
    DevStressTangent(stressdev, cmatdev, cmat, stress, cauchygreen);

    // compute volumetric stress and tangent
    stressvol.Update(-1.0, stressdev, 1.0, stress, 0.0);
    cmatvol.Update(-1.0, cmatdev, 1.0, cmat, 0.0);

    // compute nodal stress
    stress.Update(1.0 - BETA_NSTET, stressvol, 1 - ALPHA_NSTET, stressdev, 0.0);
    cmat.Update(1.0 - BETA_NSTET, cmatvol, 1 - ALPHA_NSTET, cmatdev, 0.0);
  }
#else
  {
    stress.Scale(1. - ALPHA_NSTET);
    cmat.Scale(1. - ALPHA_NSTET);
  }
#endif

  //-----------------------------------------------------------------------
  // stress output
  if (iostress != INPAR::STR::stress_none)
  {
    StressOutput(iostress, *nodalstress, stress, FnodeL, FnodeL.Determinant());
  }

  //----------------------------------------------------- internal forces
  if (force)
  {
    Epetra_SerialDenseVector stress_epetra(::View, stress.A(), stress.Rows());
    force->Multiply('T', 'N', VnodeL, bop, stress_epetra, 0.0);  // bop
  }

  //--------------------------------------------------- elastic stiffness
  if (stiff)
  {
    Epetra_SerialDenseMatrix cmat_epetra(
        ::View, cmat.A(), cmat.Rows(), cmat.Rows(), cmat.Columns());
    LINALG::SerialDenseMatrix cb(6, ndofinpatch);
    cb.Multiply('N', 'N', 1.0, cmat_epetra, bopbar, 0.0);
    stiff->Multiply('T', 'N', VnodeL, bop, cb, 0.0);  // bop
  }

  //----------------------------------------------------- geom. stiffness
  // do not use sacado for second derivative of E as it is way too expensive!
  // As long as the 2nd deriv is as easy as this, do it by hand
  if (stiff)
  {
    // loop elements in patch
    for (int ele = 0; ele < neleinpatch; ++ele)
    {
      // material deriv of element
      LINALG::Matrix<4, 3>& nxyz = adjele[ele]->Nxyz();
      // volume of element assigned to node L
      const double V = adjele[ele]->Vol() / 4;

      // loop nodes of that element
      double SmBL[3];
      for (int i = 0; i < 4; ++i)
      {
        SmBL[0] = V * (stress(0) * nxyz(i, 0) + stress(3) * nxyz(i, 1) + stress(5) * nxyz(i, 2));
        SmBL[1] = V * (stress(3) * nxyz(i, 0) + stress(1) * nxyz(i, 1) + stress(4) * nxyz(i, 2));
        SmBL[2] = V * (stress(5) * nxyz(i, 0) + stress(4) * nxyz(i, 1) + stress(2) * nxyz(i, 2));
        for (int j = 0; j < 4; ++j)
        {
          double bopstrbop = 0.0;
          for (int dim = 0; dim < 3; ++dim) bopstrbop += nxyz(j, dim) * SmBL[dim];
          (*stiff)(lmlm[ele][i * 3 + 0], lmlm[ele][j * 3 + 0]) += bopstrbop;
          (*stiff)(lmlm[ele][i * 3 + 1], lmlm[ele][j * 3 + 1]) += bopstrbop;
          (*stiff)(lmlm[ele][i * 3 + 2], lmlm[ele][j * 3 + 2]) += bopstrbop;
        }  // for (int j=0; j<4; ++j)
      }    // for (int i=0; i<4; ++i)
    }      // for (int ele=0; ele<neleinpatch; ++ele)
  }        // if (stiff)

  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 10/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::MISNodalIntegration(Epetra_SerialDenseMatrix* stiff,
    Epetra_SerialDenseVector* force, std::map<int, DRT::Node*>& adjnode,
    std::vector<DRT::ELEMENTS::NStet*>& adjele, std::vector<double>& weight, std::vector<int>& lm,
    const Epetra_Vector& disp, DRT::Discretization& dis, std::vector<double>* nodalstress,
    std::vector<double>* nodalstrain, const INPAR::STR::StressType iostress,
    const INPAR::STR::StrainType iostrain)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::MISNodalIntegration");
  typedef Sacado::Fad::DFad<double> FAD;  // for first derivs

  //-------------------------------------------------- standard quantities
  const int nnodeinpatch = (int)adjnode.size();
  const int ndofinpatch = nnodeinpatch * 3;
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
  std::vector<FAD> patchdisp(ndofinpatch);
  for (int i = 0; i < ndofinpatch; ++i)
  {
    int lid = disp.Map().LID(lm[i]);
    if (lid == -1) dserror("Cannot find degree of freedom on this proc");
    patchdisp[i] = disp[disp.Map().LID(lm[i])];
    patchdisp[i].diff(i, ndofinpatch);
  }

  //-----------------------------------------------------------------------
  // build averaged F, det(F) and volume of node (using sacado)
  double VnodeL = 0.0;
  std::vector<std::vector<int>> lmlm(neleinpatch);
  FAD fad_Jnode = 0.0;
  for (int i = 0; i < neleinpatch; ++i)
  {
    const double V = weight[i] * adjele[i]->Vol();
    VnodeL += V;

    // get the element's displacements out of the patch' displacements
    std::vector<int> elelm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    adjele[i]->LocationVector(dis, elelm, lmowner, lmstride);
    // have to find position of elelm[i] in lm
    // lmlm[i][j] : element i, degree of freedom j, lmlm[i][j] position in patchdisp[0..ndofinpatch]
    lmlm[i].resize(12);
    for (int j = 0; j < 12; ++j)
    {
      std::vector<int>::iterator k = find(lm.begin(), lm.end(), elelm[j]);
      lmlm[i][j] = k - lm.begin();  // the position of elelm[j] in lm
    }

    // copy eledisp to 4x3 format
    LINALG::TMatrix<FAD, 4, 3> eledispmat(false);
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 3; ++k) eledispmat(j, k) = patchdisp[lmlm[i][j * 3 + k]];

    // build F and det(F) of this element
    LINALG::TMatrix<FAD, 3, 3> Fele = adjele[i]->BuildF<FAD>(eledispmat, adjele[i]->nxyz_);
    FAD Jele = LINALG::Determinant3x3<FAD>(Fele);

    fad_Jnode += V * Jele;

  }  // for (int i=0; i<neleinpatch; ++i)

  // do averaging
  fad_Jnode /= VnodeL;

  FAD Jpowthird = std::pow(fad_Jnode, 1. / 3.);

  // build volumetric deformation gradient
  LINALG::TMatrix<FAD, 3, 3> fad_FnodeL(true);
  for (int i = 0; i < 3; ++i) fad_FnodeL(i, i) = Jpowthird;

  // copy to LINALG::Matrix objects for output of strain
  const double Jnode = fad_Jnode.val();
  LINALG::Matrix<3, 3> FnodeL(false);
  for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k) FnodeL(j, k) = fad_FnodeL(j, k).val();

  //-----------------------------------------------------------------------
  // build \delta B operator, this is the unmodified operator
  Epetra_SerialDenseMatrix bop(6, ndofinpatch);
  for (int ele = 0; ele < neleinpatch; ++ele)
  {
    // current element
    DRT::ELEMENTS::NStet* actele = adjele[ele];

    // volume of that element assigned to node L
    const double V = weight[ele] * actele->Vol();

    // volume ratio of volume per node of this element to
    // whole volume of node L
    const double ratio = V / VnodeL;

    // get derivatives with respect to X
    LINALG::Matrix<4, 3>& nxyz = actele->Nxyz();

    // get defgrd
    LINALG::Matrix<3, 3>& F = actele->F();

    LINALG::Matrix<6, 12> bele(false);
    for (int i = 0; i < 4; ++i)
    {
      bele(0, 3 * i + 0) = F(0, 0) * nxyz(i, 0);
      bele(0, 3 * i + 1) = F(1, 0) * nxyz(i, 0);
      bele(0, 3 * i + 2) = F(2, 0) * nxyz(i, 0);
      bele(1, 3 * i + 0) = F(0, 1) * nxyz(i, 1);
      bele(1, 3 * i + 1) = F(1, 1) * nxyz(i, 1);
      bele(1, 3 * i + 2) = F(2, 1) * nxyz(i, 1);
      bele(2, 3 * i + 0) = F(0, 2) * nxyz(i, 2);
      bele(2, 3 * i + 1) = F(1, 2) * nxyz(i, 2);
      bele(2, 3 * i + 2) = F(2, 2) * nxyz(i, 2);

      bele(3, 3 * i + 0) = F(0, 0) * nxyz(i, 1) + F(0, 1) * nxyz(i, 0);
      bele(3, 3 * i + 1) = F(1, 0) * nxyz(i, 1) + F(1, 1) * nxyz(i, 0);
      bele(3, 3 * i + 2) = F(2, 0) * nxyz(i, 1) + F(2, 1) * nxyz(i, 0);
      bele(4, 3 * i + 0) = F(0, 1) * nxyz(i, 2) + F(0, 2) * nxyz(i, 1);
      bele(4, 3 * i + 1) = F(1, 1) * nxyz(i, 2) + F(1, 2) * nxyz(i, 1);
      bele(4, 3 * i + 2) = F(2, 1) * nxyz(i, 2) + F(2, 2) * nxyz(i, 1);
      bele(5, 3 * i + 0) = F(0, 2) * nxyz(i, 0) + F(0, 0) * nxyz(i, 2);
      bele(5, 3 * i + 1) = F(1, 2) * nxyz(i, 0) + F(1, 0) * nxyz(i, 2);
      bele(5, 3 * i + 2) = F(2, 2) * nxyz(i, 0) + F(2, 0) * nxyz(i, 2);
    }

    for (int k = 0; k < 6; ++k)
      for (int j = 0; j < 12; ++j) bop(k, lmlm[ele][j]) += ratio * bele(k, j);

  }  // for (int ele=0; ele<neleinpatch; ++ele)

  //-----------------------------------------------------------------------
  // green-lagrange strains based on averaged volumetric F
  LINALG::TMatrix<FAD, 3, 3> CG(false);
  CG.MultiplyTN(fad_FnodeL, fad_FnodeL);
  std::vector<FAD> Ebar(6);
  Ebar[0] = 0.5 * (CG(0, 0) - 1.0);
  Ebar[1] = 0.5 * (CG(1, 1) - 1.0);
  Ebar[2] = 0.5 * (CG(2, 2) - 1.0);
  Ebar[3] = CG(0, 1);
  Ebar[4] = CG(1, 2);
  Ebar[5] = CG(2, 0);

  // for material law and output, copy to baci object
  LINALG::Matrix<3, 3> cauchygreen(false);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) cauchygreen(i, j) = CG(i, j).val();

  LINALG::Matrix<6, 1> glstrain(false);
  for (int i = 0; i < 6; ++i) glstrain(i) = Ebar[i].val();

  //-------------------------------------------------------- output of strain
  if (iostrain != INPAR::STR::strain_none)
    StrainOutput(iostrain, *nodalstrain, FnodeL, Jnode, BETA_NSTET, 0.0);

  //-------------------------------------------------------------------------
  // build a second B-operator from the volumetric averaged strains that are based on
  // the averaged F
  Epetra_SerialDenseMatrix bopbar(6, ndofinpatch);
  for (int i = 0; i < ndofinpatch; ++i)
    for (int k = 0; k < 6; ++k) bopbar(k, i) = Ebar[k].fastAccessDx(i);

#if 0
  Epetra_SerialDenseMatrix bopbar2(6,ndofinpatch);
  {
    double Jbar = FnodeL.Determinant();

    double fac = (1./3.)*pow(Jbar,-2./3.)/VnodeL;

    //LINALG::Matrix<3,3> Kronecker(true);
    //Kronecker(0,0) = Kronecker(1,1) = Kronecker(2,2) = 1.0;

    //LINALG::Matrix<3,3> FKron(false);
    //FKron.MultiplyTN(FnodeL,Kronecker);

    Epetra_SerialDenseMatrix sumFinvdFdm(1,ndofinpatch);
    for (int ele=0; ele<neleinpatch; ++ele)
    {
      DRT::ELEMENTS::NStet* actele = adjele[ele];
      LINALG::Matrix<4,3>& nxyz = actele->Nxyz();
      LINALG::Matrix<3,3>& F = actele->F();
      const double Ve = actele->Vol();
      const double face = Ve * F.Determinant();
      LINALG::Matrix<3,3>  Finv = F;
      Finv.Invert();
      LINALG::Matrix<9,1> FinvT_vgt(false);
      FinvT_vgt(0) = Finv(0,0);
      FinvT_vgt(1) = Finv(1,1);
      FinvT_vgt(2) = Finv(2,2);
      FinvT_vgt(3) = Finv(1,0);
      FinvT_vgt(4) = Finv(2,0);
      FinvT_vgt(5) = Finv(0,1);
      FinvT_vgt(6) = Finv(2,1);
      FinvT_vgt(7) = Finv(0,2);
      FinvT_vgt(8) = Finv(1,2);
      LINALG::Matrix<9,12> Gradde(true);
      for (int i=0; i<4; ++i)
      {
        Gradde(0,3*i+0) = nxyz(i,0);
        //Gradde(0,3*i+1) = 0.0;
        //Gradde(0,3*i+2) = 0.0;

        //Gradde(1,3*i+0) = 0.0;
        Gradde(1,3*i+1) = nxyz(i,1);
        //Gradde(1,3*i+2) = 0.0;

        //Gradde(2,3*i+0) = 0.0;
        //Gradde(2,3*i+1) = 0.0;
        Gradde(2,3*i+2) = nxyz(i,2);

        Gradde(3,3*i+0) = nxyz(i,1);
        //Gradde(3,3*i+1) = 0.0;
        //Gradde(3,3*i+2) = 0.0;
        Gradde(4,3*i+0) = nxyz(i,2);
        //Gradde(4,3*i+1) = 0.0;
        //Gradde(4,3*i+2) = 0.0;

        //Gradde(5,3*i+0) = 0.0;
        Gradde(5,3*i+1) = nxyz(i,0);
        //Gradde(5,3*i+2) = 0.0;
        //Gradde(6,3*i+0) = 0.0;
        Gradde(6,3*i+1) = nxyz(i,2);
        //Gradde(6,3*i+2) = 0.0;

        //Gradde(7,3*i+0) = 0.0;
        //Gradde(7,3*i+1) = 0.0;
        Gradde(7,3*i+2) = nxyz(i,0);
        //Gradde(8,3*i+0) = 0.0;
        //Gradde(8,3*i+1) = 0.0;
        Gradde(8,3*i+2) = nxyz(i,1);
      }
      LINALG::Matrix<1,12> FinvdFdm(true);
      for (int i=0; i<9; ++i)
        for (int j=0; j<12; ++j)
          FinvdFdm(0,j) += FinvT_vgt(i) * Gradde(i,j);
      FinvdFdm.Scale(face);

      for (int j=0; j<12; ++j)
        sumFinvdFdm(0,lmlm[ele][j]) += FinvdFdm(0,j);
    }
    //sumFinvdFdm.Scale(fac);
    //std::cout << sumFinvdFdm << std::endl;
    //std::cout << FnodeL << FKron << std::endl;

    for (int i=0; i<ndofinpatch; ++i)
    {
      bopbar2(0,i) = fac * FnodeL(0,0) * sumFinvdFdm(0,i);
      bopbar2(1,i) = fac * FnodeL(1,1) * sumFinvdFdm(0,i);
      bopbar2(2,i) = fac * FnodeL(2,2) * sumFinvdFdm(0,i);
      //bopbar2(3,i) = FKron(0,1) * sumFinvdFdm(0,i);
      //bopbar2(4,i) = FKron(1,2) * sumFinvdFdm(0,i);
      //bopbar2(5,i) = FKron(2,0) * sumFinvdFdm(0,i);
    }


    //std::cout << bopbar2 << std::endl << std::endl << std::endl << std::endl;
    bopbar = bopbar2;
  }
#endif

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
    SelectMaterial(mat, stress, cmat, density, glstrain, FnodeL, 0, -1);
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
      DRT::ELEMENTS::NStet* actele = adjele[ele];
      // volume of that element assigned to node L
      const double V = weight[ele] * actele->Vol();
      // def-gradient of the element
      Teuchos::RCP<MAT::Material> mat = actele->Material();
      // EleGID is set to -1 errorcheck is performed in
      // MAT::Evaluate. I.e if we have elementwise mat params you will catch an error
      SelectMaterial(mat, stressele, cmatele, density, glstrain, FnodeL, 0, -1);
      cmat.Update(V, cmatele, 1.0);
      stress.Update(V, stressele, 1.0);
    }  // for (int ele=0; ele<neleinpatch; ++ele)
    stress.Scale(1.0 / VnodeL);
    cmat.Scale(1.0 / VnodeL);
  }


  //-----------------------------------------------------------------------
  // stress is split as follows:
  // stress = beta * vol_misnode + (1-beta) * vol_node + (1-alpha) * dev_node + alpha * dev_ele
  {
    LINALG::Matrix<6, 1> stressdev(true);
    LINALG::Matrix<6, 6> cmatdev(true);
    LINALG::Matrix<6, 1> stressvol(false);
    LINALG::Matrix<6, 6> cmatvol(false);

    // compute deviatoric stress and tangent from total stress and tangent
    DevStressTangent(stressdev, cmatdev, cmat, stress, cauchygreen);

    // compute volumetric stress and tangent
    stressvol.Update(-1.0, stressdev, 1.0, stress, 0.0);
    cmatvol.Update(-1.0, cmatdev, 1.0, cmat, 0.0);

    // compute nodal stress
    stress.Update(BETA_NSTET, stressvol, 0.0);
    cmat.Update(BETA_NSTET, cmatvol, 0.0);
  }

  //-----------------------------------------------------------------------
  // stress output
  if (iostress != INPAR::STR::stress_none)
    StressOutput(iostress, *nodalstress, stress, FnodeL, Jnode);

  //----------------------------------------------------- internal forces
  if (force)
  {
    Epetra_SerialDenseVector stress_epetra(::View, stress.A(), stress.Rows());
    force->Multiply('T', 'N', VnodeL, bop, stress_epetra, 0.0);
  }

  //--------------------------------------------------- elastic stiffness
  if (stiff)
  {
    Epetra_SerialDenseMatrix cmat_epetra(
        ::View, cmat.A(), cmat.Rows(), cmat.Rows(), cmat.Columns());
    LINALG::SerialDenseMatrix cb(6, ndofinpatch);
    cb.Multiply('N', 'N', 1.0, cmat_epetra, bopbar, 0.0);
    stiff->Multiply('T', 'N', VnodeL, bop, cb, 0.0);
  }

  //----------------------------------------------------- geom. stiffness
  // do not use sacado for second derivative of E as it is way too expensive!
  // As long as the 2nd deriv is as easy as this, do it by hand
  if (stiff)
  {
    // loop elements in patch
    for (int ele = 0; ele < neleinpatch; ++ele)
    {
      // material deriv of that element
      LINALG::Matrix<4, 3>& nxyz = adjele[ele]->Nxyz();
      // volume of actele assigned to node L
      double V = weight[ele] * adjele[ele]->Vol();

      // loop nodes of that element
      double SmBL[3];
      for (int i = 0; i < 4; ++i)
      {
        SmBL[0] = V * (stress(0) * nxyz(i, 0) + stress(3) * nxyz(i, 1) + stress(5) * nxyz(i, 2));
        SmBL[1] = V * (stress(3) * nxyz(i, 0) + stress(1) * nxyz(i, 1) + stress(4) * nxyz(i, 2));
        SmBL[2] = V * (stress(5) * nxyz(i, 0) + stress(4) * nxyz(i, 1) + stress(2) * nxyz(i, 2));
        for (int j = 0; j < 4; ++j)
        {
          double bopstrbop = 0.0;
          for (int dim = 0; dim < 3; ++dim) bopstrbop += nxyz(j, dim) * SmBL[dim];
          (*stiff)(lmlm[ele][i * 3 + 0], lmlm[ele][j * 3 + 0]) += bopstrbop;
          (*stiff)(lmlm[ele][i * 3 + 1], lmlm[ele][j * 3 + 1]) += bopstrbop;
          (*stiff)(lmlm[ele][i * 3 + 2], lmlm[ele][j * 3 + 2]) += bopstrbop;
        }  // for (int j=0; j<4; ++j)
      }    // for (int i=0; i<4; ++i)
    }      // for (int ele=0; ele<neleinpatch; ++ele)
  }        // if (stiff)



  return;
}



/*----------------------------------------------------------------------*
 | material laws for NStet (protected)                          gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::SelectMaterial(Teuchos::RCP<MAT::Material> mat,
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
      dserror("Illegal type %d of material for element NStet tet4", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::NStet::SelectMaterial

/*----------------------------------------------------------------------*
 |  compute deviatoric tangent and stresses (private/static)   gee 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::DevStressTangent(LINALG::Matrix<6, 1>& Sdev,
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
 |                                                             gee 10/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::StrainOutput(const INPAR::STR::StrainType iostrain,
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
 |                                                             gee 10/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::StrainOutput(const INPAR::STR::StrainType iostrain,
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
 |                                                             gee 10/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::StressOutput(const INPAR::STR::StressType iostress,
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
