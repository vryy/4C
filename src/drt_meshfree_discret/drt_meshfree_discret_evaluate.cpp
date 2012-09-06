/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret_evaluate.H
 *
 * \brief evaluation of a meshfree discretisation
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "drt_meshfree_discret.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Export.h>
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensevector.H"

//#include "drt_meshfree_node.H"
//#include "drt_meshfree_cell.H"
//#include "../drt_lib/drt_utils.cpp"

/*--------------------------------------------------------------------------*
 |                                                   (public) nis Feb12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Evaluate(Teuchos::ParameterList&        params,
                                                     Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                                     Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                                     Teuchos::RCP<Epetra_Vector>    systemvector1,
                                                     Teuchos::RCP<Epetra_Vector>    systemvector2,
                                                     Teuchos::RCP<Epetra_Vector>    systemvector3)
{
  // hack ?
  systemvector2 = Teuchos::null;

  TEUCHOS_FUNC_TIME_MONITOR("DRT::MESHFREE::MeshfreeDiscretization::Evaluate");

  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  //-----------------------------------------------------------------
  // call the element's register class preevaluation method
  // for each type of element
  // for most element types, just the base class dummy is called
  // that does nothing
  {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::MESHFREE::MeshfreeDiscretization::Evaluate PreEvaluate");
    ParObjectFactory::Instance().PreEvaluate(*this,params,
                                             systemmatrix1,
                                             systemmatrix2,
                                             systemvector1,
                                             systemvector2,
                                             systemvector3);
  }

  //-----------------------------------------------------------------
  // see what we have for input
  bool assemblemat1 = (systemmatrix1!=Teuchos::null);
  bool assemblemat2 = (systemmatrix2!=Teuchos::null);
  bool assemblevec1 = (systemvector1!=Teuchos::null);
  bool assemblevec2 = (systemvector2!=Teuchos::null);
  bool assemblevec3 = (systemvector3!=Teuchos::null);
  if (assemblemat2 || (assemblevec2 || assemblevec3)) dserror("Wrong assembly expectations");

  //-----------------------------------------------------------------
  // create a temporary matrix to assemble to in a baci-unusual way
  // (across-parallel-interface assembly)
  const Epetra_Map* rmap = NULL;
  const Epetra_Map* dmap = NULL;

  RCP<Epetra_FECrsMatrix> systemmatrix_FE;
  RCP<LINALG::SparseMatrix> systemmatrix;
  if (systemmatrix1 != Teuchos::null)
  {
    rmap = &(systemmatrix1->OperatorRangeMap());
    dmap = rmap;
    systemmatrix = rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix1);
    if (systemmatrix != null && systemmatrix->Filled())
      systemmatrix_FE = rcp(new Epetra_FECrsMatrix(::Copy,systemmatrix->EpetraMatrix()->Graph()));
    else
      systemmatrix_FE = rcp(new Epetra_FECrsMatrix(::Copy,*rmap,256,false));
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
  Epetra_Vector systemvector(*this->DofColMap(),true);

  //-----------------------------------------------------------------
  // prepare location array of (single) dofset
  int row = 0;
  int col = 0;
  Element::LocationArray la(1);
  vector<int> lm;
  int rdim;
  int cdim;
  if (dofsets_.size()!=1){
    dserror("No multiple dofsets allowed for meshfree discretisations.");
  }

  //-----------------------------------------------------------------
  // prepare element matrices and vectors
  LINALG::SerialDenseMatrix elematrix1;
  LINALG::SerialDenseMatrix elematrix2;
  LINALG::SerialDenseVector elevector1;
  LINALG::SerialDenseVector elevector2;
  LINALG::SerialDenseVector elevector3;

  //-----------------------------------------------------------------
  // loop over row elements
  const int numrowele = NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = lRowElement(i);

    //-----------------------------------------------------------------
    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate LocationVector");

    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(*this,la,false);
    }

    //-----------------------------------------------------------------
    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate Resize");

    // get dimension of element matrices and vectors
    // reshape element matrices and vectors if necessary and set/init to zero
    rdim = la[row].Size(); // number of rows
    cdim = la[col].Size(); // number of cols (same as rows, since no multiple dofsets)
    lm = la[col].lm_; // relies on row==col

    if (assemblemat1)
    {
      if (elematrix1.M()!=rdim or elematrix1.N()!=cdim)
        elematrix1.Shape(rdim,cdim);
      else
        memset(elematrix1.A(),0,rdim*cdim*sizeof(double));
    }
    if (assemblevec1)
    {
      if (elevector1.Length()!=rdim)
        elevector1.Size(rdim);
      else
        memset(elevector1.Values(),0,rdim*sizeof(double));
    }
    }

    {//-----------------------------------------------------------------

    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate elements");

    // call the element evaluate method
    int err = actele->Evaluate(params,*this,la,
                               elematrix1,
                               elematrix2,
                               elevector1,
                               elevector2,
                               elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);
    }

    //-----------------------------------------------------------------
    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");

    if (assemblemat1)
    {
      vector<int> lrlm;
      vector<int> lclm;

      const Epetra_Map& dofrowmap = systemmatrix1->OperatorRangeMap();
      lrlm.resize(rdim);
      for (int i=0; i<rdim; ++i)
        lrlm[i] = dofrowmap.LID(lm[i]);
      const Epetra_Map& dofcolmap = systemmatrix->ColMap();
      lclm.resize(cdim);
      for (int i=0; i<cdim; ++i)
      lclm[i] = dofcolmap.LID(lm[i]);

      // ==============================================================
      // on-processor assembly
      for (int i=0; i<rdim; ++i)
      {
        if (lrlm[i]==-1) // off-processor row
        {
          for (int j=0; j<cdim; ++j)
          {
            int errone = systemmatrix_FE->SumIntoGlobalValues(1,&lrlm[i],1,&lclm[j],&elematrix1(i,j));
            if (errone>0)
            {
              int errtwo = systemmatrix_FE->InsertGlobalValues(1,&lrlm[i],1,&lclm[j],&elematrix1(i,j));
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
            for (int j=0; j<cdim; ++j)
            {
              int* loc = std::lower_bound(indices,indices+length,lclm[j]);
#ifdef DEBUG
              if (*loc != lclm[j]) dserror("Cannot find local column entry %d",lclm[j]);
#endif
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
                values[pos++] += elematrix1(i,j++);
                values[pos++] += elematrix1(i,j++);
                values[pos]   += elematrix1(i,j);
              }
              else
              {
                int err =  matrix.SumIntoMyValues(lrlm[i],1,&elematrix1(i,j),&lclm[j]); j++;
                err += matrix.SumIntoMyValues(lrlm[i],1,&elematrix1(i,j),&lclm[j]); j++;
                err += matrix.SumIntoMyValues(lrlm[i],1,&elematrix1(i,j),&lclm[j]);
                if (err) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d",err);
              }
            }
          }
          else // matrix not SparseMatrix (e.g. BlockMatrix) -> fall back to standard assembly
          {
            for (int j=0; j<cdim; ++j)
              systemmatrix1->Assemble(elematrix1(i,j),lrlm[i],lclm[j]);
          }
        }
      } // for (node=noderids_.begin(); node != noderids_.end(); ++node)
    } // if (assemblemat1)
    } //TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");

    //-----------------------------------------------------------------------------------
    if (assemblevec1)
    {
      for (int i=0; i<rdim; ++i)
      {
        const int rgid = lm[i];
        const int lid = systemvector.Map().LID(rgid);
        if (lid<0) dserror("global row %d does not exist in column map",rgid);
        systemvector[lid] += elevector1[i];
      }
    }
    // end on-processor assembly
    // ==============================================================

  }// for (int i=0; i<numrowele; ++i)

  // ================================================================
  // off-processor assembly

  // note that fillComplete is never called on systemmatrix_FE
  if (assemblemat1)
  {
    int err = systemmatrix_FE->GlobalAssemble(*dmap,*rmap,false);
    if (err) dserror("Epetra_FECrsMatrix::GlobalAssemble returned err=%d",err);
    const Epetra_Map& cmap = systemmatrix_FE->ColMap();
    for (int lrow=0; lrow<systemmatrix_FE->NumMyRows(); ++lrow)
    {
      int numentries;
      double* values;
      if (!systemmatrix_FE->Filled())
      {
        const int grow = systemmatrix_FE->RowMap().GID(lrow);
        int* gindices;
        int err = systemmatrix_FE->ExtractGlobalRowView(grow,numentries,values,gindices);
        if (err) dserror("Epetra_FECrsMatrix::ExtractGlobalRowView returned err=%d",err);
        for (int j=0; j<numentries; ++j)
          systemmatrix1->Assemble(values[j],grow,gindices[j]);
      }
      else
      {
        int* lindices;
        int err = systemmatrix_FE->ExtractMyRowView(lrow,numentries,values,lindices);
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
          const int grow = systemmatrix_FE->RowMap().GID(lrow);
          for (int j=0; j<numentries; ++j)
            systemmatrix1->Assemble(values[j],grow,cmap.GID(lindices[j]));
        }
      }
    }
  }

  //-----------------------------------------------------------------
  if (assemblevec1)
  {
    Epetra_Vector tmp(systemvector1->Map(),false);
    Epetra_Export exporter(systemvector.Map(),tmp.Map());
    int err = tmp.Export(systemvector,exporter,Add);
    if (err) dserror("Export using exporter returned err=%d",err);
    systemvector1->Update(1.0,tmp,1.0);
  }

  // end off-processor assembly
  // ================================================================

  return;
}
