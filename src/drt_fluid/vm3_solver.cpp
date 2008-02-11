/*!----------------------------------------------------------------------
\file vm3_solver.cpp
\brief A direct algebraic VM3 approach

Scales are separated via plain aggregation. A fine-scale 
subgrid-viscosity matrix in the form
                 M_fssv = S^T * M_sv * S
is created, using an all-scale subgrid-viscosity matrix M_sv and a
scale-separation matrix S.

To (substantially) reduce computational cost, the matrix M_fssv is
initially (i.e., in the first timestep) generated in a normalized form
(i.e., based on a virtual subgrid-viscosity "1"). In each further
timestep, a subgrid-viscosity-scaling vector is generated containing 
the "real" subgrid viscosity for each (velocity) degree of freedom, and
the normalized matrix M_fssv is left- and right-scaled by the square-root 
values of that vector.

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "vm3_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 02/08|
 *----------------------------------------------------------------------*/
VM3_Solver::VM3_Solver(RCP<Epetra_CrsMatrix> Msv,
                       const RCP<Epetra_Vector> dbctoggle,
                       ParameterList& mlparams,
                       bool compute) :
compute_(compute),
mlparams_(mlparams),
dbctoggle_(dbctoggle)
{
  if (compute_) Compute(Msv);
  return;
}

/*----------------------------------------------------------------------*
 |  routine for generating the scale-separation matrix S        vg 02/08|
 |  and precomputing the unscaled matrix S^T*M*S                        |
 |  (only called in the first timestep)                                 |
 *----------------------------------------------------------------------*/
bool VM3_Solver::Compute(RCP<Epetra_CrsMatrix> Msv)
{
  // this is important to have!!!
  MLAPI::Init();

  const Epetra_Vector& dbct = *dbctoggle_;

  // get nullspace parameters
  double* nullspace = mlparams_.get("null space: vectors",(double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int nsdim = mlparams_.get("null space: dimension",1);

  // modify nullspace to ensure that DBC are fully taken into account
  if (nullspace)
  {
    const int length = Msv->RowMap().NumMyElements();
    for (int i=0; i<nsdim; ++i)
      for (int j=0; j<length; ++j)
        if (dbct[j]!=0.0) nullspace[i*length+j] = 0.0;
  }

  // get plain aggregation Ptent
  RCP<Epetra_CrsMatrix> Ptent;
  GetPtent(*Msv,mlparams_,nullspace,Ptent);

  // compute scale-separation matrix: S = I - (Ptent*Ptent^T)
  RCP<Epetra_CrsMatrix> Sep = LINALG::Multiply(Ptent,false,Ptent,true);
  Sep->Scale(-1.0);
  RCP<Epetra_Vector> tmp = LINALG::CreateVector(Sep->RowMap(),false);
  tmp->PutScalar(1.0);
  RCP<Epetra_Vector> diag = LINALG::CreateVector(Sep->RowMap(),false);
  Sep->ExtractDiagonalCopy(*diag);
  diag->Update(1.0,*tmp,1.0);
  Sep->ReplaceDiagonalValues(*diag);

  //complete scale-separation matrix and check maps
  Sep->FillComplete(Sep->DomainMap(),Sep->RangeMap());
  if (!Sep->RowMap().SameAs(Msv->RowMap())) dserror("rowmap not equal");
  if (!Sep->RangeMap().SameAs(Msv->RangeMap())) dserror("rangemap not equal");
  if (!Sep->DomainMap().SameAs(Msv->DomainMap())) dserror("domainmap not equal");

  // precomputation of unscaled S^T*M*S:
  // pre- and post-multiply M by scale-separating operator matrix Sep
  Mnsv_ = LINALG::Multiply(Msv,false,Sep,false);
  Mnsv_ = LINALG::Multiply(Sep,true,Mnsv_,false);

  return false;
}

/*----------------------------------------------------------------------*
 |  scaling routine                                             vg 02/08|
 |  scale precomput. matrix product by subgrid-viscosity-scaling vector |
 |  (called in every timestep)                                          |
 *----------------------------------------------------------------------*/
void VM3_Solver::Scale(RCP<Epetra_CrsMatrix>& Msv,
                       RCP<Epetra_CrsMatrix>& K,
                       RCP<Epetra_Vector>& r,
                       RCP<Epetra_Vector>& rplus,
                       RCP<Epetra_Vector>& sugrvisc,
                       RCP<Epetra_Vector>& sol,
                       bool increm)
{
  // some necessary definitions
  int ierr;
  double* sgvsqrt = 0;
  int length = sugrvisc->MyLength();

  // square-root of subgrid-viscosity-scaling vector for left and right scaling
  sgvsqrt = (double*)sugrvisc->Values();
  for (int i = 0; i < length; ++i)
  {
    sgvsqrt[i] = sqrt(sgvsqrt[i]);
    sugrvisc->ReplaceMyValues(1,&sgvsqrt[i],&i);
  }

  // get unscaled S^T*M*S from Sep
  Msv = rcp(new Epetra_CrsMatrix(*Mnsv_));

  // left and right scaling of normalized fine-scale subgrid-viscosity matrix
  ierr = Msv->LeftScale(*sugrvisc);
  if (ierr) dserror("Epetra_CrsMatrix::LeftScale returned err=%d",ierr);
  ierr = Msv->RightScale(*sugrvisc);
  if (ierr) dserror("Epetra_CrsMatrix::RightScale returned err=%d",ierr);

  // compute the additional rhs and add it to existing rhs for incremental formul.
  if (increm)
  {
    // add the subgrid-viscosity-scaled fine-scale matrix to obtain complete matrix
    //LINALG::Add(Msv_,false,1.0,K,1.0);
    // multiply subgrid-viscosity-scaled fine-scale matrix by solution vector
    // (only this rhs part required for fixed-point iteration of fine-scale matrix)
    Msv->Multiply(false,*sol,*rplus);
    r->Update(-1.0,*rplus,1.0);
  }
  else
    // add the subgrid-viscosity-scaled fine-scale matrix to obtain complete matrix
    LINALG::Add(Msv,false,1.0,K,1.0);

  return;
}

#endif
