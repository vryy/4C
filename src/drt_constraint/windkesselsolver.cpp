/*!----------------------------------------------------------------------
\file windkesselsolver.cpp

\brief Class containing direct solver to solve linear system.

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/


#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <stdio.h>
#include <iostream>

#include "windkesselsolver.H"
#include "windkessel_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
UTILS::WindkesselSolver::WindkesselSolver
(
  RCP<DRT::Discretization> discr,
  LINALG::Solver& solver,
  RCP<LINALG::MapExtractor> dbcmaps,
  ParameterList params
):
actdisc_(discr),
maxIter_(params.get<int>   ("UZAWAMAXITER", 50)),
dirichtoggle_(Teuchos::null),
dbcmaps_(dbcmaps)
{
  Setup(discr,solver,dbcmaps,params);
}

/*----------------------------------------------------------------------*
 |  set-up (public)                                             tk 11/07|
 *----------------------------------------------------------------------*/
void UTILS::WindkesselSolver::Setup
(
  RCP<DRT::Discretization> discr,
  LINALG::Solver& solver,
  RCP<LINALG::MapExtractor> dbcmaps,
  ParameterList params
)
{

  solver_ = Teuchos::rcp(&solver,false);

  algochoice_ = DRT::INPUT::IntegralValue<INPAR::STR::ConSolveAlgo>(params,"UZAWAALGO");

  // different setup for #adapttol_
  isadapttol_ = true;
  isadapttol_ = (DRT::INPUT::IntegralValue<int>(params,"ADAPTCONV") == 1);
  
  // simple parameters
  adaptolbetter_ = params.get<double>("ADAPTCONV_BETTER", 0.01);
  iterationparam_ = params.get<double>("UZAWAPARAM", 1);
  minparam_ = iterationparam_*1E-3;
  iterationtol_ = params.get<double>("UZAWATOL", 1E-8);


  counter_ = 0;
  return;
}



/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear system with Windkessel boundary condition                                        |
*-----------------------------------------------------------------------*/
void UTILS::WindkesselSolver::Solve
(
  RCP<LINALG::SparseMatrix> stiff,
  RCP<LINALG::SparseMatrix> coupoffdiag_vol_d,
  RCP<LINALG::SparseMatrix> coupoffdiag_fext_p,
  RCP<LINALG::SparseMatrix> windkstiff,
  RCP<Epetra_Vector> dispinc,
  RCP<Epetra_Vector> presinc,
  const RCP<Epetra_Vector> rhsstand,
  const RCP<Epetra_Vector> rhswindk
)
{

  SolveDirect(stiff,coupoffdiag_vol_d,coupoffdiag_fext_p,windkstiff,dispinc,presinc,rhsstand,rhswindk);

  return;
}


void UTILS::WindkesselSolver::SolveDirect
(
  RCP<LINALG::SparseMatrix> stiff,
  RCP<LINALG::SparseMatrix> coupoffdiag_vol_d,
  RCP<LINALG::SparseMatrix> coupoffdiag_fext_p,
  RCP<LINALG::SparseMatrix> windkstiff,
  RCP<Epetra_Vector> dispinc,
  RCP<Epetra_Vector> presinc,
  const RCP<Epetra_Vector> rhsstand,
  const RCP<Epetra_Vector> rhswindk
)
{
  // define maps of standard dofs and additional pressures
  RCP<Epetra_Map> standrowmap = Teuchos::rcp(new Epetra_Map(stiff->RowMap()));
  RCP<Epetra_Map> windkrowmap = Teuchos::rcp(new Epetra_Map(windkstiff->RowMap()));
  //RCP<Epetra_Map> windkrowmap = Teuchos::rcp(new Epetra_Map(coupoffdiag_fext_p->DomainMap()));

  // merge maps to one large map
  RCP<Epetra_Map> mergedmap = LINALG::MergeMap(standrowmap,windkrowmap,false);
  // define MapExtractor
  LINALG::MapExtractor mapext(*mergedmap,standrowmap,windkrowmap);

  // initialize large Sparse Matrix and Epetra_Vectors
  RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap,81));
  RCP<Epetra_Vector> mergedrhs = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  RCP<Epetra_Vector> mergedsol = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = LINALG::ConvertDirichletToggleVectorToMaps(dirichtoggle_);
  // fill merged matrix using Add
  mergedmatrix -> Add(*stiff,false,1.0,1.0);
  mergedmatrix -> Add(*coupoffdiag_vol_d,true,1.0,1.0);
  mergedmatrix -> Add(*coupoffdiag_fext_p,false,1.0,1.0);
  mergedmatrix -> Add(*windkstiff,false,1.0,1.0);
  mergedmatrix -> Complete(*mergedmap,*mergedmap);

  //std::cout << "" << *coupoffdiag_vol_d << std::endl;
  //std::cout << "" << *coupoffdiag_fext_p << std::endl;
  //std::cout << "" << *windkstiff << std::endl;
  //std::cout << "" << *rhswindk << std::endl;

  // fill merged vectors using Export
  LINALG::Export(*rhswindk,*mergedrhs);
  mergedrhs -> Scale(-1.0);
  LINALG::Export(*rhsstand,*mergedrhs);

#if 0
    const int myrank=(actdisc_->Comm().MyPID());
    const double cond_number = LINALG::Condest(static_cast<LINALG::SparseMatrix&>(*mergedmatrix),Ifpack_GMRES, 100);
    // computation of significant digits might be completely bogus, so don't take it serious
    const double tmp = std::abs(std::log10(cond_number*1.11022e-16));
    const int sign_digits = (int)floor(tmp);
    if (!myrank)
      std::cout << " cond est: " << std::scientific << cond_number << ", max.sign.digits: " << sign_digits<<std::endl;
#endif

  // solve
  solver_->Solve(mergedmatrix->EpetraMatrix(),mergedsol,mergedrhs,true,counter_==0);
  solver_->ResetTolerance();
  // store results in smaller vectors
  mapext.ExtractCondVector(mergedsol,dispinc);
  mapext.ExtractOtherVector(mergedsol,presinc);

  counter_++;
  return;
}


