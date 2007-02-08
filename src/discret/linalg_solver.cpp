/*!----------------------------------------------------------------------
\file linalg_solver.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "linalg_solver.H"

extern "C" 
{

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
}
#include "dstrc.H"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(const LINALG::Solver& old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::~Solver()
{
  return;
}

/*----------------------------------------------------------------------*
 |  translate solver parameters (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::TranslateSolverParameters(ParameterList& params, 
                                          struct _SOLVAR* actsolv)
{
  // switch type of solver
  switch (actsolv->solvertyp)
  {
#ifdef PARALLEL
  case superlu://============================== superlu solver (parallel only)
    params.set("solver","superlu");
    params.set("symmetric",false);
  break;
#endif
  case amesos_klu_sym://====================================== Tim Davis' KLU
    params.set("solver","klu");
    params.set("symmetric",true);
  break;
  case amesos_klu_nonsym://=================================== Tim Davis' KLU
    params.set("solver","klu");
    params.set("symmetric",false);
  break;
  case umfpack://========================================= Tim Davis' Umfpack
    params.set("solver","umfpack");
    params.set("symmetric",false);
  break;
  case lapack_sym://================================================== Lapack
    params.set("solver","lapack");
    params.set("symmetric",true);
  break;
  case lapack_nonsym://=============================================== Lapack
    params.set("solver","lapack");
    params.set("symmetric",false);
  break;
  case aztec_msr://================================================= AztecOO
  {
    params.set("solver","aztec");
    params.set("symmetric",false);
    AZVAR* azvar = actsolv->azvar;
    ParameterList& azlist = params.sublist("Aztec Parameters");
    //--------------------------------- set scaling of linear problem
    if (azvar->azscal==1)
      azlist.set("scaling","symmetric");
    else if (azvar->azscal==2)
      azlist.set("scaling","infnorm");
    //--------------------------------------------- set type of solver
    switch (azvar->azsolvertyp)
    {
    case azsolv_CG:       azlist.set("AZ_solver",AZ_cg);       break;
    case azsolv_GMRES:    azlist.set("AZ_solver",AZ_gmres);    break;
    case azsolv_CGS:      azlist.set("AZ_solver",AZ_cgs);      break;
    case azsolv_BiCGSTAB: azlist.set("AZ_solver",AZ_bicgstab); break;
    case azsolv_LU:       azlist.set("AZ_solver",AZ_lu);       break;
    case azsolv_TFQMR:    azlist.set("AZ_solver",AZ_tfqmr);    break;
    default: dserror("Unknown solver for AztecOO");            break;
    }
    //------------------------------------- set type of preconditioner
    switch (azvar->azprectyp)
    {
    case azprec_none:
      azlist.set("AZ_precond",AZ_none);
    break;
    case azprec_ILUT:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    case azprec_ILU:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    case azprec_Jacobi:
      azlist.set("AZ_precond",AZ_Jacobi);
    break;
    case azprec_Neumann:
      azlist.set("AZ_precond",AZ_Neumann);
    break;
    case azprec_Least_Squares:
      azlist.set("AZ_precond",AZ_ls);
    break;
    case azprec_SymmGaussSeidel:
      azlist.set("AZ_precond",AZ_sym_GS);
    break;
    case azprec_LU:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    case azprec_RILU:
      azlist.set("AZ_precond",AZ_dom_decomp);
      azlist.set("AZ_subdomain_solve",AZ_rilu);
      azlist.set("AZ_graph_fill",azvar->azgfill);
    break;
    case azprec_ICC:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    case azprec_ML:
    case azprec_MLfluid:
    case azprec_MLfluid2:
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    default:
      dserror("Unknown preconditioner for AztecOO");
    break;
    }
    //------------------------------------- set other aztec parameters
    azlist.set("AZ_kspace",azvar->azsub);
    azlist.set("AZ_max_iter",azvar->aziter);
    azlist.set("AZ_overlap",0);
    azlist.set("AZ_type_overlap",AZ_symmetric);
    azlist.set("AZ_poly_ord",azvar->azpoly);
    if (!azvar->azoutput)
      azlist.set("AZ_output",AZ_none);             // AZ_none AZ_all AZ_warnings AZ_last 10
    else
      azlist.set("AZ_output",azvar->azoutput);
    azlist.set("AZ_diagnostics",AZ_none);          // AZ_none AZ_all
    azlist.set("AZ_conv",AZ_noscaled);
    azlist.set("AZ_tol",azvar->aztol);
    azlist.set("AZ_drop",azvar->azdrop);
    azlist.set("AZ_scaling",AZ_none);              
    azlist.set("AZ_keep_info",1);
    // set reuse parameters
    azlist.set("ncall",0);                         // counting number of solver calls
    azlist.set("reuse",azvar->azreuse);            // reuse info for n solver calls
    //-------------------------------- set parameters for Ifpack if used
    if (azvar->azprectyp == azprec_ILU  ||
        azvar->azprectyp == azprec_ILUT ||
        azvar->azprectyp == azprec_ICC  ||
        azvar->azprectyp == azprec_LU   )
    {
      ParameterList& ifpacklist = params.sublist("IFPACK Parameters");
      ifpacklist.set("fact: drop tolerance",azvar->azdrop);
      ifpacklist.set("fact: level-of-fill",azvar->azgfill);
      ifpacklist.set("fact: ilut level-of-fill",azvar->azfill);
      ifpacklist.set("schwarz: combine mode","Add"); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
      ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
      ifpacklist.set("amesos: solver type", "Amesos_Klu"); // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
    }
    //------------------------------------- set parameters for ML if used
    if (azvar->azprectyp == azprec_ML      ||
        azvar->azprectyp == azprec_MLfluid ||
        azvar->azprectyp == azprec_MLfluid2 )
    {
      ParameterList& mllist = params.sublist("ML Parameters");
      ML_Epetra::SetDefaults("SA",mllist);
      switch (azvar->azprectyp)
      {
      case azprec_ML: // do nothing, this is standard
      break;
      case azprec_MLfluid: // unsymmetric, unsmoothed restruction
        mllist.set("aggregation: use tentative restruction",true);
      break;
      case azprec_MLfluid2: // full Pretrov-Galerkin unsymmetric smoothed
        mllist.set("energy minimization: enable",true);
        mllist.set("energy minimization: type",2); // 1,2,3 cheap -> expensive
        mllist.set("aggregation: block scaling",false); 
      break;
      default: dserror("Unknown type of ml preconditioner");
      }
      mllist.set("output"                          ,azvar->mlprint);
      if (azvar->mlprint==10)
        mllist.set("print unused"                  ,1);
      else
        mllist.set("print unused"                  ,-2);
      mllist.set("increasing or decreasing"        ,"increasing");
      mllist.set("coarse: max size"                ,azvar->mlcsize);
      mllist.set("max levels"                      ,azvar->mlmaxlevel);
      mllist.set("smoother: pre or post"           ,"both");
      mllist.set("aggregation: threshold"          ,azvar->ml_threshold);
      mllist.set("aggregation: damping factor"     ,azvar->mldamp_prolong);
      mllist.set("aggregation: nodes per aggregate",azvar->mlaggsize);
      switch (azvar->mlcoarsentype)
      {
        case 0:  mllist.set("aggregation: type","Uncoupled");  break;
        case 1:  mllist.set("aggregation: type","METIS");      break;
        case 2:  mllist.set("aggregation: type","VBMETIS");    break;
        case 3:  mllist.set("aggregation: type","MIS");        break;
        default: dserror("Unknown type of coarsening for ML"); break;
      }
      // set ml smoothers
      for (int i=0; i<azvar->mlmaxlevel-1; ++i)
      {
        char levelstr[11];
        sprintf(levelstr,"(level %d)",i);
        int type;
        double damp;
        if (i==0)
        {
          type = azvar->mlsmotype_fine;
          damp = azvar->mldamp_fine;
        }
        else if (i < azvar->mlmaxlevel-1)
        {
          type = azvar->mlsmotype_med;
          damp = azvar->mldamp_med;
        }
        else
        {
          type = azvar->mlsmotype_coarse;
          damp = azvar->mldamp_coarse;
        }
        switch (type)
        {
        case 0:
          mllist.set("smoother: type "+(string)levelstr                    ,"symmetric Gauss-Seidel");
          mllist.set("smoother: sweeps "+(string)levelstr                  ,azvar->mlsmotimes[i]);
          mllist.set("smoother: damping factor "+(string)levelstr          ,damp);
        break;
        case 1:
          mllist.set("smoother: type "+(string)levelstr                    ,"Jacobi");
          mllist.set("smoother: sweeps "+(string)levelstr                  ,azvar->mlsmotimes[i]);
          mllist.set("smoother: damping factor "+(string)levelstr          ,damp);
        break;
        case 2:
          mllist.set("smoother: type "+(string)levelstr                    ,"MLS");
          mllist.set("smoother: MLS polynomial order "+(string)levelstr    ,azvar->mlsmotimes[i]);
        break;
        case 3:
          mllist.set("smoother: type (level 0)"                            ,"MLS");
          mllist.set("smoother: MLS polynomial order "+(string)levelstr    ,-azvar->mlsmotimes[i]);
        break;
        case 4:
          mllist.set("smoother: type "+(string)levelstr                    ,"IFPACK");
          mllist.set("smoother: ifpack type "+(string)levelstr             ,"ILU");
          mllist.set("smoother: ifpack overlap "+(string)levelstr          ,0);
          mllist.sublist("smoother: ifpack list").set("fact: level-of-fill",azvar->mlsmotimes[i]);
          mllist.sublist("smoother: ifpack list").set("schwarz: reordering type","rcm");
        break;
        case 5:
          mllist.set("smoother: type "+(string)levelstr,"Amesos-KLU");
        break;
#ifdef PARALLEL
        case 6:
          mllist.set("smoother: type "+(string)levelstr,"Amesos-Superludist");
        break;
#endif
        default: dserror("Unknown type of smoother for ML"); break;
        } // switch (type)
      } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)
      // set coarse grid solver
      const int coarse = azvar->mlmaxlevel-1;
      switch (azvar->mlsmotype_coarse)
      {
        case 0:
          mllist.set("coarse: type"          ,"symmetric Gauss-Seidel");
          mllist.set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
          mllist.set("coarse: damping factor",azvar->mldamp_coarse);
        break;
        case 1:
          mllist.set("coarse: type"          ,"Jacobi");
          mllist.set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
          mllist.set("coarse: damping factor",azvar->mldamp_coarse);
        break;
        case 2:
          mllist.set("coarse: type"                ,"MLS");
          mllist.set("coarse: MLS polynomial order",azvar->mlsmotimes[coarse]);
        break;
        case 3:
          mllist.set("coarse: type"                ,"MLS");
          mllist.set("coarse: MLS polynomial order",-azvar->mlsmotimes[coarse]);
        break;
        case 4:
          mllist.set("coarse: type"          ,"IFPACK");
          mllist.set("coarse: ifpack type"   ,"ILU");
          mllist.set("coarse: ifpack overlap",0);
          mllist.sublist("coarse: ifpack list").set("fact: level-of-fill",azvar->mlsmotimes[coarse]);
          mllist.sublist("coarse: ifpack list").set("schwarz: reordering type","rcm");
        break;
        case 5:
          mllist.set("coarse: type","Amesos-KLU");
        break;
        case 6:
          mllist.set("coarse: type","Amesos-Superludist");
        break;
        default: dserror("Unknown type of coarse solver for ML"); break;
      } // switch (azvar->mlsmotype_coarse)
      // default values for nullspace
      mllist.set("PDE equations",1);
      mllist.set("null space: dimension",1);
      mllist.set("null space: type","pre-computed");
      mllist.set("null space: add default vectors",false);
      mllist.set<double*>("null space: vectors",NULL);
    } // if ml preconditioner
  }	
  break;
#ifdef PARALLEL
  case SPOOLES_sym://================================== Spooles (parallel only)
  case SPOOLES_nonsym:
    params.set("solver","spooles");
    params.set("symmetric",false);
  break;
#endif
  default:
    dserror("Unsupported type of solver");
  break;
  }
  return;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
