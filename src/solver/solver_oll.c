/*!----------------------------------------------------------------------
\file
\brief contains functions to handle oll matrices 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"

/*! 
\addtogroup OLL 
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief controls the solver for oll matrices

<pre>                                                              mn 05/03
This function counts the number of equations on this processor
</pre>
\param *actsolv        FIELD  (i)   the active field
\param *actintra       INTRA  (i)   the active communicator
\param *oll            OLL    (i)   the matrix to solve
\param *sol            DIST_VECTOR (i)  the solution vector
\param *rhs            DIST_VECTOR (i)  the rhs vector
\param  option         INT    (i)   option for init or not

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void solver_oll(
    struct _SOLVAR         *actsolv,
    struct _INTRA          *actintra,
    struct _OLL            *oll,
    struct _DIST_VECTOR    *sol,
    struct _DIST_VECTOR    *rhs,
    INT                     option)
{

#ifdef DEBUG 
  dstrc_enter("solver_oll");
#endif
  /*----------------------------------------------------------------------*/
  switch(option)
  {
    /*----------------------------------------------------------------------*/
    /*                                                           init phase */
    /*----------------------------------------------------------------------*/
    case 1:
      /* ----------------------------allocate sysarray in solver format */
      switch(actsolv->solvertyp)
      {
        case colsol_solver:/*--------------------------------- solver is colsol */
#ifdef PARALLEL 
          /* ------------------------------------ NO COLSOL in parallel for OLL */
#ifdef SPOOLES_PACKAGE
          /* -------------------------------------------- using SPOOLES instead */
          printf("No COLSOL in parallel for OLL\nUsing SPOOLES instead!!\n");
          oll->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
          oll->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
          oll->sysarray_typ[0] = spoolmatrix;
          oll->sysarray[0].spo = (SPOOLMAT*)CCACALLOC(1,sizeof(SPOOLMAT));
          /* --------------------------------------------------- and initialize */
          oll->sysarray[0].spo->is_init    =1;
          oll->sysarray[0].spo->ncall      =0;
          oll->sysarray[0].spo->is_factored=0;
#else
          dserror("NO COLSOL in parallel for OLL");
#endif
#else
          /* -------------------allocate skyline matrix */
          oll->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
          oll->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
          oll->sysarray_typ[0] = skymatrix;
          oll->sysarray[0].sky = (SKYMATRIX*)CCACALLOC(1,sizeof(SKYMATRIX));
          /* --------------------------------------------------- and initialize */
          oll->sysarray[0].sky->is_init     = 1;
          oll->sysarray[0].sky->ncall       = 0;
          oll->sysarray[0].sky->is_factored = 0;
#endif
          break;
        case SPOOLES_nonsym:/*------------------------------- solver is spooles */
          oll->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
          oll->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
          oll->sysarray_typ[0] = spoolmatrix;
          oll->sysarray[0].spo = (SPOOLMAT*)CCACALLOC(1,sizeof(SPOOLMAT));
          /* --------------------------------------------------- and initialize */
          oll->sysarray[0].spo->is_init    =1;
          oll->sysarray[0].spo->ncall      =0;
          oll->sysarray[0].spo->is_factored=0;
          break;
        case aztec_msr:/*------------------------------------- solver is aztec */
          oll->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
          oll->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
          oll->sysarray_typ[0] = msr;
          oll->sysarray[0].msr = (AZ_ARRAY_MSR*)CCACALLOC(1,sizeof(AZ_ARRAY_MSR));
          amdef("bindx",&(oll->sysarray[0].msr->bindx),oll->nnz+1,1,"IV");
          /* --------------------------------------------------- and initialize */
          solver_az_oll(actsolv,actintra,oll->sysarray[0].msr,sol,rhs,option);
          /*    solver_control(
                actsolv,
                actintra,
                &(oll->sysarray_typ[0]),
                &(oll->sysarray[0]),
                sol,
                rhs,
                option
                );*/
          break;
        case umfpack: /*-------------------------------------- solver is umfpack */
#ifdef PARALLEL
          dserror("No UMFPACK for parallel OLL!\n");
#endif	  
          oll->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
          oll->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
          oll->sysarray_typ[0] = ccf;
	  oll->sysarray[0].ccf = (CCF*)CCACALLOC(1,sizeof(CCF));
	  oll->sysarray[0].ccf->numeq_total = oll->numeq_total;
	  oll->sysarray[0].ccf->numeq = oll->numeq;
	  oll->sysarray[0].ccf->nnz = oll->nnz;          
	  oll->sysarray[0].ccf->nnz_total = oll->nnz;          
          amdef("Ap",&(oll->sysarray[0].ccf->Ap),oll->numeq_total+1,1,"IV");
          amdef("Ai",&(oll->sysarray[0].ccf->Ai),oll->nnz     ,1,"IV");
          amdef("Ax",&(oll->sysarray[0].ccf->Ax),oll->nnz     ,1,"DV");
	  amdef("update",&(oll->sysarray[0].ccf->update),oll->numeq,1,"IV");
	  solver_umfpack(NULL,actintra,oll->sysarray[0].ccf,sol,rhs,option);
	break;
        default:
          dserror("Unknown solver typ for oll");
          break;   
      }

      /* set flag, that this matrix has been copied */   
      oll->is_copied = 0;
      break;
      /*----------------------------------------------------------------------*/
      /*                                                    end of init phase */
      /*----------------------------------------------------------------------*/
      /*----------------------------------------------------------------------*/
      /*                                                    calculation phase */
      /*----------------------------------------------------------------------*/
    case 0:
      if( oll->is_copied ==0)
      {
        /* erst die Matrix kopieren */
        switch(actsolv->solvertyp)
        {
          case colsol_solver:/*--------------------------------- solver is colsol */
#ifdef PARALLEL 
            /* ------------------------------------ NO COLSOL in parallel for OLL */
#ifdef SPOOLES_PACKAGE
            /* -------------------------------------------- using SPOOLES instead */
            oll_to_spo(oll, &(oll->sysarray[0]));
#else
            dserror("NO COLSOL in parallel for OLL");
#endif
#else
            oll_to_sky(oll, &(oll->sysarray[0]));
#endif
            break;
          case SPOOLES_nonsym:/*------------------------------- solver is spooles */
            oll_to_spo(oll, &(oll->sysarray[0]));
            break;
          case aztec_msr:/*-------------------------------------- solver is aztec */
            oll_to_msr(oll, &(oll->sysarray[0]));
            break;
          case umfpack: /*------------------------------------- solver is umfpack */
	    oll_to_ccf(oll,&(oll->sysarray[0]));
	  break; 
          default:
            dserror("Unknown solver typ for oll");
            break;   
        }
        oll->is_copied = 1;
      }

      /* print sparsity pattern to gnuplot */
      if (oll->sparsepat == 0)
      {
        oll_gnupattern(oll);
        oll->sparsepat = 1;
      }

      /* dann loesen */
      switch(oll->sysarray_typ[0])
      {
        case msr:/*-------------------------------- system matrix is msr matrix */
          solver_az_oll(actsolv,actintra,oll->sysarray[0].msr,sol,rhs,option);
          break;
        case skymatrix:/*---------------------- system matrix is skyline matrix */
          solver_colsol(actsolv,actintra,oll->sysarray[0].sky,sol,rhs,option);
          break;
        case spoolmatrix:/*-------------------- system matrix is spooles matrix */
          solver_spo_oll(actsolv,actintra,oll,oll->sysarray[0].spo,sol,rhs,option);
          break;
        case ccf:
	  solver_umfpack(actsolv,actintra,oll->sysarray[0].ccf,sol,rhs,option);
	break;
        default:
          dserror("Unknown solver typ for oll");
          break;   
      }

      break;
    default:
      dserror("Unknown option");
      break;
  }
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of solver_oll */


/*! @} (documentation module close)*/
