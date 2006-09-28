/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef TRILINOS_PACKAGE

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
//#include "Epetra_MpiComm.h"
#else
//#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
//#include "Epetra_LinearProblem.h"

//#include "Amesos_Klu.h"
//#include "Amesos_Umfpack.h"
//#include "Amesos_Lapack.h"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "../solver/solver_trilinos_ml.H"

using namespace std;
using namespace Teuchos;

static void ml_compute_nullspace(DISCRET*          actdis, 
                                 ELEMENT_TYP       etype,
                                 const int         numdf,
                                 const int         dimns,
                                 const Epetra_Map& map,
                                 double*           nullspace);

/*----------------------------------------------------------------------*
 |                                                           m.gee 9/06 |
 | create ML's parameter list from input parameters                     |
 | Also create the nullspace from the discretization                    |
 *----------------------------------------------------------------------*/
void create_ml_parameterlist(struct _SOLVAR         *actsolv,
                             ParameterList&          mllist,
                             FIELD*                  actfield,
                             int                     disnum,
                             const Epetra_CrsMatrix& matrix,
                             double**                nullspace)
{
#ifdef DEBUG
  dstrc_enter("create_ml_parameterlist");
#endif
  //-----------------------------------------------------------------------
  AZVAR* azvar  = actsolv->azvar;
  //-----------------------------------------------------------------------
  // set default values
  ML_Epetra::SetDefaults("SA",mllist);

  // set principal type of preconditioner
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
      if (matrix.Comm().NumProc()==1)
        mllist.set("aggregation: block scaling",true); // crashes in parallel
      else
        mllist.set("aggregation: block scaling",false); // crashes in parallel
    break;
    default: dserror("Unknown type of preconditioner");
  }

  // set various parameters
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
  

  // set smoothers
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
      case 6:  
        mllist.set("smoother: type "+(string)levelstr,"Amesos-Superludist");        
      break;
      default: dserror("Unknown type of smoother for ML"); break;
    }
  }
  
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
  }
  
  // set number of pde equations and nullspace dimension from fieldtyp and elementtyp
  // construct nullspace vectors
  int numdf=1; // default value
  int dimns=1; // default value
  DISCRET* actdis = &(actfield->dis[disnum]);
  if (!actdis->numele) dserror("Discretization passed to ML has no elements");
  ELEMENT* ele =  &(actdis->element[0]);
  ELEMENT_TYP etype = ele->eltyp;
  switch (etype)
  {
    case el_shell8:
      numdf = 6;
      dimns = 6;
    break;
    case el_brick1:
      numdf = 3;
      dimns = 6;
    break;
    case el_ale3:
      numdf = 3;
      dimns = 6;
    break;
    case el_wall1:
      numdf = 2;
      dimns = 3;
    break; 
    case el_ale2:
      numdf = 2;
      dimns = 3;
    break;
    case el_fluid2:
      numdf = 3;
      dimns = 3;
    break;
    case el_fluid3:
      numdf = 4;
      dimns = 4;
    break;
    default: dserror("Unknown type of element or element not compiled in");
  }
  mllist.set("PDE equations",numdf);
  mllist.set("null space: dimension",dimns);
  mllist.set("null space: type","pre-computed");
  mllist.set("null space: add default vectors",false);
  if (*nullspace) delete (*nullspace);
  *nullspace = new double[ matrix.RowMap().NumMyElements()*dimns];
  ml_compute_nullspace(actdis,etype,numdf,dimns,matrix.RowMap(),*nullspace);
  mllist.set("null space: vectors",*nullspace);
  //-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of create_ml_parameterlist */





/*----------------------------------------------------------------------*
 |                                                           m.gee 9/06 |
 | compute the nullspace from the discretization                        |
 *----------------------------------------------------------------------*/
void ml_compute_nullspace(DISCRET*          actdis, 
                          ELEMENT_TYP       etype,
                          const int         numdf,
                          const int         dimns,
                          const Epetra_Map& map,
                          double*           nullspace)
{
#ifdef DEBUG
  dstrc_enter("ml_compute_nullspace");
#endif
  //-----------------------------------------------------------------------
  const int  myrank = map.Comm().MyPID();
  const int  lrows  = map.NumMyElements();
  
  double* mode[6];
  if (dimns>6) dserror("Only upto 6 nullspace modes supported");
  for (int i=0; i<dimns; ++i)
    mode[i] = &(nullspace[i*lrows]);
  
  // nodal center of structure
  double x0[3] = {0.0,0.0,0.0};
  int    count = 0;
  for (int i=0; i<actdis->numnp; ++i)
  {
    for (int j=0; j<3; ++j) x0[j] += actdis->node[i].x[j];
    ++count;
  }
  for (int j=0; j<3; ++j) x0[j] /= count;

#ifdef D_SHELL8
  // get director field for shell8
  Epetra_SerialDenseMatrix dir(false);
  if (etype == el_shell8)
  {
    // count nodes I own
    count=0;
    for (int i=0; i<actdis->numnp; ++i)
      if (actdis->node[i].proc==myrank) ++count;
    dir.Reshape(count,3);
    int count2=0;
    for (int i=0; i<actdis->numnp; ++i)
    {
      if (actdis->node[i].proc!=myrank) continue;
      NODE*    actnode = &(actdis->node[i]);
      // get an adjacent element
      ELEMENT* actele  = actnode->element[0];
      // loop element's nodes to find our node
      int j=0;
      for (j=0; j<actele->numnp; ++j)
        if (actele->node[j]->Id == actnode->Id) break;
      dsassert(j!=actele->numnp,"Can't find matching node on element");
      // get half shell thickness
      double h2 = actele->e.s8->thick_node.a.dv[j]/2.0;
      // get director
      dir(count2,0) = actele->e.s8->a3ref.a.da[0][j]*h2;
      dir(count2,1) = actele->e.s8->a3ref.a.da[1][j]*h2;
      dir(count2,2) = actele->e.s8->a3ref.a.da[2][j]*h2;
      ++count2;
    }
    if (count!=count2) dserror("Local number of nodes wrong");
  }
#endif

  /* the rigid body modes for structures are:
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       0          a3        -a2
  dy  |    0       0       0      -a3         0          a1
  dz  |    0       0       0       a2        -a1         0
  */

  if (etype==el_shell8 ||
      etype==el_brick1 ||
      etype==el_ale3   ||
      etype==el_wall1  ||
      etype==el_ale2   )
  {
    count=0;
    for (int i=0; i<actdis->numnp; ++i)
    {
      if (actdis->node[i].proc!=myrank) continue;
      NODE* actnode = &(actdis->node[i]);
      for (int j=0; j<actnode->numdf; ++j)
      {
        int dof   = actnode->dof[j];
        int index = map.LID(dof);
        if (index==-1) dserror("Cannot find local for global dof");
        switch(j) // j is the degree of freedom
        {
        case 0:
            mode[0][index] = 1.0;
            mode[1][index] = 0.0;
          if (dimns==6) // 3D case
          {
            mode[2][index] = 0.0;
            mode[3][index] = 0.0;
            mode[4][index] = actnode->x[2] - x0[2];
            mode[5][index] = -actnode->x[1] + x0[1];
          }
          else if (dimns==3) // 2D case
            mode[2][index] = -actnode->x[1] + x0[1];
          else dserror("Something wrong with modes");
        break;
        case 1:
            mode[0][index] = 0.0;
            mode[1][index] = 1.0;
          if (dimns==6) // 3D case
          {
            mode[2][index] = 0.0;
            mode[3][index] = -actnode->x[2] + x0[2];
            mode[4][index] = 0.0;
            mode[5][index] = actnode->x[0] - x0[0];
          }
          else if (dimns==3) // 2D case
            mode[2][index] = actnode->x[0] - x0[0];
          else dserror("Something wrong with modes");
        break;
        case 2: // pure 3D case
            mode[0][index] = 0.0;
            mode[1][index] = 0.0;
            mode[2][index] = 1.0;
            mode[3][index] = actnode->x[1] - x0[1]; 
            mode[4][index] = -actnode->x[0] + x0[0];
            mode[5][index] = 0.0;
        break;
        case 3:
            mode[0][index] = 0.0;
            mode[1][index] = 0.0;
            mode[2][index] = 0.0;
#ifdef D_SHELL8
            mode[3][index] = 0.0; 
            mode[4][index] = dir(count,2);
            mode[5][index] = -dir(count,1);
#endif
        break;
        case 4:
            mode[0][index] = 0.0;
            mode[1][index] = 0.0;
            mode[2][index] = 0.0;
#ifdef D_SHELL8
            mode[3][index] = -dir(count,2); 
            mode[4][index] = 0.0;
            mode[5][index] = dir(count,0);
#endif
        break;
        case 5:
            mode[0][index] = 0.0;
            mode[1][index] = 0.0;
            mode[2][index] = 0.0;
#ifdef D_SHELL8
            mode[3][index] = dir(count,1); 
            mode[4][index] = -dir(count,0);
            mode[5][index] = 0.0;
#endif
        break;
        default: dserror("Only modes 0 to 5 supported"); break;
        }
      }
      ++count;
    }
  }
      
  /* the rigid body modes for fluids are:
        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3] 
  ----------------------------------------
  x   |    1       0       0       0      
  y   |    0       1       0       0      
  z   |    0       0       1       0      
  p   |    0       0       0       1      
  */
  if (etype==el_fluid2 ||
      etype==el_fluid3 )
  {
    for (int i=0; i<actdis->numnp; ++i)
    {
      if (actdis->node[i].proc!=myrank) continue;
      NODE* actnode = &(actdis->node[i]);
      for (int j=0; j<actnode->numdf; ++j)
      {
        int dof   = actnode->dof[j];
        int index = map.LID(dof);
        if (index==-1) dserror("Cannot find local for global dof");
        switch(j) // j is the degree of freedom
        {
        case 0:
            mode[0][index] = 1.0;
            mode[1][index] = 0.0;
            mode[2][index] = 0.0;
          if (dimns==4)
            mode[3][index] = 0.0;
        break;
        case 1:
            mode[0][index] = 0.0;
            mode[1][index] = 1.0;
            mode[2][index] = 0.0;
          if (dimns==4)
            mode[3][index] = 0.0;
        break;
        case 2: // pure 3D case
            mode[0][index] = 0.0;
            mode[1][index] = 0.0;
            mode[2][index] = 1.0;
          if (dimns==4)
            mode[3][index] = 0.0;
        break;
        case 3:
            mode[0][index] = 0.0;
            mode[1][index] = 0.0;
            mode[2][index] = 0.0;
            mode[3][index] = 1.0;
        break;
        default: dserror("Only modes 0 to 3 supported"); break;
        }
      }
    }
  }
  //-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of ml_compute_nullspace */


#endif // TRILINOS_PACKAGE
