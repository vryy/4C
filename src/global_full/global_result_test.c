/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#ifdef RESULTTEST
#include "../headers/standardtypes.h"
#include "../axishell/axishell.h"
#include "../shell9/shell9.h"
#include "../fluid_full/fluid_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;


/*
 * an array of expected results
 */
struct _RESULTDESCR      *resultdescr;


/*
 * Read a positions description. See if its name equals |name|. If so
 * read |nargs| integer arguments and store them in |args|.
 */
INT parse_position_descr(CHAR* position, CHAR* name, INT nargs, INT* args)
{
  CHAR* lp;
  INT name_length;

#ifdef DEBUG 
  dstrc_enter("parse_position_descr");
#endif
  
  lp = strpbrk(position, "(");
  if (lp == NULL) {
    dserror("Missing left parenthesis in position description");
  }
  name_length = lp - position;

  if (strncmp(position, name, name_length) == 0) {
    CHAR* pos = lp;
    CHAR* rp;
    INT i;

    for (i=0; i<nargs-1; ++i) {
      args[i] = atoi(pos+1) - 1;
      pos = strpbrk(pos+1, ",");
      if (pos == NULL) {
        dserror("Missing comma in position description");
      }
    }
    args[nargs-1] = atoi(pos+1) - 1;
    rp = strpbrk(pos+1, ")");
    if (rp == NULL) {
      dserror("Missing right parenthesis in position description");
    }
    return 1;
  }
  else {
    return 0;
  }

#ifdef DEBUG 
  dstrc_exit();
#endif
}


/*
 * return the specified value
 */
DOUBLE get_node_result_value(NODE* actnode, CHAR* position)
{
  INT args[2];
  
#ifdef DEBUG 
  dstrc_enter("get_node_result_value");
#endif

  if (parse_position_descr(position, "sol", 2, args) == 1) {
    return actnode->sol.a.da[args[0]][args[1]];
  }
  else {
    dserror("Unknown position specifier");
    return 1234567890;
  }
  
#ifdef DEBUG 
  dstrc_exit();
#endif
}


static int compare_values(FILE* err, DOUBLE actresult, DOUBLE givenresult, RESULTDESCR *res)
{
#ifdef DEBUG 
  dstrc_enter("compare_values");
#endif

  fprintf(err,"actual = %24.16f, given = %24.16f\n", actresult, givenresult);
  if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult))) {
    printf("RESULTCHECK: %s is NAN!", res->name);
    return 1;
  }
  if (FABS(actresult-givenresult) > res->tolerance) {
    printf("RESULTCHECK: %s not correct!", res->name);
    return 1;
  }

#ifdef DEBUG 
  dstrc_exit();
#endif
  return 0;
}



/*!---------------------------------------------------------------------  
\brief testing of result 

<pre>                                                         genk 10/03  

Before checking in the latest version it's necessery to check the whole
program. In this context it seems to be useful to check the numerical
results, too.
Based on the problem titel one can add verified results to this function.

</pre>  
\return void                                                              
\warning Titel must not be longer than one line!

------------------------------------------------------------------------*/
void global_result_test() 
{
  FIELD  *alefield,*structfield,*fluidfield;
  NODE   *actnode;
  DOUBLE  actresult, givenresult;
  FILE   *err = allfiles.out_err;
  INT i;
  INT nerr = 0;

#ifdef DEBUG 
  dstrc_enter("global_result_test");
#endif

  if (genprob.numff>-1) fluidfield  = &(field[genprob.numff]);
  if (genprob.numsf>-1) structfield = &(field[genprob.numsf]);
  if (genprob.numaf>-1) alefield    = &(field[genprob.numaf]);

  if (genprob.numresults>0) {
    /* let's do it in a fency style :) */
    printf("\n[37;1mChecking results ...[m\n");
  }

  for (i=0; i<genprob.numresults; ++i) {
    FIELD* current_field;
    RESULTDESCR* res = &(resultdescr[i]);
    
    switch (res->field) {
    case fluid:
      current_field = fluidfield;
      break;
    case ale:
      current_field = alefield;
      break;
    case structure:
      current_field = structfield;
      break;
    default:
      dserror("Unknown field typ");
    }

    if (res->node != -1) {
      actnode = &(current_field->dis[res->dis].node[res->node]);
      actresult = get_node_result_value(actnode, res->position);
      nerr += compare_values(err, actresult, res->value, res);
    }
    else if (res->element != -1) {
      ELEMENT* actelement = &(current_field->dis[res->dis].element[res->element]);
    
#ifdef D_AXISHELL
      if (actelement->eltyp == el_axishell) {
        INT args[3];
        if (parse_position_descr(res->position, "stress_GP", 3, args) == 1) {
          actresult = actelement->e.saxi->stress_GP.a.d3[args[0]][args[1]][args[2]];
          nerr += compare_values(err, actresult, res->value, res);
        }
        else {
          dserror("Unknown position specifier");
        }
      }
#endif
    
#ifdef D_SHELL9
      if (actelement->eltyp == el_shell9) {
        INT args[3];
        if (parse_position_descr(res->position, "stresses", 3, args) == 1) {
          actresult = actelement->e.s9->stresses.a.d3[args[0]][args[1]][args[2]];
          nerr += compare_values(err, actresult, res->value, res);
        }
        else {
          dserror("Unknown position specifier");
        }
      }
#endif

    }
    else {
      /* special cases that need further code support */
      switch (genprob.probtyp) {
      case prb_fluid:
#ifdef D_FLUID
        fluid_cal_error(fluidfield,res->dis);
        break;
#endif
      default:
        break;
      }
    }
  }

  if (nerr > 0) {
    dserror("Result check failed");
  }


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of global_result_test */
#endif
