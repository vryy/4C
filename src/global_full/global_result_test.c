/*!
\file
\brief Compare result values (after the calculation finished) with predefined ones.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Here are the functions that compare the calculated results with the
values given in the input file's result description section. It's
possible to test values that are stored at specific nodes and those
that are contained in the element structures.

Most of the functions here are internal and are concerned with
understanding the result description or comparing two given
values. You might, however, want to enhance global_result_test() if
you need to test elements that are not supported yet.

It's also possible to implement special test cases, for example if an
analytical solution is known. In this case no values are compared but
a field specific function is called that might test whatever it
likes.

One particular thing is that the testing has to honor the distribution
of elements and nodes to different processors. Every element or node
lifes on just one processor, that's why it's alright when a specified
element could not be found. It's assumed that it belongs to another
processor and will be checked there. But of course this way it might
happen that it is not checked anywhere and this failure goes
undetected.

\author uk
\date 06/04

*/

#ifdef RESULTTEST
#include "../headers/standardtypes.h"
#include "../axishell/axishell.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
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


/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;


/*!
 * \brief An array that contains the expected results.
 *
 * An array that contains the expected results.
 * Along with the real values the position in the mesh ist stored. The
 * data is read from the input file.
 *
 * \author uk
 * \date 06/04
 */
struct _RESULTDESCR      *resultdescr;


/*----------------------------------------------------------------------*/
/*!
 \brief Read a positions description.

 Read a positions description. See if its name equals \a name. If so
 read \a nargs integer arguments and store them in \a args.

 \param position   a string of the form "name(x1,...,xn)" read from the input file
 \param name       the expected name
 \param nargs      the expected number of arguments
 \param args       an array of size \a nargs that is going to be filled

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
static INT parse_position_descr(CHAR* position, CHAR* name, INT nargs, INT* args)
{
  CHAR* lp;
  INT name_length;
  INT ret = -1234567890;

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
    ret = 1;
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
 \brief return the specified value

 \param actnode    a node
 \param position   a string of the form "name(x1,...,xn)" read from the input file
                   It describes a value in one of the solution array
                   of the given node.

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
static DOUBLE get_node_result_value(NODE* actnode, CHAR* position)
{
  INT args[2];
  DOUBLE ret = -1234567890;

#ifdef DEBUG
  dstrc_enter("get_node_result_value");
#endif

  /* for debugging...
  {
    INT i, j;
    for (i=0; i<actnode->sol.fdim; ++i) {
      printf("%d ", i);
      for (j=0; j<actnode->sol.sdim; ++j) {
        printf("%f ", actnode->sol.a.da[i][j]);
      }
      printf("\n");
    }
    printf("\n");
  }
  */

  if (parse_position_descr(position, "sol", 2, args) == 1) {
    ret = actnode->sol.a.da[args[0]][args[1]];
  }
  else if (parse_position_descr(position, "sol_increment", 2, args) == 1) {
    ret = actnode->sol_increment.a.da[args[0]][args[1]];
  }
  else if (parse_position_descr(position, "sol_residual", 2, args) == 1) {
    ret = actnode->sol_residual.a.da[args[0]][args[1]];
  }
  else if (parse_position_descr(position, "sol_mf", 2, args) == 1) {
    ret = actnode->sol_mf.a.da[args[0]][args[1]];
  }
  else {
    dserror("Unknown position specifier: %s", position);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
 \brief Compare \a actresult with \a givenresult and return 0 if they are
 considered to be equal.

 Compare \a actresult with \a givenresult and return 0 if they are
 considered to be equal.

 \param err        the file where to document both values
 \param res        the describtion of the expected result including name and tolerance

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
static int compare_values(FILE* err, DOUBLE actresult, DOUBLE givenresult, RESULTDESCR *res)
{
  INT ret = 0;

#ifdef DEBUG
  dstrc_enter("compare_values");
#endif

  fprintf(err,"actual = %24.16f, given = %24.16f\n", actresult, givenresult);
  if (!(FABS(FABS(actresult-givenresult)-FABS(actresult-givenresult)) < res->tolerance)) {
    printf("RESULTCHECK: %s is NAN!\n", res->name);
    ret = 1;
  }
  else if (FABS(actresult-givenresult) > res->tolerance) {
    printf("RESULTCHECK: %s not correct. actresult=%f, givenresult=%f\n", res->name, actresult, givenresult);
    ret = 1;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}



/*----------------------------------------------------------------------*/
/*!
 \brief Find node with id \a nodenum.

 Find node with id \a nodenum. Only the given partition and
 discretization is searched.

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
static NODE* find_node(PARTITION* part, INT disnum, INT nodenum)
{
  INT i;
  NODE* res = NULL;
  PARTDISCRET* pdis;

#ifdef DEBUG
  dstrc_enter("find_node");
#endif

  pdis = &(part->pdis[disnum]);
  for (i=0; i<pdis->numnp; ++i) {
    if (pdis->node[i]->Id == nodenum) {
      res = pdis->node[i];
      break;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return res;
}


/*----------------------------------------------------------------------*/
/*!
 \brief Find element with id \a elenum.

 Find element with id \a elenum. Only the given partition and
 discretization is searched.

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
static ELEMENT* find_element(PARTITION* part, INT disnum, INT elenum)
{
  INT i;
  ELEMENT* res = 0;
  PARTDISCRET* pdis;

#ifdef DEBUG
  dstrc_enter("find_element");
#endif

  pdis = &(part->pdis[disnum]);
  for (i=0; i<pdis->numele; ++i) {
    if (pdis->element[i]->Id == elenum) {
      res = pdis->element[i];
      break;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return res;
}


/*----------------------------------------------------------------------*/
/*!
\brief testing of result

Before checking in the latest version it's necessery to check the whole
program. In this context it seems to be useful to check the numerical
results, too. This function compares all predefined values with the
calculated ones. It's called once, after the calculation finished.

\author genk
\date 10/03
*/
/*----------------------------------------------------------------------*/
void global_result_test()
{
#ifndef PARALLEL
  FIELD  *alefield     = NULL;
  FIELD  *structfield  = NULL;
  FIELD  *fluidfield   = NULL;
#endif
  PARTITION *alepart     = NULL;
  PARTITION *structpart  = NULL;
  PARTITION *fluidpart   = NULL;
  NODE   *actnode;
  DOUBLE  actresult;
  FILE   *err = allfiles.out_err;
  INT i;
  INT nerr = 0;

#ifdef DEBUG
  dstrc_enter("global_result_test");
#endif
  if (genprob.numresults<=0) goto end;

#ifndef PARALLEL
  if (genprob.numff>-1) fluidfield  = &(field[genprob.numff]);
  if (genprob.numsf>-1) structfield = &(field[genprob.numsf]);
  if (genprob.numaf>-1) alefield    = &(field[genprob.numaf]);
#endif

  if (genprob.numff>-1) fluidpart  = &(partition[genprob.numff]);
  if (genprob.numsf>-1) structpart = &(partition[genprob.numsf]);
  if (genprob.numaf>-1) alepart    = &(partition[genprob.numaf]);

  /* let's do it in a fancy style :) */
  printf("\n[35;1mChecking results ...[m\n");
  fprintf(err,"\n===========================================\n");
  fprintf(err,"Checking results ...\n");


  for (i=0; i<genprob.numresults; ++i) {
    PARTITION* actpart = NULL;
    RESULTDESCR* res = &(resultdescr[i]);

    switch (res->field) {
    case fluid:
      actpart = fluidpart;
      break;
    case ale:
      actpart = alepart;
      break;
    case structure:
      actpart = structpart;
      break;
    default:
      dserror("Unknown field typ");
    }

    if (res->node != -1) {
      /* We have a value at a node. Find the node and compare with the
       * value there. */
      actnode = find_node(actpart, res->dis, res->node);
      if (actnode != 0) {
        actresult = get_node_result_value(actnode, res->position);
        nerr += compare_values(err, actresult, res->value, res);
      }
      else
        dswarning(1,12);
    }
    else if (res->element != -1) {
      /* We have a value at an element. Find the element. If we have
       * it its depending on the element type what needs to be done. */
      ELEMENT* actelement = find_element(actpart, res->dis, res->element);
      if (actelement == 0) {
        continue;
      }

#ifdef D_AXISHELL
      if (actelement->eltyp == el_axishell) {
        INT args[3];
        if (parse_position_descr(res->position, "stress_GP", 3, args) == 1) {
          actresult = actelement->e.saxi->stress_GP.a.d3[args[0]][args[1]][args[2]];
          nerr += compare_values(err, actresult, res->value, res);
        }
        else if (parse_position_descr(res->position, "stress_ND", 3, args) == 1) {
          actresult = actelement->e.saxi->stress_ND.a.d3[args[0]][args[1]][args[2]];
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

#ifdef D_WALL1
      if (actelement->eltyp == el_wall1) {
        INT args[3];
        if (parse_position_descr(res->position, "stress_GP", 3, args) == 1) {
          actresult = actelement->e.w1->stress_GP.a.d3[args[0]][args[1]][args[2]];
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
#ifndef PARALLEL
        fluid_cal_error(fluidfield,res->dis);
#endif
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
  else
    fprintf(err,"===========================================\n");

  printf("\n[35;1mOK[m\n");


/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of global_result_test */
#endif
