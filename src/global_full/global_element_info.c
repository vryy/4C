/*!
\file
\brief Information about the elements.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is the central place where information about elements is
kept. When elements are added or changed this file needs to be
changed, too. The binary io module, the filters and maybe others
depend on it.

\author u.kue
\date 09/04

*/

#include "../headers/standardtypes.h"
#include "global_element_info.h"

#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../beam3/beam3.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
#include "../axishell/axishell.h"
#include "../interf/interf.h"
#include "../wallge/wallge.h"


/*----------------------------------------------------------------------*/
/*!
  \brief The table with all information about the elements.

  Major and minor number specify an element type. Here we have the
  table that contains everything we'd might need to know about the
  element.

  The actual information ccarat gets out of this table is the size of
  the stress entry in the binary output file. Some elements abuse
  their stress entry to store other things than stresses, so the
  numbers here are arbitrary and cannot be calculated. (If need
  araises we can add other size information here as well.)

  On the other hand the filters use this table to get the information
  that is encoded in the major/minor pair. Thus the element's shape,
  its node and gauss point numbers are contained here, too. These must
  be consistent with the functions below that figure out the minor
  numbers.

  \warning This table is accessed by index. The entries must be in the
  same order as the major element number.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
ELEMENT_INFO element_info[el_count] = {
  { "none", {
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "shell1", {
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "shell8", {
      /* The shell8 code uses MAXGAUSS to determine the size of the
       * stress array. But we cannot do this because the output file
       * must not be dependent on the arbitrarily defined MAXGAUSS
       * flag. */
      { 4, quad4,    8, 18*8*2 },
      { 8, quad8,    9, 18*9*2 },
      { 9, quad9,    9, 18*9*2 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "shell9", {
      { 4, quad4, 4, -1 },
      { 4, quad4, 9, -1 },
      { 8, quad8, 4, -1 },
      { 8, quad8, 9, -1 },
      { 9, quad9, 4, -1 },
      { 9, quad9, 9, -1 }
    } },
  { "brick1", {
      {  8, hex8,   8, 27* 8 },
      { 20, hex20, 27, 27*20 },
      { 27, hex27, 27, 27*27 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 }
    } },
  { "wall1", {
      /* For stress calculations only four gauss point version(s) are
       * supported. */
      { 3, tri3,      1, 1*4 },
      { 4, quad4, 4, 4*(7+2) },
      { 8, quad8, 9, 9*(7+2) },
      { 9, quad9, 9, 9*(7+2) },
      { 8, quad8, 4, 4*(7+2) },
      { 0, dis_none,  0,   0 }
    } },
  { "beam3", {
      { 2, line2,    1, 1*6 },
      { 2, line2,    2, 2*6 },
      { 3, line3,    2, 2*6 },
      { 3, line3,    3, 3*6 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "fluid2", {
      { 4, quad4, 4, 0 },
      { 8, quad8, 9, 0 },
      { 9, quad9, 9, 0 },
      { 3, tri3,  4, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "fluid2_pro", {
      { 4, quad4, 4, 0 },
      { 8, quad8, 9, 0 },
      { 9, quad9, 9, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "fluid2_tu", {
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "fluid3", {
      {  8, hex8,      8, 0 },
      { 20, hex20,    27, 0 },
      { 27, hex27,    27, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 }
    } },
  { "fluid3_fast", {
      {  8, hex8,      8, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 },
      {  0, dis_none,  0, 0 }
    } },
  { "ale2", {
      { 4, quad4, 1, 0 },
      { 4, quad4, 4, 0 },
      { 3, tri3,  2, 0 },
      { 3, tri3,  3, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "ale3", {
      { 8, hex8, 1, 0 },
      { 8, hex8, 8, 0 },
      { 4, tet4, 1, 0 },
      { 4, tet4, 4, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "axishell", {
      { 2, line2,    1, 5 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "interf", {
      { 4, quad4, 2, 2*5 },
      { 8, quad8, 3, 3*5 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } },
  { "wallge", {
      { 4, quad4,    4, 4*6 },
      { 8, quad8,    9, 9*6 },
      { 9, quad9,    9, 9*6 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 },
      { 0, dis_none, 0, 0 }
    } }
};


/*======================================================================*/
/* Here are the functions that calculate the minor type. There are
 * three functions per element (major) type but only one that does the
 * work. This way there's only one place to change.
 *
 * However, please make sure to update ``write_variant_group`` as well
 * if you change these functions. */
/*======================================================================*/

#ifdef D_SHELL8


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a shell8 element.

  \param numnp    (i) the number of nodes

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_shell8_minor_type(INT numnp)
{
  INT type;

#ifdef DEBUG
  dstrc_enter("get_shell8_minor_type");
#endif

  if (numnp==4) {
    type = MINOR_SHELL8_22;
  }
  else if (numnp==8) {
    type = MINOR_SHELL8_8_33;
  }
  else if (numnp==9) {
    type = MINOR_SHELL8_9_33;
  }
  else {
    dserror("unknown minor version for shell8");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 7 parameter shell
  element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_shell8_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_shell8_minor_type");
#endif

  dsassert(actele->eltyp == el_shell8, "shell8 expected");
  type = get_shell8_minor_type(actele->numnp);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 7 parameter shell element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_shell8_minor_type(MAP* group)
{
  INT type = -1;
  INT numnp;
#ifdef DEBUG
  dstrc_enter("restore_shell8_minor_type");
#endif

  numnp = map_read_int(group, "numnp");
  type = get_shell8_minor_type(numnp);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_SHELL9


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a shell9 element.

  \param numnp    (i) the number of nodes
  \param nGP0     (i) first number of gauss
  \param nGP1     (i) second number of gauss

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_shell9_minor_type(INT numnp, INT nGP0, INT nGP1)
{
  INT type;

#ifdef DEBUG
  dstrc_enter("get_shell9_minor_type");
#endif

  if ((numnp==4) && (nGP0==2) && (nGP1==2)) {
    type = MINOR_SHELL9_4_22;
  }
  else if ((numnp==4) && (nGP0==3) && (nGP1==3)) {
    type = MINOR_SHELL9_4_33;
  }
  else if ((numnp==8) && (nGP0==2) && (nGP1==2)) {
    type = MINOR_SHELL9_8_22;
  }
  else if ((numnp==8) && (nGP0==3) && (nGP1==3)) {
    type = MINOR_SHELL9_8_33;
  }
  else if ((numnp==9) && (nGP0==2) && (nGP1==2)) {
    type = MINOR_SHELL9_9_22;
  }
  else if ((numnp==9) && (nGP0==3) && (nGP1==3)) {
    type = MINOR_SHELL9_9_33;
  }
  else {
    dserror("unknown minor version for shell9");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given multi layer shell
  element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_shell9_minor_type(ELEMENT* actele)
{
  INT type = -1;
  INT numnp;
  INT nGP0;
  INT nGP1;

#ifdef DEBUG
  dstrc_enter("find_shell9_minor_type");
#endif

  dsassert(actele->eltyp == el_shell9, "shell9 expected");

  numnp = actele->numnp;
  nGP0 = actele->e.s9->nGP[0];
  nGP1 = actele->e.s9->nGP[1];

  type = get_shell9_minor_type(numnp, nGP0, nGP1);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a multi layer shell element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_shell9_minor_type(MAP* group)
{
  INT type = -1;
  INT numnp;
  INT nGP0;
  INT nGP1;

#ifdef DEBUG
  dstrc_enter("restore_shell9_minor_type");
#endif

  numnp = map_read_int(group, "numnp");
  nGP0 = map_read_int(group, "nGP0");
  nGP1 = map_read_int(group, "nGP1");

  type = get_shell9_minor_type(numnp, nGP0, nGP1);


#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_BRICK1


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a brick1 element.

  \param numnp    (i) the number of nodes

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_brick1_minor_type(INT numnp)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_brick1_minor_type");
#endif

  if (numnp==8) {
    type = MINOR_BRICK1_222;
  }
  else if (numnp==20) {
    type = MINOR_BRICK1_20_333;
  }
  else if (numnp==27) {
    type = MINOR_BRICK1_27_333;
  }
  else {
    dserror("unknown minor version for brick1");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given structural brick
  element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_brick1_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_brick1_minor_type");
#endif

  dsassert(actele->eltyp == el_brick1, "brick1 expected");
  type = get_brick1_minor_type(actele->numnp);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a structural brick element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_brick1_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_brick1_minor_type");
#endif

  type = get_brick1_minor_type(map_read_int(group, "numnp"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_WALL1


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a wall1 element.

  \param numnp    (i) the number of nodes
  \param nGP      (i) the number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_wall1_minor_type(INT numnp, INT nGP)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_wall1_minor_type");
#endif

  if ((numnp==4) && (nGP==1)) {
    type = MINOR_WALL1_11;
  }
  else if ((numnp==4) && (nGP==2)) {
    type = MINOR_WALL1_22;
  }
  else if ((numnp==8) && (nGP==3)) {
    type = MINOR_WALL1_8_33;
  }
  else if ((numnp==9) && (nGP==3)) {
    type = MINOR_WALL1_9_33;
  }
  else if ((numnp==8) && (nGP==2)) {
    type = MINOR_WALL1_8_22;
  }
  else {
    dserror("unknown minor version for wall1");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D plane stress -
  plane strain element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_wall1_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_wall1_minor_type");
#endif

  dsassert(actele->eltyp == el_wall1, "wall1 expected");
  type = get_wall1_minor_type(actele->numnp, actele->e.w1->nGP[0]);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D plane stress - plane
  strain element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_wall1_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_wall1_minor_type");
#endif

  type = get_wall1_minor_type(map_read_int(group, "numnp"),
                              map_read_int(group, "nGP0"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_BEAM3


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a beam3 element.

  \param numnp    (i) the number of nodes
  \param nGP      (i) the number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_beam3_minor_type(INT numnp, INT nGP)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_beam3_minor_type");
#endif

  if ((numnp==2) && (nGP==1)) {
    type = MINOR_BEAM3_21;
  }
  else if ((numnp==2) && (nGP==2)) {
    type = MINOR_BEAM3_22;
  }
  else if ((numnp==3) && (nGP==2)) {
    type = MINOR_BEAM3_32;
  }
  else if ((numnp==3) && (nGP==3)) {
    type = MINOR_BEAM3_33;
  }
  else {
    dserror("unknown minor version for beam3");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given structural 3D-beam
  element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_beam3_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_beam3_minor_type");
#endif

  dsassert(actele->eltyp == el_beam3, "beam3 expected");
  type = get_beam3_minor_type(actele->numnp, actele->e.b3->nGP[0]);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a structural 3D-beam element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_beam3_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_beam3_minor_type");
#endif

  type = get_beam3_minor_type(map_read_int(group, "numnp"),
                              map_read_int(group, "nGP0"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_FLUID2


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a fluid2 element.

  \param numnp    (i) the number of nodes
  \param nGP      (i) the number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_fluid2_minor_type(INT numnp, INT nGP)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_fluid2_minor_type");
#endif

  if (numnp==4) {
    type = MINOR_FLUID2_22;
  }
  else if (numnp==8) {
    type = MINOR_FLUID2_8_33;
  }
  else if (numnp==9) {
    type = MINOR_FLUID2_9_33;
  }
  else if ((numnp==3) && (nGP==4)) {
    type = MINOR_FLUID2_3_4;
  }
  else {
    dserror("unknown minor version for fluid2");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D fluid element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid2_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_fluid2_minor_type");
#endif

  dsassert(actele->eltyp == el_fluid2, "fluid2 expected");
  type = get_fluid2_minor_type(actele->numnp, actele->e.f2->nGP[0]);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D fluid element's variant
  group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid2_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_fluid2_minor_type");
#endif

  type = get_fluid2_minor_type(map_read_int(group, "numnp"),
                               map_read_int(group, "nGP0"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_FLUID2_PRO


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a fluid2_pro element.

  \param numnp    (i) the number of nodes

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_fluid2_pro_minor_type(INT numnp)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_fluid2_pro_minor_type");
#endif

  if (numnp==4) {
    type = MINOR_FLUID2_PRO_22;
  }
  else if (numnp==8) {
    type = MINOR_FLUID2_PRO_8_33;
  }
  else if (numnp==9) {
    type = MINOR_FLUID2_PRO_9_33;
  }
  else {
    dserror("unknown minor version for fluid2_pro");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D fluid pro element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid2_pro_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_fluid2_pro_minor_type");
#endif

  dsassert(actele->eltyp == el_fluid2_pro, "fluid2 pro expected");
  type = get_fluid2_pro_minor_type(actele->numnp);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D fluid pro element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid2_pro_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_fluid2_pro_minor_type");
#endif

  type = get_fluid2_pro_minor_type(map_read_int(group, "numnp"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_FLUID2TU


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a fluid2_tu element.

  \param numnp    (i) the number of nodes

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
#if 0
static INT get_fluid2_tu_minor_type(INT numnp)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_fluid2_tu_minor_type");
#endif

  dserror("fluid2_tu not yet supported");

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}
#endif


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D fluid element for
  turbulence.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid2_tu_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_fluid2_tu_minor_type");
#endif

  dsassert(actele->eltyp == el_fluid2_tu, "fluid2 for turbulence expected");
  dserror("fluid2_tu not yet supported");

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D fluid element for
  turbulence's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid2_tu_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_fluid2_tu_minor_type");
#endif

  dserror("fluid2_tu not yet supported");

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_FLUID3


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a fluid3 element.

  \param numnp    (i) the number of nodes

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_fluid3_minor_type(INT numnp)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_fluid3_minor_type");
#endif

  if (numnp==8) {
    type = MINOR_FLUID3_222;
  }
  else if (numnp==20) {
    type = MINOR_FLUID3_20_333;
  }
  else if (numnp==27) {
    type = MINOR_FLUID3_27_333;
  }
  else {
    dserror("unknown minor version for fluid3");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 3D fluid element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid3_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_fluid3_minor_type");
#endif

  dsassert(actele->eltyp == el_fluid3, "fluid3 expected");
  type = get_fluid3_minor_type(actele->numnp);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 3D fluid element's variant
  group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid3_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_fluid3_minor_type");
#endif

  type = get_fluid3_minor_type(map_read_int(group, "numnp"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_FLUID3_F


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a fluid3 element.

  \param numnp    (i) the number of nodes

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
static INT get_fluid3_fast_minor_type(INT numnp)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_fluid3_fast_minor_type");
#endif

  if (numnp==8) {
    type = MINOR_FLUID3_FAST_222;
  }
  else {
    dserror("unknown minor version for fluid3");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 3D fluid element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 11/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid3_fast_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_fluid3_fast_minor_type");
#endif

  dsassert(actele->eltyp == el_fluid3_fast, "fluid3 expected");
  type = get_fluid3_fast_minor_type(actele->numnp);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 3D fluid element's variant
  group.

  \author u.kue
  \date 11/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid3_fast_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_fluid3_fast_minor_type");
#endif

  type = get_fluid3_fast_minor_type(map_read_int(group, "numnp"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_ALE


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a ale2 element.

  \param numnp    (i) the number of nodes
  \param nGP0     (i) first number of gauss points
  \param nGP1     (i) second number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_ale2_minor_type(INT numnp, INT nGP0, INT nGP1)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_ale2_minor_type");
#endif

  if ((numnp==4) && (nGP0==1) && (nGP1==1)) {
    type = MINOR_ALE_11;
  }
  else if ((numnp==4) && (nGP0==2) && (nGP1==2)) {
    type = MINOR_ALE_22;
  }
  else if ((numnp==3) && (nGP0==1)) {
    type = MINOR_ALE_TRI_1;
  }
  else if ((numnp==3) && (nGP0==3)) {
    type = MINOR_ALE_TRI_3;
  }
  else {
    dserror("unknown minor version for ale2");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D pseudo structural
  ale element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_ale2_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_ale2_minor_type");
#endif

  dsassert(actele->eltyp == el_ale2, "ale2 expected");
  type = get_ale2_minor_type(actele->numnp,
                             actele->e.ale2->nGP[0],
                             actele->e.ale2->nGP[1]);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D pseudo structural ale
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_ale2_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_ale2_minor_type");
#endif

  type = get_ale2_minor_type(map_read_int(group, "numnp"),
                             map_read_int(group, "nGP0"),
                             map_read_int(group, "nGP1"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_ALE


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a ale3 element.

  \param numnp    (i) the number of nodes
  \param nGP0     (i) first number of gauss points
  \param nGP1     (i) second number of gauss points
  \param nGP2     (i) third number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_ale3_minor_type(INT numnp, INT nGP0, INT nGP1, INT nGP2)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_ale3_minor_type");
#endif

  if ((numnp==8) && (nGP0==1) && (nGP1==1) && (nGP2==1)) {
    type = MINOR_ALE_111;
  }
  else if ((numnp==8) && (nGP0==2) && (nGP1==2) && (nGP2==2)) {
    type = MINOR_ALE_222;
  }
  else if ((numnp==4) && (nGP0 == 1)) {
    type = MINOR_ALE_TET_1;
  }
  else if ((numnp==4) && (nGP0 == 4)) {
    type = MINOR_ALE_TET_4;
  }
  else {
    dserror("unknown minor version for ale3");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 3D pseudo structural
  ale element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_ale3_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_ale3_minor_type");
#endif

  dsassert(actele->eltyp == el_ale3, "ale3 expected");
  type = get_ale3_minor_type(actele->numnp,
                             actele->e.ale3->nGP[0],
                             actele->e.ale3->nGP[1],
                             actele->e.ale3->nGP[2]);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 3D pseudo structural ale
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_ale3_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_ale3_minor_type");
#endif

  type = get_ale3_minor_type(map_read_int(group, "numnp"),
                             map_read_int(group, "nGP0"),
                             map_read_int(group, "nGP1"),
                             map_read_int(group, "nGP2"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_AXISHELL


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a axishell element.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_axishell_minor_type()
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_axishell_minor_type");
#endif

  type = MINOR_AXISHELL;

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 1D axisymmetrical
  shell element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_axishell_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_axishell_minor_type");
#endif

  dsassert(actele->eltyp == el_axishell, "axishell expected");
  type = get_axishell_minor_type();

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 1D axisymmetrical shell
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_axishell_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_axishell_minor_type");
#endif

  type = get_axishell_minor_type();

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_INTERF


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a interf element.

  \param nGP    (i) number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_interf_minor_type(INT nGP)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_interf_minor_type");
#endif

  if (nGP==2) {
    type = MINOR_INTERF_22;
  }
  else if (nGP==3) {
    type = MINOR_INTERF_33;
  }
  else {
    dserror("unknown minor version for interf");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 1D interface element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_interf_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_interf_minor_type");
#endif

  dsassert(actele->eltyp == el_interf, "interface element expected");
  type = get_interf_minor_type(actele->e.interf->nGP);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 1D interface element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_interf_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_interf_minor_type");
#endif

  type = get_interf_minor_type(map_read_int(group, "nGP"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/

#ifdef D_WALLGE


/*----------------------------------------------------------------------*/
/*!
  \brief Determine the minor type of a wallge element.

  \param numnp    (i) the number of nodes
  \param nGP      (i) the number of gauss points

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static INT get_wallge_minor_type(INT numnp, INT nGP)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("get_wallge_minor_type");
#endif

  if (nGP==2) {
    type = MINOR_WALLGE_22;
  }
  else if ((numnp==8) && (nGP==3)) {
    type = MINOR_WALL1_8_33;
  }
  else if ((numnp==9) && (nGP==3)) {
    type = MINOR_WALL1_9_33;
  }
  else {
    dserror("unknown minor version for wallge");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given gradient enhanced wall
  element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_wallge_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_wallge_minor_type");
#endif

  dsassert(actele->eltyp == el_wallge, "wallge element expected");
  type = get_wallge_minor_type(actele->numnp, actele->e.wallge->nGP[0]);

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a gradient enhanced wall
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_wallge_minor_type(MAP* group)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("restore_wallge_minor_type");
#endif

  type = get_wallge_minor_type(map_read_int(group, "numnp"),
                               map_read_int(group, "nGP0"));

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}

#endif

/*======================================================================*/


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given element.

  This is the big one that takes elements of any type and calls the
  appropriate special function.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT find_minor_type(ELEMENT* actele)
{
  INT type = -1;
#ifdef DEBUG
  dstrc_enter("find_minor_type");
#endif

  switch (actele->eltyp) {
#ifdef D_SHELL8
  case el_shell8:               /* 7 parameter shell element */
    type = find_shell8_minor_type(actele);
    break;
#endif
#ifdef D_SHELL9
  case el_shell9:               /* multi layer shell element */
    type = find_shell9_minor_type(actele);
    break;
#endif
#ifdef D_BRICK1
  case el_brick1:               /* structural brick element */
    type = find_brick1_minor_type(actele);
    break;
#endif
#ifdef D_WALL1
  case el_wall1:                /* 2D plane stress - plane strain element */
    type = find_wall1_minor_type(actele);
    break;
#endif
#ifdef D_BEAM3
  case el_beam3:                /* structural 3D-beam element */
    type = find_beam3_minor_type(actele);
    break;
#endif
#ifdef D_FLUID2
  case el_fluid2:               /* 2D fluid element */
    type = find_fluid2_minor_type(actele);
    break;
#endif
#ifdef D_FLUID2_PRO
  case el_fluid2_pro:           /* 2D fluid element */
    type = find_fluid2_pro_minor_type(actele);
    break;
#endif
#ifdef D_FLUID2TU
  case el_fluid2_tu:            /* 2D fluid element for turbulence */
    type = find_fluid2_tu_minor_type(actele);
    break;
#endif
#ifdef D_FLUID3
  case el_fluid3:               /* 3D fluid element */
    type = find_fluid3_minor_type(actele);
    break;
#endif
#ifdef D_FLUID3_F
  case el_fluid3_fast:          /* 3D fluid element */
    type = find_fluid3_fast_minor_type(actele);
    break;
#endif
#ifdef D_ALE
  case el_ale2:                 /* 2D pseudo structural ale element */
    type = find_ale2_minor_type(actele);
    break;
#endif
#ifdef D_ALE
  case el_ale3:                 /* 3D pseudo structural ale element */
    type = find_ale3_minor_type(actele);
    break;
#endif
#ifdef D_AXISHELL
  case el_axishell:             /* 1D axisymmetrical shell element */
    type = find_axishell_minor_type(actele);
    break;
#endif
#ifdef D_INTERF
  case el_interf:               /* 1D interface element (combination only with wall) */
    type = find_interf_minor_type(actele);
    break;
#endif
#ifdef D_WALLGE
  case el_wallge:               /* gradient enhanced wall element */
    type = find_wallge_minor_type(actele);
    break;
#endif
  default:
    dserror("unknown element type %d", actele->eltyp);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return type;
}




/*----------------------------------------------------------------------*/
/*!
  \brief Find the major and minor element types belonging to the given
  criteria.

  Because result files stay around for a long time it's not guaranteed
  that the major/minor numbers used internally by the current ccarat
  will be the same as those used in old files. So a mechanism to
  translate these numbers is needed, and that's what this functions
  provides.

  There are ``element_variant`` groups in ccarat's control files that
  define the major and minor numbers used in that file. The \a group
  argument must be one such group. The criteria in this group are used
  to return the internal major and minor numbers.

  \param group    (i) group read from a control file
  \param major    (o) internal major number to this element type
  \param minor    (o) internal minor number to this element type

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restore_element_type(MAP* group, INT* major, INT* minor)
{
  INT ele;
  CHAR* element_type;
#ifdef DEBUG
  dstrc_enter("restore_element_type");
#endif

  element_type = map_read_string(group, "type");
  for (ele=0; ele<el_count; ++ele) {
    if (strcmp(element_type, element_info[ele].name)==0) {
      *major = ele;
      break;
    }
  }

  if (ele==el_count) {
    dserror("unknown element major type '%s'", element_type);
  }

  switch (ele) {
#ifdef D_SHELL8
  case el_shell8:               /* 7 parameter shell element */
    *minor = restore_shell8_minor_type(group);
    break;
#endif
#ifdef D_SHELL9
  case el_shell9:               /* multi layer shell element */
    *minor = restore_shell9_minor_type(group);
    break;
#endif
#ifdef D_BRICK1
  case el_brick1:               /* structural brick element */
    *minor = restore_brick1_minor_type(group);
    break;
#endif
#ifdef D_WALL1
  case el_wall1:                /* 2D plane stress - plane strain element */
    *minor = restore_wall1_minor_type(group);
    break;
#endif
#ifdef D_BEAM3
  case el_beam3:                /* structural 3D-beam element */
    *minor = restore_beam3_minor_type(group);
    break;
#endif
#ifdef D_FLUID2
  case el_fluid2:               /* 2D fluid element */
    *minor = restore_fluid2_minor_type(group);
    break;
#endif
#ifdef D_FLUID2_PRO
  case el_fluid2_pro:           /* 2D fluid element */
    *minor = restore_fluid2_pro_minor_type(group);
    break;
#endif
#ifdef D_FLUID2TU
  case el_fluid2_tu:            /* 2D fluid element for turbulence */
    *minor = restore_fluid2_tu_minor_type(group);
    break;
#endif
#ifdef D_FLUID3
  case el_fluid3:               /* 3D fluid element */
    *minor = restore_fluid3_minor_type(group);
    break;
#endif
#ifdef D_FLUID3_F
  case el_fluid3_fast:          /* 3D fluid element */
    *minor = restore_fluid3_fast_minor_type(group);
    break;
#endif
#ifdef D_ALE
  case el_ale2:                 /* 2D pseudo structural ale element */
    *minor = restore_ale2_minor_type(group);
    break;
#endif
#ifdef D_ALE
  case el_ale3:                 /* 3D pseudo structural ale element */
    *minor = restore_ale3_minor_type(group);
    break;
#endif
#ifdef D_AXISHELL
  case el_axishell:             /* 1D axisymmetrical shell element */
    *minor = restore_axishell_minor_type(group);
    break;
#endif
#ifdef D_INTERF
  case el_interf:               /* 1D interface element (combination only with wall) */
    *minor = restore_interf_minor_type(group);
    break;
#endif
#ifdef D_WALLGE
  case el_wallge:               /* gradient enhanced wall element */
    *minor = restore_wallge_minor_type(group);
    break;
#endif
  default:
    dserror("unknown element type %d", ele);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the element variant criterion to the given type.

  This functions writes a group that defines the criterions for the
  given element type to the control file. Using this information a
  filter can reconstruct what a certain element type means.

  \param out       (i) the file to write into
  \param ele       (i) major number
  \param variant   (i) minor number

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void write_variant_group(FILE* out, ELEMENT_TYP ele, INT variant)
{
#ifdef DEBUG
  dstrc_enter("write_variant_group");
#endif

  fprintf(out,
          "    element_variant:\n"
          "        type = \"%s\"\n"
          "        major = %d\n"
          "        minor = %d\n",
          element_info[ele].name, ele, variant);

  switch (ele) {

#ifdef D_SHELL8

  case el_shell8:               /* 7 parameter shell element */
    switch (variant) {
    case MINOR_SHELL8_22:       /* 4-noded shell8 2x2 GP */
      fprintf(out, "        numnp = 4\n");
      break;
    case MINOR_SHELL8_8_33:     /* 8-noded shell8 3x3 GP */
      fprintf(out, "        numnp = 8\n");
      break;
    case MINOR_SHELL8_9_33:     /* 9-noded shell8 3x3 GP */
      fprintf(out, "        numnp = 9\n");
      break;
    default:
      dserror("unknown variant %d for shell8", variant);
    }
    break;

#endif
#ifdef D_SHELL9

  case el_shell9:               /* multi layer shell element */
    switch (variant) {
    case MINOR_SHELL9_4_22:     /* 4-noded shell9 2x2 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 2\n");
      fprintf(out, "        nGP1 = 2\n");
      break;
    case MINOR_SHELL9_4_33:     /* 4-noded shell9 3x3 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 3\n");
      fprintf(out, "        nGP1 = 3\n");
      break;
    case MINOR_SHELL9_8_22:     /* 8-noded shell9 2x2 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 2\n");
      fprintf(out, "        nGP1 = 2\n");
      break;
    case MINOR_SHELL9_8_33:     /* 8-noded shell9 3x3 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 3\n");
      fprintf(out, "        nGP1 = 3\n");
      break;
    case MINOR_SHELL9_9_22:     /* 9-noded shell9 2x2 GP */
      fprintf(out, "        numnp = 9\n");
      fprintf(out, "        nGP0 = 2\n");
      fprintf(out, "        nGP1 = 2\n");
      break;
    case MINOR_SHELL9_9_33:     /* 9-noded shell9 3x3 GP */
      fprintf(out, "        numnp = 9\n");
      fprintf(out, "        nGP0 = 3\n");
      fprintf(out, "        nGP1 = 3\n");
      break;
    default:
      dserror("unknown variant %d for shell9", variant);
    }
    break;

#endif
#ifdef D_BRICK1

  case el_brick1:               /* structural brick element */
    switch (variant) {
    case MINOR_BRICK1_222:      /* 8-noded brick1 2x2x2 GP */
      fprintf(out, "        numnp = 8\n");
      break;
    case MINOR_BRICK1_20_333:   /* 20 noded brick1 3x3x3 GP */
      fprintf(out, "        numnp = 20\n");
      break;
    case MINOR_BRICK1_27_333:   /* 27 noded brick1 3x3x3 GP */
      fprintf(out, "        numnp = 27\n");
      break;
    default:
      dserror("unknown variant %d for brick1", variant);
    }
    break;

#endif
#ifdef D_WALL1

  case el_wall1:                /* 2D plane stress - plane strain element */
    switch (variant) {
    case MINOR_WALL1_11:        /* 3-noded wall1 1x1 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 1\n");
      break;
    case MINOR_WALL1_22:        /* 4-noded wall1 2x2 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 2\n");
      break;
    case MINOR_WALL1_8_33:      /* 8-noded wall1 3x3 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    case MINOR_WALL1_9_33:      /* 9-noded wall1 3x3 GP */
      fprintf(out, "        numnp = 9\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    case MINOR_WALL1_8_22:      /* 8-noded wall1 2x2 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 2\n");
      break;
    default:
      dserror("unknown variant %d for wall1", variant);
    }
    break;

#endif
#ifdef D_BEAM3

  case el_beam3:                /* structural 3D-beam element */
    switch (variant) {
    case MINOR_BEAM3_21:        /* 2-noded beam3 1 GP */
      fprintf(out, "        numnp = 2\n");
      fprintf(out, "        nGP0 = 1\n");
      break;
    case MINOR_BEAM3_22:        /* 2-noded beam3 2 GP */
      fprintf(out, "        numnp = 2\n");
      fprintf(out, "        nGP0 = 2\n");
      break;
    case MINOR_BEAM3_32:        /* 3-noded beam3 2 GP */
      fprintf(out, "        numnp = 3\n");
      fprintf(out, "        nGP0 = 2\n");
      break;
    case MINOR_BEAM3_33:        /* 3-noded beam3 3 GP */
      fprintf(out, "        numnp = 3\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    default:
      dserror("unknown variant %d for beam3", variant);
    }
    break;

#endif
#ifdef D_FLUID2

  case el_fluid2:               /* 2D fluid element */
    switch (variant) {
    case MINOR_FLUID2_22:       /* 4-noded fluid2 2x2 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 2\n");
      break;
    case MINOR_FLUID2_8_33:     /* 8-noded fluid2 3x3 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    case MINOR_FLUID2_9_33:     /* 9-noded fluid2 3x3 GP */
      fprintf(out, "        numnp = 9\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    case MINOR_FLUID2_3_4:      /* 3-noded fluid2 4 GP */
      fprintf(out, "        numnp = 3\n");
      fprintf(out, "        nGP0 = 4\n");
      break;
    default:
      dserror("unknown variant %d for fluid2", variant);
    }
    break;

#endif
#ifdef D_FLUID2_PRO

  case el_fluid2_pro:           /* 2D fluid element */
    switch (variant) {
    case MINOR_FLUID2_PRO_22:   /* 4-noded fluid2_pro 2x2 GP */
      fprintf(out, "        numnp = 4\n");
      break;
    case MINOR_FLUID2_PRO_8_33: /* 8-noded fluid2_pro 3x3 GP */
      fprintf(out, "        numnp = 8\n");
      break;
    case MINOR_FLUID2_PRO_9_33: /* 9-noded fluid2_pro 3x3 GP */
      fprintf(out, "        numnp = 9\n");
      break;
    default:
      dserror("unknown variant %d for fluid2 pro", variant);
    }
    break;

#endif
#ifdef D_FLUID2TU

  case el_fluid2_tu:            /* 2D fluid element for turbulence */
    switch (variant) {
    default:
      dserror("unknown variant %d for fluid2 tu", variant);
    }
    break;

#endif
#ifdef D_FLUID3_F

  case el_fluid3_fast:          /* 3D fluid element */
    switch (variant) {
    case MINOR_FLUID3_222:      /* 8-noded fluid3 2x2x2 GP */
      fprintf(out, "        numnp = 8\n");
      break;
    default:
      dserror("unknown variant %d for fluid3_fast", variant);
    }
    break;

#endif
#ifdef D_FLUID3

  case el_fluid3:               /* 3D fluid element */
    switch (variant) {
    case MINOR_FLUID3_222:      /* 8-noded fluid3 2x2x2 GP */
      fprintf(out, "        numnp = 8\n");
      break;
    case MINOR_FLUID3_20_333:   /* 20-noded fluid3 3x3x3 GP */
      fprintf(out, "        numnp = 20\n");
      break;
    case MINOR_FLUID3_27_333:   /* 27-noded fluid3 3x3x3 GP */
      fprintf(out, "        numnp = 27\n");
      break;
    default:
      dserror("unknown variant %d for fluid3", variant);
    }
    break;

#endif
#ifdef D_ALE

  case el_ale2:                 /* 2D pseudo structural ale element */
    switch (variant) {
    case MINOR_ALE_11:          /* 4-noded ale 1x1 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 1\n");
      fprintf(out, "        nGP1 = 1\n");
      break;
    case MINOR_ALE_22:          /* 4-noded ale 2x2 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 2\n");
      fprintf(out, "        nGP1 = 2\n");
      break;
    case MINOR_ALE_TRI_1:       /* 3-noded tri ale 1 GP */
      fprintf(out, "        numnp = 3\n");
      fprintf(out, "        nGP0 = 1\n");
      fprintf(out, "        nGP1 = 0\n"); /* unused, but ease reading */
      break;
    case MINOR_ALE_TRI_3:       /* 3-noded tri ale 3 GP */
      fprintf(out, "        numnp = 3\n");
      fprintf(out, "        nGP0 = 3\n");
      fprintf(out, "        nGP1 = 0\n"); /* unused, but ease reading */
      break;
    default:
      dserror("unknown variant %d for ale2", variant);
    }
    break;

#endif
#ifdef D_ALE

  case el_ale3:                 /* 3D pseudo structural ale element */
    switch (variant) {
    case MINOR_ALE_111:         /* 8-noded ale 1x1x1 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 1\n");
      fprintf(out, "        nGP1 = 1\n");
      fprintf(out, "        nGP2 = 1\n");
      break;
    case MINOR_ALE_222:         /* 8-noded ale 2x2x2 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 2\n");
      fprintf(out, "        nGP1 = 2\n");
      fprintf(out, "        nGP2 = 2\n");
      break;
    case MINOR_ALE_TET_1:       /* 4-noded tet ale 1 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 1\n");
      fprintf(out, "        nGP1 = 0\n"); /* unused, but ease reading */
      fprintf(out, "        nGP2 = 0\n"); /* unused, but ease reading */
      break;
    case MINOR_ALE_TET_4:       /* 4-noded tet ale 4 GP */
      fprintf(out, "        numnp = 4\n");
      fprintf(out, "        nGP0 = 4\n");
      fprintf(out, "        nGP1 = 0\n"); /* unused, but ease reading */
      fprintf(out, "        nGP2 = 0\n"); /* unused, but ease reading */
      break;
    default:
      dserror("unknown variant %d for ale3", variant);
    }
    break;

#endif
#ifdef D_AXISHELL

  case el_axishell:             /* 1D axisymmetrical shell element */
    switch (variant) {
    case MINOR_AXISHELL:        /* 2-noded axishell */
      /* There is just one variant -- nothing to be done. */
      break;
    default:
      dserror("unknown variant %d for axishell", variant);
    }
    break;

#endif
#ifdef D_INTERF

  case el_interf:               /* 1D interface element (combination only with wall) */
    switch (variant) {
    case MINOR_INTERF_22:       /* interface 2x2 GP */
      fprintf(out, "        nGP = 2\n");
      break;
    case MINOR_INTERF_33:       /* interface 3x3 GP */
      fprintf(out, "        nGP = 3\n");
      break;
    default:
      dserror("unknown variant %d for interf", variant);
    }
    break;

#endif
#ifdef D_WALLGE

  case el_wallge:               /* gradient enhanced wall element */
    switch (variant) {
    case MINOR_WALLGE_22:       /* gradient enhanced wall 2x2 GP */
      fprintf(out, "        numnp = 4\n"); /* ??? */
      fprintf(out, "        nGP0 = 2\n");
      break;
    case MINOR_WALLGE_8_33:     /* gradient enhanced wall 3x3 GP */
      fprintf(out, "        numnp = 8\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    case MINOR_WALLGE_9_33:     /* gradient enhanced wall 3x3 GP */
      fprintf(out, "        numnp = 9\n");
      fprintf(out, "        nGP0 = 3\n");
      break;
    default:
      dserror("unknown variant %d for wallge", variant);
    }
    break;

#endif

  default:
    dserror("unknown element type %d", ele);
  }

  fprintf(out, "\n");

#ifdef DEBUG
  dstrc_exit();
#endif
}
