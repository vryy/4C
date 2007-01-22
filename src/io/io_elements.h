/*!
\file
\brief Element version support.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Elements evolve over time. Thus the data an elements writes to the
binary files, especially the element parameters, will change. (There
is an element parameter chunk. Each element stores its parameters ---
number of nodes, number of gauss points, ... --- there.) But we
have to be able to read old element versions. This is why some support
for element versions is needed.

It works like this: Each element has an associated version number, a
plain integer, that is increased each time the structure of the
element output changes. For each element version number the meaning of
the element parameters (the meaning of each number) must be
clear. This allows to set up variables those name contains the meaning
and those value is the number's position (index) in the element
parameter entry (an entry of the element parameter chunk.) These
variables are used to access the element parameters both for writing
and reading.

This scheme brings some flexibility while it introduces only a limited
overhead and a small maintainance burden. But please be careful. You
cannot fix everything with such a crude approach. If used improper
your binary files might still become unmaintainable.

Mistakes to avoid:

- Don't ever decrease the element version number.

- Don't remove index variables that have been introduced in previous
  versions. Even if your new element version doesn't know about these
  numbers anymore.

- Don't try to introduce variable element parameter entry lengths. The
  length of each element parameter entry is given by the element
  type. Use additional chunks instead.

- Don't use a non-cvs version of ccarat to produce binary files you
  want to keep.

\author u.kue
\date 11/04

*/

#ifdef BINIO

#ifndef IO_ELEMENTS_H
#define IO_ELEMENTS_H

#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"


/* output version of the part common to all elements */
#define ELEMENT_IO_VERSION 1


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure common to all element
  output.

  These variables depend on the output version. There is a function
  that sets there variables depending of the version number. The idea
  is to allow element evolution and still be able to read old output.

  \author u.kue
  \date 11/04
  \sa setup_element_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _ELEMENT_VARIABLES
{
  /* version for which the parameters are initialized */
  INT version;

  INT ep_size_Id;
  INT ep_size_eltyp;
  INT ep_size_distyp;
  INT ep_size_numnp;

} ELEMENT_VARIABLES;


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure common to all element
  output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_element_variables
 */
/*----------------------------------------------------------------------*/
extern ELEMENT_VARIABLES element_variables;


/*======================================================================*/
/*======================================================================*/

#ifdef D_SHELL8

/* output version of shell8 element */
#define SHELL8_IO_VERSION 1


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the shell8 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_shell8_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _SHELL8_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  /* element parameter indices */
  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;
  INT ep_size_nGP_tri;
  INT ep_size_forcetyp;

  INT ep_value_sdc;

} SHELL8_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the shell8 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_shell8_variables
 */
/*----------------------------------------------------------------------*/
extern SHELL8_VARIABLES shell8_variables;


#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_SHELL9

/* output version of shell9 element */
#define SHELL9_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the shell9 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_shell9_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _SHELL9_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  /* element parameter indices */
  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;
  INT ep_size_nGP_tri;
  INT ep_size_forcetyp;

} SHELL9_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the shell9 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_shell9_variables
 */
/*----------------------------------------------------------------------*/
extern SHELL9_VARIABLES shell9_variables;


#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_BRICK1

/* output version of brick1 element */
#define BRICK1_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the brick1 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_brick1_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _BRICK1_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;
  INT ep_size_stresstyp;

} BRICK1_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the brick1 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_brick1_variables
 */
/*----------------------------------------------------------------------*/
extern BRICK1_VARIABLES brick1_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2

/* output version of fluid2 element */
#define FLUID2_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_fluid2_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID2_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;

} FLUID2_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid2_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID2_VARIABLES fluid2_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2_PRO

/* output version of fluid2_pro element */
#define FLUID2_PRO_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2_pro output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_fluid2_pro_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID2_PRO_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;

} FLUID2_PRO_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2_pro output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid2_pro_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID2_PRO_VARIABLES fluid2_pro_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3_PRO

/* output version of fluid3_pro element */
#define FLUID3_PRO_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_pro output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_pro_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID3_PRO_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;

} FLUID3_PRO_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_pro output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_pro_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID3_PRO_VARIABLES fluid3_pro_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2_IS

/* output version of fluid2_is element */
#define FLUID2_IS_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2_is output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/06
  \sa setup_fluid2_is_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID2_IS_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;

} FLUID2_IS_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2_is output.

  The only instance of this structure.

  \author u.kue
  \date 11/06
  \sa setup_fluid2_is_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID2_IS_VARIABLES fluid2_is_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3_IS

/* output version of fluid3_is element */
#define FLUID3_IS_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_is output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/06
  \sa setup_fluid3_is_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID3_IS_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;

} FLUID3_IS_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_is output.

  The only instance of this structure.

  \author u.kue
  \date 11/06
  \sa setup_fluid3_is_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID3_IS_VARIABLES fluid3_is_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2TU

/* output version of fluid2tu element */
#define FLUID2TU_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2tu output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_fluid2tu_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID2TU_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;

} FLUID2TU_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2tu output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid2tu_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID2TU_VARIABLES fluid2tu_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3

/* output version of fluid3 element */
#define FLUID3_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID3_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;

} FLUID3_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID3_VARIABLES fluid3_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3_F

/* output version of fluid3_fast element */
#define FLUID3_F_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_fast output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_fast_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _FLUID3_FAST_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;

} FLUID3_FAST_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_fast output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_fast_variables
 */
/*----------------------------------------------------------------------*/
extern FLUID3_FAST_VARIABLES fluid3_fast_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_ALE

/* output version of ale2 element */
#define ALE2_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the ale2 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_ale2_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _ALE2_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;

} ALE2_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the ale2 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_ale2_variables
 */
/*----------------------------------------------------------------------*/
extern ALE2_VARIABLES ale2_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_ALE

/* output version of ale3 element */
#define ALE3_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the ale3 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_ale3_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _ALE3_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;

} ALE3_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the ale3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_ale3_variables
 */
/*----------------------------------------------------------------------*/
extern ALE3_VARIABLES ale3_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_WALL1

/* output version of wall1 element */
#define WALL1_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the wall1 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_wall1_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _WALL1_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;
  INT ep_size_nGP3;
  INT ep_size_stresstyp;

  INT ep_value_thick;

} WALL1_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the wall1 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_wall1_variables
 */
/*----------------------------------------------------------------------*/
extern WALL1_VARIABLES wall1_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_BEAM3

/* output version of beam3 element */
#define BEAM3_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the beam3 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_beam3_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _BEAM3_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;

} BEAM3_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the beam3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_beam3_variables
 */
/*----------------------------------------------------------------------*/
extern BEAM3_VARIABLES beam3_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_AXISHELL

/* output version of axishell element */
#define AXISHELL_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the axishell output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_axishell_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _AXISHELL_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  /* Is this an element parameter? Not really... */
  INT ep_value_thick0;
  INT ep_value_thick1;

} AXISHELL_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the axishell output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_axishell_variables
 */
/*----------------------------------------------------------------------*/
extern AXISHELL_VARIABLES axishell_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_INTERF

/* output version of interf element */
#define INTERF_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the interf output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_interf_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _INTERF_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP;
  INT ep_size_stresstyp;

  INT ep_value_thick;

} INTERF_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the interf output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_interf_variables
 */
/*----------------------------------------------------------------------*/
extern INTERF_VARIABLES interf_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_WALLGE

/* output version of wallge element */
#define WALLGE_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the wallge output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_wallge_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _WALLGE_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_nGP0;
  INT ep_size_nGP1;
  INT ep_size_nGP2;
  INT ep_size_nGP3;
  INT ep_size_stresstyp;

  INT ep_value_thick;

} WALLGE_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the wallge output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_wallge_variables
 */
/*----------------------------------------------------------------------*/
extern WALLGE_VARIABLES wallge_variables;

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_SOLID3

/* output version of solid3 element */
#define SOLID3_IO_VERSION 1

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the solid3 output.

  These variables depend on the output version of the element. There
  is a function that sets there variables depending of the version
  number. The idea is to allow element evolution and still be able to
  read old output.

  \author u.kue
  \date 11/04
  \sa setup_solid3_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _SOLID3_VARIABLES
{
  /* If there is no element version in the control file we cannot
   * initialize this structure. It's illegal to use it then. */
  INT init;

  /* version for which the parameters are initialized */
  INT version;

  /* the length of the integer and the double entry */
  INT ep_size_length;
  INT ep_value_length;

  INT ep_size_gpnum0;
  INT ep_size_gpnum1;
  INT ep_size_gpnum2;
  INT ep_size_gptot;
  INT ep_size_stresstype;

} SOLID3_VARIABLES;

/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the solid3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_solid3_variables
 */
/*----------------------------------------------------------------------*/
extern SOLID3_VARIABLES solid3_variables;

#endif

/*======================================================================*/
/*======================================================================*/


/* Now here we have the same thing for nodes. But there is only one
 * type of nodes, so things are much easier here. */


/* output version to the part common to all node */
#define NODE_IO_VERSION 1


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure common to all node
  output.

  These variables depend on the output version. There is a function
  that sets there variables depending of the version number. The idea
  is to allow node evolution and still be able to read old output.

  \author u.kue
  \date 12/04
  \sa setup_node_variables
 */
/*----------------------------------------------------------------------*/
typedef struct _NODE_VARIABLES
{
  /* version for which the parameters are initialized */
  INT version;

  INT coords_size_Id;
  INT coords_size_numele;

  INT coords_size_eleid;
  /* what follows are numele element ids. */

} NODE_VARIABLES;


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure common to all node
  output.

  The only instance of this structure.

  \author u.kue
  \date 12/04
  \sa setup_node_variables
 */
/*----------------------------------------------------------------------*/
extern NODE_VARIABLES node_variables;


/*======================================================================*/
/*======================================================================*/


/*----------------------------------------------------------------------*/
/*!
  \brief Set all element variables according to the current element versions.

  This is used for writing.

  \author u.kue
  \date 11/04
 */
/*----------------------------------------------------------------------*/
void setup_element_variables_current();


/*----------------------------------------------------------------------*/
/*!
  \brief Set all element variables according to the version read from
  the map.

  \author u.kue
  \date 11/04
  \sa setup_element_variables
 */
/*----------------------------------------------------------------------*/
void setup_element_variables_map(MAP* group);


/*----------------------------------------------------------------------*/
/*!
  \brief Print the current io versions of all available elements.

  \author u.kue
  \date 11/04
 */
/*----------------------------------------------------------------------*/
void out_print_element_versions(FILE* file);


#endif /* IO_ELEMENTS_H */
#endif /* BINIO */
