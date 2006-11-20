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

#include "io_elements.h"


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
ELEMENT_VARIABLES element_variables;


static void setup_element_variables(INT version, ELEMENT_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_element_variables");
#endif

  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_Id = 0;
    vars->ep_size_eltyp = 1;
    vars->ep_size_distyp = 2;
    vars->ep_size_numnp = 3;
    break;
  default:
    dserror("element output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*======================================================================*/
/*======================================================================*/

#ifdef D_SHELL8


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the shell8 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_shell8_variables
 */
/*----------------------------------------------------------------------*/
SHELL8_VARIABLES shell8_variables;


static void setup_shell8_variables(INT version, SHELL8_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_shell8_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;
    vars->ep_size_nGP_tri = 7;
    vars->ep_size_forcetyp = 8;

    vars->ep_value_sdc = 0;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 9;
    vars->ep_value_length = 1;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("shell8 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_SHELL9


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the shell9 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_shell9_variables
 */
/*----------------------------------------------------------------------*/
SHELL9_VARIABLES shell9_variables;


static void setup_shell9_variables(INT version, SHELL9_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_shell9_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;
    vars->ep_size_nGP_tri = 7;
    vars->ep_size_forcetyp = 8;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 9;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("shell9 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_BRICK1


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the brick1 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_brick1_variables
 */
/*----------------------------------------------------------------------*/
BRICK1_VARIABLES brick1_variables;


static void setup_brick1_variables(INT version, BRICK1_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_brick1_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;
    vars->ep_size_stresstyp = 7;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 8;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("brick1 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid2_variables
 */
/*----------------------------------------------------------------------*/
FLUID2_VARIABLES fluid2_variables;


static void setup_fluid2_variables(INT version, FLUID2_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid2_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 6;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid2 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2_PRO


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2_pro output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid2_pro_variables
 */
/*----------------------------------------------------------------------*/
FLUID2_PRO_VARIABLES fluid2_pro_variables;


static void setup_fluid2_pro_variables(INT version, FLUID2_PRO_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid2_pro_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 6;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid2_pro output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3_PRO


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_pro output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_pro_variables
 */
/*----------------------------------------------------------------------*/
FLUID3_PRO_VARIABLES fluid3_pro_variables;


static void setup_fluid3_pro_variables(INT version, FLUID3_PRO_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid3_pro_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 7;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid3_pro output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2_IS


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2_is output.

  The only instance of this structure.

  \author u.kue
  \date 11/06
  \sa setup_fluid2_is_variables
 */
/*----------------------------------------------------------------------*/
FLUID2_IS_VARIABLES fluid2_is_variables;


static void setup_fluid2_is_variables(INT version, FLUID2_IS_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid2_is_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 6;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid2_is output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3_IS


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_is output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_is_variables
 */
/*----------------------------------------------------------------------*/
FLUID3_IS_VARIABLES fluid3_is_variables;


static void setup_fluid3_is_variables(INT version, FLUID3_IS_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid3_is_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 7;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid3_is output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID2TU


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid2tu output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid2tu_variables
 */
/*----------------------------------------------------------------------*/
FLUID2TU_VARIABLES fluid2tu_variables;


static void setup_fluid2tu_variables(INT version, FLUID2TU_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid2tu_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 6;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid2tu output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_variables
 */
/*----------------------------------------------------------------------*/
FLUID3_VARIABLES fluid3_variables;


static void setup_fluid3_variables(INT version, FLUID3_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid3_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 7;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid3 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_FLUID3_F


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the fluid3_fast output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_fluid3_fast_variables
 */
/*----------------------------------------------------------------------*/
FLUID3_FAST_VARIABLES fluid3_fast_variables;


static void setup_fluid3_fast_variables(INT version, FLUID3_FAST_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_fluid3_fast_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 7;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("fluid3_fast output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_ALE


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the ale2 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_ale2_variables
 */
/*----------------------------------------------------------------------*/
ALE2_VARIABLES ale2_variables;


static void setup_ale2_variables(INT version, ALE2_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_ale2_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 6;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("ale2 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_ALE


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the ale3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_ale3_variables
 */
/*----------------------------------------------------------------------*/
ALE3_VARIABLES ale3_variables;


static void setup_ale3_variables(INT version, ALE3_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_ale3_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 7;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("ale3 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_WALL1


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the wall1 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_wall1_variables
 */
/*----------------------------------------------------------------------*/
WALL1_VARIABLES wall1_variables;


static void setup_wall1_variables(INT version, WALL1_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_wall1_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;
    vars->ep_size_nGP3 = 7;
    vars->ep_size_stresstyp = 8;

    vars->ep_value_thick = 0;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 9;
    vars->ep_value_length = 1;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("wall1 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_BEAM3


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the beam3 output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_beam3_variables
 */
/*----------------------------------------------------------------------*/
BEAM3_VARIABLES beam3_variables;


static void setup_beam3_variables(INT version, BEAM3_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_beam3_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 5;
    vars->ep_value_length = 0;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("beam3 output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_AXISHELL


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the axishell output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_axishell_variables
 */
/*----------------------------------------------------------------------*/
AXISHELL_VARIABLES axishell_variables;


static void setup_axishell_variables(INT version, AXISHELL_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_axishell_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    /* Is this an element parameter? Not really... */
    vars->ep_value_thick0 = 0;
    vars->ep_value_thick1 = 1;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 4;
    vars->ep_value_length = 2;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("axishell output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_INTERF


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the interf output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_interf_variables
 */
/*----------------------------------------------------------------------*/
INTERF_VARIABLES interf_variables;


static void setup_interf_variables(INT version, INTERF_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_interf_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP = 4;
    vars->ep_size_stresstyp = 5;

    vars->ep_value_thick = 0;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 6;
    vars->ep_value_length = 1;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("interf output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/

#ifdef D_WALLGE


/*----------------------------------------------------------------------*/
/*!
  \brief Variables that describe the structure of the wallge output.

  The only instance of this structure.

  \author u.kue
  \date 11/04
  \sa setup_wallge_variables
 */
/*----------------------------------------------------------------------*/
WALLGE_VARIABLES wallge_variables;


static void setup_wallge_variables(INT version, WALLGE_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_wallge_variables");
#endif

  vars->init = 1;
  vars->version = version;

  switch (version)
  {
  case 1:
    vars->ep_size_nGP0 = 4;
    vars->ep_size_nGP1 = 5;
    vars->ep_size_nGP2 = 6;
    vars->ep_size_nGP3 = 7;
    vars->ep_size_stresstyp = 8;

    vars->ep_value_thick = 0;

    /* the length of the integer and the double entry */
    vars->ep_size_length = 9;
    vars->ep_value_length = 1;
    break;

  case 0:

    /* Oh, no version. We are not initialized. */
    vars->init = 0;
    break;

  default:
    dserror("wallge output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

/*======================================================================*/
/*======================================================================*/


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
NODE_VARIABLES node_variables;


static void setup_node_variables(INT version, NODE_VARIABLES* vars)
{
#ifdef DEBUG
  dstrc_enter("setup_node_variables");
#endif

  vars->version = version;

  switch (version)
  {
  case 1:
    vars->coords_size_Id = 0;
    vars->coords_size_numele = 1;
    vars->coords_size_eleid = 2;
    break;
  default:
    dserror("node output version %d unknown", version);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

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
void setup_element_variables_current()
{
#ifdef DEBUG
  dstrc_enter("setup_element_variables_current");
#endif

  /* we do the nodes here, too */
  setup_node_variables(NODE_IO_VERSION, &node_variables);

  setup_element_variables(ELEMENT_IO_VERSION, &element_variables);

#ifdef D_SHELL8
  setup_shell8_variables(SHELL8_IO_VERSION, &shell8_variables);
#endif
#ifdef D_SHELL9
  setup_shell9_variables(SHELL9_IO_VERSION, &shell9_variables);
#endif
#ifdef D_BRICK1
  setup_brick1_variables(BRICK1_IO_VERSION, &brick1_variables);
#endif
#ifdef D_FLUID2
  setup_fluid2_variables(FLUID2_IO_VERSION, &fluid2_variables);
#endif
#ifdef D_FLUID2_PRO
  setup_fluid2_pro_variables(FLUID2_PRO_IO_VERSION, &fluid2_pro_variables);
#endif
#ifdef D_FLUID3_PRO
  setup_fluid3_pro_variables(FLUID3_PRO_IO_VERSION, &fluid3_pro_variables);
#endif
#ifdef D_FLUID2_IS
  setup_fluid2_is_variables(FLUID2_IS_IO_VERSION, &fluid2_is_variables);
#endif
#ifdef D_FLUID3_IS
  setup_fluid3_is_variables(FLUID3_IS_IO_VERSION, &fluid3_is_variables);
#endif
#ifdef D_FLUID2TU
  setup_fluid2tu_variables(FLUID2TU_IO_VERSION, &fluid2tu_variables);
#endif
#ifdef D_FLUID3
  setup_fluid3_variables(FLUID3_IO_VERSION, &fluid3_variables);
#endif
#ifdef D_FLUID3_F
  setup_fluid3_fast_variables(FLUID3_F_IO_VERSION, &fluid3_fast_variables);
#endif
#ifdef D_ALE
  setup_ale2_variables(ALE2_IO_VERSION, &ale2_variables);
#endif
#ifdef D_ALE
  setup_ale3_variables(ALE3_IO_VERSION, &ale3_variables);
#endif
#ifdef D_WALL1
  setup_wall1_variables(WALL1_IO_VERSION, &wall1_variables);
#endif
#ifdef D_BEAM3
  setup_beam3_variables(BEAM3_IO_VERSION, &beam3_variables);
#endif
#ifdef D_AXISHELL
  setup_axishell_variables(AXISHELL_IO_VERSION, &axishell_variables);
#endif
#ifdef D_INTERF
  setup_interf_variables(INTERF_IO_VERSION, &interf_variables);
#endif
#ifdef D_WALLGE
  setup_wallge_variables(WALLGE_IO_VERSION, &wallge_variables);
#endif


#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Set all element variables according to the version read from
  the map.

  \author u.kue
  \date 11/04
  \sa setup_element_variables
 */
/*----------------------------------------------------------------------*/
void setup_element_variables_map(MAP* group)
{
#ifdef DEBUG
  dstrc_enter("setup_element_variables_map");
#endif

  setup_node_variables(map_read_int(group, "node_version"),
                       &node_variables);

  setup_element_variables(map_read_int(group, "general_element_version"),
                          &element_variables);

#ifdef D_SHELL8
  if (map_symbol_count(group, "shell8_version") > 0)
  {
    setup_shell8_variables(map_read_int(group, "shell8_version"),
                           &shell8_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_shell8_variables(0, &shell8_variables);
  }
#endif
#ifdef D_SHELL9
  if (map_symbol_count(group, "shell9_version") > 0)
  {
    setup_shell9_variables(map_read_int(group, "shell9_version"),
                           &shell9_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_shell9_variables(0, &shell9_variables);
  }
#endif
#ifdef D_BRICK1
  if (map_symbol_count(group, "brick1_version") > 0)
  {
    setup_brick1_variables(map_read_int(group, "brick1_version"),
                           &brick1_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_brick1_variables(0, &brick1_variables);
  }
#endif
#ifdef D_FLUID2
  if (map_symbol_count(group, "fluid2_version") > 0)
  {
    setup_fluid2_variables(map_read_int(group, "fluid2_version"),
                           &fluid2_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid2_variables(0, &fluid2_variables);
  }
#endif
#ifdef D_FLUID2_PRO
  if (map_symbol_count(group, "fluid2_pro_version") > 0)
  {
    setup_fluid2_pro_variables(map_read_int(group, "fluid2_pro_version"),
                           &fluid2_pro_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid2_pro_variables(0, &fluid2_pro_variables);
  }
#endif
#ifdef D_FLUID3_PRO
  if (map_symbol_count(group, "fluid3_pro_version") > 0)
  {
    setup_fluid3_pro_variables(map_read_int(group, "fluid3_pro_version"),
                           &fluid3_pro_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid3_pro_variables(0, &fluid3_pro_variables);
  }
#endif
#ifdef D_FLUID2_IS
  if (map_symbol_count(group, "fluid2_is_version") > 0)
  {
    setup_fluid2_is_variables(map_read_int(group, "fluid2_is_version"),
                           &fluid2_is_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid2_is_variables(0, &fluid2_is_variables);
  }
#endif
#ifdef D_FLUID3_IS
  if (map_symbol_count(group, "fluid3_is_version") > 0)
  {
    setup_fluid3_is_variables(map_read_int(group, "fluid3_is_version"),
                           &fluid3_is_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid3_is_variables(0, &fluid3_is_variables);
  }
#endif
#ifdef D_FLUID2TU
  if (map_symbol_count(group, "fluid2tu_version") > 0)
  {
    setup_fluid2tu_variables(map_read_int(group, "fluid2tu_version"),
                           &fluid2tu_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid2tu_variables(0, &fluid2tu_variables);
  }
#endif
#ifdef D_FLUID3
  if (map_symbol_count(group, "fluid3_version") > 0)
  {
    setup_fluid3_variables(map_read_int(group, "fluid3_version"),
                           &fluid3_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid3_variables(0, &fluid3_variables);
  }
#endif
#ifdef D_FLUID3_F
  if (map_symbol_count(group, "fluid3_fast_version") > 0)
  {
    setup_fluid3_fast_variables(map_read_int(group, "fluid3_fast_version"),
                           &fluid3_fast_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_fluid3_fast_variables(0, &fluid3_fast_variables);
  }
#endif
#ifdef D_ALE
  if (map_symbol_count(group, "ale2_version") > 0)
  {
    setup_ale2_variables(map_read_int(group, "ale2_version"),
                           &ale2_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_ale2_variables(0, &ale2_variables);
  }
#endif
#ifdef D_ALE
  if (map_symbol_count(group, "ale3_version") > 0)
  {
    setup_ale3_variables(map_read_int(group, "ale3_version"),
                           &ale3_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_ale3_variables(0, &ale3_variables);
  }
#endif
#ifdef D_WALL1
  if (map_symbol_count(group, "wall1_version") > 0)
  {
    setup_wall1_variables(map_read_int(group, "wall1_version"),
                           &wall1_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_wall1_variables(0, &wall1_variables);
  }
#endif
#ifdef D_BEAM3
  if (map_symbol_count(group, "beam3_version") > 0)
  {
    setup_beam3_variables(map_read_int(group, "beam3_version"),
                           &beam3_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_beam3_variables(0, &beam3_variables);
  }
#endif
#ifdef D_AXISHELL
  if (map_symbol_count(group, "axishell_version") > 0)
  {
    setup_axishell_variables(map_read_int(group, "axishell_version"),
                           &axishell_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_axishell_variables(0, &axishell_variables);
  }
#endif
#ifdef D_INTERF
  if (map_symbol_count(group, "interf_version") > 0)
  {
    setup_interf_variables(map_read_int(group, "interf_version"),
                           &interf_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_interf_variables(0, &interf_variables);
  }
#endif
#ifdef D_WALLGE
  if (map_symbol_count(group, "wallge_version") > 0)
  {
    setup_wallge_variables(map_read_int(group, "wallge_version"),
                           &wallge_variables);
  }
  else
  {
    /* maybe the reader known more elements than ccarat when it wrote
     * the files */
    setup_wallge_variables(0, &wallge_variables);
  }
#endif


#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Print the current io versions of all available elements.

  \author u.kue
  \date 11/04
 */
/*----------------------------------------------------------------------*/
void out_print_element_versions(FILE* file)
{
#ifdef DEBUG
  dstrc_enter("out_print_element_versions");
#endif

  fprintf(file, "# ccarat element versions used here\n");
  fprintf(file, "node_version = %d\n", NODE_IO_VERSION);
  fprintf(file, "general_element_version = %d\n", ELEMENT_IO_VERSION);
#ifdef D_SHELL8
  fprintf(file, "shell8_version = %d\n", SHELL8_IO_VERSION);
#endif
#ifdef D_SHELL9
  fprintf(file, "shell9_version = %d\n", SHELL9_IO_VERSION);
#endif
#ifdef D_BRICK1
  fprintf(file, "brick1_version = %d\n", BRICK1_IO_VERSION);
#endif
#ifdef D_FLUID2
  fprintf(file, "fluid2_version = %d\n", FLUID2_IO_VERSION);
#endif
#ifdef D_FLUID2_PRO
  fprintf(file, "fluid2_pro_version = %d\n", FLUID2_PRO_IO_VERSION);
#endif
#ifdef D_FLUID3_PRO
  fprintf(file, "fluid2_pro_version = %d\n", FLUID3_PRO_IO_VERSION);
#endif
#ifdef D_FLUID2_IS
  fprintf(file, "fluid2_is_version = %d\n", FLUID2_IS_IO_VERSION);
#endif
#ifdef D_FLUID3_IS
  fprintf(file, "fluid2_is_version = %d\n", FLUID3_IS_IO_VERSION);
#endif
#ifdef D_FLUID2TU
  fprintf(file, "fluid2tu_version = %d\n", FLUID2TU_IO_VERSION);
#endif
#ifdef D_FLUID3
  fprintf(file, "fluid3_version = %d\n", FLUID3_IO_VERSION);
#endif
#ifdef D_FLUID3_F
  fprintf(file, "fluid3_fast_version = %d\n", FLUID3_F_IO_VERSION);
#endif
#ifdef D_ALE
  fprintf(file, "ale2_version = %d\n", ALE2_IO_VERSION);
#endif
#ifdef D_ALE
  fprintf(file, "ale3_version = %d\n", ALE3_IO_VERSION);
#endif
#ifdef D_WALL1
  fprintf(file, "wall1_version = %d\n", WALL1_IO_VERSION);
#endif
#ifdef D_BEAM3
  fprintf(file, "beam3_version = %d\n", BEAM3_IO_VERSION);
#endif
#ifdef D_AXISHELL
  fprintf(file, "axishell_version = %d\n", AXISHELL_IO_VERSION);
#endif
#ifdef D_INTERF
  fprintf(file, "interf_version = %d\n", INTERF_IO_VERSION);
#endif
#ifdef D_WALLGE
  fprintf(file, "wallge_version = %d\n", WALLGE_IO_VERSION);
#endif
  fprintf(file, "\n");

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
