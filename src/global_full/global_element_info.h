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

To identify the element type there are two integers used. The first
number, also called major number, is the element type as defined by
the ELEMENT_TYP enum, the common element type. The second or minor
number identifies the number of nodes and Gauﬂ points of the element,
for brevity this is called the element variant. This scheme is
supposed to be easily adaptable. For instance you can add another
variant (minor version) without affecting any other element.

The major and minor type numbers of an element are written to the
binary files along with the mesh connectivity. This enables readers of
ccarat's output to identify the elements. To be able to write
elements we need get their major and minor types. That's what the
function ``find_minor_type`` and its subfunctions are for.

Because we cannot guarantee that the major and minor element type
numbers will never change we have to write their meaning to the
control file, too. That's done by ``write_variant_group``. Accordingly
we need to be able to convert the information from the control file to
internal major and minor numbers. The function
``restore_element_type`` does that.

Thus at the very moment major and minor numbers are written to a file
they are implicitly converted from internal quantities to external
ones. These external ones are not guaranteed (in the long run) to
equal their internal counterparts, but instead we provide a way to
convert external quantities back to internal ones. The filters are
concerned with that.

\author u.kue
\date 09/04

*/

#ifndef ELEMENT_INFO_H
#define ELEMENT_INFO_H

#include "../pss_full/pss_table.h"


/* Our elements minor type numbers. Dependent on the number of nodes
 * and number of Gauﬂ points of the element. Here we define the
 * constants.
 *
 * These names are based on the GIDSET structure that's already
 * there. However it is not sufficient. We really need to distinguish
 * elements with different numbers of nodes or gauss points. There are
 * also variants that have the same number of gauss points but
 * distribute them in different ways. These need to be distinguished,
 * too. */

#define MINOR_SHELL8_22  0      /* 4-noded shell8 2x2 GP */
#define MINOR_SHELL8_8_33  1    /* 8-noded shell8 3x3 GP */
#define MINOR_SHELL8_9_33  2    /* 9-noded shell8 3x3 GP */

#define MINOR_SHELL9_4_22  0    /* 4-noded shell9 2x2 GP */
#define MINOR_SHELL9_4_33  1    /* 4-noded shell9 3x3 GP */
#define MINOR_SHELL9_8_22  2    /* 8-noded shell9 2x2 GP */
#define MINOR_SHELL9_8_33  3    /* 8-noded shell9 3x3 GP */
#define MINOR_SHELL9_9_22  4    /* 9-noded shell9 2x2 GP */
#define MINOR_SHELL9_9_33  5    /* 9-noded shell9 3x3 GP */

#define MINOR_BRICK1_222  0     /* 8-noded brick1 2x2x2 GP */
#define MINOR_BRICK1_20_333  1  /* 20 noded brick1 3x3x3 GP */
#define MINOR_BRICK1_27_333  2  /* 27 noded brick1 3x3x3 GP */

#define MINOR_WALL1_11  0       /* 3-noded wall1 1x1 GP */
#define MINOR_WALL1_22  1       /* 4-noded wall1 2x2 GP */
#define MINOR_WALL1_8_33  2     /* 8-noded wall1 3x3 GP */
#define MINOR_WALL1_9_33  3     /* 9-noded wall1 3x3 GP */

/* a variant with insufficient gauss points */
#define MINOR_WALL1_8_22  4     /* 8-noded wall1 2x2 GP */

#define MINOR_BEAM3_21  0       /* 2-noded beam3 1 GP */
#define MINOR_BEAM3_22  1       /* 2-noded beam3 2 GP */
#define MINOR_BEAM3_32  2       /* 3-noded beam3 2 GP */
#define MINOR_BEAM3_33  3       /* 3-noded beam3 3 GP */

/* preliminary, there are many more fluid elements... */
#define MINOR_FLUID2_22  0      /* 4-noded fluid2 2x2 GP */
#define MINOR_FLUID2_8_33  1    /* 8-noded fluid2 3x3 GP */
#define MINOR_FLUID2_9_33  2    /* 9-noded fluid2 3x3 GP */
#define MINOR_FLUID2_3_4  3     /* 3-noded fluid2 4 GP */

#define MINOR_FLUID2_PRO_22   0 /* 4-noded fluid2_pro 2x2 GP */
#define MINOR_FLUID2_PRO_8_33 1 /* 8-noded fluid2_pro 3x3 GP */
#define MINOR_FLUID2_PRO_9_33 2 /* 9-noded fluid2_pro 3x3 GP */

#define MINOR_FLUID3_222  0     /* 8-noded fluid3 2x2x2 GP */
#define MINOR_FLUID3_20_333  1  /* 20-noded fluid3 3x3x3 GP */
#define MINOR_FLUID3_27_333  2  /* 27-noded fluid3 3x3x3 GP */

#define MINOR_FLUID3_FAST_222  0 /* 8-noded fluid3 2x2x2 GP */

#define MINOR_ALE_11  0         /* 4-noded ale 1x1 GP */
#define MINOR_ALE_22  1         /* 4-noded ale 2x2 GP */
#define MINOR_ALE_TRI_1  2      /* 3-noded tri ale 1 GP */
#define MINOR_ALE_TRI_3  3      /* 3-noded tri ale 3 GP */

#define MINOR_ALE_111  0        /* 8-noded ale 1x1x1 GP */
#define MINOR_ALE_222  1        /* 8-noded ale 2x2x2 GP */
#define MINOR_ALE_TET_1  2      /* 4-noded tet ale 1 GP */
#define MINOR_ALE_TET_4  3      /* 4-noded tet ale 4 GP */

#define MINOR_AXISHELL  0       /* 2-noded axishell */

#define MINOR_INTERF_22  0      /* interface 2x2 GP */
#define MINOR_INTERF_33  1      /* interface 3x3 GP */

#define MINOR_WALLGE_22  0      /* gradient enhanced wall 2x2 GP */
#define MINOR_WALLGE_8_33  1    /* gradient enhanced wall 3x3 GP */
#define MINOR_WALLGE_9_33  2    /* gradient enhanced wall 3x3 GP */


/* The maximum number of minor versions to an element. Keep this in
 * sync with the list above. */

#define MAX_EL_MINOR 6

#define ELEMENT_FLAGS_SIZE (el_count*MAX_EL_MINOR*sizeof(INT))


/*----------------------------------------------------------------------*/
/*!
  \brief The element flags.

  This is used to store information that depends both on major and
  minor number. It's mainly used to indicate what types of elements
  are used in a mesh.

  This used to be an array of chars but the HPUX version of
  MPI_Allreduce with MPI_MAX failed to handle that.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef INT ELEMENT_FLAGS[el_count][MAX_EL_MINOR];


/*----------------------------------------------------------------------*/
/*!
  \brief The structure of the information about an element type.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef struct _ELEMENT_INFO {

  /* The element's name. Depends only on the major number. */
  CHAR* name;

  /*
   * Variant information. */
  struct {

    /* The number of nodes */
    INT node_number;

    /* how the nodes are arranged */
    DIS_TYP dis_type;

    /* The number of gauss points */
    INT gauss_number;

    /*
     * The concept of stress is something common to most elements. But
     * there are differences in size.
     *
     * Really weird elements like shell9 have variable numbers of
     * dofs. In this case there's just a -1 here and special care
     * needs to be taken. */
    INT stress_matrix_size;

  } variant[MAX_EL_MINOR];

} ELEMENT_INFO;


/*----------------------------------------------------------------------*/
/*!
  \brief The table with all information about the elements.

  Major and minor number specify an element type. Here we have the
  table that contains everything we might need to know about the
  element.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
extern ELEMENT_INFO element_info[el_count];



#ifdef D_SHELL8

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
INT find_shell8_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 7 parameter shell element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_shell8_minor_type(MAP* group);


#endif
#ifdef D_SHELL9

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
INT find_shell9_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a multi layer shell element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_shell9_minor_type(MAP* group);


#endif
#ifdef D_BRICK1

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
INT find_brick1_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a structural brick element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_brick1_minor_type(MAP* group);


#endif
#ifdef D_WALL1

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
INT find_wall1_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D plane stress - plane
  strain element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_wall1_minor_type(MAP* group);


#endif
#ifdef D_BEAM3

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
INT find_beam3_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a structural 3D-beam element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_beam3_minor_type(MAP* group);


#endif
#ifdef D_FLUID2

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D fluid element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid2_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D fluid element's variant
  group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid2_minor_type(MAP* group);


#endif
#ifdef D_FLUID2_PRO

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 2D fluid pro element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid2_pro_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D fluid pro element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid2_pro_minor_type(MAP* group);


#endif
#ifdef D_FLUID2TU

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
INT find_fluid2_tu_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D fluid element for
  turbulence's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid2_tu_minor_type(MAP* group);


#endif
#ifdef D_FLUID3

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 3D fluid element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid3_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 3D fluid element's variant
  group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid3_minor_type(MAP* group);


#endif
#ifdef D_FLUID3_F

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 3D fluid element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 11/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_fluid3_fast_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 3D fluid element's variant
  group.

  \author u.kue
  \date 11/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_fluid3_fast_minor_type(MAP* group);


#endif
#ifdef D_ALE

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
INT find_ale2_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 2D pseudo structural ale
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_ale2_minor_type(MAP* group);


#endif
#ifdef D_ALE

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
INT find_ale3_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 3D pseudo structural ale
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_ale3_minor_type(MAP* group);


#endif
#ifdef D_AXISHELL

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
INT find_axishell_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 1D axisymmetrical shell
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_axishell_minor_type(MAP* group);


#endif
#ifdef D_INTERF

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a given 1D interface element.

  \param actele    (i) the element those type is requested

  \author u.kue
  \date 10/04
  \sa find_minor_type
*/
/*----------------------------------------------------------------------*/
INT find_interf_minor_type(ELEMENT* actele);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a 1D interface element's
  variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_interf_minor_type(MAP* group);


#endif
#ifdef D_WALLGE

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
INT find_wallge_minor_type(ELEMENT* actele);

/*----------------------------------------------------------------------*/
/*!
  \brief Find the element minor type to a gradient enhanced wall
  element's variant group.

  \author u.kue
  \date 10/04
  \sa restore_element_type
*/
/*----------------------------------------------------------------------*/
INT restore_wallge_minor_type(MAP* group);


#endif


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
INT find_minor_type(ELEMENT* actele);




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
void restore_element_type(MAP* group, INT* major, INT* minor);



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
void write_variant_group(FILE* out, ELEMENT_TYP ele, INT variant);

#endif
