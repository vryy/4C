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
the ELEMENT_TYP enum. The second or minor number identifies the number
of nodes and Gauß points of the element. This scheme is supposed to be
easily adaptable. For instance you can add another minor version
without affecting any other element. That's important because these
numbers are written to result files and might life there for a while.

\author u.kue
\date 09/04

*/

#ifndef ELEMENT_INFO_H
#define ELEMENT_INFO_H


/* Our elements minor type numbers. Dependent on the number of nodes
 * and number of Gauß points of the element. Here we define the
 * constants as well as macros to determine them.
 *
 * These names are based on the GIDSET structure that's already
 * there. However it is not sufficient. We really need to distinguish
 * elements with different numbers of nodes or Gauß points. */

#define MINOR_SHELL1(actele)                    \
  (dserror("shell1 no longer supported"), -1)

#define MINOR_SHELL8_22  0      /* 4-noded shell8 2x2 GP */
#define MINOR_SHELL8_8_33  1    /* 8-noded shell8 3x3 GP */
#define MINOR_SHELL8_9_33  2    /* 9-noded shell8 3x3 GP */

#ifdef D_SHELL8

#define MINOR_SHELL8(actele)                                    \
  ((actele->numnp==4) ?                                         \
   MINOR_SHELL8_22 :                                            \
   ((actele->numnp==8) ?                                        \
    MINOR_SHELL8_8_33 :                                         \
    ((actele->numnp==9) ?                                       \
     MINOR_SHELL8_9_33 :                                        \
     (dserror("unknown minor version for shell8"), -1))))

#else
#define MINOR_SHELL8(actele) (dserror("unknown shell8"), -1)
#endif


#define MINOR_SHELL9_4_22  0    /* 4-noded shell9 2x2 GP */
#define MINOR_SHELL9_4_33  1    /* 4-noded shell9 3x3 GP */
#define MINOR_SHELL9_8_22  2    /* 8-noded shell9 2x2 GP */
#define MINOR_SHELL9_8_33  3    /* 8-noded shell9 3x3 GP */
#define MINOR_SHELL9_9_22  4    /* 9-noded shell9 2x2 GP */
#define MINOR_SHELL9_9_33  5    /* 9-noded shell9 3x3 GP */

#ifdef D_SHELL9

#define MINOR_SHELL9(actele)                                    \
  ((actele->numnp==4) ?                                         \
   ((actele->e.s9->nGP[0]==2 && actele->e.s9->nGP[1]==2) ?      \
    MINOR_SHELL9_4_22 :                                         \
    ((actele->e.s9->nGP[0]==3 && actele->e.s9->nGP[1]==3) ?     \
     MINOR_SHELL9_4_33 :                                        \
     (dserror("unknown minor version for shell9"), -1))) :      \
   ((actele->numnp==8) ?                                        \
    ((actele->e.s9->nGP[0]==2 && actele->e.s9->nGP[1]==2) ?     \
     MINOR_SHELL9_8_22 :                                        \
     ((actele->e.s9->nGP[0]==3 && actele->e.s9->nGP[1]==3) ?    \
      MINOR_SHELL9_8_33 :                                       \
      (dserror("unknown minor version for shell9"), -1))) :     \
    ((actele->numnp==9) ?                                       \
     ((actele->e.s9->nGP[0]==2 && actele->e.s9->nGP[1]==2) ?    \
      MINOR_SHELL9_9_22 :                                       \
      ((actele->e.s9->nGP[0]==3 && actele->e.s9->nGP[1]==3) ?   \
       MINOR_SHELL9_9_33 :                                      \
       (dserror("unknown minor version for shell9"), -1))) :    \
     (dserror("unknown minor version for shell9"), -1))))

#else
#define MINOR_SHELL9(actele) (dserror("unknown shell9"), -1)
#endif


#define MINOR_BRICK1_222  0     /* 8-noded brick1 2x2x2 GP */
#define MINOR_BRICK1_20_333  1  /* 20 noded brick1 3x3x3 GP */
#define MINOR_BRICK1_27_333  2  /* 27 noded brick1 3x3x3 GP */

#ifdef D_BRICK1

#define MINOR_BRICK1(actele)                            \
  ((actele->numnp==8) ?                                 \
   MINOR_BRICK1_222 :                                   \
   ((actele->numnp==20) ?                               \
    MINOR_BRICK1_20_333 :                               \
    ((actele->numnp==27) ?                              \
     MINOR_BRICK1_27_333 :                              \
     (dserror("unknown minor version for brick1"), -1))))

#else
#define MINOR_BRICK1(actele) (dserror("unknown brick1"), -1)
#endif


#define MINOR_WALL1_11  0       /* 3-noded wall1 1x1 GP */
#define MINOR_WALL1_22  1       /* 4-noded wall1 2x2 GP */
#define MINOR_WALL1_8_33  2     /* 8-noded wall1 3x3 GP */
#define MINOR_WALL1_9_33  3     /* 9-noded wall1 3x3 GP */

/* a variant with insufficient gauss points */
#define MINOR_WALL1_8_22  4     /* 8-noded wall1 2x2 GP */

#ifdef D_WALL1

#define MINOR_WALL1(actele)                                     \
  (((actele->numnp==4) && (actele->e.w1->nGP[0]==1)) ?          \
   MINOR_WALL1_11 :                                             \
   (((actele->numnp==4) && (actele->e.w1->nGP[0]==2)) ?         \
    MINOR_WALL1_22 :                                            \
    (((actele->numnp==8) && (actele->e.w1->nGP[0]==3)) ?        \
     MINOR_WALL1_8_33 :                                         \
     (((actele->numnp==9) && (actele->e.w1->nGP[0]==3)) ?       \
      MINOR_WALL1_9_33 :                                        \
      (((actele->numnp==8) && (actele->e.w1->nGP[0]==2)) ?      \
       MINOR_WALL1_8_22 :                                       \
       (dserror("unknown minor version for wall1"), -1))))))

#else
#define MINOR_WALL1(actele) (dserror("unknown wall1"), -1)
#endif


#define MINOR_BEAM3_21  0       /* 2-noded beam3 1 GP */
#define MINOR_BEAM3_22  1       /* 2-noded beam3 2 GP */
#define MINOR_BEAM3_32  2       /* 3-noded beam3 2 GP */
#define MINOR_BEAM3_33  3       /* 3-noded beam3 3 GP */

#ifdef D_BEAM3

#define MINOR_BEAM3(actele)                                     \
  ((actele->numnp==2) ?                                         \
   ((actele->e.b3->nGP[0]==1) ?                                 \
    MINOR_BEAM3_21 :                                            \
    ((actele->e.b3->nGP[0]==2) ?                                \
     MINOR_BEAM3_22 :                                           \
     (dserror("unknown minor version for beam3"), -1))) :       \
   ((actele->numnp==3) ?                                        \
    ((actele->e.b3->nGP[0]==2) ?                                \
     MINOR_BEAM3_32 :                                           \
     ((actele->e.b3->nGP[0]==3) ?                               \
      MINOR_BEAM3_33 :                                          \
      (dserror("unknown minor version for beam3"), -1))) :      \
    (dserror("unknown minor version for beam3"), -1)))

#else
#define MINOR_BEAM3(actele) (dserror("unknown beam3"), -1)
#endif


/* preliminary, there are many more fluid elements... */
#define MINOR_FLUID2_22  0      /* 4-noded fluid2 2x2 GP */
#define MINOR_FLUID2_8_33  1    /* 8-noded fluid2 3x3 GP */
#define MINOR_FLUID2_9_33  2    /* 9-noded fluid2 3x3 GP */
#define MINOR_FLUID2_3_4  3     /* 3-noded fluid2 4 GP */

#ifdef D_FLUID2

#define MINOR_FLUID2(actele)                                    \
  ((actele->numnp==4) ?                                         \
   MINOR_FLUID2_22 :                                            \
   ((actele->numnp==8) ?                                        \
    MINOR_FLUID2_8_33 :                                         \
    ((actele->numnp==9) ?                                       \
     MINOR_FLUID2_9_33 :                                        \
     (((actele->numnp==3) && (actele->e.f2->nGP[0]==4)) ?       \
      MINOR_FLUID2_3_4 :                                        \
      (dserror("unknown minor version for fluid2"), -1)))))

#else
#define MINOR_FLUID2(actele) (dserror("unknown fluid2"), -1)
#endif


#define MINOR_FLUID2_PRO_22   0 /* 4-noded fluid2_pro 2x2 GP */
#define MINOR_FLUID2_PRO_8_33 1 /* 8-noded fluid2_pro 3x3 GP */
#define MINOR_FLUID2_PRO_9_33 2 /* 9-noded fluid2_pro 3x3 GP */

#ifdef D_FLUID2_PRO

#define MINOR_FLUID2_PRO(actele)                                \
  ((actele->numnp==4) ?                                         \
   MINOR_FLUID2_PRO_22 :                                        \
   ((actele->numnp==8) ?                                        \
    MINOR_FLUID2_PRO_8_33 :                                     \
    ((actele->numnp==9) ?                                       \
     MINOR_FLUID2_PRO_9_33 :                                    \
     (dserror("unknown minor version for fluid2_pro"), -1))))

#else
#define MINOR_FLUID2_PRO(actele) (dserror("unknown fluid2_pro"), -1)
#endif


#define MINOR_FLUID2_TU(actele)                 \
  (dserror("fluid2_tu not yet supported"), -1)

#define MINOR_FLUID3_222  0     /* 8-noded fluid3 2x2x2 GP */
#define MINOR_FLUID3_20_333  1  /* 20-noded fluid3 3x3x3 GP */
#define MINOR_FLUID3_27_333  2  /* 27-noded fluid3 3x3x3 GP */

#ifdef D_FLUID3

#define MINOR_FLUID3(actele)                                    \
  ((actele->numnp==8) ?                                         \
   MINOR_FLUID3_222 :                                           \
   ((actele->numnp==20) ?                                       \
    MINOR_FLUID3_20_333 :                                       \
    ((actele->numnp==27) ?                                      \
     MINOR_FLUID3_27_333 :                                      \
     (dserror("unknown minor version for fluid3"), -1))))

#else
#define MINOR_FLUID3(actele) (dserror("unknown fluid3"), -1)
#endif

#define MINOR_ALE_11  0         /* 4-noded ale 1x1 GP */
#define MINOR_ALE_22  1         /* 4-noded ale 2x2 GP */
#define MINOR_ALE_TRI_1  2      /* 3-noded tri ale 1 GP */
#define MINOR_ALE_TRI_3  3      /* 3-noded tri ale 3 GP */

#ifdef D_ALE

#define MINOR_ALE2(actele)                                              \
  ((actele->numnp==4) ?                                                 \
   ((actele->e.ale2->nGP[0]==1 && actele->e.ale2->nGP[1]==1 ) ?         \
    MINOR_ALE_11 :                                                      \
    ((actele->e.ale2->nGP[0]==2 && actele->e.ale2->nGP[1]==2 ) ?        \
     MINOR_ALE_22 :                                                     \
     (dserror("unknown minor version for ale2"), -1))) :                \
   ((actele->numnp==3) ?                                                \
    (( actele->e.ale2->nGP[0] == 1) ?                                   \
     MINOR_ALE_TRI_1 :                                                  \
     (( actele->e.ale2->nGP[0] == 3) ?                                  \
      MINOR_ALE_TRI_3 :                                                 \
      (dserror("unknown minor version for ale2"), -1))) :               \
    (dserror("unknown minor version for ale2"), -1)))

#else
#define MINOR_ALE2(actele) (dserror("unknown ale2"), -1)
#endif


#define MINOR_ALE_111  0        /* 8-noded ale 1x1x1 GP */
#define MINOR_ALE_222  1        /* 8-noded ale 2x2x2 GP */
#define MINOR_ALE_TET_1  2      /* 4-noded tet ale 1 GP */
#define MINOR_ALE_TET_4  3      /* 4-noded tet ale 4 GP */

#ifdef D_ALE

#define MINOR_ALE3(actele)                                              \
  ((actele->numnp==8) ?                                                 \
   ((actele->e.ale3->nGP[0]==1 &&                                       \
     actele->e.ale3->nGP[1]==1 &&                                       \
     actele->e.ale3->nGP[2]==1) ?                                       \
    MINOR_ALE_111 :                                                     \
    ((actele->e.ale3->nGP[0]==2 &&                                      \
      actele->e.ale3->nGP[1]==2 &&                                      \
      actele->e.ale3->nGP[2]==2) ?                                      \
     MINOR_ALE_222 :                                                    \
     (dserror("unknown minor version for ale3"), -1))) :                \
   ((actele->numnp==4) ?                                                \
    ((actele->e.ale3->nGP[0] == 1) ?                                    \
     MINOR_ALE_TET_1 :                                                  \
     ((actele->e.ale3->nGP[0] == 4) ?                                   \
      MINOR_ALE_TET_4 :                                                 \
      (dserror("unknown minor version for ale3"), -1))) :               \
    (dserror("unknown minor version for ale3"), -1)))

#else
#define MINOR_ALE3(actele) (dserror("unknown ale3"), -1)
#endif


#define MINOR_AXISHELL  0       /* 2-noded axishell */


#define MINOR_INTERF_22  0      /* interface 2x2 GP */
#define MINOR_INTERF_33  1      /* interface 3x3 GP */

#ifdef D_INTERF

#define MINOR_INTERF(actele)                            \
  ((actele->e.interf->nGP==2) ?                         \
   MINOR_INTERF_22 :                                    \
   ((actele->e.interf->nGP==3) ?                        \
    MINOR_INTERF_33 :                                   \
    (dserror("unknown minor version for interf"), -1)))

#else
#define MINOR_INTERF(actele) (dserror("unknown interf"), -1)
#endif


#define MINOR_WALLGE_22  0      /* gradient enhanced wall 2x2 GP */
#define MINOR_WALLGE_8_33  1    /* gradient enhanced wall 3x3 GP */
#define MINOR_WALLGE_9_33  2    /* gradient enhanced wall 3x3 GP */

#ifdef D_WALLGE

#define MINOR_WALLGE(actele)                                    \
  ((actele->e.wallge->nGP[0]==2) ?                              \
   MINOR_WALLGE_22 :                                            \
    (((actele->numnp==8) && (actele->e.wallge->nGP[0]==3)) ?    \
     MINOR_WALL1_8_33 :                                         \
     (((actele->numnp==9) && (actele->e.wallge->nGP[0]==3)) ?   \
      MINOR_WALL1_9_33 :                                        \
      (dserror("unknown minor version for wallge"), -1))))

#else
#define MINOR_WALLGE(actele) (dserror("unknown wallge"), -1)
#endif


#define MINOR_LS2_33  0         /* 3-noded ls 3 GP */
#define MINOR_LS2_22  1         /* 4-noded ls 2x2 GP */

#ifdef D_LS

#define MINOR_LS2(actele)                               \
  ((actele->e.ls2->nGP[0]==2) ?                         \
   MINOR_LS2_22 :                                       \
   ((actele->e.ls2->nGP[1]==3) ?                        \
    MINOR_LS2_33 :                                      \
    (dserror("unknown minor version for ls2"), -1)))

#else
#define MINOR_LS2(actele) (dserror("unknown minor version for ls2"), -1)
#endif


#define MINOR_FLUID2_XFEM_22  0 /* 4-noded fluid2_xfem 2x2 GP */
#define MINOR_FLUID2_XFEM_33  1 /* 3-noded fluid2_xfem 3 GP */

#ifdef D_XFEM

#define MINOR_FLUID2_XFEM(actele)                               \
  ((actele->numnp==4) ?                                         \
   MINOR_FLUID2_XFEM_22 :                                       \
   ((actele->numnp==3) ?                                        \
    MINOR_FLUID2_XFEM_33 :                                      \
    (dserror("unknown minor version for fluid2_xfem"), -1)))

#else
#define MINOR_FLUID2_XFEM(actele) (dserror("unknown minor version for fluid2_xfem"), -1)
#endif


/* The quick way to find the minor element version. */
#define FIND_MINOR(actele)                                              \
  ((actele->eltyp == el_shell1) ? MINOR_SHELL1(actele) :                \
   ((actele->eltyp == el_shell8) ? MINOR_SHELL8(actele) :               \
    ((actele->eltyp == el_shell9) ? MINOR_SHELL9(actele) :              \
     ((actele->eltyp == el_brick1) ? MINOR_BRICK1(actele) :             \
      ((actele->eltyp == el_wall1) ? MINOR_WALL1(actele) :              \
       ((actele->eltyp == el_beam3) ? MINOR_BEAM3(actele) :             \
        ((actele->eltyp == el_fluid2) ? MINOR_FLUID2(actele) :          \
         ((actele->eltyp == el_fluid2_pro) ? MINOR_FLUID2_PRO(actele) : \
          ((actele->eltyp == el_fluid2_tu) ? MINOR_FLUID2_TU(actele) :  \
           ((actele->eltyp == el_fluid3) ? MINOR_FLUID3(actele) :       \
            ((actele->eltyp == el_ale2) ? MINOR_ALE2(actele) :          \
             ((actele->eltyp == el_ale3) ? MINOR_ALE3(actele) :         \
              ((actele->eltyp == el_axishell) ? MINOR_AXISHELL :        \
               ((actele->eltyp == el_interf) ? MINOR_INTERF(actele) :   \
                ((actele->eltyp == el_wallge) ? MINOR_WALLGE(actele) :  \
                 ((actele->eltyp == el_ls2) ? MINOR_LS2(actele) :       \
                  ((actele->eltyp == el_fluid2_xfem) ? MINOR_FLUID2_XFEM(actele) : \
                   (dserror("unknown element type %d", actele->eltyp), -1))))))))))))))))))


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
  \brief post element types.

  These are the types of elements GiD knows. It's useful to have them
  here because we need to maintain the mapping information.

  Be careful! The values are written to GiD's binary file. For this
  reason this enum must not be changed.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef enum _GiD_ElementType {
  GiD_NoElement = 0,
  GiD_Point,
  GiD_Linear,
  GiD_Triangle,
  GiD_Quadrilateral,
  GiD_Tetrahedra,
  GiD_Hexahedra
} GiD_ElementType;


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
   * Variant information. The entries in this group do not necessarily
   * tell what the element does internally but instead what the
   * element pretends to do. For example the four noded shell8 element
   * used to pretend to be hexahedra. This way the thickness of the
   * element was shown in GiD.
   *
   * The clean way is, of course, to give the real values here and not
   * to lie to the output routines. It's up to the filters to tell
   * lies to the postprocessor. But it's going to take time to make
   * this transition. */
  struct {

    /* The number of nodes */
    INT node_number;

    /* The geometry type of this element according to GiD */
    GiD_ElementType gid_type;

    /* The number of gauss points */
    INT gauss_number;

    /*
     * The concept of stress is something common to most elements. But
     * there are differences in size.
     *
     * We could avoid to store these numbers if stress arrays were
     * always of the size dimension times dof times gauss point number
     * (or symmetric part of it). But often stresses are extrapolated
     * to the node. In this case stress is a nodal value. We need to
     * know both, the element and the node based stress size. (I'd
     * like to remove the node stresses. It's due to the filter to
     * extrapolate the stresses if needed.)
     *
     * Really weird elements like shell9 have variable sized dofs. In
     * this case there's just a -1 here and special care needs to be
     * taken. */
    INT el_stress_matrix_size;
    INT nd_stress_matrix_size;

  } variant[MAX_EL_MINOR];

} ELEMENT_INFO;


/*----------------------------------------------------------------------*/
/*!
  \brief The table with all information about the elements.

  Major and minor number specify an element type. Here we have the
  table that contains everything we'd might need to know about the
  element.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
extern ELEMENT_INFO element_info[el_count];

#endif
