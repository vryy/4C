/*!----------------------------------------------------------------------
\file
\brief chimera.h

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_CHIMERA

#include "../fluid_full/fluid_prototypes.h"

/*!
\addtogroup CHIMERA
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief struct CHIMERA_DATA

<pre>                                                            irhan 09/04
This structure contains data related to chimera module
</pre>

*----------------------------------------------------------------------*/
typedef struct _CHIMERA_DATA
{
  enum _CHIMERA_DEFINE_BACKGROUND_BOUNDARY chimera_define_back_bound;
  INT                            overlapping_length;
  INT                            n_Hintergrunddiskretisierung;
  INT                            nmax_Gebietsiteration;
  INT                            Konvergenz_Gebietsiteration;
  DOUBLE                         rel_TOL;
  DOUBLE                         abs_TOL;
  DOUBLE                        *rel_err;
  DOUBLE                        *abs_err;
  INT                            search_action;       /* 0 brute-force       1 quadtree*/
  INT                            conti_interpolation; /* 1 mit               0 ohne    */
  DOUBLE                         search_TOL;
  DOUBLE                         quadtree_TOL;
  INT                            nmax_pro_Blatt_in_quadtree;
  struct quadtree_node_struct  **quad_root;
  INT                            write_out_every;

  INT                            interact_dis_on_off;
} CHIMERA_DATA;



/*!----------------------------------------------------------------------
\brief quadtree_node_struct

<pre>                                                            irhan 09/04
This structure is used for quadtree search
</pre>

*----------------------------------------------------------------------*/
typedef struct quadtree_node_struct
{
  int      n_enthaltene_Elemente;
  double **node_Ecken;
  int     *enthaltene_Elemente;
  struct quadtree_node_struct *Nachfolger_lu;
  struct quadtree_node_struct *Nachfolger_ru;
  struct quadtree_node_struct *Nachfolger_ro;
  struct quadtree_node_struct *Nachfolger_lo;
} quadtree_node_type;

#include "chimera_prototypes.h"

/*! @} (documentation module close)*/
#endif
