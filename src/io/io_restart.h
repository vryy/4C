/*!
\file
\brief Restart functions based on the binary io.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Here the restart functions for all algorithms are collected. These
functions are supposed to be rather high level and easy to
understand. Feel free to add your own.

\author u.kue
\date 08/04

*/

#ifdef BINIO

#ifndef BIN_RESTART_H
#define BIN_RESTART_H

#include "../headers/standardtypes.h"

#include "io_packing.h"
#include "io_singlefile.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Write structure restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_nlnstructdyn(struct _BIN_OUT_FIELD *context,
                                    STRUCT_DYNAMIC  *sdyn,
                                    STRUCT_DYN_CALC *dynvar,
                                    INT nrhs,  DIST_VECTOR *rhs,
                                    INT nsol,  DIST_VECTOR *sol,
                                    INT ndis,  DIST_VECTOR *dispi,
                                    INT nvel,  DIST_VECTOR *vel,
                                    INT nacc,  DIST_VECTOR *acc,
                                    INT nfie,  DIST_VECTOR *fie,
                                    INT nwork, DIST_VECTOR *work);


/*----------------------------------------------------------------------*/
/*!
  \brief Read structure restart data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_nlnstructdyn(STRUCT_DYNAMIC  *sdyn,
                                   STRUCT_DYN_CALC *dynvar,
                                   SPARSE_TYP      *sysarray_typ,
                                   SPARSE_ARRAY    *sysarray,
                                   FIELD           *actfield,
                                   PARTITION       *actpart,
                                   INT             disnum,
                                   INTRA           *actintra,
                                   INT nrhs,  DIST_VECTOR *rhs,
                                   INT nsol,  DIST_VECTOR *sol,
                                   INT ndis,  DIST_VECTOR *dispi,
                                   INT nvel,  DIST_VECTOR *vel,
                                   INT nacc,  DIST_VECTOR *acc,
                                   INT nfie,  DIST_VECTOR *fie,
                                   INT nwork, DIST_VECTOR *work,
                                   INT step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write static structure restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_nlnstructstat(struct _BIN_OUT_FIELD *context,
                                     STATIC_VAR      *statvar,
                                     STANLN          *nln_data,
                                     INT kstep,
                                     INT nrhs,  DIST_VECTOR *rhs,
                                     INT nsol,  DIST_VECTOR *sol,
                                     INT ndis,  DIST_VECTOR *dispi);


/*----------------------------------------------------------------------*/
/*!
  \brief Read static structure restart data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_nlnstructstat(STATIC_VAR      *statvar,
                                    STANLN          *nln_data,
                                    SPARSE_TYP      *sysarray_typ,
                                    SPARSE_ARRAY    *sysarray,
                                    FIELD           *actfield,
                                    PARTITION       *actpart,
                                    INT              disnum,
                                    INTRA           *actintra,
                                    INT nrhs,  DIST_VECTOR *rhs,
                                    INT nsol,  DIST_VECTOR *sol,
                                    INT ndis,  DIST_VECTOR *dispi,
                                    INT step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write fluid restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_fluiddyn(struct _BIN_OUT_FIELD   *context,
                                FLUID_DYNAMIC   *fdyn);


/*----------------------------------------------------------------------*/
/*!
  \brief Read fluid restart data.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_fluiddyn(FLUID_DYNAMIC   *fdyn,
                               SPARSE_TYP      *sysarray_typ,
                               SPARSE_ARRAY    *sysarray,
                               FIELD           *actfield,
                               PARTITION       *actpart,
                               INT              disnum,
                               INTRA           *actintra,
                               INT              step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write ale restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_aledyn(struct _BIN_OUT_FIELD *context,
                              ALE_DYNAMIC *adyn);


/*----------------------------------------------------------------------*/
/*!
  \brief Read ale restart data.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_aledyn(ALE_DYNAMIC     *adyn,
                             SPARSE_TYP      *sysarray_typ,
                             SPARSE_ARRAY    *sysarray,
                             FIELD	     *actfield,
                             PARTITION       *actpart,
                             INT              disnum,
                             INTRA	     *actintra,
                             INT              step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write fsi restart data.

  This is additionally to the individual fields data. So only very
  little is stored here. No field at all.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_fsidyn(FSI_DYNAMIC *fsidyn);


/*----------------------------------------------------------------------*/
/*!
  \brief Read fsi restart data.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_fsidyn(FSI_DYNAMIC *fsidyn, INT step);


#endif
#endif
