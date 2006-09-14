/*!
\file
\brief Postprocessing utility that shows the results using visual2.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

For visual we indeed do load all the results at once.

\author u.kue
\date 10/04

*/

#ifndef POST_OUT_H
#define POST_OUT_H

/* Include these first. */

#include "../post_common/post_common.h"
#include "../post_common/post_design.h"
#include "../post_common/post_octtree.h"

/* The following stems from the visual3 header file. */

#ifdef mips

#define PREFIX extern void
#define UNDERSCORE
#define MAINPROG MAIN__
#ifdef MIPSEL
#define BYTE_SWAP 1
#else
#define BYTE_SWAP 0
#endif /* MIPSEL */

#endif /* SGI & DECStation */

#ifdef __mips

#define PREFIX extern void
#define UNDERSCORE
#define MAINPROG MAIN__
#define BYTE_SWAP 0

#endif /* SGI */


#ifdef __MACH__

#define PREFIX extern void
#define UNDERSCORE
#define MAINPROG MAIN__
#define BYTE_SWAP 0

#endif /* MAC OS X */


#ifdef __alpha

#define PREFIX extern void
#define UNDERSCORE
#define MAINPROG MAIN__
#define BYTE_SWAP 1

#endif /* DEC Alpha */


#ifdef __linux

#define PREFIX extern void
#ifdef PGF
#define UNDERSCORE
#define MAINPROG MAIN_
#define BYTE_SWAP 1
#else
#ifdef INTEL
#define UNDERSCORE
#define MAINPROG main
#define BYTE_SWAP 1
#else
#ifdef IFORT
#define UNDERSCORE
#define MAINPROG MAIN__
#define BYTE_SWAP 1
#else
#ifdef ABS
#define UNDERSCORE
#define MAINPROG main
#define BYTE_SWAP 1
#else   /* assume g77 */
#define UNDERSCORE
#define MAINPROG MAIN__
#define BYTE_SWAP 1
#endif
#endif
#endif
#endif

#endif /* x86 Linux */


#ifdef __hpux

#define PREFIX extern void
#undef  UNDERSCORE
#define MAINPROG main
#define BYTE_SWAP 0

#endif /* HP 9000/700 */


#ifdef _IBMR2

#define PREFIX extern void
#undef  UNDERSCORE
#define MAINPROG main
#define BYTE_SWAP 0

#endif /* IBM RS6000 */


#ifdef __sparc

#define PREFIX extern void
#define UNDERSCORE
#define MAINPROG main
#define BYTE_SWAP 0

#endif /* SUN */


#ifdef WIN32 /* windows */

#ifdef ABSOFT
#undef WIN32
#define PREFIX extern void __cdecl
#undef  UNDERSCORE
#define MAINPROG main
#define BYTE_SWAP 1
#else
#define PREFIX extern void __stdcall
#define MAINPROG MAIN__
#define BYTE_SWAP 1
#endif

#endif /* windows */


#ifndef WIN32
#ifdef UNDERSCORE

#if defined(__linux) || defined(__MACH__)

#if defined(PGF) || defined(ABS) || defined(INTEL) || defined(IFORT)
#define V3_INIT       v3_init__
#define V3_CURSOR     v3_cursor_
#define V3_STAT       v3_stat_
#define V3_OBJECT3D   v3_object3d_
#define V3_OBJECT2D   v3_object2d_
#define V3_LINE       v3_line_
#define V3_GETSTATE   v3_getstate_
#define V3_SETSTATE   v3_setstate_
#define V3_GETSTRUC   v3_getstruc_
#define V3_INT        v3_int_
#define V3_FILL2DSCAL v3_fill2dscal_
#define V3_FILL2DVECT v3_fill2dvect_
#else
#define V3_INIT       v3_init__
#define V3_CURSOR     v3_cursor__
#define V3_STAT       v3_stat__
#define V3_OBJECT3D   v3_object3d__
#define V3_OBJECT2D   v3_object2d__
#define V3_LINE       v3_line__
#define V3_GETSTATE   v3_getstate__
#define V3_SETSTATE   v3_setstate__
#define V3_GETSTRUC   v3_getstruc__
#define V3_INT        v3_int__
#define V3_FILL2DSCAL v3_fill2dscal__
#define V3_FILL2DVECT v3_fill2dvect__
#endif

#else

#define V3_INIT       v3_init__
#define V3_CURSOR     v3_cursor_
#define V3_STAT       v3_stat_
#define V3_OBJECT3D   v3_object3d_
#define V3_OBJECT2D   v3_object2d_
#define V3_LINE       v3_line_
#define V3_GETSTATE   v3_getstate_
#define V3_SETSTATE   v3_setstate_
#define V3_GETSTRUC   v3_getstruc_
#define V3_INT        v3_int_
#define V3_FILL2DSCAL v3_fill2dscal_
#define V3_FILL2DVECT v3_fill2dvect_

#endif /* x86 Linux */

#define V3CELL      v3cell_
#define V3SURFACE   v3surface_
#define V3GRID      v3grid_
#define V3SCAL      v3scal_
#define V3VECT      v3vect_
#define V3THRES     v3thres_
#define V3EQUIV     v3equiv_
#define V3UPDATE    v3update_
#define V3ZPRIME    v3zprime_
#define V3XYPRIME   v3xyprime_
#define V3SURF      v3surf_
#define V3XYSURF    v3xysurf_
#define V3SSURF     v3ssurf_
#define V3VSURF     v3vsurf_
#define V3EVENT     v3event_
#define V3EVENTHELP v3eventhelp_
#define V3STRING    v3string_
#define V3PEVENTS   v3pevents_
#define V3DRAW3D    v3draw3d_
#define V3DRAW2D    v3draw2d_
#define V3PROBE     v3probe_
#define V3INIT      v3init_
#define V3DIALOG    v3dialog_
#define XFTNSTRING  xftnstring_

#else

#define V3_INIT     v3_init
#define V3_CURSOR   v3_cursor
#define V3_STAT     v3_stat
#define V3CELL      v3cell
#define V3SURFACE   v3surface
#define V3GRID      v3grid
#define V3SCAL      v3scal
#define V3VECT      v3vect
#define V3THRES     v3thres
#define V3EQUIV     v3equiv
#define V3UPDATE    v3update
#define V3ZPRIME    v3zprime
#define V3XYPRIME   v3xyprime
#define V3SURF      v3surf
#define V3XYSURF    v3xysurf
#define V3SSURF     v3ssurf
#define V3VSURF     v3vsurf
#define V3EVENT     v3event
#define V3EVENTHELP v3eventhelp
#define V3STRING    v3string
#define V3PEVENTS   v3pevents
#define V3_GETSTATE v3_getstate
#define V3_SETSTATE v3_setstate
#define V3_GETSTRUC v3_getstruc
#define V3DRAW3D    v3draw3d
#define V3_OBJECT3D v3_object3d
#define V3DRAW2D    v3draw2d
#define V3_OBJECT2D v3_object2d
#define V3PROBE     v3probe
#define V3_LINE     v3_line
#define V3INIT      v3init
#define V3_INT      v3_int
#define V3DIALOG    v3dialog
#define XFTNSTRING  xftnstring

#endif
#endif

#ifdef __ProtoGlarp__
#undef __ProtoGlarp__
#endif
#if defined(__STDC__) || defined(__cplusplus)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
PREFIX	V3_INIT __ProtoGlarp__(( char *title, int title_len, int *iopt,
                                 char *cmfile, int cmfile_len, int *unit,
                                 int *win3d, int *nkeys, int *ikeys,
                                 char *tkeys, int tkeys_len, float *fkeys,
                                 float *flims, int *mirror, int *knode,
                                 int *kequiv, int *kcel1, int *kcel2,
                                 int *kcel3, int *kcel4, int *knptet,
                                 int *kptet, int *knblock, int *blocks,
                                 int *ksurf, int *knsurf ));
#else
PREFIX	V3_INIT __ProtoGlarp__(( char *title, int *iopt, char *cmfile,
				 int *unit, int *win3d, int *nkeys, int *ikeys,
				 char *tkeys, int *fkeys, float *flims,
				 int *mirror, int *knode, int *kequiv,
				 int *kcel1, int *kcel2, int *kcel3,
				 int *kcel4, int *knptet, int *kptet,
				 int *knblock, int *blocks, int *ksurf,
				 int *knsurf, int title_len, int cmfile_len,
				 int tkeys_len ));
#endif
PREFIX	V3_STAT  __ProtoGlarp__(( int *istat, int *ostat, float *time ));
PREFIX	V3_CURSOR __ProtoGlarp__(( int *flag ));
PREFIX	V3_GETSTATE __ProtoGlarp__(( int *opt, long *ivec,  float *rvec,
				     char *str, int str_len ));
PREFIX	V3_SETSTATE __ProtoGlarp__(( int *opt, long *ivec,  float *rvec,
				     char *str, int str_len ));
PREFIX	V3_GETSTRUC __ProtoGlarp__(( int *opt, void **ptr, int *len ));
PREFIX	V3_OBJECT3D __ProtoGlarp__(( int *itype, int *icolor, float *xyzp,
				     float *radii, float *col, int *nc,
				     int *np));
PREFIX	V3_OBJECT2D __ProtoGlarp__(( int *itype, int *icolor, float *xyw,
				     float *col, int *np));
PREFIX	V3_LINE __ProtoGlarp__(( float *x, float *y, int *np));
PREFIX	V3_INT __ProtoGlarp__(( int *kc, float *xyz, int *nnodes, int *nodes,
	                        float *weights, int *nsufg ));
PREFIX	V3_FILL2DSCAL __ProtoGlarp__(( float *s ));
PREFIX	V3_FILL2DVECT __ProtoGlarp__(( float *v ));
PREFIX	XFTNSTRING __ProtoGlarp__(( long *iwin, int *type, float *color,
                                    float *pos, char *str, int str_len ));

#ifdef __cplusplus
}
#endif


#endif

