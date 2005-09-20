/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

/*--------------------------------------- check plausibility of defines */
/* spooles only in parallel */
#if defined(SPOOLES_PACKAGE) && !defined(PARALLEL)
#error "SPOOLES only possible in  PARALLEL."
#error "Use -DPARALLEL in the makefile."
#endif

/* wallcontact */
#if defined(WALLCONTACT) && !defined(D_CONTACT)
#error "WALLCONTACT needs D_CONTACT."
#error "Use -DD_CONTACT in the makefile."
#endif

#if defined(S8CONTACT) && defined(WALLCONTACT)
#error "WALLCONTACT may not come together with S8CONTACT"
#endif

/* check shell contact definitions */
#if defined(S8CONTACT) && !defined(D_CONTACT)
#error "S8CONTACT needs D_CONTACT."
#error "Use -DD_CONTACT in the makefile."
#endif

/* no s8_hex8 without shell8 */
#if defined(S8_HEX8) && !defined(D_SHELL8)
#error "S8_HEX8 needs D_SHELL8."
#error "Use -DD_SHELL8 in the makefile."
#endif

/* no shell9 without mat */
#if defined(D_SHELL9) && !defined(D_MAT)
#error "D_SHELL9 needs D_MAT."
#error "Use -DD_MAT in the makefile."
#endif

/* no fsi without ale */
#if defined(D_FSI) && !defined(D_ALE)
#error "D_FSI needs D_ALE."
#error "Use -DD_ALE in the makefile."
#endif

/* no fsi without fluid */
#if defined(D_FSI) && !defined(D_FLUID)
#error "D_FSI needs D_FLUID."
#error "Use -DD_FLUID in the makefile."
#endif

/* no fluid2 without fluid */
#if defined(D_FLUID2) && !defined(D_FLUID)
#error "D_FLUID2 needs D_FLUID."
#error "Use -DD_FLUID in the makefile."
#endif

/* no fluid3 without fluid */
#if defined(D_FLUID3) && !defined(D_FLUID)
#error "D_FLUID3 needs D_FLUID."
#error "Use -DD_FLUID in the makefile."
#endif

/* no fluid2pro without fluid */
#if defined(D_FLUID2_PRO) && !defined(D_FLUID)
#error "D_FLUID2_PRO needs D_FLUID."
#error "Use -DD_FLUID in the makefile."
#endif

/* no fsi with gemm */
#if defined(D_FSI) && defined(GEMM)
#error "FSI not possible with GEMM."
#endif

/* no fluid2_ml without fluid2 */
#if defined(FLUID2_ML) && !defined(D_FLUID2)
#error "FLUID2_ML not possible without FLUID2."
#endif

/* no fluid3_ml without fluid3 */
#if defined(FLUID3_ML) && !defined(D_FLUID3)
#error "FLUID3_ML not possible without FLUID3."
#endif

/* no Visual2 with Visual3 */
#if defined(VISUAL2_PACKAGE) && defined(VISUAL3_PACKAGE)
#error "VISUAL2 not possible with VISUAL3."
#endif


/* solver and assembling */
/* fast_ass together with fast_ass2 */
#if defined(FAST_ASS) && defined(FAST_ASS2)
#error "FAST_ASS and FAST_ASS2 are not possible together."
#endif



/* check the plausability of some defines */
/*
#if MAXNOD < MAXNOD_AXISHELL
#error "MAXNOD must be larger then or equal to MAXNOD_AXISHELL!!"
#endif

#if MAXNOD < MAXNOD_BEAM3
#error "MAXNOD must be larger then or equal to MAXNOD_BEAM3!!"
#endif

#if MAXNOD < MAXNOD_BRICK1
#error "MAXNOD must be larger then or equal to MAXNOD_BRICK1!!"
#endif

#if MAXNOD < MAXNOD_F2
#error "MAXNOD must be larger then or equal to MAXNOD_F2!!"
#endif

#if MAXNOD < MAXNOD_F3
#error "MAXNOD must be larger then or equal to MAXNOD_F3!!"
#endif

#if MAXNOD < MAXNOD_ALE3
#error "MAXNOD must be larger then or equal to MAXNOD_ALE3!!"
#endif

#if MAXNOD < MAXNOD_SHELL8
#error "MAXNOD must be larger then or equal to MAXNOD_SHELL8!!"
#endif

#if MAXNOD < MAXNOD_SHELL9
#error "MAXNOD must be larger then or equal to MAXNOD_SHELL9!!"
#endif

#if MAXNOD < MAXNOD_WALL1
#error "MAXNOD must be larger then or equal to MAXNOD_WALL1!!"
#endif
*/


/*----------------------------------------------------------------------*/
