/*!---------------------------------------------------------------------
\file
\brief 

---------------------------------------------------------------------*/

/*--------------------------------------- check plausibility of defines */
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
#error "s*CONTACT needs D_CONTACT."
#error "Use -DD_CONTACT in the makefile."
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
/*----------------------------------------------------------------------*/
