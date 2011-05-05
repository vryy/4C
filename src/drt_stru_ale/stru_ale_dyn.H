/*!----------------------------------------------------------------------
\file stru_ale_dyn.H
\brief Control routine for structure with ale problems.


<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef STRU_ALE_DYN_H
#define STRU_ALE_DYN_H

#ifdef __cplusplus
extern "C"
{
#endif
/* headers and prototypes of functions also accessible from c */

/* entry point for the solution of structural ale problems */
void stru_ale_dyn_drt(int disnumsf,int disnumaf,int restart);

#ifdef __cplusplus
}
#endif

/*----------------------------------------------------------------------*/
#endif /* #ifndef STRU_ALE_DYN_H */
#endif /* #ifdef CCADISCRET */
