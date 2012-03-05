/*!----------------------------------------------------------------------
\file global_inp_control2.H
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef GLOBAL_INP_CONTROL2_H
#define GLOBAL_INP_CONTROL2_H


/* headers and prototypes of functions also accessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

void ntainp_ccadiscret(std::string& inputfile_name,std::string& outputfile_kenner,std::string& restartfile_kenner);
void dyn_nlnstructural_drt();
void caldyn_drt();

#ifdef __cplusplus
}
#endif



#endif  /* #ifndef GLOBAL_INP_CONTROL2_H */

#endif  /* #ifdef CCADISCRET */
