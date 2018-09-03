/*!----------------------------------------------------------------------
\file loma_dyn.H

\brief Control routine for low-Mach-number flow module.

\level 2

\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245

*----------------------------------------------------------------------*/

#ifndef LOMA_DYN_H
#define LOMA_DYN_H

/*! entry point for the solution of low-Mach-number flow problems */
void loma_dyn(int restart /* do we have to perform a restart?  */
);


#endif /* #ifndef LOMA_DYN_H */
