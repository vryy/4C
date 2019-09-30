/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for fluid/xfluid and ScaTra coupled routines.
\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236

*----------------------------------------------------------------------*/

#ifndef TWO_PHASE_DYN_H
#define TWO_PHASE_DYN_H

/*! entry point for the solution of two phase flow problems */
void two_phase_dyn(int restart /* do we have to perform a restart?  */
);

/*! entry point for the solution of two phase flow problems utilizing XFEM */
void fluid_xfem_ls_drt(int restart /* do we have to perform a restart?  */
);

void DoSafetyCheck();

#endif /* #ifndef TWO_PHASE_DYN_H */
