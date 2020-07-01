/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport control algorithm
\level 1


*----------------------------------------------------------------------*/

#ifndef SCATRA_DYN_H
#define SCATRA_DYN_H

//! enumeration of coupling type
enum FieldCouplingType
{
  coupl_none,       //< unknown coupling type
  coupl_volmortar,  //< volmortar coupling
  coupl_match       //< matching coupling (clone scatra from fluid)
};

/*! entry point for the solution of scalar transport problems */
void scatra_dyn(int restart /* do we have to perform a restart?  */
);

#endif /* #ifndef SCATRA_DYN_H */
