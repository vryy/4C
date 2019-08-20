/*----------------------------------------------------------------------------*/
/** \file
\brief inequality level-set approach a.k.a. xcontact or extended contact

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------------*/


#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_DYN_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_DYN_H_

/*! \brief Entry point for the solution of a eXtended contact problem
 *
 *  Inequality level-set approach to handle contact problems.
 *
 *  \param[in] restart Do we have to perform a restart?
 *
 *  \date 06/16 */
void xcontact_dyn(const int& restart);


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_DYN_H_ */
