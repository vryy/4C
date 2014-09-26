/*----------------------------------------------------------------------*/
/*!
\file drt_uq_dyn.H
\brief  Uncertainty Quantification

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/


#ifndef DRT_UQ_DYN_H_
#define DRT_UQ_DYN_H_


void dyn_uq();

//! brief setup and run uq with redairways
void redairways_uq();

//! brief setup and run uq with structure problem
void structure_uq();

#endif /* DRT_UQ_DYN_H_ */
