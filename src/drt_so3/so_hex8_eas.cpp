/*!----------------------------------------------------------------------
\file so_hex8_eas.cpp
\brief Everything concernign EAS technology for so_hex8

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"


extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                              maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_easinit()
{
  DSTraceHelper dst("So_hex8::soh8_easinit");
  
  // EAS enhanced strain parameters
  Epetra_SerialDenseVector alpha(neas_);
  // EAS portion of internal forces, also called enhacement vector s or Rtilde
  Epetra_SerialDenseVector feas(neas_);
  // EAS matrix K_{alpha alpha}, also called Dtilde
  Epetra_SerialDenseMatrix Kaa(neas_,neas_);
  // EAS matrix K_{d alpha}
  Epetra_SerialDenseMatrix Kda(neas_,NUMDOF_SOH8);
  
  // save EAS data into element container easdata_
  easdata_.Add("alpha",alpha);
  easdata_.Add("feas",feas);
  easdata_.Add("Kaa",Kaa);
  easdata_.Add("Kda",Kda);
  
  return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
