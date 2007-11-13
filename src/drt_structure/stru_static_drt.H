/*!----------------------------------------------------------------------
\file stru_dyn_nln_drt.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef STRU_STA_NLN_DRT_H
#define STRU_STA_NLN_DRT_H


#ifdef __cplusplus /* stuff which is c++ specific */

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"

using namespace std;
using namespace Teuchos;



#endif /* end of #ifdef __cplusplus */





/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C" 
{
#endif

void stru_static_drt();

void synchronize_VolConstraint(const Epetra_Comm& Comm, 
								ParameterList& params,
								Epetra_SerialDenseVector& vect,
								int& mi, int& ma, int & numID);

void setupVolDofrowmaps(RefCountPtr<DRT::Discretization> actdis,
		  map<int, RCP<Epetra_Map> >& voldofrowmaps,
		  const int numconstID,
		  const string& condstring);

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef STRU_STA_NLN_DRT_H */
#endif  /* #ifdef CCADISCRET */
