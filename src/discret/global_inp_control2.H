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
#ifdef TRILINOS_PACKAGE

#ifndef GLOBAL_INP_CONTROL2_H
#define GLOBAL_INP_CONTROL2_H


#ifdef __cplusplus // stuff which is c++ specific

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_RefCountPtr.hpp"

#include "discret.H"
#ifdef D_SHELL8
#include "shell8.H"
#endif

using namespace std;
using namespace Teuchos;

#endif // end of #ifdef __cplusplus



// headers and prototypes of functions also acessible from c
#ifdef __cplusplus
extern "C" 
{
#endif

#include "../headers/standardtypes.h"
#include "dstrc.H"


void ntainp_ccadiscret();
void inpfield_ccadiscret();
void inpnodes_ccadiscret(Epetra_SerialDenseMatrix& tmpnodes);
void input_structural_field(FIELD *structfield, RefCountPtr<Epetra_Comm> comm);
void input_assign_nodes(CCADISCRETIZATION::Discretization& actdis,Epetra_SerialDenseMatrix& tmpnodes);


#ifdef __cplusplus
}
#endif






#endif  // #ifndef GLOBAL_INP_CONTROL2_H

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
