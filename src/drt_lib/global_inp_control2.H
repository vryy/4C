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


#ifdef __cplusplus /* stuff which is c++ specific */

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "drt_inputreader.H"

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_condition.H"
#include "drt_utils.H"
#include "linalg_utils.H"
#ifdef D_SHELL8
#include "../drt_s8/shell8.H"
#endif
#ifdef D_FLUID2
#include "../drt_f2/fluid2.H"
#include "../drt_f2/condif2.H"
#endif
#ifdef D_FLUID3
#include "../drt_f3/fluid3.H"
#endif
#ifdef D_FLUID3_XFEM
#include "../drt_f3_xfem/fluid3_xfem.H"
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
#include "dstrc.H" /* needs standardtypes.h */

using namespace std;
using namespace Teuchos;

// headers and prototypes of functions acessible from c++
void inpnodes_ccadiscret_jumbo(RefCountPtr<DRT::Discretization>& dis,
                               RefCountPtr<Epetra_Map>& rownodes,
                               RefCountPtr<Epetra_Map>& colnodes,
                               const string& searchword);
void input_design_dpoint_fenode_read(vector<vector<int> >& dnode_fenode,
                                     vector<int>& ndnode_fenode);
void input_design_dline_fenode_read(vector<vector<int> >& dline_fenode,
                                     vector<int>& ndline_fenode);
void input_design_dsurf_fenode_read(vector<vector<int> >& dsurf_fenode,
                                     vector<int>& ndsurf_fenode);
void input_design_dvol_fenode_read(vector<vector<int> >& dvol_fenode,
                                   vector<int>& ndvol_fenode);


#endif /* end of #ifdef __cplusplus */



/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif


void ntainp_ccadiscret();
void inpfield_ccadiscret_jumbo();
void input_conditions();
void distribute_drt_grids();
void dyn_nlnstructural_drt();

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef GLOBAL_INP_CONTROL2_H */

#endif  /* #ifdef TRILINOS_PACKAGE */
#endif  /* #ifdef CCADISCRET */
