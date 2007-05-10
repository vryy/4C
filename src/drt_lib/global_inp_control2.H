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

#include "drt_discret.H"
#include "drt_condition.H"
#include "drt_utils.H"
#ifdef D_SHELL8
#include "../drt_s8/shell8.H"
#endif
#ifdef D_FLUID2
#include "fluid2.H"
#endif
#ifdef D_FLUID3
#include "fluid3.H"
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
#include "dstrc.H" /* needs standardtypes.h */

using namespace std;
using namespace Teuchos;

// headers and prototypes of functions acessible from c++
void inpnodes_ccadiscret(Epetra_SerialDenseMatrix& tmpnodes);
void inpnodes_ccadiscret_jumbo(RefCountPtr<DRT::Discretization>& dis,
                               RefCountPtr<Epetra_Map>& rownodes,
                               RefCountPtr<Epetra_Map>& colnodes);
void input_structural_field(FIELD *structfield, RefCountPtr<Epetra_Comm> comm);
void input_fluid_field(FIELD *structfield, RefCountPtr<Epetra_Comm> comm);
void input_field_jumbo(RefCountPtr<DRT::Discretization>& dis,
                                  RefCountPtr<Epetra_Map>& rownodes,
                                  RefCountPtr<Epetra_Map>& colnodes,
                                  RefCountPtr<Epetra_Map>& roweles,
                                  RefCountPtr<Epetra_Map>& coleles,
                                  const string& searchword);
void input_field_graph(int& nele,
                                  RefCountPtr<Epetra_Comm> comm,
                                  RefCountPtr<Epetra_Map>& rownodes,
                                  RefCountPtr<Epetra_Map>& colnodes,
                                  RefCountPtr<Epetra_CrsGraph>& graph,
                                  const string& searchword);
void input_assign_nodes(DRT::Discretization& actdis,Epetra_SerialDenseMatrix* tmpnodes);
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
void inpfield_ccadiscret();
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
