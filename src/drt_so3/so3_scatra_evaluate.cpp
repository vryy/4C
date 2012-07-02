/*!----------------------------------------------------------------------
\file So3_scatra_evaluate.cpp

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "so3_scatra.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/robinson.H"
#include "../drt_mat/micromaterial.H"
#include <iterator>

//#include "../drt_mat/fluidporo.H"
//#include "../drt_mat/structporo.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_globalproblem.H"

//#include "Sacado.hpp"

/*----------------------------------------------------------------------*
 |  preevaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<distype>::PreEvaluate(ParameterList& params,
                                        DRT::Discretization&      discretization,
                                        DRT::Element::LocationArray& la)
{
  // TODO: read
  if(la.Size()>1)
  {
    //  dofs per node of second dofset
    const int numdofpernode = NumDofPerNode(1,*(Nodes()[0]));

    if (la[1].Size() != numnod_*numdofpernode)
      dserror("calc_struct_nlnstiff: Location vector length for velocities does not match!");

    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> tempnp
        = discretization.GetState(1,"temperature");

      if (tempnp==Teuchos::null)
        dserror("calc_struct_nlnstiff: Cannot get state vector 'fluidvel' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double> >mytemp = rcp(new std::vector<double>(la[1].lm_.size()) );
      DRT::UTILS::ExtractMyValues(*tempnp,*mytemp,la[1].lm_);

      params.set<Teuchos::RCP<vector<double> > >("scalar",mytemp);
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Scatra<distype>::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{

  return 0;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<distype>::InitJacobianMapping()
{
  //const static vector<LINALG::Matrix<numdim_,numnod_> > derivs;// = soh8_derivs();
  LINALG::Matrix<numdim_,numnod_> deriv ;
  LINALG::Matrix<numnod_,numdim_> xrefe;
  for (int i=0; i<numnod_; ++i)
  {
    Node** nodes=Nodes();
    if(!nodes) dserror("Nodes() returned null pointer");
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  for (int gp=0; gp<numgpt_; ++gp)
  {
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim=0;idim<numdim_;idim++)
    {
       xsi_[gp](idim) = gpcoord[idim];
    }

    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    //invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(deriv,xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0",detJ_[gp]);
  }

  return;
}

#include "so3_scatra_fwd.hpp"

