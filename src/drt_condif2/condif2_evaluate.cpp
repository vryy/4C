/*!----------------------------------------------------------------------
\file condif2_evaluate.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "condif2.H"
#include "condif2_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include <Epetra_SerialDenseSolver.h>

#include "../drt_lib/drt_globalproblem.H"

using namespace DRT::UTILS;

/*---------------------------------------------------------------------*
|  converts a string into an action for this element(private) gjb 06/08|
 *---------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2::ActionType DRT::ELEMENTS::Condif2::convertStringToActionType(
  const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::Condif2::ActionType act = Condif2::none;
  if (action == "calc_condif_systemmat_and_residual")
    act = Condif2::calc_condif_systemmat_and_residual;
  else if (action == "initialize_one_step_theta")
    act = Condif2::initialize_one_step_theta;
  else if (action == "calc_condif_flux")
    act = Condif2::calc_condif_flux;
  else
    dserror("Unknown type of action for Condif2");
  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                               vg 05/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2::Evaluate(ParameterList& params,
                                     DRT::Discretization&      discretization,
                                     vector<int>&              lm,
                                     Epetra_SerialDenseMatrix& elemat1,
                                     Epetra_SerialDenseMatrix& elemat2,
                                     Epetra_SerialDenseVector& elevec1,
                                     Epetra_SerialDenseVector& elevec2,
                                     Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const string action = params.get<string>("action","none");
  const DRT::ELEMENTS::Condif2::ActionType act = convertStringToActionType(action);

  // get the material
  RefCountPtr<MAT::Material> mat = Material();

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
  else if (mat->MaterialType()== m_matlist)
    actmat = static_cast<MAT::MatList*>(mat.get())->MaterialData();
  else
    dserror("condif material expected but got type %d", mat->MaterialType());

  switch(act)
  {
    case calc_condif_systemmat_and_residual:
    {
      // need current history vector and
      // density*specific heat capacity at constant pressure vector
      RefCountPtr<const Epetra_Vector> hist = discretization.GetState("hist");
      RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
      if (hist==null || densnp==null)
        dserror("Cannot get state vector 'hist' and/or 'densnp'");

      // extract local values from the global vector
      vector<double> myhist(lm.size());
      vector<double> mydensnp(lm.size());
      DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
      DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);

      // get control parameter
      const bool is_stationary = params.get<bool>("using stationary formulation");
      const double time = params.get<double>("total time");

      // One-step-Theta: timefac = theta*dt
      // BDF2:           timefac = 2/3 * dt
      double timefac = 0.0;
      if (not is_stationary)
      {
        timefac = params.get<double>("thsl");
        if (timefac < 0.0) dserror("No thsl supplied");
      }

      // get (weighted) velocity values at the nodes (3rd component of velocity field is ignored!)
      // compare also with DRT::UTILS::ExtractMyValues()
      const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
      const int iel = NumNode();
      const int nsd=2;
      Epetra_SerialDenseVector evel(nsd*iel);
      ExtractMyNodeBasedValues2D(evel,velocity);

      // get flag for fine-scale subgrid diffusivity
      string fssgd = params.get<string>("fs subgrid diffusivity","No");

      // set flag for type of scalar whether it is temperature or not
      string scaltypestr=params.get<string>("problem type");
      bool temperature = false;
      if(scaltypestr =="loma") temperature = true;

      // calculate element coefficient matrix and rhs
      DRT::ELEMENTS::Condif2Impl::Impl(this)->Sysmat(
          this,
          myhist,
          mydensnp,
          &elemat1,
          &elemat2,
          &elevec1,
          elevec2,
          actmat,
          time,
          timefac,
          evel,
          temperature,
          fssgd,
          is_stationary);
    }
    break;
    // calculate time derivative for time value t_0
    case initialize_one_step_theta:
    {
      const double time = params.get<double>("total time");
      const double timefac = params.get<double>("thsl");

      // need initial field
      RefCountPtr<const Epetra_Vector> phi0 = discretization.GetState("phi0");
      RefCountPtr<const Epetra_Vector> dens0 = discretization.GetState("dens0");
      if (phi0==null || dens0==null)
        dserror("Cannot get state vector 'phi0' and/or 'densnp'");

      // extract local values from the global vector
      vector<double> myphi0(lm.size());
      vector<double> mydens0(lm.size());
      DRT::UTILS::ExtractMyValues(*phi0,myphi0,lm);
      DRT::UTILS::ExtractMyValues(*dens0,mydens0,lm);

      // get initial velocity values at the nodes
      // compare also with DRT::UTILS::ExtractMyValues()
      const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
      const int iel = NumNode();
      const int nsd=2;
      Epetra_SerialDenseVector evel(nsd*iel);
      DRT::UTILS::ExtractMyNodeBasedValues(this,evel,velocity);

      // get flag for fine-scale subgrid diffusivity
      string fssgd = params.get<string>("fs subgrid diffusivity","No");

      // set flag for type of scalar whether it is temperature or not
      string scaltypestr=params.get<string>("problem type");
      bool temperature = false;
      if(scaltypestr =="loma") temperature = true;

      // calculate mass matrix and rhs
      DRT::ELEMENTS::Condif2Impl::Impl(this)->InitializeOST(
           this,
           myphi0,
           mydens0,
           elemat1,
           elevec1,
           elevec2,
           actmat,
           time,
           timefac,
           evel,
           temperature,
           fssgd);
    }
    break;
    case calc_condif_flux:
      // do nothing here instead of throwing a dserror 
      // this keeps the result test on fluxes alive
    break;
    default:
      dserror("Unknown type of action for Condif2");
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Condif2::Evaluate



/*----------------------------------------------------------------------*
 |  extract my velocity dof from global vector (private)       gjb 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2::ExtractMyNodeBasedValues2D(
     Epetra_SerialDenseVector&      local,
     const RCP<Epetra_MultiVector>& global)
{
  if (global==null) dserror("received a TEUCHOS::null pointer");
  //const int nsd = global->NumVectors(); // get dimension
  const int nsd = 2; // get dimension
  const int iel = NumNode(); // number of nodes
  if (local.Length()!=(iel*nsd)) dserror("vector size mismatch.");

  for (int i=0; i<nsd; i++)
  {
    // get actual component column of velocity multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      local(i+(nsd*j))=globalcolumn[lid];
    }
  }
  return;
} //Condif2::ExtractMyNodeBasedValues2D


/*----------------------------------------------------------------------*
 |  do nothing (public)                                         vg 08/07|
 |                                                                      |
 |  The function is just a dummy. For the condif2 elements, the         |
 |  integration of the surface neumann loads takes place in the element.|
 |  We need it there for the stabilization terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2::EvaluateNeumann(ParameterList&            params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            vector<int>&              lm,
                                            Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
