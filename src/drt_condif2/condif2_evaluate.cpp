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
#include "condif2_utils.H"
#include "condif2_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // for CalculateFlux()
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
  else if (action == "calc_subgrid_diffusivity_matrix")
    act = Condif2::calc_subgrid_diffusivity_matrix;
  else if (action == "calc_condif_flux")
    act = Condif2::calc_condif_flux;
  else if (action == "calc_temp_and_dens")
    act = Condif2::calc_temp_and_dens;
  else
    dserror("Unknown type of action for Condif2");
  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                               vg 05/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2::Evaluate(ParameterList&            params,
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
      RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (hist==null || densnp==null || phinp==null)
        dserror("Cannot get state vector 'hist', 'densnp' and/or 'phinp'");

      // extract local values from the global vector
      vector<double> myhist(lm.size());
      vector<double> mydensnp(lm.size());
      vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
      DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

      // get control parameter
      const bool is_stationary = params.get<bool>("using stationary formulation");
      const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");
      const double time = params.get<double>("total time");

      // get time-step length
      const double dt = params.get<double>("time-step length");

      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
      double timefac = 1.0;
      double alphaF  = 1.0;
      if (not is_stationary)
      {
        timefac = params.get<double>("time factor");
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
        if (timefac < 0.0) dserror("time factor is negative.");
      }

      // set parameters for stabilization
      ParameterList& stablist = params.sublist("STABILIZATION");

      // select tau definition
      Condif2::TauType whichtau = Condif2::tau_not_defined;
      {
        const string taudef = stablist.get<string>("DEFINITION_TAU");

        if(taudef == "Franca_Valentin") whichtau = Condif2::franca_valentin;
        else if(taudef == "Bazilevs")   whichtau = Condif2::bazilevs;
      }

      // get (weighted) velocity values at the nodes (3rd component of velocity field is ignored!)
      // compare also with DRT::UTILS::ExtractMyValues()
      const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
      const int iel = NumNode();
      const int nsd=2;
      Epetra_SerialDenseVector evelnp(nsd*iel);
      DRT::UTILS::ExtractMyNodeBasedValues2D(this,evelnp,velocity);

      // get flag for fine-scale subgrid diffusivity
      string fssgd = params.get<string>("fs subgrid diffusivity","No");

      // check for non-existing subgrid-diffusivity models
      if (fssgd == "artificial_small" || fssgd == "Smagorinsky_all" || fssgd == "Smagorinsky_small")
        dserror("only all-scale artficial diffusivity for convection-diffusion problems possible so far!\n");

      // set flag for type of scalar whether it is temperature or not
      string scaltypestr=params.get<string>("problem type");
      bool temperature = false;
      if(scaltypestr =="loma") temperature = true;

      // calculate element coefficient matrix and rhs
      DRT::ELEMENTS::Condif2Impl::Impl(this)->Sysmat(
          this,
          myphinp,
          myhist,
          mydensnp,
          &elemat1,
          &elevec1,
          elevec2,
          actmat,
          time,
          dt,
          timefac,
          alphaF,
          evelnp,
          temperature,
          whichtau,
          fssgd,
          is_stationary,
          is_genalpha);
    }
    break;
    // calculate time derivative for time value t_0
    case initialize_one_step_theta:
    {
      const double time    = params.get<double>("total time");
      const double dt      = params.get<double>("time-step length");

      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
      double timefac = 1.0;
      double alphaF  = 1.0;
      timefac = params.get<double>("time factor");
      alphaF = params.get<double>("alpha_F");
      timefac *= alphaF;
      if (timefac < 0.0) dserror("time factor is negative.");

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

      // set parameters for stabilization
      ParameterList& stablist = params.sublist("STABILIZATION");

      // select tau definition
      Condif2::TauType whichtau = Condif2::tau_not_defined;
      {
        const string taudef = stablist.get<string>("DEFINITION_TAU");

        if(taudef == "Franca_Valentin") whichtau = Condif2::franca_valentin;
        else if(taudef == "Bazilevs")   whichtau = Condif2::bazilevs;
      }

      // get initial velocity values at the nodes
      // compare also with DRT::UTILS::ExtractMyValues()
      const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
      const int iel = NumNode();
      const int nsd=2;
      Epetra_SerialDenseVector evel0(nsd*iel);
      DRT::UTILS::ExtractMyNodeBasedValues2D(this,evel0,velocity);

      // get flag for fine-scale subgrid diffusivity
      string fssgd = params.get<string>("fs subgrid diffusivity","No");

      // check for non-existing subgrid-diffusivity models
      if (fssgd == "artificial_small" || fssgd == "Smagorinsky_all" || fssgd == "Smagorinsky_small")
        dserror("only all-scale artficial diffusivity for convection-diffusion problems possible so far!\n");

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
           dt,
           timefac,
           evel0,
           temperature,
           whichtau,
           fssgd);
    }
    break;
    // calculate normalized subgrid-diffusivity matrix
    case calc_subgrid_diffusivity_matrix:
    {
      // get control parameter
      const bool is_stationary = params.get<bool>("using stationary formulation");

      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
      double timefac = 1.0;
      double alphaF  = 1.0;
      if (not is_stationary)
      {
        timefac = params.get<double>("time factor");
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
        if (timefac < 0.0) dserror("time factor is negative.");
      }

      // calculate mass matrix and rhs
      DRT::ELEMENTS::Condif2Impl::Impl(this)->CalcSubgridDiffMatrix(
           this,
           elemat1,
           timefac,
           is_stationary);
    }
    break;
    case calc_condif_flux:
    {
      // get velocity values at the nodes
      // compare also with DRT::UTILS::ExtractMyValues()
      const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
      const int iel = NumNode();
      const int nsd=2;
      Epetra_SerialDenseVector evel(nsd*iel);
      DRT::UTILS::ExtractMyNodeBasedValues2D(this,evel,velocity);

      // need current values of transported scalar
      RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp==null) dserror("Cannot get state vector 'phinp'");

      // extract local values from the global vectors
      vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

      // assure, that the values are in the same order as the element nodes
      for(int k=0;k<iel;++k)
      {
        Node* node = (Nodes())[k];
        vector<int> dof = discretization.Dof(node);
        int numdof = dof.size();
        // up to now, there's only one dof per node
        for (int i=0;i<numdof;++i)
        {
          if (dof[i]!=lm[k*numdof+i])
          {
            cout<<"dof[i]= "<<dof[i]<<"  lm[k*numdof+i]="<<lm[k*numdof+i]<<endl;
            dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
          }
        }
      }

      // access control parameter
      Condif2::FluxType fluxtype;
      string fluxtypestring = params.get<string>("fluxtype","noflux");
      if (fluxtypestring == "totalflux")
        fluxtype = Condif2::totalflux;
      else if (fluxtypestring == "diffusiveflux")
        fluxtype = Condif2::diffusiveflux;
      else
        fluxtype=Condif2::noflux;  //default value

      // set flag for type of scalar
      string scaltypestr=params.get<string>("problem type");
      int numscal = numdofpernode_;
      bool temperature = false;
      if (scaltypestr =="loma") temperature = true;
      if (scaltypestr =="elch") numscal -= 1; // ELCH case: last dof is for el. potential

      // do a loop for systems of transported scalars
      for (int i = 0; i<numscal; ++i)
      {
        Epetra_SerialDenseMatrix eflux = CalculateFlux(myphinp,actmat,temperature,evel,fluxtype,i);

        for (int k=0;k<iel;k++)
        { // form arithmetic mean of assembled nodal flux vectors
          // => factor is the number of adjacent elements for each node
          double factor = (Nodes()[k])->NumElement();
          elevec1[k*numdofpernode_+i]+=eflux(0,k)/factor;
          elevec2[k*numdofpernode_+i]+=eflux(1,k)/factor;
          elevec3[k*numdofpernode_+i]+=eflux(2,k)/factor;
        }
      } // loop over numdofpernode
    }
    break;
    // calculate mean temperature and density
    case calc_temp_and_dens:
    {
      // need current scalar and density vector
      RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
      RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
      if (phinp==null || densnp==null)
        dserror("Cannot get state vector 'phinp' and/or 'densnp'");

      // extract local values from the global vectors
      vector<double> myphinp(lm.size());
      vector<double> mydensnp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
      DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);

      // calculate temperature, density and domain integral
      CalculateTempAndDens(params,myphinp,mydensnp);
    }
    break;
    default:
      dserror("Unknown type of action for Condif2");
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Condif2::Evaluate


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


/*----------------------------------------------------------------------*
 |  calculate mass flux                               (private) vg 09/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::ELEMENTS::Condif2::CalculateFlux(
    vector<double>&           ephinp,
    struct _MATERIAL*         material,
    bool                      temperature,
    Epetra_SerialDenseVector& evel,
    Condif2::FluxType         fluxtype,
    const int&                dofindex
)
{
  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  const int nsd = 2;

  Epetra_SerialDenseMatrix xyze(nsd,iel);
  Epetra_SerialDenseMatrix flux(nsd,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  // get diffusivity
  double diffus = 0;

  if (material->mattyp == m_matlist)
  {
    const int matid = material->m.matlist->matids[dofindex];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_condif)
      diffus = singlemat.m.condif->diffusivity;
    else if (singlemat.mattyp == m_ion)
      diffus = singlemat.m.ion->diffusivity;
    else
      dserror("type of material found in material list is not supported.");
  }
  else if (material->mattyp == m_condif)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for condif material"); // paranoia?
    if (temperature)
      diffus = material->m.condif->diffusivity/material->m.condif->shc;
    else
      diffus = material->m.condif->diffusivity;
  }
  else
    dserror("Material type is not supported");

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector        funct(iel);
  Epetra_SerialDenseMatrix        deriv(nsd,iel);
  static Epetra_SerialDenseMatrix xjm(nsd,nsd);
  Epetra_SerialDenseMatrix        derxy(nsd,iel);

  vector< LINALG::Matrix<3,1> > nodecoords;
  nodecoords = DRT::UTILS::getEleNodeNumbering_nodes_reference(distype);

  if ((int) nodecoords.size() != iel) dserror("number of nodes does not match");

  // loop over all nodes
  for (int iquad=0; iquad<iel; ++iquad)
  {
    // coordiantes of the current integration point
    const double e1 = nodecoords[iquad](0);
    const double e2 = nodecoords[iquad](1);

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

    /*----------------------------------------- compute Jacobian matrix */
    double dum;
    /*-------------------------------- determine jacobian at point r,s ---*/
    for (int i=0; i<nsd; i++)
    {
      for (int j=0; j<nsd; j++)
      {
        dum=0.0;
        for (int l=0; l<iel; l++)
        {
          dum += deriv(i,l)*xyze(j,l);
        }
        xjm(i,j)=dum;
      } /* end of loop j */
    } /* end of loop i */

    // The determinant is computed using Sarrus's rule:
    const double det = xjm(0,0)*xjm(1,1)-xjm(0,1)*xjm(1,0);

    if (det < 0.0)
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", Id(), det);
    if (abs(det) < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", Id(), det);

    /*------------------------------------------------------------------*/
    /*                                         compute global derivates */
    /*------------------------------------------------------------------*/

    /*------------------------------------------------- initialization */
    for(int k=0;k<iel;k++)
    {
      derxy(0,k)=0.0;
      derxy(1,k)=0.0;
    } /* end of loop over k */

    // ---------------------------------------inverse of transposed jacobian
    static Epetra_SerialDenseMatrix       xij(nsd,nsd);
    xij(0,0) =  xjm(1,1)/det;
    xij(1,0) = -xjm(1,0)/det;
    xij(0,1) = -xjm(0,1)/det;
    xij(1,1) =  xjm(0,0)/det;

    // ---------------------------------------- calculate global derivatives
    for (int k=0;k<iel;k++)
    {
      derxy(0,k) +=  xij(0,0) * deriv(0,k) + xij(0,1) * deriv(1,k) ;
      derxy(1,k) +=  xij(1,0) * deriv(0,k) + xij(1,1) * deriv(1,k) ;
    }

    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
      case Condif2::totalflux:
        //convective flux terms
        flux(0,iquad)+=evel[iquad*nsd]*ephinp[iquad*numdofpernode_+dofindex];
        flux(1,iquad)+=evel[1+iquad*nsd]*ephinp[iquad*numdofpernode_+dofindex];
        // no break statement here!
      case Condif2::diffusiveflux:
        //diffusive flux terms
        for (int k=0;k<iel;k++)
        {
          flux(0,iquad)+=-diffus*derxy(0,k)*ephinp[k*numdofpernode_+dofindex];
          flux(1,iquad)+=-diffus*derxy(1,k)*ephinp[k*numdofpernode_+dofindex];
        }
        break;
      case Condif2::noflux:
        dserror("received noflux flag inside CONDIF2 flux evaluation");
    };

  } // loop over nodes

  return flux;
} // Condif2::CalculateFlux


/*----------------------------------------------------------------------*
 |  calculate temperature, density and domain integral          vg 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2::CalculateTempAndDens(
    ParameterList&            params,
    vector<double>&           ephinp,
    vector<double>&           edensnp
)
{
  // get variables for integrals
  double tempint = params.get<double>("temperature integral");
  double densint = params.get<double>("density integral");
  double domint  = params.get<double>("domain integral");

  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  const int nsd = 2;

  // get node coordinates
  Epetra_SerialDenseMatrix xyze(nsd,iel);
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector        funct(iel);
  Epetra_SerialDenseMatrix        deriv(nsd,iel);
  static Epetra_SerialDenseMatrix xjm(nsd,nsd);

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(SCATRA::get2DOptimalGaussrule(distype));

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

    /*----------------------------------------- compute Jacobian matrix */
    double dum;
    /*-------------------------------- determine jacobian at point r,s ---*/
    for (int i=0; i<nsd; i++)
    {
      for (int j=0; j<nsd; j++)
      {
        dum=0.0;
        for (int l=0; l<iel; l++)
        {
          dum += deriv(i,l)*xyze(j,l);
        }
        xjm(i,j)=dum;
      } /* end of loop j */
    } /* end of loop i */

    // The determinant is computed using Sarrus's rule:
    const double det = xjm(0,0)*xjm(1,1)-xjm(0,1)*xjm(1,0);

    if (det < 0.0)
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", Id(), det);
    if (abs(det) < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", Id(), det);

    const double fac = intpoints.qwgt[iquad]*det; // Gauss weight * det(J)

    // calculate integrals of temperature, density and domain
    for (int i=0; i<iel; i++)
    {
      tempint += fac*funct[i]*ephinp[i];
      densint += fac*funct[i]*edensnp[i];
      domint  += fac*funct[i];
    }
  } // loop over nodes

  // return variables for integrals
  params.set<double>("temperature integral",tempint);
  params.set<double>("density integral",densint);
  params.set<double>("domain integral",domint);

  return;
} // Condif2::CalculateTempAndDens


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
