/*!----------------------------------------------------------------------
\file fluid3_boundary_impl.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_boundary_impl.H"
#include "fluid3_weak_dbc.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

#include "fluid3_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_geometry/position_array.H"


#include <blitz/array.h>

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3BoundaryImplInterface* DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(const DRT::Element* ele)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));

  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    static Fluid3BoundaryImpl<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new Fluid3BoundaryImpl<DRT::Element::quad4>(numdofpernode);
      return cp4;
  }
  case DRT::Element::quad8:
  {
    static Fluid3BoundaryImpl<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new Fluid3BoundaryImpl<DRT::Element::quad8>(numdofpernode);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static Fluid3BoundaryImpl<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new Fluid3BoundaryImpl<DRT::Element::quad9>(numdofpernode);
    return cp9;
  }
  case DRT::Element::tri3:
  {
    static Fluid3BoundaryImpl<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new Fluid3BoundaryImpl<DRT::Element::tri3>(numdofpernode);
      return cp3;
  }
  /*  case DRT::Element::tri6:
  {
    static ScaTraBoundaryImpl<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new ScaTraBoundaryImpl<DRT::Element::tri6>(numdofpernode,numscal);
    return cp6;
  }*/
  case DRT::Element::line2:
  {
    static Fluid3BoundaryImpl<DRT::Element::line2>* cl2;
    if (cl2==NULL)
      cl2 = new Fluid3BoundaryImpl<DRT::Element::line2>(numdofpernode);
      return cl2;
  }
  /*
  case DRT::Element::line3:
  {
    static ScaTraBoundaryImpl<DRT::Element::line3>* cl3;
    if (cl3==NULL)
      cl3 = new ScaTraBoundaryImpl<DRT::Element::line3>(numdofpernode,numscal);
    return cl3;
  }*/
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    static Fluid3BoundaryImpl<DRT::Element::nurbs2>* cn2;
    if (cn2==NULL)
      cn2 = new Fluid3BoundaryImpl<DRT::Element::nurbs2>(numdofpernode);
      return cn2;
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    static Fluid3BoundaryImpl<DRT::Element::nurbs3>* cn3;
    if (cn3==NULL)
      cn3 = new Fluid3BoundaryImpl<DRT::Element::nurbs3>(numdofpernode);
      return cn3;
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    static Fluid3BoundaryImpl<DRT::Element::nurbs4>* cn4;
    if (cn4==NULL)
      cn4 = new Fluid3BoundaryImpl<DRT::Element::nurbs4>(numdofpernode);
      return cn4;
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    static Fluid3BoundaryImpl<DRT::Element::nurbs9>* cn9;
    if (cn9==NULL)
      cn9 = new Fluid3BoundaryImpl<DRT::Element::nurbs9>(numdofpernode);
      return cn9;
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::Fluid3BoundaryImpl(int numdofpernode)
  : numdofpernode_(numdofpernode)
    {
      return;
    }


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::Evaluate(DRT::ELEMENTS::Fluid3Boundary* ele,
                                                ParameterList&            params,
                                                DRT::Discretization&      discretization,
                                                vector<int>&              lm,
                                                Epetra_SerialDenseMatrix& elemat1,
                                                Epetra_SerialDenseMatrix& elemat2,
                                                Epetra_SerialDenseVector& elevec1,
                                                Epetra_SerialDenseVector& elevec2,
                                                Epetra_SerialDenseVector& elevec3)
{
    DRT::ELEMENTS::Fluid3Boundary::ActionType act = Fluid3Boundary::none;
    string action = params.get<string>("action","none");
    if (action == "none") dserror("No action supplied");
    else if (action == "integrate_Shapefunction")
        act = Fluid3Boundary::integrate_Shapefunction;
    else if (action == "area calculation")
        act = Fluid3Boundary::areacalc;
    else if (action == "flowrate calculation")
        act = Fluid3Boundary::flowratecalc;
    // general action to calculate the flow rate (replaces flowratecalc soon)
    else if (action == "calc_line_flowrate")
      act = Fluid3Boundary::calc_flowrate;
    else if (action == "flowrate_deriv")
        act = Fluid3Boundary::flowratederiv;
    else if (action == "Outlet impedance")
        act = Fluid3Boundary::Outletimpedance;
    else if (action == "calc_node_normal")
        act = Fluid3Boundary::calc_node_normal;
    else if (action == "calc_node_curvature")
        act = Fluid3Boundary::calc_node_curvature;
    else if (action == "calc_surface_tension")
        act = Fluid3Boundary::calc_surface_tension;
    //TODO: weak dirichlet boundary condition
    else if (action == "enforce_weak_dbc")
        act = Fluid3Boundary::enforce_weak_dbc;
    else if (action == "conservative_outflow_bc")
        act = Fluid3Boundary::conservative_outflow_bc;
    else if (action == "calc_Neumann_inflow")
        act = Fluid3Boundary::calc_Neumann_inflow;
    else if (action == "calculate integrated pressure")
        act = Fluid3Boundary::integ_pressure_calc;
    else dserror("Unknown type of action for Fluid3_Boundary");

    // get status of Ale
    const bool isale = ele->ParentElement()->IsAle();

    switch(act)
    {
    case integrate_Shapefunction:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      if (isale)
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp!=null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
      }

      IntegrateShapeFunction(ele,params,discretization,lm,elevec1,mydispnp);
      break;
    }
    case areacalc:
    {
      AreaCaculation(ele, params, discretization,lm);
        break;
    }
    case flowratecalc:
    {
        FlowRateParameterCalculation(ele, params,discretization,lm);
        break;
    }
    case integ_pressure_calc:
    {
      IntegratedPressureParameterCalculation(ele, params,discretization,lm);
      break;
    }
    // general action to calculate the flow rate (replaces flowratecalc soon)
    case calc_flowrate:
    {
        ComputeFlowRate(ele, params,discretization,lm, elevec1);
        break;
    }
    case flowratederiv:
    {
      FlowRateDeriv(ele,params,discretization,lm,elemat1,elemat2,elevec1,elevec2,elevec3);
      break;
    }
    case Outletimpedance:
    {
        ImpedanceIntegration(ele,params,discretization,lm,elevec1);
        break;
    }
    case calc_node_normal:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      if (isale)
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp!=null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
      }
      ElementNodeNormal(ele,params,discretization,lm,elevec1,mydispnp);
      break;
    }
    case calc_node_curvature:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      if (isale)
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp!=null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
      }

      RefCountPtr<const Epetra_Vector> normals;
      vector<double> mynormals;

      normals = discretization.GetState("normals");
      if (normals!=null)
      {
        mynormals.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*normals,mynormals,lm);
      }

      // what happens, if the mynormals vector is empty? (ehrl)
      dserror("the action calc_node_curvature has not been called by now. What happens, if the mynormal vector is empty");

      ElementMeanCurvature(ele, params,discretization,lm,elevec1,mydispnp,mynormals);
      break;
    }
    case enforce_weak_dbc:
    {
      return DRT::ELEMENTS::Fluid3BoundaryWeakDBCInterface::Impl(ele)->EvaluateWeakDBC(
        ele,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;
    }
    case conservative_outflow_bc:
    {
      ConservativeOutflowConsistency(
        ele,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;
    }
    case calc_Neumann_inflow:
    {
      NeumannInflow(
        ele,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;
    }
    case calc_surface_tension:
    {
      // employs the divergence theorem acc. to Saksono eq. (24) and does not
      // require second derivatives.

      RCP<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }

      // mynormals and mycurvature are not used in the function
      vector<double> mynormals;
      vector<double> mycurvature;

      ElementSurfaceTension(ele,params,discretization,lm,elevec1,mydispnp,mynormals,mycurvature);
      break;
    }
    default:
        dserror("Unknown type of action for Fluid3_Surface");
    } // end of switch(act)

    return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::EvaluateNeumann(
                              DRT::ELEMENTS::Fluid3Boundary*         ele,
                              ParameterList&              params,
                              DRT::Discretization&        discretization,
                              DRT::Condition&             condition,
                              vector<int>&                lm,
                              Epetra_SerialDenseVector&   elevec1_epetra,
                              Epetra_SerialDenseMatrix*   elemat1_epetra)
{
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // get time-curve factor/ n = - grad phi / |grad phi|
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

  // get time parameter
  const double thsl = params.get("thsl",0.0);

  // get flag for type of fluid flow
  const INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  // allocate vector for shape functions and matrix for derivatives
  LINALG::Matrix<iel,1>   funct(true);
  LINALG::Matrix<bdrynsd_,iel>   deriv(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get local node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  LINALG::Matrix<nsd_,iel>   xyze(true);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  // get scalar vector
  RCP<const Epetra_Vector> scanp = discretization.GetState("scanp");
  if (scanp==null) dserror("Cannot get state vector 'scanp'");

  // extract local values from global vector
  vector<double> myscanp(lm.size());
  DRT::UTILS::ExtractMyValues(*scanp,myscanp,lm);

  LINALG::Matrix<iel,1>   escanp(true);

  // insert scalar into element array
  // the scalar is stored to the pressure dof
  for (int inode=0;inode<iel;++inode)
  {
    escanp(inode) = myscanp[(nsd_)+(inode*numdofpernode_)];
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is miultiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(iel);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(IsNurbs<distype>::isnurbs)
  {
     bool zero_size = GetKnotVectorAndWeightsForNurbs(ele, discretization, mypknots, myknots, weights, normalfac);
     if(zero_size)
     {
       return 0;
     }
  }
  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor, the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.) and density

    // integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, &myknots, &weights, xsi, funct, deriv);

    // compute temperature and density at the gauss point for low-Mach-number flow
    double dens = 1.0;
    if (physicaltype == INPAR::FLUID::loma)
    {
      // This is a hack for low-Mach-number flow with temperature
      // equation until material data will be available here
      // get thermodynamic pressure and its time derivative or history
      double thermpress = params.get<double>("thermpress at n+1",0.0);
      double gasconstant = 287.0;

      double temp = 0.0;
      temp = funct.Dot(escanp);
      dens = thermpress/(gasconstant*temp);
    }
    // compute density at the gauss point for varying density and Boussinesq
    else if(physicaltype == INPAR::FLUID::varying_density)
    {
      dens = funct.Dot(escanp);
    }

    const double fac_curvefac_thsl_dens =  fac *(curvefac * thsl * dens);

    // factor given by spatial function
    double functfac = 1.0;

    // global coordinates of gausspoint
    LINALG::Matrix<(nsd_),1>  coordgp(0.0);

    // determine coordinates of current Gauss point
    coordgp.Multiply(xyze,funct);

    int functnum = -1;
    const double* coordgpref = &coordgp(0); // needed for function evaluation

    for(int idim=0; idim<(nsd_); ++idim)
    {
      if((*onoff)[idim])  // Is this dof activated
      {
        if (func) functnum = (*func)[idim];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,coordgpref,time,NULL);
          }
          else
            functfac = 1.0;
        }
        const double val_fac_curvefac_thsl_dens_functfac = (*val)[idim]*fac_curvefac_thsl_dens*functfac;

        for(int inode=0; inode < iel; ++inode )
        {
        elevec1_epetra[inode*numdofpernode_+idim] += funct(inode)*val_fac_curvefac_thsl_dens_functfac;
        }
      }  // if (*onoff)
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 | apply outflow boundary condition which is necessary for the          |
 | conservative element formulation (since the convective term was      |
 | partially integrated)                                                |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::ConservativeOutflowConsistency(
    DRT::ELEMENTS::Fluid3Boundary*  ele,
    ParameterList&                  params,
    DRT::Discretization&            discretization,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       elemat1_epetra,
    Epetra_SerialDenseVector&       elevec1_epetra)
{
  // ------------------------------------
  //     GET TIME INTEGRATION DATA
  // ------------------------------------
  // we use two timefacs for matrix and right hand side to be able to
  // use the method for both time integrations
  const double timefac_mat = params.get<double>("timefac_mat");
  const double timefac_rhs = params.get<double>("timefac_rhs");

  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // ------------------------------------
  //     GET GENERAL ELEMENT DATA
  // ------------------------------------

  // vector for shape functions and matrix for derivatives
  LINALG::Matrix<iel,1>     funct(true);
  LINALG::Matrix<bdrynsd_,iel>  deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel>    xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // the element's normal vector
  LINALG::Matrix<nsd_,1>    norm(true);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  // ------------------------------------
  // get statevectors from discretisation
  // ------------------------------------

  // extract local displacements from the global vectors
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int inode=0;inode<iel;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // extract local velocities from the global vectors
  LINALG::Matrix<nsd_, iel>   evel(true);

  RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
  if (vel==null) dserror("Cannot get state vector 'u and p (trial)'");

  // extract local values from the global vectors
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

  for (int inode=0;inode<iel;++inode)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      evel(idim,inode) = myvel[numdofpernode_*inode+idim];
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is miultiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(iel);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights

  if(IsNurbs<distype>::isnurbs)
  {
     bool zero_size = GetKnotVectorAndWeightsForNurbs(ele, discretization, mypknots, myknots, weights, normalfac);
     if(zero_size)
     {
       return;
     }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, &myknots, &weights, xsi, funct, deriv, &norm);

    // Multiply the normal vector with the integration factor
    norm.Scale(fac);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs)
      norm.Scale(normalfac);

    // get velocity at the gauss point
    // velocity at gausspoint
    LINALG::Matrix<nsd_,1>    velint(true);
    velint.Multiply(evel,funct);

    // compute normal flux
    const double u_o_n = velint.Dot(norm);

    // rescaled flux (according to time integration)
    const double timefac_mat_u_o_n = timefac_mat*u_o_n;

   // dyadic product of element's normal vector and velocity
   LINALG::Matrix<nsd_,nsd_>  n_x_u(true);

   // dyadic product of u and n
   n_x_u.MultiplyNT(timefac_mat,velint,norm);

    /*
              /                \
             |                  |
           + |  Du o n , u o v  |
             |                  |
              \                /
    */

    // fill all velocity elements of the matrix
    for (int ui=0; ui<iel; ++ui) // loop columns
    {
      //Epetra_SerialDenseMatrix  temp(nsd_,nsd_) = n_x_u (copy);
      LINALG::Matrix<nsd_,nsd_>   temp(n_x_u, false);

      // temp(nsd_,nsd) = n_x_u(nsd_,nsd_)*funct(ui)
      temp.Scale(funct(ui));

      for (int idimcol=0; idimcol < (nsd_); ++idimcol) // loop over dimensions for the columns
      {
        const int fui   = numdofpernode_*ui+idimcol;

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          // temp(nsd_,nsd) *= funct(vi)
          temp.Scale(funct(vi));

          for (int idimrow = 0; idimrow < nsd_; ++idimrow) // loop over dimensions for the rows
          {
            const int fvi = numdofpernode_*vi+idimrow;
            elemat1_epetra(fvi  ,fui  ) = temp(fvi, fui);
          }  // end loop over dimensions for the rows
        } // end loop over rows (vi)
      } // end oop over dimensions for the columns
    } // end loop over columns (ui)

    /*
              /                \
             |                  |
           + |  u o n , Du o v  |
             |                  |
              \                /
    */

   // fill only diagonal velocity elements of the matrix
   for (int idim=0; idim < (nsd_); ++idim) // loop dimensions
   {
     for (int ui=0; ui<iel; ++ui) // loop columns
     {
       const int fui   = numdofpernode_*ui+idim;
       const double timefac_mat_u_o_n_funct_ui = timefac_mat_u_o_n*funct(ui);

       for (int vi=0; vi<iel; ++vi)  // loop rows
       {
         const int fvi = numdofpernode_*vi + idim;
         const double timefac_mat_u_o_n_funct_ui_funct_vi
                   =
                   timefac_mat_u_o_n_funct_ui*funct(vi);

         elemat1_epetra(fvi  ,fui  ) += timefac_mat_u_o_n_funct_ui_funct_vi;
       }  // loop rows
     }  // loop columns
   }  //loop over dimensions

  // rhs
  {
    // 3 temp vector
    LINALG::Matrix<nsd_,1>    temp(velint, false);

    // temp(nsd, nsd_) *= timefac_rhs * u_o_n
    temp.Scale(timefac_rhs*u_o_n);

    for (int vi=0; vi<iel; ++vi) // loop rows  (test functions)
    {
      for (int idim = 0; idim<(nsd_); ++idim) // loop over dimensions
      {
        int fvi=numdofpernode_*vi + idim;

    /*


                /               \
               |                 |
             + |  u o n , u o v  |
               |                 |
                \               /
    */

        elevec1_epetra(fvi) -= temp(fvi)*funct(vi);
      } // end loop over dimensions
    } // ui
  } // end rhs
  } // end gaussloop

  return;
}// DRT::ELEMENTS::Fluid3Surface::SurfaceConservativeOutflowConsistency


/*----------------------------------------------------------------------*
 | compute potential Neumann inflow                            vg 03/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::NeumannInflow(
    DRT::ELEMENTS::Fluid3Boundary*  ele,
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1,
    Epetra_SerialDenseVector&  elevec1)
{
  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  // check whether we have a generalized-alpha time-integration scheme
  const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");

  // get timefactor for left hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // get flag for type of fluid flow
  const INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // (density-weighted) shape functions and first derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<iel,1> densfunct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // node coordinate
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // the element's normal vector
  LINALG::Matrix<nsd_,1> normal(true);


  LINALG::Matrix<nsd_,1> momint(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  // determine outward-pointing normal to this element
  //normal(0) = (xyze(1,1)-xyze(1,0))*(xyze(2,2)-xyze(2,0))-(xyze(2,1)-xyze(2,0))*(xyze(1,2)-xyze(1,0));
  //normal(1) = (xyze(2,1)-xyze(2,0))*(xyze(0,2)-xyze(0,0))-(xyze(0,1)-xyze(0,0))*(xyze(2,2)-xyze(2,0));
  //normal(2) = (xyze(0,1)-xyze(0,0))*(xyze(1,2)-xyze(1,0))-(xyze(1,1)-xyze(1,0))*(xyze(0,2)-xyze(0,0));

  // Calculation of the normalized normal vector and the metrictensor
  //DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze,deriv,metrictensor,drs, &normal);

  // length of normal
/*  double length = 0.0;
  length = sqrt(normal(0)*normal(0)+normal(1)*normal(1)+normal(2)*normal(2));

  // outward-pointing normal of unit length
  for(int inode=0;inode<3;inode++)
  {
    normal(inode) = normal(inode)/length;
  }*/

  // get velocity and scalar vector at time n+alpha_F/n+1
  RCP<const Epetra_Vector> velaf = discretization.GetState("velaf");
  RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  if (velaf==null or scaaf==null)
    dserror("Cannot get state vector 'velaf' and/or 'scaaf'");

  // extract local values from global vector
  vector<double> myvelaf(lm.size());
  vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  // create Epetra objects for scalar array and velocities
  LINALG::Matrix<nsd_,iel> evelaf(true);
  LINALG::Matrix<iel,1> escaaf(true);

  // insert velocity and scalar into element array
  for (int inode=0;inode<iel;++inode)
  {
    for (int idim=0; idim<(nsd_);++idim)
    {
      evelaf(idim,inode) = myvelaf[idim+(inode*numdofpernode_)];
    }
    escaaf(inode) = myscaaf[(nsd_)+(inode*numdofpernode_)];
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // TODO: check with Volker, if normal at the Gauss point is similar to normal at node
    // Alternative: Implement a methode to calulate the normal at the node
    // Check the requirements for all normals

    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &normal);

    // compute velocity and normal velocity
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    // velocity at gausspoint
    LINALG::Matrix<nsd_,1> velint(true);

    double normvel = 0.0;
    velint.Multiply(evelaf,funct);
    normvel = velint.Dot(normal);

    // computation only required for negative normal velocity
    if (normvel<-0.0001)
    {
      // compute temperature and density for low-Mach-number flow
      double dens = 1.0;
      if(physicaltype == INPAR::FLUID::loma)
      {
        // This is a hack for low-Mach-number flow with temperature
        // equation until material data will be available here
        // get thermodynamic pressure and its time derivative or history
        double thermpress = params.get<double>("thermpress at n+alpha_F/n+1",0.0);
        double gasconstant = 287.0;

        double temp = 0.0;
        temp = funct.Dot(escaaf);
        dens = thermpress/(gasconstant*temp);
      }
      // TODO: varyingDensity??? Imnplement a method to get the density

      // integration factor for left- and right-hand side
      const double lhsfac = dens*normvel*timefac*fac;
      double rhsfac = dens*normvel*fac;

      // genalpha does not use a time factor for the rhs
      if (not is_genalpha) rhsfac *= timefac;

      // matrix
      // fill diagonal elements
      for (int idim = 0; idim < nsd_; ++idim) // loop over dimensions
      {
        for (int vi=0; vi<iel; ++vi) // loop over rows
        {
          const double vlhs = lhsfac*funct(vi);

          const int fvi = numdofpernode_*vi+idim;

          for (int ui=0; ui<iel; ++ui) // loop over columns
          {
            const int fui = numdofpernode_*ui+idim;

            // negative contribution to the matrix
            elemat1(fvi  ,fui  ) -= vlhs*funct(ui);
          } // end loop over columns
        }  // end loop over rows
      }  // end loop over dimensions

      // rhs
      LINALG::Matrix<nsd_,1> vrhs(velint, false);
      vrhs.Scale(rhsfac);

      for (int vi=0; vi<iel; ++vi) //loop over rows
      {
        for (int idim = 0; idim < nsd_; ++idim)  //loop over dimensions
        {
          const int fvi   = numdofpernode_*vi+idim;
          elevec1(fvi) += funct(vi)*vrhs(idim);
        } // end loop over dimensions
      }  // end loop over rows
    }
  }

  return;
}// DRT::ELEMENTS::Fluid3Surface::NeumannInflow


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)            gjb 07/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::IntegrateShapeFunction(
                  DRT::ELEMENTS::Fluid3Boundary* ele,
                  ParameterList& params,
                  DRT::Discretization&       discretization,
                  vector<int>&               lm,
                  Epetra_SerialDenseVector&  elevec1,
                  const std::vector<double>& edispnp)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // allocate vector for shape functions and matrix for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,iel> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<iel;++inode)
    {
      for (int idim=0;idim<(nsd_);++idim)
      {
        xyze(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv);

    for (int inode=0;inode<iel;++inode)
    {
      for(int idim=0;idim<(nsd_);idim++)
      {
        elevec1(inode*numdofpernode_+idim)+= funct(inode) * fac;
      }
    }

  } /* end of loop over integration points gpid */


return;
} // DRT::ELEMENTS::Fluid3Surface::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::ElementNodeNormal(
                                                     DRT::ELEMENTS::Fluid3Boundary*   ele,
                                                     ParameterList&                   params,
                                                     DRT::Discretization&             discretization,
                                                     vector<int>&                     lm,
                                                     Epetra_SerialDenseVector&        elevec1,
                                                     const std::vector<double>&       edispnp)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  //get gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // allocate vector for shape functions and matrix for derivatives
  LINALG::Matrix<iel,1> funct(0.0);
  LINALG::Matrix<bdrynsd_,iel> deriv(0.0);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(0.0);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(0.0);

  // normal vector in gausspoint
  LINALG::Matrix<nsd_,1> norm(0.0);


  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<iel; ++inode)
    {
      for (int idim=0;idim<(nsd_); ++idim)
      {
        xyze(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  //const IntegrationPoints2D  intpoints(gaussrule);

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &norm);

    for (int inode=0; inode<iel; ++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
      {
        elevec1(inode*numdofpernode_+idim) += norm(idim) * funct(inode) * fac;
      }
      // pressure dof is set to zero
      elevec1(inode*numdofpernode_+(nsd_)) = 0.0;
    }
  } /* end of loop over integration points gpid */

  return;

} // DRT::ELEMENTS::Fluid3Surface::ElementNodeNormal


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::ElementMeanCurvature(
                                                        DRT::ELEMENTS::Fluid3Boundary*    ele,
                                                        ParameterList&                    params,
                                                        DRT::Discretization&              discretization,
                                                        vector<int>&                      lm,
                                                        Epetra_SerialDenseVector&         elevec1,
                                                        const std::vector<double>&        edispnp,
                                                        std::vector<double>&              enormals)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // node normals &
  LINALG::Matrix<nsd_,iel> norm_elem(true);
  LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);

  // coordinates of current node in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi_node(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<iel; ++inode)
    {
      for (int idim=0;idim<nsd_; ++idim)
      {
        xyze(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // set normal vectors to length = 1.0
  // normal vector is coming from outside
  for (int inode=0; inode<iel; ++inode)
  {
    //double length = 0.0;
    for (int idim=0; idim < nsd_; ++idim)
    {
      norm_elem(idim,inode) = enormals[numdofpernode_*inode+idim];
    }
  }
  // compute normalized normal vector
  norm_elem.Scale(1/norm_elem.Norm2());

  // get local node coordinates of the element
  // function gives back a matrix with the local node coordinates of the element (nsd_,iel)
  // the function gives back an Epetra_SerialDenseMatrix!!!
  Epetra_SerialDenseMatrix xsi_ele = getEleNodeNumbering_nodes_paramspace(distype);

  // ============================== loop over nodes ==========================
  for (int inode=0;inode<iel; ++inode)
  {
    // the local node coordinates matrix is split to a vector containing the local coordinates of the actual node
    for (int idim = 0; idim < bdrynsd_; idim++)
    {
      xsi_node(idim) = xsi_ele(idim,inode);
    }

    // get shape derivatives at this node
    // shape_function_2D_deriv1(deriv, e0, e1, distype);
    DRT::UTILS::shape_function<distype>(xsi_node,funct);

    // the metric tensor and its determinant
    //Epetra_SerialDenseMatrix      metrictensor(nsd_,nsd_);
    LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);
    double drs;

    // Addionally, compute metric tensor
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze,deriv,metrictensor,drs);

    dxyzdrs.MultiplyNT(deriv,xyze);

    // calculate mean curvature H at node.
    double H = 0.0;
    LINALG::Matrix<bdrynsd_,nsd_> dn123drs(0.0);

    dn123drs.MultiplyNT(deriv,norm_elem);

    //Acc. to Bronstein ..."mittlere Kruemmung":
    // calculation of the mean curvature for a surface element
    if (bdrynsd_==2)
    {
      double L = 0.0, twoM = 0.0, N = 0.0;
      for (int i=0;i<3;i++)
      {
        L += (-1.0) * dxyzdrs(0,i) * dn123drs(0,i);
        twoM += (-1.0) * dxyzdrs(0,i) * dn123drs(1,i) - dxyzdrs(1,i) * dn123drs(0,i);
        N += (-1.0) * dxyzdrs(1,i) * dn123drs(1,i);
      }
      //mean curvature: H = 0.5*(k_1+k_2)
      H = 0.5 *
          (metrictensor(0,0)*N - twoM*metrictensor(0,1) + metrictensor(1,1)*L)
          / (drs*drs);
    }
    else
     dserror("Calcualtion of the mean curvature is only implemented for a 2D surface element");


    // get the number of elements adjacent to this node. Find out how many
    // will contribute to the interpolated mean curvature value.
    int contr_elements = 0;
    DRT::Node* thisNode = (ele->Nodes())[inode];
#ifdef DEBUG
    if (thisNode == NULL) dserror("No node!\n");
#endif
    int NumElement = thisNode->NumElement();
    DRT::Element** ElementsPtr = thisNode->Elements();

    // loop over adjacent Fluid3 elements
    for (int ele=0;ele<NumElement;ele++)
    {
      DRT::Element* Element = ElementsPtr[ele];

      // get surfaces
      vector< RCP< DRT::Element > > surfaces = Element->Surfaces();

      // loop over surfaces: how many free surfaces with this node on it?
      for (unsigned int surf=0; surf<surfaces.size(); ++surf)
      {
        RCP< DRT::Element > surface = surfaces[surf];
        DRT::Node** NodesPtr = surface->Nodes();
        int numfsnodes = 0;
        bool hasthisnode = false;

        for (int surfnode = 0; surfnode < surface->NumNode(); ++surfnode)
        {
          DRT::Node* checkNode = NodesPtr[surfnode];
          // check whether a free surface condition is active on this node
          if (checkNode->GetCondition("FREESURFCoupling") != NULL)
          {
            numfsnodes++;
          }
          if (checkNode->Id() == thisNode->Id())
          {
            hasthisnode = true;
          }
        }

        if (numfsnodes == surface->NumNode() and hasthisnode)
        {
          // this is a free surface adjacent to this node.
          contr_elements++;
        }

      }

    }
#ifdef DEBUG
    if (!contr_elements) dserror("No contributing elements found!\n");
#endif

    for(int idim=0; idim<nsd_; ++idim)
    {
      elevec1[inode*numdofpernode_+idim] = H / contr_elements;
    }
    elevec1[inode*numdofpernode_+(numdofpernode_-1)] = 0.0;
  } // END: loop over nodes

} // DRT::ELEMENTS::Fluid3Surface::ElementMeanCurvature



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::ElementSurfaceTension(
                                                         DRT::ELEMENTS::Fluid3Boundary*   ele,
                                                         ParameterList&                   params,
                                                         DRT::Discretization&             discretization,
                                                         vector<int>&                     lm,
                                                         Epetra_SerialDenseVector&        elevec1,
                                                         const std::vector<double>&       edispnp,
                                                         std::vector<double>&             enormals,
                                                         std::vector<double>&             ecurvature)
                                                         // Attention: mynormals and mycurvature are not used in the function
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get time integration parameters
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // isotropic and isothermal surface tension coefficient
  double SFgamma = 0.0;
  // get material data
  RCP<MAT::Material> mat = ele->ParentElement()->Material();
  if (mat==null)
    dserror("no mat from parent!");
  else if (mat->MaterialType()==INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    SFgamma = actmat->Gamma();
  }
  else
    dserror("Newtonian fluid material expected but got type %d", mat->MaterialType());

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1>       funct(true);
  LINALG::Matrix<bdrynsd_,iel>    deriv(true);

  // node coordinates
  LINALG::Matrix<nsd_, iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_, 1> xsi(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<iel; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector

    double drs;
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, NULL, &drs);

    // fac multiplied by the timefac
    const double fac_timefac = fac * timefac;

    // Compute dxyzdrs
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);
    dxyzdrs.MultiplyNT(deriv,xyze);

    if (bdrynsd_==2)
    {
      double abs_dxyzdr = 0.0;
      double abs_dxyzds = 0.0;
      double pointproduct_rs = 0.0;

      for (int dim=0;dim<3;dim++)
      {
        abs_dxyzdr += dxyzdrs(0,dim) * dxyzdrs(0,dim);
        abs_dxyzds += dxyzdrs(1,dim) * dxyzdrs(1,dim);
        pointproduct_rs += dxyzdrs(0,dim) * dxyzdrs(1,dim);
      }
      abs_dxyzdr = sqrt(abs_dxyzdr);
      abs_dxyzds = sqrt(abs_dxyzds);

      for (int node=0;node<iel;++node)
      {
        for (int dim=0;dim<3;dim++)
        {
          // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
          // function) on Gamma_FS
          // See Saksono eq. (26)
          // discretized as surface gradient * ( Shapefunction-Matrix
          // transformed )

          // This uses a surface_gradient extracted from gauss general
          // formula for 2H...
          // this gives convincing results with TET elements, but HEX
          // elements seem more difficult -> due to edge problems?
          // too many nonlinear iterations
          elevec1[node*numdofpernode_+dim] += SFgamma *
                                     (-1.0) / (
                                       drs * drs //= abs_dxyzdr * abs_dxyzdr * abs_dxyzds * abs_dxyzds - pointproduct_rs * pointproduct_rs
                                       )
                                     *
                                     (
                                       abs_dxyzds * abs_dxyzds * deriv(0,node) * dxyzdrs(0,dim)
                                       - pointproduct_rs * deriv(0,node) * dxyzdrs(1,dim)
                                       - pointproduct_rs * deriv(1,node) * dxyzdrs(0,dim)
                                       + abs_dxyzdr * abs_dxyzdr * deriv(1,node) * dxyzdrs(1,dim)
                                       )
                                     * fac_timefac;

        }
        elevec1[node*numdofpernode_+3] = 0.0;
      }
    } // end if (nsd_==2)
    else if (bdrynsd_==1)
    {
      for (int inode=0;inode<iel;++inode)
      {
         for(int idim=0;idim<2;idim++)
         {
            // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
            // function) on Gamma_FS
            // See Saksono eq. (26)
            // discretized as surface gradient * ( Shapefunction-Matrix
            // transformed )
            // 2D: See Slikkerveer ep. (17)
            elevec1[inode*numdofpernode_+idim]+= SFgamma / drs / drs *
                                      (-1.0) * deriv(0, inode) * dxyzdrs(0,idim)
                                      * fac_timefac;
         }
      }
    } // end if else (nsd_=1)
    else
      dserror("There are no 3D boundary elements implemented");
  } /* end of loop over integration points gpid */
} // DRT::ELEMENTS::Fluid3Surface::ElementSurfaceTension

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::AreaCaculation(
                                                  DRT::ELEMENTS::Fluid3Boundary*  ele,
                                                  ParameterList&                  params,
                                                  DRT::Discretization&            discretization,
                                                  vector<int>&                    lm)
{

  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // get the required material information
  RCP<MAT::Material> mat = ele->ParentElement()->Material();
  double density=0.0, viscosity=0.0;

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    density = actmat->Density();
    viscosity =  actmat->Viscosity();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    density = actmat->Density();
    dserror("How to extract viscosity from Carreau Yasuda material law for artery tree??");
  }
  else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
    density = actmat->Density();
    dserror("How to extract viscosity from modified power law material for artery tree??");
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());

  // set required material data for communication
  params.set<double>("density", density);
  params.set<double>("viscosity", viscosity);

  // get gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  double area        = params.get<double>("Area calculation");

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);


#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<iel;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

#endif // D_ALE_BFLOW

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv);

    area += fac;
  }

  params.set<double>("Area calculation", area);
}//DRT::ELEMENTS::Fluid3Surface::AreaCaculation


// TODO: Replace function by ComputeFlowRate
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::FlowRateParameterCalculation(
                                                                DRT::ELEMENTS::Fluid3Boundary*    ele,
                                                                ParameterList&                    params,
                                                                DRT::Discretization&              discretization,
                                                                vector<int>&                      lm)
{
  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  //get gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");

  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // alocate local velocity vector
  LINALG::Matrix<nsd_,iel> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0; inode<iel; ++inode)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
    }
  }

  // get  actual outflowrate
  double flowrate    = params.get<double>("Outlet flowrate");

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<iel;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  // TODO: Normal at Gauss point or at node

  //const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector at gausspoint
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &normal);

    //Compute elment flowrate (add to actual frow rate obtained before
    for (int inode=0;inode<iel;++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
      {
        flowrate += funct(inode) * evelnp(idim,inode)*normal(idim) *fac;
      }
    }
  }  // end Gauss loop
  // set new flow rate
  params.set<double>("Outlet flowrate", flowrate);
}//DRT::ELEMENTS::Fluid3Surface::FlowRateParameterCalculation


/*----------------------------------------------------------------------*
 |                                                        ismail 04/2010|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::IntegratedPressureParameterCalculation(
  DRT::ELEMENTS::Fluid3Boundary*    ele,
  ParameterList&                    params,
  DRT::Discretization&              discretization,
  vector<int>&                      lm)
{

  //------------------------------------------------------------------
  // This calculates the integrated the pressure from the 
  // the actual pressure values
  //------------------------------------------------------------------
#if 1
  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = ele->ParentElement()->Material();

  double density = 0.0;

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
    density = actmat->Density();
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());

  //get gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");

  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate local velocity vector
  LINALG::Matrix<nsd_,iel> evelnp(true);
  // allocate local pressure vector
  LINALG::Matrix<1,iel>   eprenp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0; inode<iel; ++inode)
  {
    eprenp(inode) = density*myvelnp[nsd_+inode*numdofpernode_];
  }

  // get  actual outflowrate
  double pressure    = params.get<double>("Inlet integrated pressure");

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<iel;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  // TODO: Normal at Gauss point or at node

  //const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector at gausspoint
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &normal);

    //Compute elment flowrate (add to actual frow rate obtained before
    for (int inode=0;inode<iel;++inode)
    {
      pressure += funct(inode) * eprenp(inode) *fac;
    }
  }  // end Gauss loop
  // set new flow rate

  params.set<double>("Inlet integrated pressure", pressure);
#endif
}//DRT::ELEMENTS::Fluid3Surface::IntegratedPressureParameterCalculation


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::ComputeFlowRate(
                                                                DRT::ELEMENTS::Fluid3Boundary*    ele,
                                                                ParameterList&                    params,
                                                                DRT::Discretization&              discretization,
                                                                vector<int>&                      lm,
                                                                Epetra_SerialDenseVector&         elevec1)
{
  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");

  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate velocity vector
  LINALG::Matrix<nsd_,iel> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<iel;inode++)
  {
    for (int idim=0; idim< nsd_; idim++)
    {
      evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
    }
  }

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<iel;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW


  //const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector at gausspoint
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &normal);

    //compute flowrate at gauss point
    LINALG::Matrix<nsd_,1> velint(true);
    velint.Multiply(evelnp,funct);

    // flowrate = uint o normal
    const double flowrate = velint.Dot(normal);

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode=0;inode<iel;++inode)
    {
      elevec1[inode*numdofpernode_] -= funct(inode)* fac * flowrate; //
    }
  }
}//DRT::ELEMENTS::Fluid3Surface::ComputeFlowRate


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::FlowRateDeriv(
                                                 DRT::ELEMENTS::Fluid3Boundary*   ele,
                                                 ParameterList&                   params,
                                                 DRT::Discretization&             discretization,
                                                 vector<int>&                     lm,
                                                 Epetra_SerialDenseMatrix&        elemat1,
                                                 Epetra_SerialDenseMatrix&        elemat2,
                                                 Epetra_SerialDenseVector&        elevec1,
                                                 Epetra_SerialDenseVector&        elevec2,
                                                 Epetra_SerialDenseVector&        elevec3)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  RefCountPtr<const Epetra_Vector> dispnp;
  vector<double> edispnp;

  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    edispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,edispnp,lm);
  }

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // order of accuracy of grid velocity determination
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  int gridvel = fdyn.get<int>("GRIDVEL");

  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<iel; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get nodal velocities and pressures
  RefCountPtr<const Epetra_Vector> convelnp = discretization.GetState("convectivevel");

  if (convelnp==null)
    dserror("Cannot get state vector 'convectivevel'");

  // extract local values from the global vectors
  vector<double> myconvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*convelnp,myconvelnp,lm);

  // allocate velocities vector
  LINALG::Matrix<nsd_,iel> evelnp(true);

  for (int inode=0; inode<iel; ++inode)
  {
    for (int idim=0;idim<nsd_; ++idim)
    {
      evelnp(idim,inode) = myconvelnp[(numdofpernode_*inode)+idim];
    }
  }


  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    //TODO: normal computed at gausspoint or node??

    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normalized normal vector at gausspoint
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &normal);

    // dxyzdrs vector
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(deriv,xyze);

    if(bdrynsd_==2)
    {
      //-------------------------------------------------------------------
      //  Q
      blitz::Array<double,1> u(3);
      u = 0.;
      for (int dim=0;dim<3;++dim)
          for (int node=0;node<iel;++node)
            u(dim) += funct(node) * evelnp(dim,node);

      for(int dim=0;dim<3;++dim)
        elevec3[0] += u(dim) * normal(dim) * fac;

      if (params.get<bool>("flowrateonly", false)==false)
      {
        //-------------------------------------------------------------------
        // dQ/du
        for (int node=0;node<iel;++node)
        {
          for (int dim=0;dim<3;++dim)
            elevec1[node*numdofpernode_+dim] += funct(node) * normal(dim) * fac;
          elevec1[node*numdofpernode_+3] = 0.0;
        }

        //-------------------------------------------------------------------
        // dQ/dd

        // determine derivatives of surface normals wrt mesh displacements
        blitz::Array<double,2> normalderiv(3,iel*3);

        for (int node=0;node<iel;++node)
        {
          normalderiv(0,3*node)   = 0.;
          normalderiv(0,3*node+1) = deriv(0,node)*dxyzdrs(1,2)-deriv(1,node)*dxyzdrs(0,2);
          normalderiv(0,3*node+2) = deriv(1,node)*dxyzdrs(0,1)-deriv(0,node)*dxyzdrs(1,1);

          normalderiv(1,3*node)   = deriv(1,node)*dxyzdrs(0,2)-deriv(0,node)*dxyzdrs(1,2);
          normalderiv(1,3*node+1) = 0.;
          normalderiv(1,3*node+2) = deriv(0,node)*dxyzdrs(1,0)-deriv(1,node)*dxyzdrs(0,0);

          normalderiv(2,3*node)   = deriv(0,node)*dxyzdrs(1,1)-deriv(1,node)*dxyzdrs(0,1);
          normalderiv(2,3*node+1) = deriv(1,node)*dxyzdrs(0,0)-deriv(0,node)*dxyzdrs(1,0);
          normalderiv(2,3*node+2) = 0.;
        }

        for (int node=0;node<iel;++node)
        {
          for (int dim=0;dim<3;++dim)
            for (int iterdim=0;iterdim<3;++iterdim)
              elevec2[node*numdofpernode_+dim] += u(iterdim) * normalderiv(iterdim,3*node+dim) * fac;
          elevec2[node*numdofpernode_+3] = 0.0;
        }

        // consideration of grid velocity
        if (isale)
        {
          // get time step size
          const double dt = params.get<double>("dt", -1.0);
          if (dt < 0.) dserror("invalid time step size");

          if (gridvel == 1)  // BE time discretization
          {
            for (int node=0;node<iel;++node)
            {
              for (int dim=0;dim<3;++dim)
                elevec2[node*numdofpernode_+dim] -= 1.0/dt * funct(node) * normal(dim) * fac;
            }
          }
          else
            dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
        }

        //-------------------------------------------------------------------
        // (d^2 Q)/(du dd)

        for (int unode=0;unode<iel;++unode)
        {
          for (int udim=0;udim<numdofpernode_;++udim)
          {
            for (int nnode=0;nnode<iel;++nnode)
            {
              for (int ndim=0;ndim<numdofpernode_;++ndim)
              {
                if (udim == 3 or ndim == 3)
                  elemat1(unode*numdofpernode_+udim,nnode*numdofpernode_+ndim) = 0.0;
                else
                  elemat1(unode*numdofpernode_+udim,nnode*numdofpernode_+ndim) = funct(unode) * normalderiv(udim,3*nnode+ndim) * fac;
              }
            }
          }
        }

        //-------------------------------------------------------------------
        // (d^2 Q)/(dd)^2

        // determine second derivatives of surface normals wrt mesh displacements
        blitz::Array<double,3> normalderiv2(3,iel*3,iel*3);
        normalderiv2 = 0.;

        for (int node1=0;node1<iel;++node1)
        {
          for (int node2=0;node2<iel;++node2)
          {
            double temp = deriv(0,node1)*deriv(1,node2)-deriv(1,node1)*deriv(0,node2);

            normalderiv2(0,node1*3+1,node2*3+2) = temp;
            normalderiv2(0,node1*3+2,node2*3+1) = - temp;

            normalderiv2(1,node1*3  ,node2*3+2) = - temp;
            normalderiv2(1,node1*3+2,node2*3  ) = temp;

            normalderiv2(2,node1*3  ,node2*3+1) = temp;
            normalderiv2(2,node1*3+1,node2*3  ) = - temp;
          }
        }

        for (int node1=0;node1<iel;++node1)
        {
          for (int dim1=0;dim1<numdofpernode_;++dim1)
          {
            for (int node2=0;node2<iel;++node2)
            {
              for (int dim2=0;dim2<numdofpernode_;++dim2)
              {
                if (dim1 == 3 or dim2 == 3)
                  elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) = 0.0;
                else
                {
                  for (int iterdim=0;iterdim<3;++iterdim)
                    elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) +=
                      u(iterdim) * normalderiv2(iterdim,node1*3+dim1,node2*3+dim2) * fac;
                }
              }
            }
          }
        }

        // consideration of grid velocity
        if (isale)
        {
          // get time step size
          const double dt = params.get<double>("dt", -1.0);
          if (dt < 0.) dserror("invalid time step size");

          if (gridvel == 1)
          {
            for (int node1=0;node1<iel;++node1)
            {
              for (int dim1=0;dim1<3;++dim1)
              {
                for (int node2=0;node2<iel;++node2)
                {
                  for (int dim2=0;dim2<3;++dim2)
                  {
                    elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) -= (1.0/dt * funct(node1) * normalderiv(dim1, 3*node2+dim2)
                                                                   + 1.0/dt * funct(node2) * normalderiv(dim2, 3*node1+dim1))
                                                                  * fac;
                  }
                }
              }
            }
          }
          else
            dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
        }


        //-------------------------------------------------------------------
      }
    }
    else
      dserror("FlowRateDeriv is only implemented for 3D -> 2D with 1D boundary elements is still missing!! ");
  }
}//DRT::ELEMENTS::Fluid3Surface::FlowRateDeriv


 /*----------------------------------------------------------------------*
  |  Impedance related parameters on boundary elements          AC 03/08  |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::ImpedanceIntegration(
                  DRT::ELEMENTS::Fluid3Boundary*    ele,
                  ParameterList&                    params,
                  DRT::Discretization&              discretization,
                  vector<int>&                      lm,
                  Epetra_SerialDenseVector&         elevec1)
{
  const double thsl = params.get("thsl",0.0);

  double invdensity=0.0; // inverse density of my parent element

  // allocate vector for shape functions and for derivatives
  LINALG::Matrix<iel,1> funct(true);
  LINALG::Matrix<bdrynsd_,iel> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,iel> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = ele->ParentElement()->Material();

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());

  double pressure = params.get<double>("ConvolutedPressure");

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze);

#ifdef D_ALE_BFLOW
  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<iel;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // integration factor * Gauss weights & local Gauss point coordinates
    // shape function at the Gauss point & derivative of the shape function at the Gauss point
    // normal vector of unit length at gausspoint (outward pointing)
    const double fac = EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, &normal);

    const double fac_thsl_pres_inve = fac * thsl * pressure * invdensity;

    for (int inode=0;inode<iel;++inode)
      for(int idim=0;idim<nsd_;++idim)
        // inward pointing normal of unit length
        elevec1[inode*numdofpernode_+idim] += funct(inode) * fac_thsl_pres_inve * (-normal(idim));
  }
  return;
} //DRT::ELEMENTS::Fluid3Surface::ImpedanceIntegration


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<bdrynsd_>&  intpoints,
    const int                                         gpid,
    const LINALG::Matrix<nsd_,iel>&                   xyze,
    const std::vector<Epetra_SerialDenseVector>*      myknots,
    const Epetra_SerialDenseVector*                   weights,
    LINALG::Matrix<bdrynsd_,1>&                       xsi,
    LINALG::Matrix<iel,1>&                            funct,
    LINALG::Matrix<bdrynsd_,iel>&                     deriv,
    LINALG::Matrix<nsd_,1>*                           normal,
    double*                                           drs)
{
  // local coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[gpid];
  for (int idim=0;idim<bdrynsd_;++idim)
  {
    xsi(idim) = gpcoord[idim];
  }

  // get shape functions and derivatives in the plane of the element
  if(not IsNurbs<distype>::isnurbs)
  {
    // shape function derivs of boundary element at gausspoint
    DRT::UTILS::shape_function<distype>(xsi,funct);
    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
  }
  else
  {
    if (bdrynsd_==2)  // TODO: Nurbs for 2D and 3D
    {
      // this is just a temporary work-around
      Epetra_SerialDenseVector gp(2);
      gp(0)=xsi(0);
      gp(1)=xsi(1);

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct  ,
         deriv  ,
         gp     ,
         (*myknots),
         (*weights),
         distype);
    }
    else if(bdrynsd_==1)
    {
      //const double gp = xsi_(0);
      dserror("1d Fluid3Boundary nurbs elements not yet implemented");
      //DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv(funct,deriv,gp,myknots,weights,distype);
    }
    else dserror("Discretisation type %s not yet implemented",DRT::DistypeToString(distype).c_str());
  }

  // compute measure tensor for surface element and the infinitesimal
  // area element drs for the integration
  double sqrtdetg(0.0);
  LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);

  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze,deriv,metrictensor, sqrtdetg, normal);

  // return drs
  if (drs != NULL)
    (*drs)=sqrtdetg;

  return intpoints.IP().qwgt[gpid]*sqrtdetg;  // return integration factor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::GetKnotVectorAndWeightsForNurbs(
    DRT::ELEMENTS::Fluid3Boundary*              ele,
    DRT::Discretization&                        discretization,
    std::vector<Epetra_SerialDenseVector>&      mypknots,
    std::vector<Epetra_SerialDenseVector>&      myknots,
    Epetra_SerialDenseVector&                   weights,
    double&                                     normalfac)
{
  // TODO: Check function 1D / 2D
  // ehrl
  if (bdrynsd_ == 1)
    dserror("1D line element -> It is not check if it is working.");

  // get pointer to parent element
  DRT::ELEMENTS::Fluid3* parent_ele = ele->ParentElement();

  // local surface id
  const int surfaceid = ele->SurfaceNumber();

  // --------------------------------------------------
  // get knotvector

  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

  bool zero_size = knots->GetBoundaryEleAndParentKnots(mypknots     ,
                                                     myknots      ,
                                                     normalfac    ,
                                                     parent_ele->Id(),
                                                     surfaceid    );

  // --------------------------------------------------
  // get node weights for nurbs elements
  for (int inode=0; inode<iel; ++inode)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (ele->Nodes()[inode]);

    weights(inode) = cp->W();
  }
  return zero_size;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
