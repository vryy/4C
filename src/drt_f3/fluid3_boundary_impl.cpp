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
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
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
#include "../drt_lib/drt_globalproblem.H"


#include <blitz/array.h>

#include <cstdlib>

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3BoundaryImplInterface* DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(const DRT::Element* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return Fluid3BoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return Fluid3BoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return Fluid3BoundaryImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return Fluid3BoundaryImpl<DRT::Element::tri3>::Instance();
  }
  /*  case DRT::Element::tri6:
  {
    return Fluid3BoundaryImpl<DRT::Element::tri6>::Instance();
  }*/
  case DRT::Element::line2:
  {
    return Fluid3BoundaryImpl<DRT::Element::line2>::Instance();
  }
  /*
  case DRT::Element::line3:
  {
    return Fluid3BoundaryImpl<DRT::Element::line3>::Instance();
  }*/
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return Fluid3BoundaryImpl<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return Fluid3BoundaryImpl<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return Fluid3BoundaryImpl<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return Fluid3BoundaryImpl<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3BoundaryImpl<distype> * DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::Instance()
{
  static Fluid3BoundaryImpl<distype> * instance;
  if ( instance==NULL )
    instance = new Fluid3BoundaryImpl<distype>();
  return instance;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::Fluid3BoundaryImpl()
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
    // TODO: remove this action -> talk with Mahmoud
    else if (action == "flowrate calculation")
        act = Fluid3Boundary::flowratecalc;
    else if (action == "calc_flowrate")
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
    //TODO: weak Dirichlet boundary condition:
    //There is still a 2D and a 3D implementation
    else if (action == "enforce_weak_dbc")
        act = Fluid3Boundary::enforce_weak_dbc;
    else if (action == "MixedHybridDirichlet")
        act = Fluid3Boundary::mixed_hybrid_dbc;
    else if (action == "conservative_outflow_bc")
        act = Fluid3Boundary::conservative_outflow_bc;
    else if (action == "calc_Neumann_inflow")
        act = Fluid3Boundary::calc_Neumann_inflow;
    else if (action == "calculate integrated pressure")
        act = Fluid3Boundary::integ_pressure_calc;
    else if (action == "center of mass calculation")
        act = Fluid3Boundary::center_of_mass_calc;
    else dserror("Unknown type of action for Fluid3_Boundary: %s",action.c_str());

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
      if (ele->Owner() == discretization.Comm().MyPID())
        AreaCaculation(ele, params, discretization,lm);
      break;
    }
    case flowratecalc:
    {
      if (ele->Owner() == discretization.Comm().MyPID())
         FlowRateParameterCalculation(ele, params,discretization,lm);
      break;
    }
    case integ_pressure_calc:
    {
      if(ele->Owner() == discretization.Comm().MyPID())
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
    case mixed_hybrid_dbc:
    {
      switch (distype)
      {
      // 2D:
      case DRT::Element::line2:
      {
        if(ele->ParentElement()->Shape()==DRT::Element::quad4)
        {
          MixHybDirichlet<DRT::Element::line2,DRT::Element::quad4>(ele,
                                                                   params,
                                                                   discretization,
                                                                   lm,
                                                                   elemat1,
                                                                   elevec1);
        }
        else
        {
          dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
        }
        break;
      }
      // 3D:
      case DRT::Element::quad4:
      {
        if(ele->ParentElement()->Shape()==DRT::Element::hex8)
        {
          MixHybDirichlet<DRT::Element::quad4,DRT::Element::hex8>(ele,
                                                                  params,
                                                                  discretization,
                                                                  lm,
                                                                  elemat1,
                                                                  elevec1);
        }
        else
        {
          dserror("expected combination quad4/hex8 for surface/parent pair");
        }
        break;
      }
      default:
      {
        dserror("not implemented yet\n");
      }

      }

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
    case center_of_mass_calc:
    {
      // evaluate center of mass
      if(ele->Owner() == discretization.Comm().MyPID())
        CenterOfMassCalculation(ele, params,discretization,lm);
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

  // get time factor for Neumann term
  const double timefac = params.get("thsl",0.0);

  // get flag for type of fluid flow
  const INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  // allocate vector for shape functions and matrix for derivatives
  LINALG::Matrix<bdrynen_,1>   funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_>   deriv(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get local node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  LINALG::Matrix<nsd_,bdrynen_>   xyze(true);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

  // get scalar vector
  RCP<const Epetra_Vector> scanp = discretization.GetState("scanp");
  if (scanp==null) dserror("Cannot get state vector 'scanp'");

  // extract local values from global vector
  vector<double> myscanp(lm.size());
  DRT::UTILS::ExtractMyValues(*scanp,myscanp,lm);

  LINALG::Matrix<bdrynen_,1>   escanp(true);

  // insert scalar into element array
  // the scalar is stored to the pressure dof
  for (int inode=0;inode<bdrynen_;++inode)
  {
    escanp(inode) = myscanp[(nsd_)+(inode*numdofpernode_)];
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

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

    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurbs specific stuff
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, &myknots, &weights, xsi, funct, deriv, drs, NULL);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    // compute temperature and density at the gauss point for low-Mach-number flow
    double dens = 1.0;
    if (physicaltype == INPAR::FLUID::loma)
    {
      // This is a hack for low-Mach-number flow with temperature
      // equation until material data will be available here
      // get thermodynamic pressure and its time derivative or history
      double thermpress  = params.get<double>("thermodynamic pressure",0.0);
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

    const double fac_drs_curvefac_timefac_dens = fac_drs*curvefac*timefac*dens;

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
        const double valfac = (*val)[idim]*fac_drs_curvefac_timefac_dens*functfac;

        for(int inode=0; inode < bdrynen_; ++inode )
        {
          elevec1_epetra[inode*numdofpernode_+idim] += funct(inode)*valfac;
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
  LINALG::Matrix<bdrynen_,1>     funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_>  deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_>    xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // the element's normal vector
  LINALG::Matrix<nsd_,1>    unitnormal(true);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // extract local velocities from the global vectors
  LINALG::Matrix<nsd_, bdrynen_>   evel(true);

  RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
  if (vel==null) dserror("Cannot get state vector 'u and p (trial)'");

  // extract local values from the global vectors
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

  for (int inode=0;inode<bdrynen_;++inode)
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
  Epetra_SerialDenseVector weights(bdrynen_);

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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurbs specific stuff
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, &myknots, &weights, xsi, funct, deriv, drs, &unitnormal);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    // Multiply the normal vector with the integration factor
    unitnormal.Scale(fac_drs);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs)
      unitnormal.Scale(normalfac);

    // get velocity at the gauss point
    // velocity at gausspoint
    LINALG::Matrix<nsd_,1>    velint(true);
    velint.Multiply(evel,funct);

    // compute normal flux
    const double u_o_n = velint.Dot(unitnormal);

    // rescaled flux (according to time integration)
    const double timefac_mat_u_o_n = timefac_mat*u_o_n;

   // dyadic product of element's normal vector and velocity
   LINALG::Matrix<nsd_,nsd_>  n_x_u(true);

   // dyadic product of u and n
   n_x_u.MultiplyNT(timefac_mat,velint,unitnormal);

    /*
              /                \
             |                  |
           + |  Du o n , u o v  |
             |                  |
              \                /
    */

    // fill all velocity elements of the matrix
    for (int ui=0; ui<bdrynen_; ++ui) // loop columns
    {
      //Epetra_SerialDenseMatrix  temp(nsd_,nsd_) = n_x_u (copy);
      LINALG::Matrix<nsd_,nsd_>   temp(n_x_u, false);

      // temp(nsd_,nsd) = n_x_u(nsd_,nsd_)*funct(ui)
      temp.Scale(funct(ui));

      for (int idimcol=0; idimcol < (nsd_); ++idimcol) // loop over dimensions for the columns
      {
        const int fui   = numdofpernode_*ui+idimcol;

        for (int vi=0; vi<bdrynen_; ++vi)  // loop rows
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
     for (int ui=0; ui<bdrynen_; ++ui) // loop columns
     {
       const int fui   = numdofpernode_*ui+idim;
       const double timefac_mat_u_o_n_funct_ui = timefac_mat_u_o_n*funct(ui);

       for (int vi=0; vi<bdrynen_; ++vi)  // loop rows
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

    for (int vi=0; vi<bdrynen_; ++vi) // loop rows  (test functions)
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

  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // get timefactor for right-hand side
  // One-step-Theta:            timefacrhs = theta*dt
  // BDF2:                      timefacrhs = 2/3 * dt
  // af-generalized-alpha:      timefacrhs = (alpha_F/alpha_M) * gamma * dt
  // Peter's generalized-alpha: timefacrhs = 1.0
  double timefacrhs;
  if (is_genalpha)
  {
     timefacrhs = params.get<double>("rhs time factor",-1.0);
     if (timefacrhs < 0.0) dserror("incorrect time factor for right-hand side!");
  }
  else timefacrhs = timefac;

  // get flag for type of fluid flow
  const INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // (density-weighted) shape functions and first derivatives
  LINALG::Matrix<bdrynen_,1>        funct(true);
  LINALG::Matrix<bdrynen_,1>        densfunct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // node coordinate
  LINALG::Matrix<nsd_,bdrynen_>     xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1>   xsi(true);

  // the element's normal vector
  LINALG::Matrix<nsd_,1>       unitnormal(true);

  // velocity and momentum at gausspoint
  LINALG::Matrix<nsd_,1>       momint(true);
  LINALG::Matrix<nsd_,1>       velint(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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
  LINALG::Matrix<nsd_,bdrynen_> evelaf(true);
  LINALG::Matrix<bdrynen_,1> escaaf(true);

  // insert velocity and scalar into element array
  for (int inode=0;inode<bdrynen_;++inode)
  {
    for (int idim=0; idim<(nsd_);++idim)
    {
      evelaf(idim,inode) = myvelaf[idim+(inode*numdofpernode_)];
    }
    escaaf(inode) = myscaaf[(nsd_)+(inode*numdofpernode_)];
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is miultiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;

    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, &myknots, &weights, xsi, funct, deriv, drs, &unitnormal);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs)
      unitnormal.Scale(normalfac);

    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    // compute velocity and normal velocity
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    // velocity at gausspoint

    double normvel = 0.0;
    velint.Multiply(evelaf,funct);
    normvel = velint.Dot(unitnormal);

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
      else if(physicaltype == INPAR::FLUID::varying_density)
        dserror("The actual density is not calculated for the BC NeumannInflow in the case of VaryingDensity");

      // integration factor for left- and right-hand side
      const double lhsfac = dens*normvel*timefac*fac_drs;
      const double rhsfac = dens*normvel*timefacrhs*fac_drs;

      // matrix
      // fill diagonal elements
      for (int idim = 0; idim < nsd_; ++idim) // loop over dimensions
      {
        for (int vi=0; vi<bdrynen_; ++vi) // loop over rows
        {
          const double vlhs = lhsfac*funct(vi);

          const int fvi = numdofpernode_*vi+idim;

          for (int ui=0; ui<bdrynen_; ++ui) // loop over columns
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

      for (int vi=0; vi<bdrynen_; ++vi) //loop over rows
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<bdrynen_;++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, NULL);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);idim++)
      {
        elevec1(inode*numdofpernode_+idim)+= funct(inode) * fac_drs;
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
  LINALG::Matrix<bdrynen_,1> funct(0.0);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(0.0);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(0.0);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(0.0);

  // normal vector in gausspoint
  LINALG::Matrix<nsd_,1> unitnormal(0.0);


  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<bdrynen_; ++inode)
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

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, &unitnormal);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
      {
        elevec1(inode*numdofpernode_+idim) += unitnormal(idim) * funct(inode) * fac_drs;
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // node normals &
  LINALG::Matrix<nsd_,bdrynen_> norm_elem(true);
  LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);

  // coordinates of current node in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi_node(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0;idim<nsd_; ++idim)
      {
        xyze(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // set normal vectors to length = 1.0
  // normal vector is coming from outside
  for (int inode=0; inode<bdrynen_; ++inode)
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
  // function gives back a matrix with the local node coordinates of the element (nsd_,bdrynen_)
  // the function gives back an Epetra_SerialDenseMatrix!!!
  Epetra_SerialDenseMatrix xsi_ele = getEleNodeNumbering_nodes_paramspace(distype);

  // ============================== loop over nodes ==========================
  for (int inode=0;inode<bdrynen_; ++inode)
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
  LINALG::Matrix<bdrynen_,1>       funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_>    deriv(true);

  // node coordinates
  LINALG::Matrix<nsd_, bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_, 1> xsi(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, NULL);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    // fac multiplied by the timefac
    const double fac_drs_timefac = fac_drs * timefac;

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

      for (int node=0;node<bdrynen_;++node)
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
                                     * fac_drs_timefac;

        }
        elevec1[node*numdofpernode_+3] = 0.0;
      }
    } // end if (nsd_==2)
    else if (bdrynsd_==1)
    {
      for (int inode=0;inode<bdrynen_;++inode)
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
                                      * fac_drs_timefac;
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

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
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);


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
    for (int inode=0;inode<bdrynen_;++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, NULL);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    area += fac_drs;
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> unitnormal(true);

  //get gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");

  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // alocate local velocity vector
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0; inode<bdrynen_; ++inode)
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
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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
    for (int inode=0;inode<bdrynen_;++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
     EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, &unitnormal);
     const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    //Compute elment flowrate (add to actual frow rate obtained before
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0; idim<nsd_; ++idim)
      {
        flowrate += funct(inode) * evelnp(idim,inode)*unitnormal(idim) *fac_drs;
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> unitnormal(true);

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
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);
  // allocate local pressure vector
  LINALG::Matrix<1,bdrynen_>   eprenp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0; inode<bdrynen_; ++inode)
  {
    eprenp(inode) = density*myvelnp[nsd_+inode*numdofpernode_];
  }

  // get  actual outflowrate
  double pressure    = params.get<double>("Inlet integrated pressure");

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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
    for (int inode=0;inode<bdrynen_;++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, &unitnormal);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    //Compute elment flowrate (add to actual frow rate obtained before
    for (int inode=0;inode<bdrynen_;++inode)
    {
      pressure += funct(inode) * eprenp(inode) *fac_drs;
    }
  }  // end Gauss loop
  // set new flow rate

  params.set<double>("Inlet integrated pressure", pressure);
#endif
}//DRT::ELEMENTS::Fluid3Surface::IntegratedPressureParameterCalculation


/*----------------------------------------------------------------------*
 |                                                        ismail 10/2010|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::CenterOfMassCalculation(
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> unitnormal(true);


  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }
#endif // D_ALE_BFLOW

  // first evaluate the area of the surface element
  params.set<double>("Area calculation",0.0);
  this->AreaCaculation(ele, params, discretization,lm);

  // get the surface element area
  const double elem_area = params.get<double>("Area calculation");

  LINALG::Matrix<(nsd_),1>  xyzGe(true);

  for (int i = 0; i< nsd_;i++)
  {
    //const IntegrationPoints2D  intpoints(gaussrule);
    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
      // Computation of the unit normal vector at the Gauss points
      // Computation of nurb specific stuff is not activated here
      double drs = 0.0;
      EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, &unitnormal);

      // global coordinates of gausspoint
      LINALG::Matrix<(nsd_),1>  coordgp(true);

      // determine coordinates of current Gauss point
      coordgp.Multiply(xyze,funct);

      //Compute elment center of gravity
      xyzGe(i) += intpoints.IP().qwgt[gpid]*coordgp(i)*drs;

    }  // end Gauss loop
    xyzGe(i) /= elem_area;
  }

  // Get the center of mass of the already calculate surface elements
  RCP<std::vector<double> > xyzG  = params.get<RCP<std::vector<double> > >("center of mass");

  RCP<std::vector<double> > normal  = params.get<RCP<std::vector<double> > >("normal");

  // Get the area of the of the already calculate surface elements
  double area = params.get<double>("total area");

  for (int i = 0; i<nsd_;i++)
  {
    (*xyzG)  [i] = ((*xyzG)[i]*area   + xyzGe(i)     *elem_area)/(area+elem_area);
    (*normal)[i] = ((*normal)[i]*area + unitnormal(i)*elem_area)/(area+elem_area);
  }

  // set new center of mass
  params.set("total area", area+elem_area);

#endif
}//DRT::ELEMENTS::Fluid3Surface::CenterOfMassCalculation



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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> unitnormal(true);

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");

  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  // allocate velocity vector
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode=0;inode<bdrynen_;inode++)
  {
    for (int idim=0; idim< nsd_; idim++)
    {
      evelnp(idim,inode) = myvelnp[idim+(inode*numdofpernode_)];
    }
  }

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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
    for (int inode=0;inode<bdrynen_;++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, &unitnormal);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    //compute flowrate at gauss point
    LINALG::Matrix<nsd_,1> velint(true);
    velint.Multiply(evelnp,funct);

    // flowrate = uint o normal
    const double flowrate = velint.Dot(unitnormal);

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode=0;inode<bdrynen_;++inode)
    {
      elevec1[inode*numdofpernode_] -= funct(inode)* fac_drs * flowrate; //
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
  // This function is only implemented for 3D
  if(bdrynsd_!=2)
    dserror("FlowRateDeriv is only implemented for 3D!");

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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of Fluid3Boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
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
  LINALG::Matrix<nsd_,bdrynen_> evelnp(true);

  for (int inode=0; inode<bdrynen_; ++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points is not activated here
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, NULL);
    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(deriv,xyze);
    normal(0,0) = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    normal(1,0) = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    normal(2,0) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

    //-------------------------------------------------------------------
    //  Q
    blitz::Array<double,1> u(3);
    u = 0.;
    for (int dim=0;dim<3;++dim)
      for (int node=0;node<bdrynen_;++node)
        u(dim) += funct(node) * evelnp(dim,node);

    for(int dim=0;dim<3;++dim)
      elevec3[0] += u(dim) * normal(dim,0) * fac;

    if (params.get<bool>("flowrateonly", false)==false)
    {
      //-------------------------------------------------------------------
      // dQ/du
      for (int node=0;node<bdrynen_;++node)
      {
        for (int dim=0;dim<3;++dim)
          elevec1[node*numdofpernode_+dim] += funct(node) * normal(dim,0) * fac;
        elevec1[node*numdofpernode_+3] = 0.0;
      }

      //-------------------------------------------------------------------
      // dQ/dd

      // determine derivatives of surface normals wrt mesh displacements
      blitz::Array<double,2> normalderiv(3,bdrynen_*3);

      for (int node=0;node<bdrynen_;++node)
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

      for (int node=0;node<bdrynen_;++node)
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
          for (int node=0;node<bdrynen_;++node)
          {
            for (int dim=0;dim<3;++dim)
              elevec2[node*numdofpernode_+dim] -= 1.0/dt * funct(node) * normal(dim,0) * fac;
          }
        }
        else
          dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(du dd)

      for (int unode=0;unode<bdrynen_;++unode)
      {
        for (int udim=0;udim<numdofpernode_;++udim)
        {
          for (int nnode=0;nnode<bdrynen_;++nnode)
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
      blitz::Array<double,3> normalderiv2(3,bdrynen_*3,bdrynen_*3);
      normalderiv2 = 0.;

      for (int node1=0;node1<bdrynen_;++node1)
      {
        for (int node2=0;node2<bdrynen_;++node2)
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

      for (int node1=0;node1<bdrynen_;++node1)
      {
        for (int dim1=0;dim1<numdofpernode_;++dim1)
        {
          for (int node2=0;node2<bdrynen_;++node2)
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
          for (int node1=0;node1<bdrynen_;++node1)
          {
            for (int dim1=0;dim1<3;++dim1)
            {
              for (int node2=0;node2<bdrynen_;++node2)
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
  LINALG::Matrix<bdrynen_,1> funct(true);
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv(true);

  // global node coordinates
  LINALG::Matrix<nsd_,bdrynen_> xyze(true);

  // coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi(true);

  // normal vector
  LINALG::Matrix<nsd_,1> unitnormal(true);

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
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze);

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
    for (int inode=0;inode<bdrynen_;++inode)
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
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    double drs = 0.0;
    EvalShapeFuncAndIntFac(intpoints, gpid, xyze, NULL, NULL, xsi, funct, deriv, drs, &unitnormal);
    const double fac_drs = intpoints.IP().qwgt[gpid]*drs;

    const double fac_drs_thsl_pres_inve = fac_drs * thsl * pressure * invdensity;

    for (int inode=0;inode<bdrynen_;++inode)
      for(int idim=0;idim<nsd_;++idim)
        // inward pointing normal of unit length
        elevec1[inode*numdofpernode_+idim] += funct(inode) * fac_drs_thsl_pres_inve * (-unitnormal(idim));
  }
  return;
} //DRT::ELEMENTS::Fluid3Surface::ImpedanceIntegration


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<bdrynsd_>&  intpoints,
    const int                                         gpid,
    const LINALG::Matrix<nsd_,bdrynen_>&                   xyze,
    const std::vector<Epetra_SerialDenseVector>*      myknots,
    const Epetra_SerialDenseVector*                   weights,
    LINALG::Matrix<bdrynsd_,1>&                       xsi,
    LINALG::Matrix<bdrynen_,1>&                            funct,
    LINALG::Matrix<bdrynsd_,bdrynen_>&                     deriv,
    double&                                           drs,
    LINALG::Matrix<nsd_,1>*                           unitnormal
)
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
  // only for NURBS!!!
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
  LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);

  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze,deriv,metrictensor, drs, unitnormal);

  return;
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
  // TODO: Check function 1D / 2D for Nurbs
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
  for (int inode=0; inode<bdrynen_; ++inode)
  {
    DRT::NURBS::ControlPoint* cp
      =
      dynamic_cast<DRT::NURBS::ControlPoint* > (ele->Nodes()[inode]);

    weights(inode) = cp->W();
  }
  return zero_size;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bndydistype,
          DRT::Element::DiscretizationType pdistype>
   void  DRT::ELEMENTS::Fluid3BoundaryImpl<distype>::MixHybDirichlet(
     DRT::ELEMENTS::Fluid3Boundary*  surfele,
     ParameterList&                  params,
     DRT::Discretization&            discretization,
     vector<int>&                    lm,
     Epetra_SerialDenseMatrix&       elemat_epetra,
     Epetra_SerialDenseVector&       elevec_epetra)
{
  // check whether we have a generalized-alpha time-integration scheme
  // and give a warning
  const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");
  if (is_genalpha) dserror("Mixed-hybrid formulation not yet available for af-gen-alpha!");

  //--------------------------------------------------
  // time integration related parameters
  const double timefac= params.get<double>("timefac");

  //--------------------------------------------------
  // get my parent element
  DRT::Element* parent=surfele->ParentElement();

  // get the required material information
  RCP<MAT::Material> mat = parent->Material();
  double density  =0.0;
  double viscosity=0.0;

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    density   = actmat->Density();
    viscosity = actmat->Viscosity();
  }
  else
    dserror("other material expected but got type %d", mat->MaterialType());

  /// number of parentnodes
  static const int piel    = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  /// number of surfacenodes
  static const int siel    = DRT::UTILS::DisTypeToNumNodePerEle<bndydistype>::numNodePerElement;

  /// number of spatial dimensions
  static const int nsd     = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  static const int bndynsd = DRT::UTILS::DisTypeToDim<bndydistype>::dim;

  // number of internal stress dofs is equivalent to number of second derivatives
  static const int numstressdof_= DRT::UTILS::DisTypeToNumDeriv2<pdistype>::numderiv2;

  // --------------------------------------------------
  // Reshape element matrices and vectors and init to zero, construct views
  const int peledim = (nsd +1)*piel;

  elemat_epetra.Shape(peledim,peledim);
  elevec_epetra.Size (peledim);

  LINALG::Matrix<peledim,peledim> elemat(elemat_epetra.A(),true);
  LINALG::Matrix<peledim,      1> elevec(elevec_epetra.A(),true);

  //--------------------------------------------------
  // scaling for constitutive law
  const double scaling=1.0/(2.0*viscosity);
  //  double rhsscaling=1./(2.0*viscosity);

  //  const double scaling=1.0;
  const double rhsscaling=1.0;

  //--------------------------------------------------
  // get the condition information
  RefCountPtr<DRT::Condition> hixhybdbc_cond
    =
    params.get<RefCountPtr<DRT::Condition> >("condition");

  // get value for boundary condition
  const vector<double>* val = (*hixhybdbc_cond).Get<vector<double> >("val");

  // find out whether we will use a time curve
  const double time = params.get<double>("total time");

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = (*hixhybdbc_cond).Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and switches from the condition
  // (assumed to be constant on element boundary)
  const vector<int>* functions = (*hixhybdbc_cond).Get<vector<int> >   ("funct");


  LINALG::Matrix<nsd,1> u_dirich(true);

  for(int rr=0;rr<nsd;++rr)
  {
    u_dirich(rr)=(*val)[rr]*curvefac;
  }

  // --------------------------------------------------
  // Extra matrices

  // for r / sigma: indices ordered according to
  //
  //      0    1    2    (2D)
  //     11 , 22 , 12
  //

  //      0    1    2    3    4   5  (3D)
  //     11 , 22 , 33 , 12 , 13 ,23

  // for volume integrals

  LINALG::Matrix<numstressdof_*piel,              piel> mat_r_p(true);
  LINALG::Matrix<numstressdof_*piel,numstressdof_*piel> mat_r_sigma(true);
  LINALG::Matrix<numstressdof_*piel,nsd          *piel> mat_r_epsu(true);

  // for boundary integrals

  LINALG::Matrix<nsd          *piel,numstressdof_*piel> mat_v_sigma_o_n(true);
  LINALG::Matrix<numstressdof_*piel,nsd          *piel> mat_r_o_n_u(true);

  // rearranging and computational arrays

  LINALG::Matrix<numstressdof_*piel,(nsd+1)*piel>       mat_r_up_block(true);
  LINALG::Matrix<numstressdof_*piel,numstressdof_*piel> inv_r_sigma(true);


  // --------------------------------------------------
  // Extra vectors

  // for volume integrals

  LINALG::Matrix<numstressdof_*piel,                 1> vec_r_p(true);
  LINALG::Matrix<numstressdof_*piel,                 1> vec_r_epsu(true);

  // for boundary integrals
  LINALG::Matrix<numstressdof_*piel,                 1> vec_r_o_n_u_minus_g(true);

  //--------------------------------------------------
  // get parent elements location vector and ownerships

  // the vectors have been allocated outside in
  // EvaluateConditionUsingParentData
  RCP<vector<int> > plm
    =
    params.get<RCP<vector<int> > >("plm");
  RCP<vector<int> > plmowner
    =
    params.get<RCP<vector<int> > >("plmowner");
  RCP<vector<int> > plmstride
    =
    params.get<RCP<vector<int> > >("plmstride");

  parent->LocationVector(discretization,*plm,*plmowner,*plmstride);

  // extract local velocities and pressure from the global vectors
  LINALG::Matrix<nsd ,piel>    pevel (true);
  LINALG::Matrix<piel,   1>    pepres(true);

  RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
  if (vel==null) dserror("Cannot get state vector 'u and p (trial)'");

  // extract local node values for pressure and velocities from global vectors
  vector<double> mypvel((*plm).size());
  DRT::UTILS::ExtractMyValues(*vel,mypvel,*plm);

  for (int inode=0;inode<piel;++inode)
  {
    for (int idim=0; idim<nsd ; ++idim)
    {
      pevel(idim,inode) = mypvel[(nsd +1)*inode+idim];
    }
    pepres(inode) = mypvel[(nsd +1)*inode+nsd ];
  }


  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 1: Gaussloop for volume integrals of parent element
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  {
    // allocate vector for shape functions and matrix for derivatives
    LINALG::Matrix<piel,1>       pfunct(true);
    LINALG::Matrix<nsd ,piel>    pderiv(true);

    // get local node coordinates
    LINALG::Matrix<nsd ,piel>    pxyze(true);
    GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

    //--------------------------------------------------
    // Gaussian integration points
    const DRT::UTILS::IntPointsAndWeights<nsd>
      pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

    //--------------------------------------------------
    // vectors/scalars for Gausspoint values

    // velocity at gausspoint
    LINALG::Matrix<nsd ,1>       pvelint(true);
    // velocity derivatives at gausspoint
    LINALG::Matrix<nsd ,nsd >    pvderxy(true);
    // pressure at gausspoint
    double                       ppressure=0.0;

    // global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd ,piel>    pderxy(true);
    // transposed jacobian "dx/ds"
    LINALG::Matrix<nsd ,nsd >    pxjm(true);
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd ,nsd >    pxji(true);

    LINALG::Matrix<nsd ,   1>    pxsi(true);

    //--------------------------------------------------
    // the actual loop
    for (int iquad=0; iquad<pintpoints.IP().nquad; ++iquad)
    {
      // coordinates of the current integration point
      const double* gpcoord = (pintpoints.IP().qxg)[iquad];
      for (int idim=0;idim<nsd ;idim++)
      {
        pxsi(idim) = gpcoord[idim];
      }

      // get parent elements shape functions
      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

      // get Jacobian matrix and determinant
      // actually compute its transpose....
      /*
        +-            -+ T      +-            -+
        | dx   dx   dx |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dr   dr   dr |
        |              |        |              |
        | dy   dy   dy |        | dx   dy   dz |
        | --   --   -- |   =    | --   --   -- |
        | dr   ds   dt |        | ds   ds   ds |
        |              |        |              |
        | dz   dz   dz |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dt   dt   dt |
        +-            -+        +-            -+
      */
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);

      if (det < 1E-16)
        dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", parent->Id(), det);

      // compute integration factor
      const double fac = pintpoints.IP().qwgt[iquad]*det*timefac;

      // compute global first derivates
      pderxy.Multiply(pxji,pderiv);

      // interpolate to gausspoint
      pvelint.Multiply(pevel,pfunct);

      // get velocity derivatives at integration point
      pvderxy.MultiplyNT(pevel,pderxy);

      // interpolate pressure to gausspoint
      ppressure = pfunct.Dot(pepres);

      /*
                     /          \
                    |  h       h |
                  - | r : sigma  |
                    |            |
                     \          / Omega
      */
      for(int A=0;A<piel;++A)
      {
        for(int B=0;B<piel;++B)
        {
          for(int i=0;i<nsd ;++i)
          {
            mat_r_sigma(A*numstressdof_+i,B*numstressdof_+i)-=fac*pfunct(A)*pfunct(B)*scaling;
          }
          for(int i=nsd;i<numstressdof_;++i)
          {
            mat_r_sigma(A*numstressdof_+i,B*numstressdof_+i)-=fac*pfunct(A)*pfunct(B)*2.0*scaling;
          }
        }
      }

      /*
                     /         \
                    |  h   h    |
                  - | r : p * I |
                    |           |
                     \         / Omega
      */
      for(int A=0;A<piel;++A)
      {
        for(int B=0;B<piel;++B)
        {
          for(int i=0;i<nsd ;++i)
          {
            mat_r_p(A*numstressdof_+i,B)-=fac*pfunct(A)*pfunct(B)*scaling;
          }
        }
      }

      for(int A=0;A<piel;++A)
      {
        for(int i=0;i<nsd ;++i)
        {
          vec_r_p(A*numstressdof_+i)-=fac*ppressure*pfunct(A)*scaling;
        }
      }


      /*
                              /              \
                             |  h       / h\  |
                  + 2 * nu * | r : eps | u  | |
                             |          \  /  |
                              \              / Omega
      */
      const double viscfac=fac*2.0*viscosity*scaling;
      if(nsd ==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_epsu(A*numstressdof_  ,B*nsd   )+=viscfac*pfunct(A)*pderxy(0,B);
            mat_r_epsu(A*numstressdof_+1,B*nsd +1)+=viscfac*pfunct(A)*pderxy(1,B);

            mat_r_epsu(A*numstressdof_+2,B*nsd   )+=viscfac*pfunct(A)*pderxy(1,B);
            mat_r_epsu(A*numstressdof_+2,B*nsd +1)+=viscfac*pfunct(A)*pderxy(0,B);
          }
        }
      }
      else if(nsd ==3)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_epsu(A*numstressdof_  ,B*nsd   )+=viscfac*pfunct(A)*pderxy(0,B);
            mat_r_epsu(A*numstressdof_+1,B*nsd +1)+=viscfac*pfunct(A)*pderxy(1,B);
            mat_r_epsu(A*numstressdof_+2,B*nsd +2)+=viscfac*pfunct(A)*pderxy(2,B);

            mat_r_epsu(A*numstressdof_+3,B*nsd   )+=viscfac*pfunct(A)*pderxy(1,B);
            mat_r_epsu(A*numstressdof_+3,B*nsd +1)+=viscfac*pfunct(A)*pderxy(0,B);

            mat_r_epsu(A*numstressdof_+4,B*nsd   )+=viscfac*pfunct(A)*pderxy(2,B);
            mat_r_epsu(A*numstressdof_+4,B*nsd +2)+=viscfac*pfunct(A)*pderxy(0,B);

            mat_r_epsu(A*numstressdof_+5,B*nsd +1)+=viscfac*pfunct(A)*pderxy(2,B);
            mat_r_epsu(A*numstressdof_+5,B*nsd +2)+=viscfac*pfunct(A)*pderxy(1,B);
          }
        }
      }

      if(nsd ==2)
      {
        for(int A=0;A<piel;++A)
        {
          vec_r_epsu(A*numstressdof_  )+=fac*pfunct(A)*pvderxy(0,0)*2.0*viscosity*scaling;
          vec_r_epsu(A*numstressdof_+1)+=fac*pfunct(A)*pvderxy(1,1)*2.0*viscosity*scaling;

          vec_r_epsu(A*numstressdof_+2)+=fac*pfunct(A)*(pvderxy(0,1)+pvderxy(1,0))*2.0*viscosity*scaling;
        }
      }
      else if(nsd ==3)
      {
        for(int A=0;A<piel;++A)
        {
          vec_r_epsu(A*numstressdof_  )+=fac*pfunct(A)*pvderxy(0,0)*2.0*viscosity*scaling;
          vec_r_epsu(A*numstressdof_+1)+=fac*pfunct(A)*pvderxy(1,1)*2.0*viscosity*scaling;
          vec_r_epsu(A*numstressdof_+2)+=fac*pfunct(A)*pvderxy(2,2)*2.0*viscosity*scaling;

          vec_r_epsu(A*numstressdof_+3)+=fac*pfunct(A)*(pvderxy(0,1)+pvderxy(1,0))*2.0*viscosity*scaling;
          vec_r_epsu(A*numstressdof_+4)+=fac*pfunct(A)*(pvderxy(0,2)+pvderxy(2,0))*2.0*viscosity*scaling;
          vec_r_epsu(A*numstressdof_+5)+=fac*pfunct(A)*(pvderxy(1,2)+pvderxy(2,1))*2.0*viscosity*scaling;
        }
      }
    }
  }

  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 2: Gaussloop for line integrals of boundary element
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  {
    // allocate vector/matrix for shape functions and derivatives
    LINALG::Matrix<siel   ,1>     funct(true);
    LINALG::Matrix<bndynsd,siel>  deriv(true);

    // allocate vector for parents shape functions and matrix for derivatives
    LINALG::Matrix<piel,1>        pfunct(true);
    LINALG::Matrix<nsd ,piel>     pderiv(true);


    // get local node coordinates
    LINALG::Matrix<nsd ,siel>    xyze(true);
    GEO::fillInitialPositionArray<bndydistype,nsd ,LINALG::Matrix<nsd ,siel> >(surfele,xyze);

    // get local node coordinates
    LINALG::Matrix<nsd ,piel>    pxyze(true);
    GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

    //--------------------------------------------------
    // Gaussian integration points
    const DRT::UTILS::IntPointsAndWeights<bndynsd>
      intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<bndydistype>::rule);

    const DRT::UTILS::IntPointsAndWeights<nsd>
      pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

    // coordinates of current integration point in reference coordinates
    LINALG::Matrix<bndynsd,1>    xsi(true);
    LINALG::Matrix<nsd    ,1>    pxsi(true);


    Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,nsd);

    {
      Epetra_SerialDenseMatrix gps(intpoints.IP().nquad,bndynsd);


      // coordinates of the current integration point
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        const double* gpcoord = (intpoints.IP().qxg)[iquad];

        for (int idim=0;idim<bndynsd ;idim++)
        {
          gps(iquad,idim) = gpcoord[idim];
        }
      }
      BoundaryGPToParentGP<nsd>(pqxg     ,
                                gps,
                                pdistype   ,
                                bndydistype,
                                surfele->SurfaceNumber());
    }


    //--------------------------------------------------
    // vectors/scalars for Gausspoint values

    // the element's normal vector
    LINALG::Matrix<nsd ,1>       unitnormal(true);
    // velocity at gausspoint
    LINALG::Matrix<nsd ,1>       velint(true);

    // transposed jacobian "dx/ds"
    LINALG::Matrix<nsd ,nsd >    xjm(true);
    // inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd ,nsd >    xji(true);

    // transposed jacobian "dx/ds" for parent
    LINALG::Matrix<nsd ,nsd >    pxjm(true);
    // inverse of transposed jacobian "ds/dx" for parent
    LINALG::Matrix<nsd ,nsd >    pxji(true);


    //--------------------------------------------------
    // the actual loop
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {

       // coordinates of the current integration point
      const double* gpcoord = (intpoints.IP().qxg)[iquad];
      for (int idim=0;idim<bndynsd ;idim++)
      {
        xsi(idim) = gpcoord[idim];
      }

      DRT::UTILS::shape_function       <bndydistype>(xsi,funct);
      DRT::UTILS::shape_function_deriv1<bndydistype>(xsi,deriv);

      for (int idim=0;idim<nsd ;idim++)
      {
        pxsi(idim) = pqxg(iquad,idim);
      }

      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

      double drs=0.0;

      // compute measure tensor for surface element and the infinitesimal
      // area element drs for the integration
      LINALG::Matrix<bndynsd,bndynsd> metrictensor(true);

      DRT::UTILS::ComputeMetricTensorForBoundaryEle<bndydistype>(xyze,deriv,metrictensor,drs,&unitnormal);

      // compute integration factor
      const double fac = intpoints.IP().qwgt[iquad]*drs*timefac;
#if 0

      // get Jacobian matrix and determinant
      // actually compute its transpose....
      /*
        +-            -+ T      +-            -+
        | dx   dx   dx |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dr   dr   dr |
        |              |        |              |
        | dy   dy   dy |        | dx   dy   dz |
        | --   --   -- |   =    | --   --   -- |
        | dr   ds   dt |        | ds   ds   ds |
        |              |        |              |
        | dz   dz   dz |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dt   dt   dt |
        +-            -+        +-            -+
      */
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);

      if (det < 1E-16)
        dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", parent->Id(), det);

      //-----------------------------------------------------
      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      LINALG::Matrix<nsd,nsd> G;

      for (int nn=0;nn<nsd;++nn)
      {
        for (int rr=0;rr<nsd;++rr)
        {
          G(nn,rr) = pxji(nn,0)*pxji(rr,0);
          for (int mm=1;mm<nsd;++mm)
          {
            G(nn,rr) += pxji(nn,mm)*pxji(rr,mm);
          }
        }
      }

      //
      //                           2.0
      //             h  = ---------------------
      //              b        +-------------+
      //                      / /  T       \ |
      //                   \ / |  n * G * n |
      //                    +   \          /
      //

      double nGn=0;
      for (int nn=0;nn<nsd;++nn)
      {
        for (int rr=0;rr<nsd;++rr)
        {
          nGn+=unitnormal(rr)*G(rr,nn)*unitnormal(nn);
        }
      }
      const double h =2.0/sqrt(nGn);

      //      rhsscaling=1./h;
#endif

      // interpolate to gausspoint
      velint.Multiply(pevel,pfunct);

      /*
                              /              \
                             |  h       h     |
                           - | v , sigma  o n |
                             |                |
                              \              / Gamma
      */
      if(nsd ==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_  )-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+1)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);

            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_+2)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+2)-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
          }
        }
      }
      else if(nsd ==3)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_  )-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_+3)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_v_sigma_o_n(A*nsd  ,B*numstressdof_+4)-=fac*pfunct(A)*pfunct(B)*unitnormal(2);

            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+3)-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+1)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_v_sigma_o_n(A*nsd+1,B*numstressdof_+5)-=fac*pfunct(A)*pfunct(B)*unitnormal(2);

            mat_v_sigma_o_n(A*nsd+2,B*numstressdof_+4)-=fac*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_v_sigma_o_n(A*nsd+2,B*numstressdof_+5)-=fac*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_v_sigma_o_n(A*nsd+2,B*numstressdof_+2)-=fac*pfunct(A)*pfunct(B)*unitnormal(2);
          }
        }
      }

      /*
                     /          \
                    |  h       h |
                  - | r o n , u  |
                    |            |
                     \          / Gamma
      */
      const double invert=1.;

      if(nsd==2)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_o_n_u(A*numstressdof_  ,B*nsd  )-=invert*fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+1)-=invert*fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(1);

            mat_r_o_n_u(A*numstressdof_+2,B*nsd  )-=invert*fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+1)-=invert*fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(0);
          }
        }
      }
      else if(nsd==3)
      {
        for(int A=0;A<piel;++A)
        {
          for(int B=0;B<piel;++B)
          {
            mat_r_o_n_u(A*numstressdof_  ,B*nsd  )-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(0);
            mat_r_o_n_u(A*numstressdof_+1,B*nsd+1)-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_r_o_n_u(A*numstressdof_+2,B*nsd+2)-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(2);

            mat_r_o_n_u(A*numstressdof_+3,B*nsd  )-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(1);
            mat_r_o_n_u(A*numstressdof_+3,B*nsd+1)-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(0);

            mat_r_o_n_u(A*numstressdof_+4,B*nsd  )-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(2);
            mat_r_o_n_u(A*numstressdof_+4,B*nsd+2)-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(0);

            mat_r_o_n_u(A*numstressdof_+5,B*nsd+1)-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(2);
            mat_r_o_n_u(A*numstressdof_+5,B*nsd+2)-=fac*rhsscaling*pfunct(A)*pfunct(B)*unitnormal(1);
          }
        }
      }

      // ------------------------------------------------
      // factor given by spatial function
      LINALG::Matrix<nsd,1> functionfac(true);
      for(int i=0;i<nsd;++i)
      {
        functionfac(i)= 1.0;
      }

      // determine coordinates of current Gauss point
      LINALG::Matrix<3,1> coordgp(true);

      for (int A=0;A<siel;++A)
      {
        for(int j=0;j<nsd;++j)
        {
          coordgp(j)+=xyze(j,A)*funct(A);
        }
      }


#if 0

      // determine coordinates of current Gauss point
      LINALG::Matrix<3,1> check(true);
      LINALG::Matrix<3,1> diff(true);

      for (int A=0;A<piel;++A)
      {
        for(int j=0;j<nsd;++j)
        {
          check(j)+=pxyze(j,A)*pfunct(A);
        }
      }

      diff=check;
      diff-=coordgp;

      const double norm=diff.Norm2();

      if(norm>1e-9)
      {
        for(int j=0;j<nsd;++j)
        {
          printf("%12.5e %12.5e\n",check(j),coordgp(j));
        }

        dserror("Gausspoint matching error %12.5e\n",norm);
      }

#endif

      int functnum = -1;

      for(int dim=0;dim<nsd;++dim)
      {
        // factor given by spatial function
        if (functions)
        {
          functnum = (*functions)[dim];
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functionfac(dim) = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgp.A(),0.0,NULL);
          }
          else
          {
            functionfac(dim) = 1.0;
          }
        }
      }

      LINALG::Matrix<nsd,1> delta_vel(true);

      for(int rr=0;rr<nsd;++rr)
      {
        delta_vel(rr)=velint(rr)-u_dirich(rr)*functionfac(rr);
      }

      if(nsd==2)
      {
        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g(A*numstressdof_  )-=invert*fac*rhsscaling*pfunct(A)*unitnormal(0)*delta_vel(0);
          vec_r_o_n_u_minus_g(A*numstressdof_+1)-=invert*fac*rhsscaling*pfunct(A)*unitnormal(1)*delta_vel(1);

          vec_r_o_n_u_minus_g(A*numstressdof_+2)-=invert*fac*rhsscaling*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
        }
      }
      else if(nsd==3)
      {
        for(int A=0;A<piel;++A)
        {
          vec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*rhsscaling*pfunct(A)*unitnormal(0)*delta_vel(0);
          vec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*rhsscaling*pfunct(A)*unitnormal(1)*delta_vel(1);
          vec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*rhsscaling*pfunct(A)*unitnormal(2)*delta_vel(2);

          vec_r_o_n_u_minus_g(A*numstressdof_+3)-=fac*rhsscaling*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
          vec_r_o_n_u_minus_g(A*numstressdof_+4)-=fac*rhsscaling*pfunct(A)*(unitnormal(2)*delta_vel(0)+unitnormal(0)*delta_vel(2));
          vec_r_o_n_u_minus_g(A*numstressdof_+5)-=fac*rhsscaling*pfunct(A)*(unitnormal(2)*delta_vel(1)+unitnormal(1)*delta_vel(2));
        }
      }



#if 1

      //--------------------------------------------------
      // adjoint consistency term, pressure/continuity part
      /*
      //
      //
      //             /              \
      //            |                |
      //          - |  q , Dacc * n  |
      //            |                |
      //             \              / boundaryele
      //
      */
      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          for(int i=0;i<nsd;++i)
          {
            elemat(vi*(nsd+1)+nsd,ui*(nsd+1)+i) -= fac*timefac*pfunct(vi)*pfunct(ui)*unitnormal(i);
          }
        }
      }

      /*
      // factor: 1.0
      //
      //             /                       \
      //            |       / n+1     \       |
      //          + |  q , | u   - u   | * n  |
      //            |       \ (i)   B /       |
      //             \                       / boundaryele
      //
      */
      for (int vi=0; vi<piel; ++vi)
      {
        for(int i=0;i<nsd;++i)
        {
          elevec(vi*(nsd+1)+nsd) += fac*timefac*pfunct(vi)*delta_vel(i)*unitnormal(i);
        }
      }
#endif





#if 0
      const double penalty=4*1000*viscosity/h;

      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          const double temp=fac*penalty*timefac*pfunct(ui)*pfunct(vi);
          for(int i=0;i<nsd;++i)
          {
            elemat(vi*(nsd+1) + i,ui*(nsd+1) + i) +=temp;
          }
        }
      }

      /*
      // factor: nu*Cb/h
      //
      //    /                \
      //   |        n+af      |
      // + |   w , u    - u   |
      //   |               b  |
      //    \                / boundaryele
      //
      */

      for (int vi=0; vi<piel; ++vi)
      {
        for(int i=0;i<nsd;++i)
        {
          elevec(vi*(nsd+1) + i) -= fac*penalty*timefac*pfunct(vi)*delta_vel(i);
        }
      }
#endif

    }
  }
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 3: Local condensation (Matrix inversion etc)
    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/

  // rearrange to pattern uvwp uvwp ...
  for(int A=0;A<piel;++A)
  {
    for(int B=0;B<piel;++B)
    {
      for(int i=0;i<numstressdof_;++i)
      {
        for(int j=0;j<nsd;++j)
        {
          mat_r_up_block(A*numstressdof_+i,B*(nsd+1)+j)+=mat_r_epsu (A*numstressdof_+i,B*nsd+j);
          mat_r_up_block(A*numstressdof_+i,B*(nsd+1)+j)+=mat_r_o_n_u(A*numstressdof_+i,B*nsd+j);
        }
        mat_r_up_block(A*numstressdof_+i,B*(nsd+1)+nsd)+=mat_r_p(A*numstressdof_+i,B);
      }
    }
  }

  // matrix inversion of stress-stress block
  inv_r_sigma=mat_r_sigma;

  LINALG::FixedSizeSerialDenseSolver<numstressdof_*piel,numstressdof_*piel> solver;

  solver.SetMatrix(inv_r_sigma);
  solver.Invert();

  // computation of matrix-matrix and matrix vector products, local assembly
  for(int A=0;A<piel;++A)
  {
    for(int i=0;i<nsd;++i)
    {
      for(int B=0;B<piel;++B)
      {
        for(int rr=0;rr<numstressdof_*piel;++rr)
        {
          for(int mm=0;mm<numstressdof_*piel;++mm)
          {
            for(int j=0;j<nsd+1;++j)
            {
              elemat(A*(nsd+1)+i,B*(nsd+1)+j)-=mat_v_sigma_o_n(A*nsd+i,rr)*inv_r_sigma(rr,mm)*mat_r_up_block(mm,B*(nsd+1)+j);
            }
          }
        }
      }
    }
  }

  for(int A=0;A<piel;++A)
  {
    for(int i=0;i<nsd;++i)
    {
      for(int rr=0;rr<numstressdof_*piel;++rr)
      {
        for(int mm=0;mm<numstressdof_*piel;++mm)
        {
          elevec(A*(nsd+1)+i)-=mat_v_sigma_o_n(A*nsd+i,rr)*inv_r_sigma(rr,mm)*(-vec_r_o_n_u_minus_g(mm)-vec_r_epsu(mm)-vec_r_p(mm));
        }
      }
    }
  }

#if 0
  // --------------------------------------------------
  //
  //                       FDCHECK
  //
  // --------------------------------------------------
  for(int fd=0;fd<nsd+1;fd++)
  {
    // Extra vectors

    // for volume integrals

    LINALG::Matrix<numstressdof_*piel,                 1> FDvec_r_p(true);
    LINALG::Matrix<numstressdof_*piel,                 1> FDvec_r_epsu(true);

    // for boundary integrals
    LINALG::Matrix<numstressdof_*piel,                 1> FDvec_r_o_n_u_minus_g(true);

    for (int inode=0;inode<piel;++inode)
    {
      for (int idim=0; idim<nsd ; ++idim)
      {
        pevel(idim,inode) = mypvel[(nsd +1)*inode+idim];
      }
      pepres(inode) = mypvel[(nsd +1)*inode+nsd ];
    }

    if(fd==nsd)
    {
      pepres(0)+=1;
    }
    else
    {
      pevel(fd,0)+=1;
    }


    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
         PART 1: Gaussloop for volume integrals of parent element
      <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    {
      // allocate vector for shape functions and matrix for derivatives
      LINALG::Matrix<piel,1>       pfunct(true);
      LINALG::Matrix<nsd ,piel>    pderiv(true);

      // get local node coordinates
      LINALG::Matrix<nsd ,piel>    pxyze(true);
      GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

      //--------------------------------------------------
      // Gaussian integration points
      const DRT::UTILS::IntPointsAndWeights<nsd>
        pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

      //--------------------------------------------------
      // vectors/scalars for Gausspoint values

      // velocity at gausspoint
      LINALG::Matrix<nsd ,1>       pvelint(true);
      // velocity derivatives at gausspoint
      LINALG::Matrix<nsd ,nsd >    pvderxy(true);
      // pressure at gausspoint
      double                       ppressure=0.0;

      // global derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<nsd ,piel>    pderxy(true);
      // transposed jacobian "dx/ds"
      LINALG::Matrix<nsd ,nsd >    pxjm(true);
      // inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd ,nsd >    pxji(true);

      LINALG::Matrix<nsd ,   1>    pxsi(true);

      //--------------------------------------------------
      // the actual loop
      for (int iquad=0; iquad<pintpoints.IP().nquad; ++iquad)
      {
        // coordinates of the current integration point
        const double* gpcoord = (pintpoints.IP().qxg)[iquad];
        for (int idim=0;idim<nsd ;idim++)
        {
          pxsi(idim) = gpcoord[idim];
        }

        // get parent elements shape functions
        DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
        DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

        // get Jacobian matrix and determinant
        // actually compute its transpose....
        /*
        +-            -+ T      +-            -+
        | dx   dx   dx |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dr   dr   dr |
        |              |        |              |
        | dy   dy   dy |        | dx   dy   dz |
        | --   --   -- |   =    | --   --   -- |
        | dr   ds   dt |        | ds   ds   ds |
        |              |        |              |
        | dz   dz   dz |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dt   dt   dt |
        +-            -+        +-            -+
        */
        pxjm.MultiplyNT(pderiv,pxyze);
        const double det = pxji.Invert(pxjm);

        if (det < 1E-16)
        dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", parent->Id(), det);

        // compute integration factor
        const double fac = pintpoints.IP().qwgt[iquad]*det*timefac;

        // compute global first derivates
        pderxy.Multiply(pxji,pderiv);

        // interpolate to gausspoint
        pvelint.Multiply(pevel,pfunct);

        // get velocity derivatives at integration point
        pvderxy.MultiplyNT(pevel,pderxy);

        // interpolate pressure to gausspoint
        ppressure = pfunct.Dot(pepres);

        /*
                     /         \
                    |  h   h    |
                  - | r : p * I |
                    |           |
                     \         / Omega
        */

        for(int A=0;A<piel;++A)
        {
          for(int i=0;i<nsd ;++i)
          {
            FDvec_r_p(A*numstressdof_+i)-=fac*ppressure*pfunct(A)/(2.0*viscosity);
          }
        }


        /*
                              /              \
                             |  h       / h\  |
                  + 2 * nu * | r : eps | u  | |
                             |          \  /  |
                              \              / Omega
        */

        if(nsd ==2)
        {
          for(int A=0;A<piel;++A)
          {
            FDvec_r_epsu(A*numstressdof_  )+=fac*pfunct(A)*pvderxy(0,0);
            FDvec_r_epsu(A*numstressdof_+1)+=fac*pfunct(A)*pvderxy(1,1);

            FDvec_r_epsu(A*numstressdof_+2)+=fac*pfunct(A)*(pvderxy(0,1)+pvderxy(1,0));
          }
        }
        else if(nsd ==3)
        {
          for(int A=0;A<piel;++A)
          {
            FDvec_r_epsu(A*numstressdof_  )+=fac*pfunct(A)*pvderxy(0,0);
            FDvec_r_epsu(A*numstressdof_+1)+=fac*pfunct(A)*pvderxy(1,1);
            FDvec_r_epsu(A*numstressdof_+2)+=fac*pfunct(A)*pvderxy(2,2);

            FDvec_r_epsu(A*numstressdof_+3)+=fac*pfunct(A)*(pvderxy(0,1)+pvderxy(1,0));
            FDvec_r_epsu(A*numstressdof_+4)+=fac*pfunct(A)*(pvderxy(0,2)+pvderxy(2,0));
            FDvec_r_epsu(A*numstressdof_+5)+=fac*pfunct(A)*(pvderxy(1,2)+pvderxy(2,1));
          }
        }
      }
    }

    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      PART 2: Gaussloop for line integrals of boundary element
      <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    {
      // allocate vector/matrix for shape functions and derivatives
      LINALG::Matrix<siel   ,1>     funct(true);
      LINALG::Matrix<bndynsd,siel>  deriv(true);

      // allocate vector for parents shape functions and matrix for derivatives
      LINALG::Matrix<piel,1>        pfunct(true);
      LINALG::Matrix<nsd ,piel>     pderiv(true);


      // get local node coordinates
      LINALG::Matrix<nsd ,siel> xyze(true);
      GEO::fillInitialPositionArray<bndydistype,nsd ,LINALG::Matrix<nsd ,siel> >(surfele,xyze);

      // get local node coordinates
      LINALG::Matrix<nsd ,piel>    pxyze(true);
      GEO::fillInitialPositionArray<pdistype,nsd ,LINALG::Matrix<nsd ,piel> >(parent,pxyze);

      //--------------------------------------------------
      // Gaussian integration points
      const DRT::UTILS::IntPointsAndWeights<bndynsd>
        intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<bndydistype>::rule);

      const DRT::UTILS::IntPointsAndWeights<nsd>
        pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);



      // coordinates of current integration point in reference coordinates
      LINALG::Matrix<bndynsd,1>   xsi(true);
      LINALG::Matrix<nsd    ,1>   pxsi(true);

      Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,nsd);

      BoundaryGPToParentGP(pqxg     ,
                           intpoints,
                           pdistype   ,
                           bndydistype,
                           surfele->SurfaceNumber());

      //--------------------------------------------------
      // vectors/scalars for Gausspoint values

      // the element's normal vector
      LINALG::Matrix<nsd ,1>       unitnormal(true);
      // velocity at gausspoint
      LINALG::Matrix<nsd ,1>       velint(true);

      // transposed jacobian "dx/ds"
      LINALG::Matrix<nsd ,nsd >    xjm(true);
      // inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd ,nsd >    xji(true);

      //--------------------------------------------------
      // the actual loop
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {

        // coordinates of the current integration point
        const double* gpcoord = (intpoints.IP().qxg)[iquad];
        for (int idim=0;idim<bndynsd ;idim++)
        {
          xsi(idim) = gpcoord[idim];
        }

        DRT::UTILS::shape_function       <bndydistype>(xsi,funct);
        DRT::UTILS::shape_function_deriv1<bndydistype>(xsi,deriv);

        for (int idim=0;idim<nsd ;idim++)
        {
          pxsi(idim) = pqxg(iquad,idim);
        }

        DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);

        double drs=0.0;

        // compute measure tensor for surface element and the infinitesimal
        // area element drs for the integration
        LINALG::Matrix<bndynsd,bndynsd> metrictensor(true);

        DRT::UTILS::ComputeMetricTensorForBoundaryEle<bndydistype>(xyze,deriv,metrictensor,drs,&unitnormal);

        // compute integration factor
        const double fac = intpoints.IP().qwgt[iquad]*drs*timefac;

        // interpolate to gausspoint
        velint.Multiply(pevel,pfunct);

        /*
                     /          \
                    |  h       h |
                  - | r o n , u  |
                    |            |
                     \          / Gamma
        */

        // ------------------------------------------------
        // factor given by spatial function
        LINALG::Matrix<nsd,1> functionfac(true);
        for(int i=0;i<nsd;++i)
        {
          functionfac(i)= 1.0;
        }

        // determine coordinates of current Gauss point
        LINALG::Matrix<3,1> coordgp(true);

        for (int A=0;A<siel;++A)
        {
          for(int j=0;j<nsd;++j)
          {
            coordgp(j)+=xyze(j,A)*funct(A);
          }
        }

        int functnum = -1;

        for(int dim=0;dim<nsd;++dim)
        {
          // factor given by spatial function
          if (functions)
          {
            functnum = (*functions)[dim];
            if (functnum>0)
            {
              // evaluate function at current gauss point
              functionfac(dim) = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgp.A(),0.0,NULL);
            }
            else
            {
              functionfac(dim) = 1.0;
            }
          }
        }

        LINALG::Matrix<nsd,1> delta_vel(true);

        for(int rr=0;rr<nsd;++rr)
        {
          delta_vel(rr)=velint(rr)-u_dirich(rr)*functionfac(rr);
        }

        if(nsd==2)
        {
          for(int A=0;A<piel;++A)
          {
            FDvec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*pfunct(A)*unitnormal(0)*delta_vel(0);
            FDvec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*pfunct(A)*unitnormal(1)*delta_vel(1);

            FDvec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
          }
        }
        else if(nsd==3)
        {
          for(int A=0;A<piel;++A)
          {
            FDvec_r_o_n_u_minus_g(A*numstressdof_  )-=fac*pfunct(A)*unitnormal(0)*delta_vel(0);
            FDvec_r_o_n_u_minus_g(A*numstressdof_+1)-=fac*pfunct(A)*unitnormal(1)*delta_vel(1);
            FDvec_r_o_n_u_minus_g(A*numstressdof_+2)-=fac*pfunct(A)*unitnormal(2)*delta_vel(2);

            FDvec_r_o_n_u_minus_g(A*numstressdof_+3)-=fac*pfunct(A)*(unitnormal(1)*delta_vel(0)+unitnormal(0)*delta_vel(1));
            FDvec_r_o_n_u_minus_g(A*numstressdof_+4)-=fac*pfunct(A)*(unitnormal(2)*delta_vel(0)+unitnormal(0)*delta_vel(2));
            FDvec_r_o_n_u_minus_g(A*numstressdof_+5)-=fac*pfunct(A)*(unitnormal(2)*delta_vel(1)+unitnormal(1)*delta_vel(2));
          }
        }
      }
    }

    if(fd==nsd)
    {
      for(int rr=0;rr<numstressdof_;++rr)
      {
        printf("FDvec_r_p(%2d,%2d)             %12.5e  %12.5e %12.5e  %12.5e\n",rr,0,mat_r_p(rr,0)  ,FDvec_r_p(rr)-vec_r_p(rr),FDvec_r_p(rr),vec_r_p(rr));
      }
    }
    else
    {
      for(int rr=0;rr<numstressdof_;++rr)
      {
        printf("FDvec_r_epsu(%2d,%2d)          %12.5e  %12.5e %12.5e  %12.5e\n",rr,fd,mat_r_epsu (rr,fd),FDvec_r_epsu(rr)-vec_r_epsu(rr),FDvec_r_epsu(rr),vec_r_epsu(rr));
        printf("FDvec_r_o_n_u_minus_g(%2d,%2d) %12.5e  %12.5e %12.5e  %12.5e\n",rr,fd,mat_r_o_n_u(rr,fd),FDvec_r_o_n_u_minus_g(rr)-vec_r_o_n_u_minus_g(rr),FDvec_r_o_n_u_minus_g(rr),vec_r_o_n_u_minus_g(rr));
      }
    }
  }
#endif
  return;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
