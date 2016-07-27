/*!----------------------------------------------------------------------
\file fluid_ele_boundary_calc.cpp

\brief evaluation of fluid terms at integration points

\level 1

<pre>
\maintainer Martin Kronbichler & Volker Gravemeier
             {kronbichler,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*----------------------------------------------------------------------*/

#include "fluid_ele_boundary_calc.H"
#include "fluid_ele.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_action.H"

#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_material.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/herschelbulkley.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/permeablefluid.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/cavitationfluid.H"
#include "../drt_mat/yoghurt.H"
#include "../drt_mat/fluidporo.H"
#include "../drt_mat/matlist.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
  : //DRT::ELEMENTS::FluidBoundaryInterface(),
    xyze_(true),
    funct_(true),
    deriv_(true),
    unitnormal_(true),
    velint_(true),
    drs_(0.0),
    fac_(0.0),
    visc_(0.0),
    densaf_(1.0)
{
  // pointer to class FluidImplParameterTimInt for time integration
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
  // initialize also general parameter list, also it will be overwritten in derived subclasses
  fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::EvaluateAction(DRT::ELEMENTS::FluidBoundary*  ele1,
                                                          Teuchos::ParameterList&         params,
                                                          DRT::Discretization&            discretization,
                                                          std::vector<int>&               lm,
                                                          Epetra_SerialDenseMatrix&       elemat1,
                                                          Epetra_SerialDenseMatrix&       elemat2,
                                                          Epetra_SerialDenseVector&       elevec1,
                                                          Epetra_SerialDenseVector&       elevec2,
                                                          Epetra_SerialDenseVector&       elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params,"action");

  // get status of Ale
  const bool isale = ele1->ParentElement()->IsAle();

  switch(act)
  {
  case FLD::integrate_Shapefunction:
  {
    Teuchos::RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=Teuchos::null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }

    IntegrateShapeFunction(
        ele1,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp);
    break;
  }
  case FLD::calc_area:
  {
    if (ele1->Owner() == discretization.Comm().MyPID())
      AreaCalculation(
          ele1,
          params,
          discretization,
          lm);
    break;
  }
  case FLD::calc_pressure_bou_int:
  {
    if(ele1->Owner() == discretization.Comm().MyPID())
      PressureBoundaryIntegral(
          ele1,
          params,
          discretization,
          lm);
    break;
  }
  // general action to calculate the flow rate
  case FLD::calc_flowrate:
  {
    ComputeFlowRate(
        ele1,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FLD::flowratederiv:
  {
    FlowRateDeriv(
        ele1,
        params,
        discretization,
        lm,
        elemat1,
        elemat2,
        elevec1,
        elevec2,
        elevec3);
    break;
  }
  case FLD::Outletimpedance:
  {
    ImpedanceIntegration(
        ele1,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }

  case FLD::dQdu:
  {
    dQdu(
        ele1,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FLD::ba_calc_node_normal:
  {
    Teuchos::RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=Teuchos::null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }

    DRT::UTILS::ElementNodeNormal<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                           ele1,discretization,elevec1,mydispnp,
                                           IsNurbs<distype>::isnurbs, ele1->ParentElement()->IsAle());
    break;
  }
  case FLD::calc_node_curvature:
  {
    Teuchos::RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=Teuchos::null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }

    Teuchos::RCP<const Epetra_Vector> normals;
    std::vector<double> mynormals;

    normals = discretization.GetState("normals");
    if (normals!=Teuchos::null)
    {
      mynormals.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*normals,mynormals,lm);
    }

    // what happens, if the mynormals vector is empty? (ehrl)
    dserror("the action calc_node_curvature has not been called by now. What happens, if the mynormal vector is empty");

    ElementMeanCurvature(
        ele1,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp,
        mynormals);
    break;
  }
  case FLD::conservative_outflow_bc:
  {
    ConservativeOutflowConsistency(
        ele1,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::calc_Neumann_inflow:
  {
    NeumannInflow(
        ele1,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::calc_surface_tension:
  {
    // employs the divergence theorem acc. to Saksono eq. (24) and does not
        // require second derivatives.

        Teuchos::RCP<const Epetra_Vector> dispnp;
        std::vector<double> mydispnp;

        dispnp = discretization.GetState("dispnp");
        if (dispnp!=Teuchos::null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }

        // mynormals and mycurvature are not used in the function
        std::vector<double> mynormals;
        std::vector<double> mycurvature;

        ElementSurfaceTension(
            ele1,
            params,
            discretization,
            lm,
            elevec1,
            mydispnp,
            mynormals,
            mycurvature);
        break;
  }
  case FLD::center_of_mass_calc:
  {
    // evaluate center of mass
    if(ele1->Owner() == discretization.Comm().MyPID())
      CenterOfMassCalculation(
          ele1,
          params,
          discretization,
          lm);
    break;
  }
  case FLD::traction_velocity_component:
  {
    CalcTractionVelocityComponent(
        ele1,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FLD::traction_Uv_integral_component:
  {
    ComputeNeumannUvIntegral(
        ele1,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  default:
  {
    dserror("Unknown type of action for FluidBoundaryImpl!");
    break;
  }
  } // end of switch(act)


}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidBoundaryImpl<distype>::EvaluateNeumann(
                              DRT::ELEMENTS::FluidBoundary*  ele,
                              Teuchos::ParameterList&        params,
                              DRT::Discretization&           discretization,
                              DRT::Condition&                condition,
                              std::vector<int>&              lm,
                              Epetra_SerialDenseVector&      elevec1_epetra,
                              Epetra_SerialDenseMatrix*      elemat1_epetra)
{
  // find out whether we will use a time curve
  const double time = fldparatimint_->Time();
  const bool usetime = (time<0.0) ? false : true;

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");
  const std::vector<int>*    curve = condition.Get<std::vector<int> >   ("curve");
  const std::string*         type  = condition.Get<std::string>("type");

  // get time factor for Neumann term
  const double timefac = fldparatimint_->TimeFacRhs();
  const double timefacn = (1.0-fldparatimint_->Theta())*fldparatimint_->Dt();

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get local node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);


  // THESE ARE NEEDED IF DENSITY OR A SCALAR IS USED in the BC (which normally is NOT the case)
  //========================================================
  // get scalar vector
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  if (scaaf==Teuchos::null) dserror("Cannot get state vector 'scaaf'");

  // extract local values from global vector
  std::vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  LINALG::Matrix<bdrynen_,1> escaaf(true);

  // insert scalar into element array
  // the scalar is stored to the pressure dof
  for (int inode=0;inode<bdrynen_;++inode)
  {
    escaaf(inode) = myscaaf[(nsd_)+(inode*numdofpernode_)];
  }

  // get thermodynamic pressure at n+1/n+alpha_F
  const double thermpressaf = params.get<double>("thermpress at n+alpha_F/n+1",0.0);
  //========================================================


  // add potential ALE displacements
  if (ele->ParentElement()->IsAle())
  {
    Teuchos::RCP<const Epetra_Vector>  dispnp;
    std::vector<double>                mydispnp;
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(IsNurbs<distype>::isnurbs)
  {
    std::vector<Epetra_SerialDenseVector> mypknots(nsd_);

    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
        ele, ele->SurfaceNumber(), ele->ParentElement()->Id(), discretization, mypknots, myknots, weights, normalfac);
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
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,&myknots,&weights,
                                                    IsNurbs<distype>::isnurbs);

    // get the required material information
    Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

    // get density
    // (evaluation always at integration point, in contrast to parent element)
    GetDensity(material,escaaf,thermpressaf);

    const double tol=1e-8;
    if(densfac_<(1.0-tol) or densfac_>(1.0+tol))
    {
      // If you got this warning have a look in GetDensity(...) for clarification!
      std::cout << "                                                                    " << std::endl;
      std::cout << "                                                                    " << std::endl;
      std::cout << "          WARNING:                                                  " << std::endl;
      std::cout << "                  BACI scales your NEUMANN BC with the DENSITY!!!   " << std::endl;
      std::cout << "                                                                    " << std::endl;
      std::cout << "                  Do you really want this?                          " << std::endl;
      std::cout << "                  Like, super sure about this?                      " << std::endl;
      std::cout << "                  Might want to think this through one more time...." << std::endl;
      std::cout << "                                                                    " << std::endl;
      std::cout << "                                                                    " << std::endl;
      std::cout << "        Fine... Do what you want...                                 " << std::endl;
      std::cout << "        densfac_=               " << densfac_ << std::endl;
    }

    // factor given by temporal curve
    double curvefac = 1.0;
    double curvefacn = 1.0;
    // number of temporal curve to be evaluated (from input file)
    int curvenum = -1;

    const double fac_time_dens = fac_*timefac*densfac_;
    const double fac_time_densn = fac_*timefacn;

    // factor given by spatial function
    double functfac = 1.0;
    double functfacn = 1.0;

    // global coordinates of gausspoint
    LINALG::Matrix<(nsd_),1>  coordgp(0.0);

    // determine coordinates of current Gauss point
    coordgp.Multiply(xyze_,funct_);

    // we need a 3D position vector for function evaluation!
    double coordgp3D[3];
    coordgp3D[0]=0.0;
    coordgp3D[1]=0.0;
    coordgp3D[2]=0.0;
    for (int i=0;i<nsd_;i++)
      coordgp3D[i]=coordgp(i);

    int functnum = -1;
    const double* coordgpref = &coordgp3D[0]; // needed for function evaluation

    for(int idim=0; idim<(nsd_); ++idim)
    {
      if (*type == "neum_live")
      {
        if((*onoff)[idim])  // Is this dof activated
        {
          if (func)
            functnum = (*func)[idim];
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,coordgpref,time,NULL);
            functfacn = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,coordgpref,time-fldparatimint_->Dt(),NULL);
          }
          else
          {
            functfac = 1.0;
            functfacn = 1.0;
          }

          // get time-curve factor/ n = - grad phi / |grad phi|
          if (curve)
            curvenum = (*curve)[idim];
          if (curvenum>=0 and usetime)
          {
            curvefac  = DRT::Problem::Instance()->Curve(curvenum).f(time);
            curvefacn = DRT::Problem::Instance()->Curve(curvenum).f(time-fldparatimint_->Dt());
          }
          else
          {
            curvefac  = 1.0;
            curvefacn = 1.0;
          }

          const double valfac = (*val)[idim]*fac_time_dens*functfac*curvefac;
          const double valfacn = (*val)[idim]*fac_time_densn*functfacn*curvefacn;
          for(int inode=0; inode < bdrynen_; ++inode )
          {
            elevec1_epetra[inode*numdofpernode_+idim] += funct_(inode)*valfac;
            if(fldparatimint_->IsNewOSTImplementation())
            {
              if(fldparatimint_->IsOneStepTheta())
                elevec1_epetra[inode*numdofpernode_+idim] += funct_(inode)*valfacn;
            }
          } //end IsNewOSTImplementation
        }  // if (*onoff)
      }
      else if (*type == "neum_pseudo_orthopressure")
      {
        if (idim != 0 and (*onoff)[idim])
          dserror("If you apply a pseudo_orthopressure load on the fluid, only a load in\n"
              "the first component (which corresponds to the normal direction) is valid!");

        if((*onoff)[0])  // Do we have a load in normal direction?
        {
          if (func)
            functnum = (*func)[0];
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,coordgpref,time,NULL);
            functfacn = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,coordgpref,time-fldparatimint_->Dt(),NULL);
          }
          else
          {
            functfac = 1.0;
            functfacn = 1.0;
          }

          // get time-curve factor/ n = - grad phi / |grad phi|
          if (curve)
            curvenum = (*curve)[0];
          if (curvenum>=0 and usetime)
          {
            curvefac  = DRT::Problem::Instance()->Curve(curvenum).f(time);
            curvefacn = DRT::Problem::Instance()->Curve(curvenum).f(time-fldparatimint_->Dt());
          }
          else
          {
            curvefac  = 1.0;
            curvefacn = 1.0;
          }

          const double valfac = (*val)[0]*fac_time_dens*functfac*curvefac;
          const double valfacn = (*val)[0]*fac_time_densn*functfacn*curvefacn;
          for(int inode=0; inode < bdrynen_; ++inode )
          {
            elevec1_epetra[inode*numdofpernode_+idim] += funct_(inode)*valfac*(-unitnormal_(idim));

            if(fldparatimint_->IsNewOSTImplementation())
            {
              if(fldparatimint_->IsOneStepTheta())
                elevec1_epetra[inode*numdofpernode_+idim] += funct_(inode)*valfacn*(-unitnormal_(idim));
            }
          } //end IsNewOSTImplementation
        }  // if (*onoff)
      }
      else
        dserror("The type '%s' is not supported in the fluid neumann condition!",type->c_str());

    } //for(int idim=0; idim<(nsd_); ++idim)
  }

  return 0;
}

/*----------------------------------------------------------------------*
 | apply outflow boundary condition which is necessary for the          |
 | conservative element formulation (since the convective term was      |
 | partially integrated)                                                |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ConservativeOutflowConsistency(
    DRT::ELEMENTS::FluidBoundary*  ele,
    Teuchos::ParameterList&         params,
    DRT::Discretization&            discretization,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       elemat1_epetra,
    Epetra_SerialDenseVector&       elevec1_epetra)
{
  if(fldparatimint_->TimeAlgo()== INPAR::FLUID::timeint_afgenalpha or
      fldparatimint_->TimeAlgo()== INPAR::FLUID::timeint_npgenalpha or
      fldparatimint_->TimeAlgo()== INPAR::FLUID::timeint_one_step_theta)
       dserror("The boundary condition ConservativeOutflowConsistency is not supported by ost/afgenalpha/npgenalpha!!\n"
               "the convective term is not partially integrated!");

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
  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // ------------------------------------
  // get statevectors from discretisation
  // ------------------------------------

  // extract local displacements from the global vectors
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;

  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // extract local velocities from the global vectors
  LINALG::Matrix<nsd_, bdrynen_>   evel(true);

  Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
  if (vel==Teuchos::null) dserror("Cannot get state vector 'u and p (trial)'");

  // extract local values from the global vectors
  std::vector<double> myvel(lm.size());
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
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
        ele, ele->SurfaceNumber(), ele->ParentElement()->Id(), discretization, mypknots, myknots, weights, normalfac);
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
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,&myknots,&weights,
                                                    IsNurbs<distype>::isnurbs);

    // Multiply the normal vector with the integration factor
    unitnormal_.Scale(fac_);

    // in the case of nurbs the normal vector must be scaled with a special factor
    if (IsNurbs<distype>::isnurbs)
      unitnormal_.Scale(normalfac);

    // get velocity at integration point
    velint_.Multiply(evel,funct_);

    // compute normal flux
    const double u_o_n = velint_.Dot(unitnormal_);

    // rescaled flux (according to time integration)
    const double timefac_mat_u_o_n = timefac_mat*u_o_n;

   // dyadic product of element's normal vector and velocity
   LINALG::Matrix<nsd_,nsd_>  n_x_u(true);

   // dyadic product of u and n
   n_x_u.MultiplyNT(timefac_mat,velint_,unitnormal_);

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

      // temp(nsd_,nsd) = n_x_u(nsd_,nsd_)*funct_(ui)
      temp.Scale(funct_(ui));

      for (int idimcol=0; idimcol < (nsd_); ++idimcol) // loop over dimensions for the columns
      {
        const int fui   = numdofpernode_*ui+idimcol;

        for (int vi=0; vi<bdrynen_; ++vi)  // loop rows
        {
          // temp(nsd_,nsd) *= funct_(vi)
          temp.Scale(funct_(vi));

          for (int idimrow = 0; idimrow < nsd_; ++idimrow) // loop over dimensions for the rows
          {
            const int fvi = numdofpernode_*vi+idimrow;
            elemat1_epetra(fvi  ,fui  ) += temp(fvi, fui);
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
       const double timefac_mat_u_o_n_funct_ui = timefac_mat_u_o_n*funct_(ui);

       for (int vi=0; vi<bdrynen_; ++vi)  // loop rows
       {
         const int fvi = numdofpernode_*vi + idim;
         const double timefac_mat_u_o_n_funct_ui_funct_vi
                   =
                   timefac_mat_u_o_n_funct_ui*funct_(vi);

         elemat1_epetra(fvi  ,fui  ) += timefac_mat_u_o_n_funct_ui_funct_vi;
       }  // loop rows
     }  // loop columns
   }  //loop over dimensions

  // rhs
  {
    // 3 temp vector
    LINALG::Matrix<nsd_,1>    temp(velint_, false);

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

        elevec1_epetra(fvi) -= temp(fvi)*funct_(vi);
      } // end loop over dimensions
    } // ui
  } // end rhs
  } // end gaussloop

  return;
}// DRT::ELEMENTS::FluidSurface::SurfaceConservativeOutflowConsistency


/*----------------------------------------------------------------------*
 | compute additional term at Neumann inflow boundary          vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::NeumannInflow(
    DRT::ELEMENTS::FluidBoundary*  ele,
    Teuchos::ParameterList&        params,
    DRT::Discretization&           discretization,
    std::vector<int>&              lm,
    Epetra_SerialDenseMatrix&      elemat1,
    Epetra_SerialDenseVector&      elevec1)
{

  if(fldparatimint_->IsNewOSTImplementation())
  {
    if(fldparatimint_->IsOneStepTheta())
    {
      dserror("NEUMANN INFLOW IS NOT IMPLEMENTED FOR NEW OST AS OF YET!");
    }
  }//end IsNewOSTImplementation

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // af-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // np-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // genalpha:    timefac =  alpha_F * gamma * dt
  const double timefac = fldparatimint_->TimeFac();

  // get timefactor for right-hand side
  // One-step-Theta:            timefacrhs = theta*dt
  // BDF2:                      timefacrhs = 2/3 * dt
  // af-genalpha:               timefacrhs = (1/alpha_M) * gamma * dt
  // np-genalpha:               timefacrhs = (1/alpha_M) * gamma * dt
  // genalpha:                  timefacrhs = 1.0
  double timefacrhs = fldparatimint_->TimeFacRhs();

  // check ALE status
  const bool isale = ele->ParentElement()->IsAle();

  // set flag for type of linearization to default value (fixed-point-like)
  bool is_newton = fldpara_->IsNewton();

  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates for nsd_-dimensional domain
  // (nsd_: number of spatial dimensions of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;
  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get velocity and scalar vector at time n+alpha_F/n+1
  Teuchos::RCP<const Epetra_Vector> velaf = discretization.GetState("velaf");
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  if (velaf==Teuchos::null or scaaf==Teuchos::null)
    dserror("Cannot get state vector 'velaf' and/or 'scaaf'");

  // extract local values from global vector
  std::vector<double> myvelaf(lm.size());
  std::vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  // create Epetra objects for scalar array and velocities
  LINALG::Matrix<nsd_,bdrynen_> evelaf(true);
  LINALG::Matrix<bdrynen_,1>    escaaf(true);

  // insert velocity and scalar into element array
  for (int inode=0;inode<bdrynen_;++inode)
  {
    for (int idim=0; idim<(nsd_);++idim)
    {
      evelaf(idim,inode) = myvelaf[idim+(inode*numdofpernode_)];
    }
    escaaf(inode) = myscaaf[(nsd_)+(inode*numdofpernode_)];
  }

  // get thermodynamic pressure at n+1/n+alpha_F
  const double thermpressaf = params.get<double>("thermpress at n+alpha_F/n+1",1.0);

  // --------------------------------------------------
  // nurbs-specific stuff
  // --------------------------------------------------
  // normal vector multiplied by normalfac for nurbs
  double normalfac = 0.0;
  std::vector<Epetra_SerialDenseVector> mypknots(nsd_);
  std::vector<Epetra_SerialDenseVector> myknots (bdrynsd_);
  Epetra_SerialDenseVector weights(bdrynen_);

  // get knotvectors for parent element and surface element as well as weights
  // for isogeometric elements
  if(IsNurbs<distype>::isnurbs)
  {
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
        ele, ele->SurfaceNumber(), ele->ParentElement()->Id(), discretization, mypknots, myknots, weights, normalfac);
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
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,&myknots,&weights,
                                                    IsNurbs<distype>::isnurbs);

    // normal vector scaled by special factor in case of nurbs
    if (IsNurbs<distype>::isnurbs) unitnormal_.Scale(normalfac);

    // compute velocity vector and normal velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double normvel = 0.0;
    velint_.Multiply(evelaf,funct_);
    normvel = velint_.Dot(unitnormal_);

    // check normal velocity -> further computation only required for
    // negative normal velocity, that is, inflow at this Neumann boundary
    if (normvel<-0.0001)
    {
      // get the required material information
      Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

      // get density
      // (evaluation always at integration point, in contrast to parent element)
      GetDensity(material,escaaf,thermpressaf);

      // extended integration factors for left- and right-hand side, respectively
      const double lhsfac = densaf_*normvel*timefac*fac_;
      const double rhsfac = densaf_*normvel*timefacrhs*fac_;

      // compute matrix contribution (fill diagonal elements)
      /*
              /                        \
             |                          |
           - |  v , rho * Du ( u o n )  |
             |                          |
              \                        /
      */
      for (int idim = 0; idim < nsd_; ++idim) // loop over dimensions
      {
        for (int vi=0; vi<bdrynen_; ++vi) // loop over rows
        {
          const double vlhs = lhsfac*funct_(vi);

          const int fvi = numdofpernode_*vi+idim;

          for (int ui=0; ui<bdrynen_; ++ui) // loop over columns
          {
            const int fui = numdofpernode_*ui+idim;

            elemat1(fvi,fui) -= vlhs*funct_(ui);
          } // end loop over columns
        }  // end loop over rows
      }  // end loop over dimensions

      // compute additional matrix contribution for Newton linearization
      if (is_newton)
      {
        // integration factor
        const double lhsnewtonfac = densaf_*timefac*fac_;

        // dyadic product of unit normal vector and velocity vector
        LINALG::Matrix<nsd_,nsd_>  n_x_u(true);
        n_x_u.MultiplyNT(velint_,unitnormal_);

        /*
                /                        \
               |                          |
             - |  v , rho * u ( Du o n )  |
               |                          |
                \                        /

               rho * v_i * u_i * Du_j * n_j

        */
        for (int vi=0; vi<bdrynen_; ++vi) // loop rows
        {
          const double dens_dt_v = lhsnewtonfac*funct_(vi);

          for (int idimrow=0; idimrow < nsd_; ++idimrow) // loop row dim.
          {
            const int fvi = numdofpernode_*vi+idimrow;

            for (int ui=0; ui<bdrynen_; ++ui) // loop columns
            {
              const double dens_dt_v_Du = dens_dt_v * funct_(ui);

              for (int idimcol = 0; idimcol < nsd_; ++idimcol) // loop column dim.
              {
                const int fui = numdofpernode_*ui+idimcol;

                elemat1(fvi,fui) -= dens_dt_v_Du*n_x_u(idimrow,idimcol);
              } // end loop row dimensions
            } // end loop rows
          } // end loop column dimensions
        } // end loop columns
      } // end of Newton loop

      // compute rhs contribution
      LINALG::Matrix<nsd_,1> vrhs(velint_, false);
      vrhs.Scale(rhsfac);

      for (int vi=0; vi<bdrynen_; ++vi) // loop over rows
      {
        for (int idim = 0; idim < nsd_; ++idim)  // loop over dimensions
        {
          const int fvi = numdofpernode_*vi+idim;

          elevec1(fvi) += funct_(vi)*vrhs(idim);
        } // end loop over dimensions
      }  // end loop over rows
    }
  }

  return;
}// DRT::ELEMENTS::FluidSurface::NeumannInflow


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)            gjb 07/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::IntegrateShapeFunction(
                  DRT::ELEMENTS::FluidBoundary* ele,
                  Teuchos::ParameterList& params,
                  DRT::Discretization&       discretization,
                  std::vector<int>&          lm,
                  Epetra_SerialDenseVector&  elevec1,
                  const std::vector<double>& edispnp)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0;idim<(nsd_);++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
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
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<(nsd_);idim++)
      {
        elevec1(inode*numdofpernode_+idim)+= funct_(inode) * fac_;
      }
    }

  } /* end of loop over integration points gpid */


return;
} // DRT::ELEMENTS::FluidSurface::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ElementMeanCurvature(
                                                        DRT::ELEMENTS::FluidBoundary*    ele,
                                                        Teuchos::ParameterList&           params,
                                                        DRT::Discretization&              discretization,
                                                        std::vector<int>&                 lm,
                                                        Epetra_SerialDenseVector&         elevec1,
                                                        const std::vector<double>&        edispnp,
                                                        std::vector<double>&              enormals)
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // node normals &
  LINALG::Matrix<nsd_,bdrynen_> norm_elem(true);
  LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);

  // coordinates of current node in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi_node(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0;idim<nsd_; ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
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
  Epetra_SerialDenseMatrix xsi_ele = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  // ============================== loop over nodes ==========================
  for (int inode=0;inode<bdrynen_; ++inode)
  {
    // the local node coordinates matrix is split to a vector containing the local coordinates of the actual node
    for (int idim = 0; idim < bdrynsd_; idim++)
    {
      xsi_node(idim) = xsi_ele(idim,inode);
    }

    // get shape derivatives at this node
    // shape_function_2D_deriv1(deriv_, e0, e1, distype);
    DRT::UTILS::shape_function<distype>(xsi_node,funct_);

    // the metric tensor and its determinant
    //Epetra_SerialDenseMatrix      metrictensor(nsd_,nsd_);
    LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);

    // Addionally, compute metric tensor
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor,drs_);

    dxyzdrs.MultiplyNT(deriv_,xyze_);

    // calculate mean curvature H at node.
    double H = 0.0;
    LINALG::Matrix<bdrynsd_,nsd_> dn123drs(0.0);

    dn123drs.MultiplyNT(deriv_,norm_elem);

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
          / (drs_*drs_);
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

    // loop over adjacent Fluid elements
    for (int ele=0;ele<NumElement;ele++)
    {
      DRT::Element* Element = ElementsPtr[ele];

      // get surfaces
      std::vector< Teuchos::RCP< DRT::Element > > surfaces = Element->Surfaces();

      // loop over surfaces: how many free surfaces with this node on it?
      for (unsigned int surf=0; surf<surfaces.size(); ++surf)
      {
        Teuchos::RCP< DRT::Element > surface = surfaces[surf];
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

} // DRT::ELEMENTS::FluidSurface::ElementMeanCurvature



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ElementSurfaceTension(
                                                         DRT::ELEMENTS::FluidBoundary*   ele,
                                                         Teuchos::ParameterList&          params,
                                                         DRT::Discretization&             discretization,
                                                         std::vector<int>&                lm,
                                                         Epetra_SerialDenseVector&        elevec1,
                                                         const std::vector<double>&       edispnp,
                                                         std::vector<double>&             enormals,
                                                         std::vector<double>&             ecurvature)
                                                         // Attention: mynormals and mycurvature are not used in the function
{
  // get status of Ale
  const bool isale = ele->ParentElement()->IsAle();

  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // af-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // np-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // genalpha:    timefac =  alpha_F * gamma * dt
  const double timefac = fldparatimint_->TimeFac();

  // isotropic and isothermal surface tension coefficient
  double SFgamma = 0.0;
  // get material data
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();
  if (mat==Teuchos::null)
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

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    // fac multiplied by the timefac
    const double fac_timefac = fac_ * timefac;

    // Compute dxyzdrs
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(true);
    dxyzdrs.MultiplyNT(deriv_,xyze_);

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
                                       drs_ * drs_ //= abs_dxyzdr * abs_dxyzdr * abs_dxyzds * abs_dxyzds - pointproduct_rs * pointproduct_rs
                                       )
                                     *
                                     (
                                       abs_dxyzds * abs_dxyzds * deriv_(0,node) * dxyzdrs(0,dim)
                                       - pointproduct_rs * deriv_(0,node) * dxyzdrs(1,dim)
                                       - pointproduct_rs * deriv_(1,node) * dxyzdrs(0,dim)
                                       + abs_dxyzdr * abs_dxyzdr * deriv_(1,node) * dxyzdrs(1,dim)
                                       )
                                     * fac_timefac;

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
            elevec1[inode*numdofpernode_+idim]+= SFgamma / drs_ / drs_ *
                                      (-1.0) * deriv_(0, inode) * dxyzdrs(0,idim)
                                      * fac_timefac;
         }
      }
    } // end if else (nsd_=1)
    else
      dserror("There are no 3D boundary elements implemented");
  } /* end of loop over integration points gpid */
} // DRT::ELEMENTS::FluidSurface::ElementSurfaceTension

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::AreaCalculation(
  DRT::ELEMENTS::FluidBoundary*  ele,
  Teuchos::ParameterList&        params,
  DRT::Discretization&           discretization,
  std::vector<int>&              lm)
{
  //------------------------------------------------------------------
  // get and set density and viscosity (still required for following routines:
  // FluidImpedanceBc/FluidVolumetricSurfaceFlowBc/Fluid_couplingBc::Area)
  //------------------------------------------------------------------
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();
  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    densaf_ = actmat->Density();
    visc_   = actmat->Viscosity();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_permeable_fluid)
  {
    const MAT::PermeableFluid* actmat = static_cast<const MAT::PermeableFluid*>(mat.get());
    densaf_ = actmat->Density();
    visc_   = actmat->SetViscosity();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    densaf_ = actmat->Density();

    const double nu_inf = actmat->NuInf();  // parameter for infinite-shear viscosity

    // dynamic viscosity = kinematic viscosity * density
    visc_ = nu_inf *densaf_;
  }

  params.set<double>("density",   densaf_);
  params.set<double>("viscosity", visc_);
  //------------------------------------------------------------------
  // end of get and set density and viscosity
  //------------------------------------------------------------------

  //------------------------------------------------------------------
  // start of actual area calculation
  //------------------------------------------------------------------
  // get node coordinates (nsd_: dimension of boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;
  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get initial value for area
  double area = params.get<double>("area");

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    // add to area integral
    area += fac_;
  }

  // set final value for area
  params.set<double>("area",area);
  //------------------------------------------------------------------
  // end of actual area calculation
  //------------------------------------------------------------------

}//DRT::ELEMENTS::FluidSurface::AreaCalculation


/*----------------------------------------------------------------------*
 |                                                       ismail 04/2010 |
 |                                                           vg 06/2013 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::PressureBoundaryIntegral(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&          params,
  DRT::Discretization&             discretization,
  std::vector<int>&                lm)
{
  // extract pressure values from global velocity/pressure vector
  //renamed to "velaf" to be consistent in fluidimplicitintegration.cpp (krank 12/13)
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velaf");
  if (velnp == Teuchos::null) dserror("Cannot get state vector 'velaf'");

  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  LINALG::Matrix<1,bdrynen_> eprenp(true);
  for (int inode=0;inode<bdrynen_;inode++)
  {
    eprenp(inode) = myvelnp[nsd_+inode*numdofpernode_];
  }

  // get node coordinates (nsd_: dimension of boundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector>  dispnp;
  std::vector<double>                mydispnp;
  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get initial value for pressure boundary integral
  double press_int = params.get<double>("pressure boundary integral");

  // get Gauss rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    // add to pressure boundary integral
    for (int inode=0;inode<bdrynen_;++inode)
    {
      press_int += funct_(inode) * eprenp(inode) *fac_;
    }
  }

  // set final value for pressure boundary integral
  params.set<double>("pressure boundary integral",press_int);

}//DRT::ELEMENTS::FluidSurface::PressureBoundaryIntegral


/*----------------------------------------------------------------------*
 |                                                        ismail 10/2010|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::CenterOfMassCalculation(
  DRT::ELEMENTS::FluidBoundary*     ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm)
{

  //------------------------------------------------------------------
  // This calculates the integrated the pressure from the
  // the actual pressure values
  //------------------------------------------------------------------

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze_);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // first evaluate the area of the surface element
  params.set<double>("area",0.0);
  this->AreaCalculation(ele, params, discretization,lm);

  // get the surface element area
  const double elem_area = params.get<double>("area");

  LINALG::Matrix<(nsd_),1>  xyzGe(true);

  for (int i = 0; i< nsd_;i++)
  {
    //const IntegrationPoints2D  intpoints(gaussrule);
    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
      // Computation of the unit normal vector at the Gauss points
      // Computation of nurb specific stuff is not activated here
      DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                      intpoints,gpid,NULL,NULL,
                                                      IsNurbs<distype>::isnurbs);

      // global coordinates of gausspoint
      LINALG::Matrix<(nsd_),1>  coordgp(true);

      // determine coordinates of current Gauss point
      coordgp.Multiply(xyze_,funct_);

      //Compute elment center of gravity
      xyzGe(i) += intpoints.IP().qwgt[gpid]*coordgp(i)*drs_;

    }  // end Gauss loop
    xyzGe(i) /= elem_area;
  }

  // Get the center of mass of the already calculate surface elements
  Teuchos::RCP<std::vector<double> > xyzG  = params.get<Teuchos::RCP<std::vector<double> > >("center of mass");

  Teuchos::RCP<std::vector<double> > normal  = params.get<Teuchos::RCP<std::vector<double> > >("normal");

  // Get the area of the of the already calculate surface elements
  double area = params.get<double>("total area");

  for (int i = 0; i<nsd_;i++)
  {
    (*xyzG)  [i] = ((*xyzG)[i]*area   + xyzGe(i)     *elem_area)/(area+elem_area);
    (*normal)[i] = ((*normal)[i]*area + unitnormal_(i)*elem_area)/(area+elem_area);
  }

  // set new center of mass
  params.set("total area", area+elem_area);

}//DRT::ELEMENTS::FluidSurface::CenterOfMassCalculation



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ComputeFlowRate(
                                                                DRT::ELEMENTS::FluidBoundary*     ele,
                                                                Teuchos::ParameterList&           params,
                                                                DRT::Discretization&              discretization,
                                                                std::vector<int>&                 lm,
                                                                Epetra_SerialDenseVector&         elevec1)
{
  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  //renamed to "velaf" to be consistent in fluidimplicitintegration.cpp (krank 12/13)
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velaf");

  if (velnp==Teuchos::null)
    dserror("Cannot get state vector 'velaf'");

  std::vector<double> myvelnp(lm.size());
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
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  //GEO::fillInitialPositionArray<distype,nsd_,Epetra_SerialDenseMatrix>(ele,xyze_);
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  //const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    //compute flowrate at gauss point
    velint_.Multiply(evelnp,funct_);

    // flowrate = uint o normal
    const double flowrate = velint_.Dot(unitnormal_);

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode=0;inode<bdrynen_;++inode)
    {
      // see "A better consistency for low order stabilized finite element methods"
      // Jansen, Collis, Whiting, Shakib
      //
      // Here the principle is used to bring the flow rate to the outside world!!
      //
      // funct_ *  velint * n * fac
      //   |      |________________|
      //   |              |
      //   |         flow rate * fac  -> integral over Gamma
      //   |
      // flow rate is distributed to the single nodes of the element
      // = flow rate per node
      //
      // adding up all nodes (ghost elements are handled by the assembling strategy)
      // -> total flow rate at the desired boundary
      //
      // it can be interpreted as a rhs term
      //
      //  ( v , u o n)
      //               Gamma
      //
      elevec1[inode*numdofpernode_] += funct_(inode)* fac_ * flowrate;

      // alternative way:
      //
      //  velint * n * fac
      // |________________|
      //         |
      //    flow rate * fac  -> integral over Gamma
      //     = flow rate per element
      //
      //  adding up all elements (be aware of ghost elements!!)
      //  -> total flow rate at the desired boundary
      //     (is identical to the total flow rate computed above)
    }
  }
}//DRT::ELEMENTS::FluidSurface::ComputeFlowRate


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::FlowRateDeriv(
                                                 DRT::ELEMENTS::FluidBoundary*    ele,
                                                 Teuchos::ParameterList&          params,
                                                 DRT::Discretization&             discretization,
                                                 std::vector<int>&                lm,
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

  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> edispnp;

  if (isale)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==Teuchos::null) dserror("Cannot get state vectors 'dispnp'");
    edispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,edispnp,lm);
  }

  // get integration rule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // order of accuracy of grid velocity determination
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  const int gridvel = DRT::INPUT::IntegralValue<INPAR::FLUID::Gridvel>(fdyn, "GRIDVEL");

  // normal vector
  LINALG::Matrix<nsd_,1> normal(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  if (isale)
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int inode=0; inode<bdrynen_; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode) += edispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // get nodal velocities and pressures
  Teuchos::RCP<const Epetra_Vector> convelnp = discretization.GetState("convectivevel");

  if (convelnp==Teuchos::null)
    dserror("Cannot get state vector 'convectivevel'");

  // extract local values from the global vectors
  std::vector<double> myconvelnp(lm.size());
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
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal
    // Therefore it cancels out!!
    const double fac = intpoints.IP().qwgt[gpid];

    // dxyzdrs vector -> normal which is not normalized
    LINALG::Matrix<bdrynsd_,nsd_> dxyzdrs(0.0);
    dxyzdrs.MultiplyNT(deriv_,xyze_);
    normal(0,0) = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    normal(1,0) = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    normal(2,0) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

    //-------------------------------------------------------------------
    //  Q
    LINALG::Matrix<3,1> u(true);
    for (int dim=0;dim<3;++dim)
      for (int node=0;node<bdrynen_;++node)
        u(dim) += funct_(node) * evelnp(dim,node);

    for(int dim=0;dim<3;++dim)
      elevec3[0] += u(dim) * normal(dim,0) * fac;

    if (params.get<bool>("flowrateonly", false)==false)
    {
      //-------------------------------------------------------------------
      // dQ/du
      for (int node=0;node<bdrynen_;++node)
      {
        for (int dim=0;dim<3;++dim)
          elevec1[node*numdofpernode_+dim] += funct_(node) * normal(dim,0) * fac;
        elevec1[node*numdofpernode_+3] = 0.0;
      }

      //-------------------------------------------------------------------
      // dQ/dd

      // determine derivatives of surface normals wrt mesh displacements
      LINALG::Matrix<3,bdrynen_*3> normalderiv(true);

      for (int node=0;node<bdrynen_;++node)
      {
        normalderiv(0,3*node)   = 0.;
        normalderiv(0,3*node+1) = deriv_(0,node)*dxyzdrs(1,2)-deriv_(1,node)*dxyzdrs(0,2);
        normalderiv(0,3*node+2) = deriv_(1,node)*dxyzdrs(0,1)-deriv_(0,node)*dxyzdrs(1,1);

        normalderiv(1,3*node)   = deriv_(1,node)*dxyzdrs(0,2)-deriv_(0,node)*dxyzdrs(1,2);
        normalderiv(1,3*node+1) = 0.;
        normalderiv(1,3*node+2) = deriv_(0,node)*dxyzdrs(1,0)-deriv_(1,node)*dxyzdrs(0,0);

        normalderiv(2,3*node)   = deriv_(0,node)*dxyzdrs(1,1)-deriv_(1,node)*dxyzdrs(0,1);
        normalderiv(2,3*node+1) = deriv_(1,node)*dxyzdrs(0,0)-deriv_(0,node)*dxyzdrs(1,0);
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

        if (gridvel == INPAR::FLUID::BE)  // BE time discretization
        {
          for (int node=0;node<bdrynen_;++node)
          {
            for (int dim=0;dim<3;++dim)
              elevec2[node*numdofpernode_+dim] -= 1.0/dt * funct_(node) * normal(dim,0) * fac;
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
                elemat1(unode*numdofpernode_+udim,nnode*numdofpernode_+ndim) = funct_(unode) * normalderiv(udim,3*nnode+ndim) * fac;
            }
          }
        }
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(dd)^2

      // determine second derivatives of surface normals wrt mesh displacements
      std::vector<LINALG::Matrix<bdrynen_*3,bdrynen_*3> > normalderiv2(3);

      for (int node1=0;node1<bdrynen_;++node1)
      {
        for (int node2=0;node2<bdrynen_;++node2)
        {
          double temp = deriv_(0,node1)*deriv_(1,node2)-deriv_(1,node1)*deriv_(0,node2);

          normalderiv2[0](node1*3+1,node2*3+2) = temp;
          normalderiv2[0](node1*3+2,node2*3+1) = - temp;

          normalderiv2[1](node1*3  ,node2*3+2) = - temp;
          normalderiv2[1](node1*3+2,node2*3  ) = temp;

          normalderiv2[2](node1*3  ,node2*3+1) = temp;
          normalderiv2[2](node1*3+1,node2*3  ) = - temp;
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
                    u(iterdim) * normalderiv2[iterdim](node1*3+dim1,node2*3+dim2) * fac;
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

        if (gridvel == INPAR::FLUID::BE)
        {
          for (int node1=0;node1<bdrynen_;++node1)
          {
            for (int dim1=0;dim1<3;++dim1)
            {
              for (int node2=0;node2<bdrynen_;++node2)
              {
                for (int dim2=0;dim2<3;++dim2)
                {
                  elemat2(node1*numdofpernode_+dim1,node2*numdofpernode_+dim2) -= (1.0/dt * funct_(node1) * normalderiv(dim1, 3*node2+dim2)
                                                                                   + 1.0/dt * funct_(node2) * normalderiv(dim2, 3*node1+dim1))
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
}//DRT::ELEMENTS::FluidSurface::FlowRateDeriv


 /*----------------------------------------------------------------------*
  |  Impedance related parameters on boundary elements          AC 03/08  |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ImpedanceIntegration(
                  DRT::ELEMENTS::FluidBoundary*     ele,
                  Teuchos::ParameterList&           params,
                  DRT::Discretization&              discretization,
                  std::vector<int>&                 lm,
                  Epetra_SerialDenseVector&         elevec1)
{
  const double tfacrhs = fldparatimint_->TimeFacRhs();

  const double pressure = params.get<double>("WindkesselPressure");

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  if (ele->ParentElement()->IsAle())
  {
    Teuchos::RCP<const Epetra_Vector>  dispnp;
    std::vector<double>                mydispnp;

    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // add traction in the inward normal direction with norm 'pressure'
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    const double fac_facrhs_pres = fac_ * tfacrhs * pressure;

    for (int inode=0;inode<bdrynen_;++inode)
      for(int idim=0;idim<nsd_;++idim)
        elevec1[inode*numdofpernode_+idim] += fac_facrhs_pres * funct_(inode) * (-unitnormal_(idim));
  }

  return;
} //DRT::ELEMENTS::FluidSurface::ImpedanceIntegration


/*---------------------------------------------------------------------------*
 |  linearization of flux w.r.t velocities on boundary elements  Thon 10/15  |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::dQdu(
                  DRT::ELEMENTS::FluidBoundary*     ele,
                  Teuchos::ParameterList&           params,
                  DRT::Discretization&              discretization,
                  std::vector<int>&                 lm,
                  Epetra_SerialDenseVector&         elevec1)
{
  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  if (ele->ParentElement()->IsAle())
  {
    Teuchos::RCP<const Epetra_Vector>  dispnp;
    std::vector<double>                mydispnp;

    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  // compute dQ/du were (dQ/du)_i= (\phi_i \dot n)_Gamma
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    for (int node=0;node<bdrynen_;++node)
    {
      for (int dim=0;dim<nsd_;++dim)
      {
        elevec1[node*numdofpernode_+dim] += funct_(node) * unitnormal_(dim,0) * fac_;
      }
      elevec1[node*numdofpernode_+nsd_] = 0.0;
    }
  }

  params.set<double>("tfaclhs",fldparatimint_->TimeFac());
}


/*----------------------------------------------------------------------*
 |  get density                                                vg 06/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::GetDensity(
  Teuchos::RCP<const MAT::Material>    material,
  const LINALG::Matrix<bdrynen_,1>&    escaaf,
  const double                         thermpressaf
)
{
  // initially set density and density factor for Neumann boundary conditions to 1.0
  // (the latter only changed for low-Mach-number flow/combustion problems)
  // (This is due to the nature of the Neumann BC's for these problems as they are usually given as h_N=\rho*h.)
  // (Thus in the input file for these problems, h is set which is then scaled by \rho here.)
  //
  //          See Gravemeier, Wall 2011 IJNMF example 4.1.2
  //
  densaf_  = 1.0;
  densfac_ = 1.0;

  if (material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

    // varying density
    if (fldpara_->PhysicalType() == INPAR::FLUID::varying_density)
      densaf_ = funct_.Dot(escaaf);
    // Boussinesq approximation: Calculation of delta rho
    else if (fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
      dserror("Boussinesq approximation not yet supported for boundary terms!");
    else
      densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(material.get());

    densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(material.get());

    densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_herschelbulkley)
  {
    const MAT::HerschelBulkley* actmat = static_cast<const MAT::HerschelBulkley*>(material.get());

    densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_yoghurt)
  {
    const MAT::Yoghurt* actmat = static_cast<const MAT::Yoghurt*>(material.get());

    // get constant density
    densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
  {
    const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

    // compute mixture fraction at n+alpha_F or n+1
    const double mixfracaf = funct_.Dot(escaaf);

    // compute density at n+alpha_F or n+1 based on mixture fraction
    densaf_ = actmat->ComputeDensity(mixfracaf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
  {
    const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

    // compute temperature at n+alpha_F or n+1
    const double tempaf = funct_.Dot(escaaf);

    // compute density at n+alpha_F or n+1 based on temperature
    // and thermodynamic pressure

    densaf_ = actmat->ComputeDensity(tempaf,thermpressaf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->MaterialType() == INPAR::MAT::m_cavitation)
  {
    const MAT::CavitationFluid* actmat = static_cast<const MAT::CavitationFluid*>(material.get());

    // get constant base density
    const double density_0 = actmat->Density();

    // get fluid fraction at at n+alpha_F
    // const double fluidfracaf = funct_.Dot(escaaf);

    // compute density at at n+alpha_F; no density scaling necessary here due
    // to the special choice of forces applied to the fluid
    //  densaf_ = fluidfracaf*density_0;
    densaf_ = density_0;
  }
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
  {
    const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

    // get progress variable at n+alpha_F or n+1
    const double provaraf = funct_.Dot(escaaf);

    // compute density at n+alpha_F or n+1 based on progress variable
    densaf_ = actmat->ComputeDensity(provaraf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
  {
    const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

    // get progress variable at n+alpha_F or n+1
    const double provaraf = funct_.Dot(escaaf);

    // compute density at n+alpha_F or n+1 based on progress variable
    densaf_ = actmat->ComputeDensity(provaraf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->MaterialType() == INPAR::MAT::m_permeable_fluid)
  {
    const MAT::PermeableFluid* actmat = static_cast<const MAT::PermeableFluid*>(material.get());

    densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_fluidporo)
  {
    const MAT::FluidPoro* actmat = static_cast<const MAT::FluidPoro*>(material.get());

    densaf_ = actmat->Density();
  }
  else if (material->MaterialType() == INPAR::MAT::m_matlist)
  {

    densaf_=1.0;
    densfac_=1.0;

    //This is only necessary if the BC is dependant on the density!

  }// end else if m_matlist
  else dserror("Material type is not supported for density evaluation for boundary element!");

//  // check whether there is zero or negative density
  if (densaf_ < EPS15) dserror("zero or negative density!");




  return;
} // FluidBoundaryImpl::GetDensity

/*----------------------------------------------------------------------*
 |  Evaluating the velocity component of the traction      ismail 05/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::CalcTractionVelocityComponent(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm,
  Epetra_SerialDenseVector&         elevec1)
{

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velaf");

  if (velnp==Teuchos::null)
    dserror("Cannot get state vector 'velaf'");

  std::vector<double> myvelnp(lm.size());
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


  Teuchos::RCP<Epetra_Vector> cond_velocities = params.get<Teuchos::RCP<Epetra_Vector> > ("condition velocities");
  Teuchos::RCP<Epetra_Map>    cond_dofrowmap  = params.get<Teuchos::RCP<Epetra_Map> > ("condition dofrowmap");

  double density=0.0; // inverse density of my parent element

  // get material of volume element this surface belongs to
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();

  if( mat->MaterialType() != INPAR::MAT::m_carreauyasuda
   && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
   && mat->MaterialType() != INPAR::MAT::m_herschelbulkley
   && mat->MaterialType() != INPAR::MAT::m_fluid
   && mat->MaterialType() != INPAR::MAT::m_permeable_fluid)
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
  else if(mat->MaterialType()== INPAR::MAT::m_herschelbulkley)
  {
    const MAT::HerschelBulkley* actmat = static_cast<const MAT::HerschelBulkley*>(mat.get());
    density = actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_permeable_fluid)
  {
    const MAT::PermeableFluid* actmat = static_cast<const MAT::PermeableFluid*>(mat.get());
    density = actmat->Density();
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());

  //-------------------------------------------------------------------
  // get the tractions velocity component
  //-------------------------------------------------------------------

  // get Gaussrule
  //  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToGaussRuleForExactSol<distype>::rule);
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,bdrynen_> >(ele,xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector>      dispnp;
  std::vector<double>                mydispnp;

  if (ele->ParentElement()->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }
    dsassert(mydispnp.size()!=0,"paranoid");
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        xyze_(idim,inode)+=mydispnp[numdofpernode_*inode+idim];
      }
    }
  }

  const double timefac = fldparatimint_->TimeFacRhs();

  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the shape function at the Gauss point
    // Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(funct_,deriv_,fac_,unitnormal_,drs_,xsi_,xyze_,
                                                    intpoints,gpid,NULL,NULL,
                                                    IsNurbs<distype>::isnurbs);

    // Get the velocity value at the corresponding Gauss point.
    std::vector<double> vel_gps(nsd_,0.0);
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<nsd_;++idim)
      {
        vel_gps[idim] += myvelnp[inode*numdofpernode_+idim]*funct_(inode);
      }
    }

    // Evaluate the normal velocity at the corresponding Gauss point
    double n_vel = 0.0;
    for(int idim = 0 ;idim<nsd_;++idim)
    {
      n_vel += vel_gps[idim]*(unitnormal_(idim));
    }
    // loop over all node and add the corresponding effect of the Neumann-Inflow condition
    for (int inode=0;inode<bdrynen_;++inode)
    {
      for(int idim=0;idim<nsd_;++idim)
      {
        // evaluate the value of the Un.U at the corresponding Gauss point
        const double  uV = n_vel*vel_gps[idim] * density;
        const double fac_thsl_pres_inve = fac_ * timefac  * uV;

        // remove the Neumann-inflow contribution only if the normal velocity is an inflow velocity
        // i.e n_vel < 0
        if (n_vel<0.0)
        {
          elevec1[inode*numdofpernode_+idim] -= fac_thsl_pres_inve*funct_(inode);
        }
      }
      //      double radius = sqrt(pow(xyze_(0,inode),2.0)+pow(xyze_(1,inode),2.0));
      //      cout<<"n_vel("<<n_vel<<") vel: "<<n_vel<<" rad: "<<radius<<endl;
    }
  }
  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidBoundaryImpl<distype>::ComputeNeumannUvIntegral(
  DRT::ELEMENTS::FluidBoundary*    ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm,
  Epetra_SerialDenseVector&         elevec1)
{
}
//DRT::ELEMENTS::FluidSurface::ComputeNeumannUvIntegral

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::line2>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::line3>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::nurbs2>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::nurbs4>;
template class DRT::ELEMENTS::FluidBoundaryImpl<DRT::Element::nurbs9>;
