/*----------------------------------------------------------------------*/
/*!
\file thermo_ele_boundary_impl.cpp

\brief Internal implementation of thermo boundary elements (ThermoBoundary)

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
 */

/*----------------------------------------------------------------------*
 |  definitions                                              dano 09/09 |
 *----------------------------------------------------------------------*/
#ifdef D_THERMO

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/
#include <cstdlib>
#include "thermo_ele_boundary_impl.H"
#include "thermo_ele_impl.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_function.H"

#include "../drt_inpar/inpar_thermo.H"

#include "../drt_tsi/tsi_defines.H"

/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TemperBoundaryImplInterface* DRT::ELEMENTS::TemperBoundaryImplInterface::Impl(DRT::Element* ele)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));

  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    static TemperBoundaryImpl<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new TemperBoundaryImpl<DRT::Element::quad4>(numdofpernode);
      return cp4;
  }
  /*  case DRT::Element::quad8:
  {
    static TemperBoundaryImpl<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new TemperImpl<DRT::Element::quad8>(numdofpernode);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static TemperBoundaryImpl<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new TemperImpl<DRT::Element::quad9>(numdofpernode);
    return cp9;
  }*/
  case DRT::Element::tri3:
  {
    static TemperBoundaryImpl<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new TemperBoundaryImpl<DRT::Element::tri3>(numdofpernode);
      return cp3;
  }
  /*  case DRT::Element::tri6:
  {
    static TemperBoundaryImpl<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new TemperBoundaryImpl<DRT::Element::tri6>(numdofpernode);
    return cp6;
  }*/
  case DRT::Element::line2:
  {
    static TemperBoundaryImpl<DRT::Element::line2>* cl2;
    if (cl2==NULL)
      cl2 = new TemperBoundaryImpl<DRT::Element::line2>(numdofpernode);
      return cl2;
  }/*
  case DRT::Element::line3:
  {
    static TemperBoundaryImpl<DRT::Element::line3>* cl3;
    if (cl3==NULL)
      cl3 = new TemperBoundaryImpl<DRT::Element::line3>(numdofpernode);
    return cl3;
  }*/
  default:
    dserror("Shape %d (%d nodes) not supported", ele->Shape(), ele->NumNode());
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TemperBoundaryImpl<distype>::TemperBoundaryImpl(
  int numdofpernode
  )
: numdofpernode_(numdofpernode),
    xyze_(true),  // initialize to zero
    xsi_(true),
    funct_(true),
    deriv_(true),
    derxy_(true),
    normal_(true),
    fac_(0.0)
{
  return;
}

/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperBoundaryImpl<distype>::Evaluate(
  DRT::ELEMENTS::ThermoBoundary* ele,
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  std::vector<int>& lm,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{
  // First, do the things that are needed for all actions:
  // get the material (of the parent element)
  DRT::ELEMENTS::Thermo* parentele = ele->ParentElement();
  Teuchos::RCP<MAT::Material> mat = parentele->Material();

  // Now, check for the action parameter
  const std::string action = params.get<std::string>("action","none");
  if (action == "calc_normal_vectors")
  {
    // access the global vector
    const Teuchos::RCP<Epetra_MultiVector> normals
      = params.get< Teuchos::RCP<Epetra_MultiVector> >("normal vectors", Teuchos::null);
    if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

    // determine constant normal to this element
    GetConstNormal(normal_,xyze_);

    // loop over the element nodes
    for (int j=0;j<nen_;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      if (normals->Map().MyGID(nodegid) )
      { // OK, the node belongs to this processor

        // scaling to a unit vector is performed on the global level after
        // assembly of nodal contributions since we have no reliable information
        // about the number of boundary elements adjacent to a node
        for (int dim=0; dim<(nsd_+1); dim++)
        {
          normals->SumIntoGlobalValue(nodegid,dim,normal_(dim));
        }
      }
      // else: the node belongs to another processor; the ghosted
      //       element will contribute the right value on that proc
    }
  }
  else if (action == "integrate_shape_functions")
  {
    // NOTE: add area value only for elements which are NOT ghosted!
    const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
    IntegrateShapeFunctions(ele,params,elevec1_epetra,addarea);
  }
  // surface heat transfer boundary condition
  else if (action == "calc_thermo_fextconvection")
  {
    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

    // set views
    LINALG::Matrix<nen_,nen_> etang(elemat1_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_,1> efext(elevec1_epetra.A(),true);  // view only!

    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null)
      dserror("Cannot access condition 'ThermoConvections'");

    // access parameters of the condition
    const std::string* tempstate = cond->Get<std::string>("temperature state");
    double coeff = cond->GetDouble("coeff");
    const int curvenum = cond->GetInt("curve");
    const double time = params.get<double>("total time");

    double surtemp = cond->GetDouble("surtemp");
    // increase the surrounding temperature step by step (scale with a time curve)
    const int surtempcurvenum = cond->GetInt("surtempcurve");

    // find out whether we shall use a time curve and get the factor
    double curvefac = 1.0;
    if (curvenum>=0)
    {
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    }
    // multiply heat convection coefficient with the timecurve factor
    coeff *= curvefac;

    // we can increase or decrease the surrounding (fluid) temperature T_oo
    // enabling for instance a load cycle due to combustion of a fluid
    double surtempcurvefac = 1.0;
    // find out whether we shall use a time curve for T_oo and get the factor
    if (surtempcurvenum>=0)
    {
      surtempcurvefac = DRT::Problem::Instance()->Curve(surtempcurvenum).f(time);
    }
    // complete surrounding temperatures T_oo: multiply with the timecurve factor
    surtemp *= surtempcurvefac;

    // use current temperature T_{n+1} of current time step for heat convection
    // boundary condition
    if (*tempstate == "Tempnp")
    {
      // disassemble temperature
      if (discretization.HasState("temperature"))
      {
        // get actual values of temperature from global location vector
        std::vector<double> mytempnp(lm.size());
        Teuchos::RCP<const Epetra_Vector> tempnp
          = discretization.GetState("temperature");
        if (tempnp == Teuchos::null) dserror("Cannot get state vector 'tempnp'");

        DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,lm);
        // build the element temperature
        LINALG::Matrix<nen_,1> etemp(&(mytempnp[0]),true);  // view only!
        etemp_.Update(etemp);  // copy
      } // discretization.HasState("temperature")
      else
        dserror("No old temperature T_n+1 available");
    }
    // use temperature T_n of last known time step t_n
    else if (*tempstate == "Tempn")
    {
      if (discretization.HasState("old temperature"))
      {
        // get actual values of temperature from global location vector
        std::vector<double> mytempn(lm.size());
        Teuchos::RCP<const Epetra_Vector> tempn
          = discretization.GetState("old temperature");
        if (tempn == Teuchos::null) dserror("Cannot get state vector 'tempn'");

        DRT::UTILS::ExtractMyValues(*tempn,mytempn,lm);
        // build the element temperature
        LINALG::Matrix<nen_,1> etemp(&(mytempn[0]),true);  // view only!
        etemp_.Update(etemp);  // copy
      } // discretization.HasState("old temperature")
      else
        dserror("No old temperature T_n available");
    }
    else
      dserror("Unknown type of convection boundary condition");

#ifdef THRASOUTPUT
    if (ele->Id()==0)
    {
      cout << "ele Id= " << ele->Id() << endl;
      // print all parameters read from the current condition
      cout<<"type of boundary condition  = "<<*tempstate<<endl;
      cout<<"heat convection coefficient = "<<coeff<<endl;
      cout<<"surrounding temperature     = "<<surtemp<<endl;
      cout<<"time curve                  = "<<curvenum<<endl;
      cout<<"total time                  = "<<time<<endl;
    }
#endif // THRASOUTPUT

    // and now check if there is a convection heat transfer boundary condition
    EvaluateThermoConvection(
      ele, // current boundary element
      &etang, // element-matrix
      &efext, // element-rhs
      coeff,
      surtemp,
      *tempstate
      );

    // BUILD EFFECTIVE TANGENT AND RESIDUAL ACC TO TIME INTEGRATOR
    // check the time integrator
    const INPAR::THR::DynamicType timint
      = DRT::INPUT::get<INPAR::THR::DynamicType>(params, "time integrator",INPAR::THR::dyna_undefined);
    switch (timint)
    {
      case INPAR::THR::dyna_statics :
      {
        if (*tempstate == "Tempn")
          dserror("Old temperature T_n is not allowed with static time integrator");
        // continue
        break;
      }
      case INPAR::THR::dyna_onesteptheta :
      {
        // Note: efext is scaled with theta in thrtimint_ost.cpp. Because the
        // convective boundary condition is nonlinear and produces a term in the
        // tangent, consider the factor theta here, too
        const double theta = params.get<double>("theta");
        // combined tangent and conductivity matrix to one global matrix
        etang.Scale(theta);
        break;
      }
      case INPAR::THR::dyna_genalpha :
      {
        dserror("Genalpha not yet implemented");
        break;
      }
      case INPAR::THR::dyna_undefined :
      default :
      {
        dserror("Don't know what to do...");
        break;
      }
    }  // end of switch(timint)
  }  // calc_thermo_fextconvection
  else
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());

  return 0;
} // Evaluate


/*----------------------------------------------------------------------*
 | Evaluate the element in case of volume coupling           dano 02/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperBoundaryImpl<distype>::Evaluate(
  DRT::ELEMENTS::ThermoBoundary* ele,
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{
  //if (la.Size()==1)
  {
    return Evaluate(
      ele,
      params,
      discretization,
      la[0].lm_, // location vector is build by the first column of la
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }

  // First, do the things that are needed for all actions:
  // get the material (of the parent element)
  DRT::ELEMENTS::Thermo* parentele = ele->ParentElement();
  Teuchos::RCP<MAT::Material> mat = parentele->Material();

  // Now, check for the action parameter
  const std::string action = params.get<std::string>("action","none");
  if (action == "calc_normal_vectors")
  {
    // access the global vector
    const Teuchos::RCP<Epetra_MultiVector> normals
      = params.get< Teuchos::RCP<Epetra_MultiVector> >("normal vectors", Teuchos::null);
    if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

    // determine constant normal to this element
    GetConstNormal(normal_,xyze_);

    // loop over the element nodes
    for (int j=0;j<nen_;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      if (normals->Map().MyGID(nodegid) )
      { // OK, the node belongs to this processor

        // scaling to a unit vector is performed on the global level after
        // assembly of nodal contributions since we have no reliable information
        // about the number of boundary elements adjacent to a node
        for (int dim=0; dim<(nsd_+1); dim++)
        {
          normals->SumIntoGlobalValue(nodegid,dim,normal_(dim));
        }
      }
      // else: the node belongs to another processor; the ghosted
      //       element will contribute the right value on that proc
    }
  }
  else if (action == "integrate_shape_functions")
  {
    // NOTE: add area value only for elements which are NOT ghosted!
    const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
    IntegrateShapeFunctions(ele,params,elevec1_epetra,addarea);
  }
  // surface heat convection boundary condition
  else if (action == "calc_thermo_fextconvection")
  {
    dserror("EvaluateCondition is used with the location vector lm");
  }
  else
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());

  return 0;
} // Evaluate


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperBoundaryImpl<distype>::EvaluateNeumann(
  DRT::Element* ele,
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Condition& condition,
  std::vector<int>& lm,
  Epetra_SerialDenseVector& elevec1
  )
{
  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());

    // multiply integration factor with the timecurve factor
    // fac_ = fac_ * curvefac
    fac_ *= curvefac;

    // factor given by spatial function
    double functfac = 1.0;
    // determine global coordinates of current Gauss point
    double coordgp[3]; // coordinate has always to be given in 3D!
    for (int i = 0; i< 3; i++)
      coordgp[i] = 0.0;

    for (int i = 0; i< nsd_; i++)
    {
      for (int j = 0; j < nen_; j++)
      {
        // node coordinate * shape function
        coordgp[i] += xyze_(i,j) * funct_(j);
      }
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    for(int dof=0;dof<numdofpernode_;dof++)
    {
      if ((*onoff)[dof]) // is this dof activated?
      {
        // factor given by spatial function
        if (func) functnum = (*func)[dof];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,coordgpref,time,NULL);
          }
          else
            functfac = 1.0;
        }

        const double val_fac_functfac = (*val)[dof]*fac_*functfac;

        for (int node=0;node<nen_;++node)
        {
          // fext  = fext +  N^T  * q * detJ * w(gp) * spatial_fac * timecurve_fac
          elevec1[node*numdofpernode_+dof] += funct_(node)*val_fac_functfac;
        }
      } // if ((*onoff)[dof])
    }
  } //end of loop over integration points

  return 0;
}

/*----------------------------------------------------------------------*
 | evaluate an convective thermo boundary condition (private) dano 12/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::EvaluateThermoConvection(
  const DRT::Element* ele,
  LINALG::Matrix<nen_,nen_>* etang,
  LINALG::Matrix<nen_,1>* efext,
  const double coeff,
  const double surtemp,
  const std::string tempstate
  )
{
  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // integration loop
  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int iquad=0; iquad<intpoints.IP().nquad; iquad++)
  {
    EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());
    // fac_ = Gauss weight * det(J) is calculated in EvalShapeFuncAndIntFac()

    // ------------right-hand-side
    // q . n = h ( T - T_sur )

    // multiply fac_ * coeff
    // --> must be insert in balance equation es positive term,
    // but fext is included as negative --> scale with (-1)
    double coefffac_ = fac_ * coeff;

    // get the current temperature
    // Theta = Ntemp = N . T
    // caution: funct_ implemented as (4,1)--> use transposed in code for
    // theoretic part
    // funct_ describes a 2D area, for hex8: 4 nodes
    // (1x1)= (1x4)(4x1) = (nen_*numdofpernode_ x 1)^T(nen_*numdofpernode_ x 1)
    LINALG::Matrix<1,1> Ntemp(true);
    Ntemp.MultiplyTN(funct_,etemp_);

    // - T_surf * I
    LINALG::Matrix<1,1> Tsurf(true);
    for (int i=0; i<1; ++i)
    {
      Tsurf(i) = (surtemp);
    }

    // Ntemp = T - T_surf
    Ntemp.Update(-1.0,Tsurf,1.0);

    // ---------------------------------------------- right-hand-side
    if (efext != NULL)
    {
      // efext = efext - N^T . coeff . ( N . T - T_surf) . detJ * w(gp)
      // in energy balance: q_c positive, but fext = r + q^bar  q_c^bar
      // we want to define q_c = - q . n
      // q_c^bar = k . Grad T = h (T - T_oo), vgl. Farhat(1992)
      efext->Multiply(coefffac_,funct_,Ntemp,1.0);
    }

    // ------------------------------------------------------ tangent
    // if current temperature is used in boundary condition, additional term in
    // thermal tangent occurs
    if (tempstate == "Tempnp")
    {
      // if the boundary condition shall be dependent on the current temperature
      // solution T_n+1 --> linearisation must be considered
      // ---------------------matrix
      if (etang != NULL)
      {
        // ke = ke + (N^T . coeff . N) * detJ * w(gp)
        etang->MultiplyNT((-1.0)*coefffac_,funct_,funct_,1.0);
      }
    }

   /* =======================================================================*/
  }/* ================================================== end of Loop over GP */
   /* =======================================================================*/

  return;

} // TemperBoundaryImpl<distype>::EvaluateThermoConvection()



/*----------------------------------------------------------------------*
 |  evaluate shape functions and int. factor at int. point    gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::EvalShapeFuncAndIntFac(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
  const int& iquad,  ///< id of current Gauss point
  const int& eleid  ///< the element id
  )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // the metric tensor and the area of an infinitesimal surface/line element
  double drs(0.0);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs);

  // set the integration factor
  fac_ = intpoints.IP().qwgt[iquad] * drs;

  // say goodbye
  return;
}


/*----------------------------------------------------------------------*
 |  get constant normal                                       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::GetConstNormal(
  LINALG::Matrix<nsd_+1,1>& normal,
  const LINALG::Matrix<nsd_+1,nen_>& xyze
  )
{
  // determine normal to this element
  switch(nsd_)
  {
    case 2:
    {
      LINALG::Matrix<3,1> dist1(true), dist2(true);
      for (int i=0; i<3; i++)
      {
        dist1(i) = xyze(i,1)-xyze(i,0);
        dist2(i) = xyze(i,2)-xyze(i,0);
      }

      normal(0) = dist1(1)*dist2(2) - dist1(2)*dist2(1);
      normal(1) = dist1(2)*dist2(0) - dist1(0)*dist2(2);
      normal(2) = dist1(0)*dist2(1) - dist1(1)*dist2(0);
    }
    break;
    case 1:
    {
      normal(0) = xyze(1,1) - xyze(1,0);
      normal(1) = (-1.0)*(xyze(0,1) - xyze(0,0));
    }
    break;
    default:
      dserror("Illegal number of space dimensions: %d",nsd_);
  } // switch(nsd)

  // length of normal to this element
  const double length = normal.Norm2();
  // outward-pointing normal of length 1.0
  normal.Scale(1/length);

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)           gjb 02/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::IntegrateShapeFunctions(
  const DRT::Element* ele,
  Teuchos::ParameterList& params,
  Epetra_SerialDenseVector& elevec1,
  const bool addarea
  )
{
  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("boundaryint");

  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad=0; iquad<intpoints.IP().nquad; iquad++)
  {
    EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());

    // compute integral of shape functions
    for (int node=0; node<nen_; ++node)
    {
      for (int k=0; k< numdofpernode_; k++)
      {
        elevec1[node*numdofpernode_+k] += funct_(node) * fac_;
      }
    }

    if (addarea)
    {
      // area calculation
      boundaryint += fac_;
    }

  } //loop over integration points

  // add contribution to the global value
  params.set<double>("boundaryint",boundaryint);

  return;

} //TemperBoundaryImpl<distype>::IntegrateShapeFunction


/*----------------------------------------------------------------------*/
#endif // D_THERMO
