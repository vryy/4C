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
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef D_THERMO

#include <cstdlib>
#include "thermo_ele_boundary_impl.H"
#include "thermo_ele_impl.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/linalg_serialdensematrix.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_lib/drt_function.H"


/*----------------------------------------------------------------------*
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
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TemperBoundaryImpl<distype>::TemperBoundaryImpl
(int numdofpernode)
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
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperBoundaryImpl<distype>::Evaluate
(
  DRT::ELEMENTS::ThermoBoundary*    ele,
  Teuchos::ParameterList&           params,
  DRT::Discretization&              discretization,
  std::vector<int>&                 lm,
  Epetra_SerialDenseMatrix&         elemat1_epetra,
  Epetra_SerialDenseMatrix&         elemat2_epetra,
  Epetra_SerialDenseVector&         elevec1_epetra,
  Epetra_SerialDenseVector&         elevec2_epetra,
  Epetra_SerialDenseVector&         elevec3_epetra
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
    const Teuchos::RCP<Epetra_MultiVector> normals = params.get< Teuchos::RCP<Epetra_MultiVector> >("normal vectors", Teuchos::null);
    if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

    // get node coordinates (we have a nsd_+1 dimensional domain!)
    GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

    // determine constant normal to this element
    GetConstNormal(normal_,xyze_);

    // loop over the element nodes
    for (int j=0;j<iel;j++)
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
  else if (action =="integrate_shape_functions")
  {
    // NOTE: add area value only for elements which are NOT ghosted!
    const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
    IntegrateShapeFunctions(ele,params,elevec1_epetra,addarea);
  }
  else
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());

  return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperBoundaryImpl<distype>::EvaluateNeumann
(
  DRT::Element*             ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  DRT::Condition&           condition,
  std::vector<int>&         lm,
  Epetra_SerialDenseVector& elevec1
)
{
  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

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
    fac_ *= curvefac;

    // factor given by spatial function
    double functfac = 1.0;
    // determine global coordinates of current Gauss point
    double coordgp[nsd_];
    for (int i = 0; i< nsd_; i++)
    {
      coordgp[i] = 0.0;
      for (int j = 0; j < iel; j++)
      {
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
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,coordgpref,0.0,NULL);
          }
          else
            functfac = 1.0;
        }

        const double val_fac_functfac = (*val)[dof]*fac_*functfac;

        for (int node=0;node<iel;++node)
        {
          elevec1[node*numdofpernode_+dof] += funct_(node)*val_fac_functfac;
        }
      } // if ((*onoff)[dof])
    }
  } //end of loop over integration points

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate shape functions and int. factor at int. point    gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::EvalShapeFuncAndIntFac
(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
  const int&                                   iquad,      ///< id of current Gauss point
  const int&                                   eleid       ///< the element id
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
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::GetConstNormal
(
  LINALG::Matrix<nsd_+1,1>&          normal,
  const LINALG::Matrix<nsd_+1,iel>&  xyze
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
} // TemperBoundaryImpl<distype>::

/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)           gjb 02/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperBoundaryImpl<distype>::IntegrateShapeFunctions
(
  const DRT::Element*        ele,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseVector&  elevec1,
  const bool                 addarea
)
{
  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("boundaryint");

  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,iel> >(ele,xyze_);

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of shape functions
    for (int node=0;node<iel;++node)
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


#endif // D_THERMO
#endif // CCADISCRET
