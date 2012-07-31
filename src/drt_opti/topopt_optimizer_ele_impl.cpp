/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_impl.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer_ele_impl.H"
#include "topopt_optimizer_ele_parameter.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_mat/optimization_density.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptImplInterface* DRT::ELEMENTS::TopOptImplInterface::Impl(
  const DRT::Element* ele
)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    return TopOptImpl<DRT::Element::hex8>::Instance();
  }
  case DRT::Element::hex20:
  {
    return TopOptImpl<DRT::Element::hex20>::Instance();
  }
  case DRT::Element::hex27:
  {
    return TopOptImpl<DRT::Element::hex27>::Instance();
  }
  case DRT::Element::tet4:
  {
    return TopOptImpl<DRT::Element::tet4>::Instance();
  }
  case DRT::Element::tet10:
  {
    return TopOptImpl<DRT::Element::tet10>::Instance();
  }
  case DRT::Element::quad4:
  {
    return TopOptImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return TopOptImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return TopOptImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return TopOptImpl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return TopOptImpl<DRT::Element::tri6>::Instance();

  case DRT::Element::line2:
  {
    return TopOptImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return TopOptImpl<DRT::Element::line3>::Instance();
  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(ele->Shape()).c_str());
  }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptImpl<distype> * DRT::ELEMENTS::TopOptImpl<distype>::Instance(
    bool create)
{
  static TopOptImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new TopOptImpl<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptImpl<distype>::TopOptImpl()
: xyze_(true),
funct_(true),
deriv_(true),
derxy_(true),
efluidvel_(true),
eadjointvel_(true),
fluidvelint_(true),
adjointvelint_(true),
fluidvelxy_(true),
poroint_(0.0),
intpoints_( distype ),
xsi_(true),
det_(0.0),
fac_(0.0),
visc_(0.0),
reacoeff_(0.0),
dens_(0.0),
is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateObjective(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  vector<int>&               lm
)
{
  return EvaluateObjective(
      ele,
      params,
      optidis,
      mat,
      lm,
      intpoints_
  );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateObjective(
  DRT::Element*                 ele,
  ParameterList&                params,
  DRT::Discretization&          optidis,
  RCP<MAT::Material>            mat,
  vector<int>&                  lm,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  double& objective = params.get<double>("objective_value");

  RCP<DRT::Discretization> fluiddis = params.get<RCP<DRT::Discretization> >("fluiddis");

  RCP<map<int,RCP<Epetra_Vector> > > fluidvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("fluidvel");

  map<int,LINALG::Matrix<nsd_,nen_> > efluidvels;

  LINALG::Matrix<nsd_,nen_> efluidvel;

  // extract element data of all time steps from fluid solution
  vector<int> fluidlm;
  {
    vector<int> lmowner; // dummy for function call
    vector<int> lmstride; // dummy for function call
    ele->LocationVector(*fluiddis,fluidlm,lmowner,lmstride);
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=fluidvels->begin();
      i!=fluidvels->end();i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis,fluidlm,&efluidvel,NULL,i->second);
    efluidvels.insert(pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,efluidvel));
  }

  RCP<const Epetra_Vector> dens = optidis.GetState("density");
  LINALG::Matrix<nen_,1> edens(true);

  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*dens,mymatrix,lm);
  for (int inode=0; inode<nen_; ++inode) edens(inode,0) = mymatrix[inode];

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);


  Objective(
      ele->Id(),
      efluidvels,
      edens,
      objective,
      mat,
      intpoints
  );

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Objective(
  const int eid,
  map<int,LINALG::Matrix<nsd_,nen_> >& efluidvel,
  LINALG::Matrix<nen_,1>& edens,
  double& objective,
  RCP<MAT::Material> mat,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  double value = 0.0; // total unscaled entries at one gauss point

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    EvalPorosityAtIntPoint(edens);

    // dissipation in objective present?
    if (optiparams_->ObjDissipationTerm())
    {
      switch (optiparams_->TimeIntScheme())
      {
      case INPAR::FLUID::timeint_stationary:
      case INPAR::FLUID::timeint_one_step_theta: // handle these two cases together
      {
        for (int timestep=0;timestep<=optiparams_->NumTimesteps();timestep++)
        {
          if (optiparams_->IsStationary())
            timestep=1; // just one stationary time step 1

          efluidvel_ = efluidvel.find(timestep)->second;

          fluidvelint_.Multiply(efluidvel_,funct_);
          fluidvelxy_.MultiplyNT(efluidvel_,derxy_);

          // weighting of the timesteps
          // if stationary, we have the second case with theta = 1, so all is ok
          double timefac = 0.0;
          if (timestep==0) // first time step -> factor 1-theta (old sol at first time step)
            timefac = 1.0 - optiparams_->Theta();
          else if (timestep==optiparams_->NumTimesteps()) // last time step -> factor theta (new sol at last time step)
            timefac = optiparams_->Theta();
          else // all other time steps -> factor 1-theta as old sol, factor theta as new sol -> overall factor 1
            timefac = 1.0;

          for (int idim=0;idim<nsd_;idim++)
          {
            value += poroint_*fluidvelint_.Dot(fluidvelint_);

            for (int jdim=0;jdim<nsd_;jdim++)
            {
              value += visc_*fluidvelxy_(idim,jdim)*(fluidvelxy_(idim,jdim)+fluidvelxy_(jdim,idim));
            }
          }
        }
      }
      break;
      default:
        dserror("unknown time integration scheme while evaluating objective gradient");
      }

      objective += optiparams_->Dt()*fac_*optiparams_->ObjDissipationFac()*value;
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradient(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  vector<int>&               lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{
  return EvaluateGradient(
      ele,
      params,
      optidis,
      mat,
      lm,
      elevec1_epetra,
      intpoints_
  );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradient(
  DRT::Element*                 ele,
  ParameterList&                params,
  DRT::Discretization&          optidis,
  RCP<MAT::Material>            mat,
  vector<int>&                  lm,
  Epetra_SerialDenseVector&     elevec1_epetra,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  LINALG::Matrix<nen_,1> egrad(elevec1_epetra,true);

  RCP<DRT::Discretization> fluiddis = params.get<RCP<DRT::Discretization> >("fluiddis");

  RCP<map<int,RCP<Epetra_Vector> > > fluidvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("fluidvel");
  RCP<map<int,RCP<Epetra_Vector> > > adjointvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("adjointvel");

  map<int,LINALG::Matrix<nsd_,nen_> > efluidvels;
  map<int,LINALG::Matrix<nsd_,nen_> > eadjointvels;

  LINALG::Matrix<nsd_,nen_> efluidvel;
  LINALG::Matrix<nsd_,nen_> eadjointvel;

  // extract element data of all time steps from fluid solution
  vector<int> fluidlm;
  {
    vector<int> lmowner; // dummy for function call
    vector<int> lmstride; // dummy for function call
    ele->LocationVector(*fluiddis,fluidlm,lmowner,lmstride);
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=fluidvels->begin();
      i!=fluidvels->end();i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis,fluidlm,&efluidvel,NULL,i->second);
    efluidvels.insert(pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,efluidvel));
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=adjointvels->begin();
      i!=adjointvels->end();i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis,fluidlm,&eadjointvel,NULL,i->second);
    eadjointvels.insert(pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,eadjointvel));
  }

  RCP<const Epetra_Vector> dens = optidis.GetState("density");
  LINALG::Matrix<nen_,1> edens(true);

  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*dens,mymatrix,lm);
  for (int inode=0; inode<nen_; ++inode) edens(inode,0) = mymatrix[inode];

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);


  Gradient(
      ele->Id(),
      efluidvels,
      eadjointvels,
      edens,
      egrad,
      mat,
      intpoints
  );

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Gradient(
  const int eid,
  map<int,LINALG::Matrix<nsd_,nen_> >& efluidvel,
  map<int,LINALG::Matrix<nsd_,nen_> >& eadjointvel,
  LINALG::Matrix<nen_,1>& edens,
  LINALG::Matrix<nen_,1>& egrad,
  RCP<MAT::Material> mat,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  double value = 0.0;

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    EvalPorosityAtIntPoint(edens);

    // dissipation in objective present?
    double dissipation_fac = 0.0;
    if (optiparams_->ObjDissipationTerm())
      dissipation_fac = optiparams_->ObjDissipationFac();


    switch (optiparams_->TimeIntScheme())
    {
    case INPAR::FLUID::timeint_stationary:
    case INPAR::FLUID::timeint_one_step_theta: // handle these two cases together
    {
      for (int timestep=0;timestep<=optiparams_->NumTimesteps();timestep++)
      {
        if (optiparams_->IsStationary())
          timestep=1; // just one stationary time step 1

        efluidvel_ = efluidvel.find(timestep)->second;
        eadjointvel_ = eadjointvel.find(timestep)->second;

        fluidvelint_.Multiply(efluidvel_,funct_);
        adjointvelint_.Multiply(eadjointvel_,funct_);

        double timefac = 0.0;
        if (timestep==0)                                timefac = 1.0 - optiparams_->Theta();
        else if (timestep==optiparams_->NumTimesteps()) timefac = optiparams_->Theta();
        else                                            timefac = 1.0;

        for (int idim=0;idim<nsd_;idim++)
        {
          value += timefac*(
              dissipation_fac*fluidvelint_.Dot(fluidvelint_) // dissipation part
              -fluidvelint_.Dot(adjointvelint_)); // adjoint part
        }
      }
    }
    break;
    default:
      dserror("unknown time integration scheme while evaluating objective gradient");
    }

    value*= fac_*optiparams_->Dt()*poroderdens_; // scale with integration factor, time step size and poro-derivative

    for (int vi=0;vi<nen_;vi++)
    {
      egrad(vi,0) += value*funct_(vi);
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  return;
}



/*----------------------------------------------------------------------*
  | evaluate shape functions and derivatives at int. point     gjb 08/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
  DRT::UTILS::GaussIntegration::iterator&      iquad,      ///< actual integration point
  const int                                    eleid       ///< the element id
  )
{
  // coordinates of the current integration point
  const double* gpcoord = iquad.Point();
  for (int idim=0;idim<nsd_;idim++)
     xsi_(idim) = gpcoord[idim];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // compute Jacobian matrix and determinant
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
  LINALG::Matrix<nsd_,nsd_> xjm(true); // jacobian dx/ds
  LINALG::Matrix<nsd_,nsd_> xij(true); // invers transposed jacobian ds/dx

  xjm.MultiplyNT(deriv_,xyze_);
  const double det = xij.Invert(xjm);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);

  // compute integration factor
  fac_ = iquad.Weight()*det_;

  // compute global derivatives
  derxy_.Multiply(xij,deriv_);

} //TopOptImpl::CalcSubgrVelocity



/*---------------------------------------------------------------------------------*
 | evaluate porosity at gauss point                               winklmaier 05/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::EvalPorosityAtIntPoint(
    LINALG::Matrix<nen_,1> edens
)
{
  double densint = edens.Dot(funct_);
  poroint_ = optiparams_->MaxPoro() + (optiparams_->MinPoro()-optiparams_->MaxPoro())*densint
                                      *(1+optiparams_->SmearFac())/(densint+optiparams_->SmearFac());

  poroderdens_ = (optiparams_->MinPoro()-optiparams_->MaxPoro())
                *(optiparams_->SmearFac()+optiparams_->SmearFac()*optiparams_->SmearFac())
                /(densint+optiparams_->SmearFac());
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::ExtractValuesFromGlobalVector(
    DRT::Discretization&         discretization, ///< discretization
    const vector<int>&           lm,                  ///<
    LINALG::Matrix<nsd_,nen_> *  matrixtofill,        ///< vector field
    LINALG::Matrix<nen_,1> *     vectortofill,        ///< scalar field
    RCP<Epetra_Vector>&          globalvector         ///< global vector
) const
{
  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*globalvector,mymatrix,lm);

  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*(nsd_+1))];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*(nsd_+1))];
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptBoundaryImplInterface* DRT::ELEMENTS::TopOptBoundaryImplInterface::Impl(
    const DRT::Element* ele
)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return TopOptBoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return TopOptBoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return TopOptBoundaryImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return TopOptBoundaryImpl<DRT::Element::tri3>::Instance();
  }
  /*  case DRT::Element::tri6:
  {
    return TopOptBoundaryImpl<DRT::Element::tri6>::Instance();
  }*/
  case DRT::Element::line2:
  {
    return TopOptBoundaryImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return TopOptBoundaryImpl<DRT::Element::line3>::Instance();
  }
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptBoundaryImpl<distype> * DRT::ELEMENTS::TopOptBoundaryImpl<distype>::Instance(
    bool create
)
{
  static TopOptBoundaryImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new TopOptBoundaryImpl<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptBoundaryImpl<distype>::TopOptBoundaryImpl
()
: intpoints_( distype ),
xsi_(true),
det_(0.0),
fac_(0.0),
visc_(0.0),
reacoeff_(0.0),
dens_(0.0),
is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryObjective(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  vector<int>&               lm
  )
{
  // TODO coming...

  // work is done
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryGradient(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  vector<int>&               lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{
  // TODO coming...

  // work is done
  return 0;
}





///*----------------------------------------------------------------------*
// | evaluate shape functions and int. factor at int. point     gjb 01/09 |
// *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//double DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvalShapeFuncAndIntFac(
//    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
//    const int                                    iquad,      ///< id of current Gauss point
//    const int                                    eleid,      ///< the element id
//    LINALG::Matrix<1 + nsd_,1>*         normalvec ///< normal vector at Gauss point(optional)
//)
//{
//  // coordinates of the current integration point
//  const double* gpcoord = (intpoints.IP().qxg)[iquad];
//  for (int idim=0;idim<nsd_;idim++)
//  {xsi_(idim) = gpcoord[idim];}
//
//  if(not DRT::NURBS::IsNurbs(distype))
//  {
//    // shape functions and their first derivatives
//    DRT::UTILS::shape_function<distype>(xsi_,funct_);
//    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
//  }
//  else // nurbs elements are always somewhat special...
//  {
//    DRT::NURBS::UTILS::nurbs_get_funct_deriv(
//        funct_  ,
//        deriv_  ,
//        xsi_    ,
//        myknots_,
//        weights_,
//        distype );
//  }
//
//  // the metric tensor and the area of an infinitesimal surface/line element
//  // optional: get normal at integration point as well
//  // Note: this is NOT yet a unit normal. Its norm corresponds to the area/length of the element
//  double drs(0.0);
//  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs,normalvec);
//
//  // for nurbs elements the normal vector must be scaled with a special orientation factor!!
//  if(DRT::NURBS::IsNurbs(distype))
//  {
//    if (normalvec != NULL)
//      normal_.Scale(normalfac_);
//  }
//
//  // return the integration factor
//  return intpoints.IP().qwgt[iquad] * drs;
//}
//
//
//
///*----------------------------------------------------------------------*
// |  Integrate shapefunctions over surface (private)           gjb 02/09 |
// *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::IntegrateShapeFunctions(
//    const DRT::Element*        ele,
//    ParameterList&             params,
//    Epetra_SerialDenseVector&  elevec1,
//    const bool                 addarea
//)
//{
//  // access boundary area variable with its actual value
//  double boundaryint = params.get<double>("boundaryint");
//
//  // integrations points and weights
//  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
//
//  // loop over integration points
//  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
//  {
//    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());
//
//    // compute integral of shape functions
//    for (int node=0;node<nen_;++node)
//    {
//      for (int k=0; k< numscal_; k++)
//      {
//        elevec1[node*numdofpernode_+k] += funct_(node) * fac;
//      }
//    }
//
//    if (addarea)
//    {
//      // area calculation
//      boundaryint += fac;
//    }
//
//  } //loop over integration points
//
//  // add contribution to the global value
//  params.set<double>("boundaryint",boundaryint);
//
//  return;
//
//} //TopOptBoundaryImpl<distype>::IntegrateShapeFunction



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::ExtractValuesFromGlobalVector(
    DRT::Discretization&         discretization, ///< discretization
    const vector<int>&           lm,                  ///<
    LINALG::Matrix<nsd_,nen_> *  matrixtofill,        ///< vector field
    LINALG::Matrix<nen_,1> *     vectortofill,        ///< scalar field
    RCP<Epetra_Vector>&          globalvector         ///< global vector
) const
{
  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*globalvector,mymatrix,lm);

  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*(nsd_+1))];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*(nsd_+1))];
  }
}


