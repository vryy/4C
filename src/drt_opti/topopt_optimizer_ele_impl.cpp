/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_impl.cpp

\brief element routines of the topology optimization element

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
#include "../linalg/linalg_utils.H"


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
fac_(0.0),
is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateValues(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  std::vector<int>&          lm
)
{
  return EvaluateValues(
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
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateValues(
  DRT::Element*                 ele,
  Teuchos::ParameterList&       params,
  DRT::Discretization&          optidis,
  RCP<MAT::Material>            mat,
  std::vector<int>&             lm,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  double& objective = params.get<double>("objective_value");

  Teuchos::RCP<Epetra_SerialDenseVector> constraints = params.get<Teuchos::RCP<Epetra_SerialDenseVector> >("constraint_values");

  RCP<DRT::Discretization> fluiddis = params.get<RCP<DRT::Discretization> >("fluiddis");

  RCP<map<int,RCP<Epetra_Vector> > > fluidvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("fluidvel");

  map<int,LINALG::Matrix<nsd_,nen_> > efluidvels;

  LINALG::Matrix<nsd_,nen_> efluidvel;

  // extract element data of all time steps from fluid solution
  vector<int> fluidlm;
  DRT::UTILS::DisBasedLocationVector(*fluiddis,*ele,fluidlm,nsd_+1);

  for (map<int,RCP<Epetra_Vector> >::iterator i=fluidvels->begin();
      i!=fluidvels->end();i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis,fluidlm,&efluidvel,NULL,i->second);
    efluidvels.insert(std::pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,efluidvel));
  }

  RCP<const Epetra_Vector> dens = optidis.GetState("density");
  LINALG::Matrix<nen_,1> edens(true);

  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*dens,mymatrix,lm);
  for (int inode=0; inode<nen_; ++inode) edens(inode,0) = mymatrix[inode];

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);


  Values(
      ele->Id(),
      efluidvels,
      edens,
      objective,
      constraints,
      mat,
      intpoints
  );

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Values(
  const int eid,
  map<int,LINALG::Matrix<nsd_,nen_> >& efluidvel,
  LINALG::Matrix<nen_,1>& edens,
  double& objective,
  Teuchos::RCP<Epetra_SerialDenseVector> constraints,
  RCP<MAT::Material> mat,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  double& densint = (*constraints)[0]; // integrated density

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    double value = 0.0;

    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    EvalPorosityAtIntPoint(edens);

    // volume constraint
    densint += fac_*(edens.Dot(funct_) - optiparams_->VolBd()); // TODO c missing

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

          value += timefac*poroint_*fluidvelint_.Dot(fluidvelint_);
          for (int idim=0;idim<nsd_;idim++)
          {
            for (int jdim=0;jdim<nsd_;jdim++)
            {
              value += timefac*optiparams_->Viscosity()*fluidvelxy_(idim,jdim)*(fluidvelxy_(idim,jdim)+fluidvelxy_(jdim,idim));
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
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradients(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  std::vector<int>&          lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{
  return EvaluateGradients(
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
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradients(
  DRT::Element*                 ele,
  Teuchos::ParameterList&       params,
  DRT::Discretization&          optidis,
  RCP<MAT::Material>            mat,
  std::vector<int>&             lm,
  Epetra_SerialDenseVector&     elevec1_epetra,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  LINALG::Matrix<nen_,1> egrad(elevec1_epetra,true);

  Epetra_SerialDenseVector constr_deriv(nen_);
  LINALG::Matrix<nen_,1> econstr_der(constr_deriv,true); // volume constraint

  RCP<DRT::Discretization> fluiddis = params.get<RCP<DRT::Discretization> >("fluiddis");

  RCP<map<int,RCP<Epetra_Vector> > > fluidvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("fluidvel");
  RCP<map<int,RCP<Epetra_Vector> > > adjointvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("adjointvel");

  map<int,LINALG::Matrix<nsd_,nen_> > efluidvels;
  map<int,LINALG::Matrix<nsd_,nen_> > eadjointvels;

  LINALG::Matrix<nsd_,nen_> efluidvel;
  LINALG::Matrix<nsd_,nen_> eadjointvel;

  // extract element data of all time steps from fluid solution
  vector<int> fluidlm;
  DRT::UTILS::DisBasedLocationVector(*fluiddis,*ele,fluidlm,nsd_+1);

  for (map<int,RCP<Epetra_Vector> >::iterator i=fluidvels->begin();
      i!=fluidvels->end();i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis,fluidlm,&efluidvel,NULL,i->second);
    efluidvels.insert(std::pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,efluidvel));
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=adjointvels->begin();
      i!=adjointvels->end();i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis,fluidlm,&eadjointvel,NULL,i->second);
    eadjointvels.insert(std::pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,eadjointvel));
  }

  RCP<const Epetra_Vector> dens = optidis.GetState("density");
  LINALG::Matrix<nen_,1> edens(true);

  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*dens,mymatrix,lm);
  for (int inode=0; inode<nen_; ++inode) edens(inode,0) = mymatrix[inode];

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);


  Gradients(
      ele->Id(),
      efluidvels,
      eadjointvels,
      edens,
      egrad,
      econstr_der,
      mat,
      intpoints
  );

  // the derivation of the constraints is not handled by the standard assembly process
  // since it is a MultiVector. Thus it is handled manually here
  RCP<Epetra_MultiVector> constr_der = params.get<RCP<Epetra_MultiVector> >("constraints_derivations");

  vector<int> dummylm; // the same as lm
  vector<int> dummylmstride; // not required
  vector<int> lmowner; // owners of the lm-dofs
  ele->LocationVector(optidis,dummylm,lmowner,dummylmstride);
  if (dummylm!=lm) dserror("non fitting local maps which shall be identical");

  LINALG::Assemble(*constr_der,0,constr_deriv,lm,lmowner);

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Gradients(
  const int eid,
  map<int,LINALG::Matrix<nsd_,nen_> >& efluidvel,
  map<int,LINALG::Matrix<nsd_,nen_> >& eadjointvel,
  LINALG::Matrix<nen_,1>& edens,
  LINALG::Matrix<nen_,1>& egrad,
  LINALG::Matrix<nen_,1>& econstr_der,
  RCP<MAT::Material> mat,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    double value = 0.0;

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
      econstr_der(vi,0) += fac_*funct_(vi);
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
  fac_ = iquad.Weight()*det;

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
  // poro = poro_max + (poro_min - poro_max) * rho * (1+fac) / (rho+fac)
  double densint = edens.Dot(funct_);
  poroint_ = optiparams_->MaxPoro() + (optiparams_->MinPoro()-optiparams_->MaxPoro())*densint
                                      *(1+optiparams_->SmearFac())/(densint+optiparams_->SmearFac());

  // dporo/ddens = (poro_min - poro_max) * (fac+fac*fac) / (rho+fac)
  poroderdens_ = (optiparams_->MinPoro()-optiparams_->MaxPoro())
                *(optiparams_->SmearFac()+optiparams_->SmearFac()*optiparams_->SmearFac())
                /pow(densint+optiparams_->SmearFac(),2);
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::ExtractValuesFromGlobalVector(
    DRT::Discretization&         discretization, ///< discretization
    const std::vector<int>&      lm,                  ///<
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
  if ( create ) // create instance if not present
  {
    if ( instance==NULL )
    {
      instance = new TopOptBoundaryImpl<distype>();
    }
  }
  else // delete instance if present
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
fac_(0.0),
is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryValues(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  std::vector<int>&          lm
  )
{
  // TODO coming...

  // work is done
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryGradients(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       optidis,
  RCP<MAT::Material>         mat,
  std::vector<int>&          lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{
  // TODO coming...

  // work is done
  return 0;
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::ExtractValuesFromGlobalVector(
    DRT::Discretization&         discretization, ///< discretization
    const std::vector<int>&      lm,                  ///<
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


