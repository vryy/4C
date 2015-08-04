/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_poro.cpp

\brief evaluation of scatra boundary terms at integration points

 <pre>
   Maintainer: Moritz Thon & Andre Hemmler
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 */
/*----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_poro.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_action.H"
#include "scatra_ele.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"

/*----------------------------------------------------------------------*
 |  Singleton access method                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype> * DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    bool create
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcPoro<distype>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcPoro<distype>(numdofpernode,numscal,disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcPoro<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
     {
      delete i->second;
      i->second = NULL;
     }

    instances.clear();
    return NULL;
  }

  return instances[disname];
}

/*----------------------------------------------------------------------*
 |  Clean up                                              hemmler 07/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( 0, 0, "", false );
}


/*----------------------------------------------------------------------*
 |  Private constructor                                   hemmler 07/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::ScaTraEleBoundaryCalcPoro(const int numdofpernode, const int numscal,const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ScaTraEleBoundaryCalc(numdofpernode,numscal,disname),
    xyze0_(true),
    eporosity_(true),
    eprenp_(true),
    isnodalporosity_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::EvaluateAction(
    DRT::FaceElement*                   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    SCATRA::BoundaryAction              action,
    DRT::Element::LocationArray&        la,
    Epetra_SerialDenseMatrix&           elemat1_epetra,
    Epetra_SerialDenseMatrix&           elemat2_epetra,
    Epetra_SerialDenseVector&           elevec1_epetra,
    Epetra_SerialDenseVector&           elevec2_epetra,
    Epetra_SerialDenseVector&           elevec3_epetra
    )
{
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::bd_calc_fps3i_surface_permeability:
  case SCATRA::bd_calc_fs3i_surface_permeability:
  case SCATRA::bd_calc_Neumann:
  case SCATRA::bd_calc_normal_vectors:
  case SCATRA::bd_integrate_shape_functions:
  {
    my::EvaluateAction(
        ele,
        params,
        discretization,
        action,
        la,
        elemat1_epetra,
        elemat2_epetra,
        elevec1_epetra,
        elevec2_epetra,
        elevec3_epetra
        );
    break;
  }
  case SCATRA::bd_add_convective_mass_flux:
  {
    //calculate integral of convective mass/heat flux
    // NOTE: since results are added to a global vector via normal assembly
    //       it would be wrong to suppress results for a ghosted boundary!

    // get actual values of transported scalars
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> ephinp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,la[0].lm_);

    // get velocity values at nodes
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field",Teuchos::null);

    // we deal with a (nsd_+1)-dimensional flow field
    LINALG::Matrix<my::nsd_+1,my::nen_>  evel(true);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,my::nsd_+1);

    {
      const Teuchos::RCP<Epetra_MultiVector> pre = params.get< Teuchos::RCP<Epetra_MultiVector> >("pressure field");
      LINALG::Matrix<1,my::nen_> eprenp;
      DRT::UTILS::ExtractMyNodeBasedValues(ele,eprenp,pre,1);

      //pressure values
      for (int i=0;i<my::nen_;++i)
      {
        eprenp_(i) = eprenp(0,i);
      }

      // this is a hack. Check if the structure (assumed to be the dofset 1) has more DOFs than dimension. If so,
      // we assume that this is the porosity
      if( discretization.NumDof(1,ele->Nodes()[0])==my::nsd_+2 )
      {
        isnodalporosity_=true;

        Teuchos::RCP<const Epetra_Vector> disp= discretization.GetState(1,"displacement");

        if(disp!=Teuchos::null)
        {
          std::vector<double> mydisp(la[1].lm_.size());
          DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

          for (int inode=0; inode<my::nen_; ++inode)  // number of nodes
            eporosity_(inode,0) = mydisp[my::nsd_+1+(inode*(my::nsd_+2))];
        }
        else
          dserror("Cannot get state vector displacement");
      }
      else
        isnodalporosity_=false;
    }

    // for the moment we ignore the return values of this method
    CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
    //vector<double> locfluxintegral = CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
    //std::cout<<"locfluxintegral[0] = "<<locfluxintegral[0]<<std::endl;

    break;
  }
  default:
  {
    dserror("Invalid action parameter nr. %i!",action);
    break;
  }
  } // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
 | calculate integral of convective flux across boundary    vuong 07/15 |
 | (overwrites method in ScaTraEleBoundaryCalc)                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double> DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::CalcConvectiveFlux(
    const DRT::FaceElement*                     ele,
    const std::vector<double>&                  ephinp,
    const LINALG::Matrix<my::nsd_+1,my::nen_>&  evelnp,
    Epetra_SerialDenseVector&                   erhs
)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<my::nen_,1> phinod(true);

  std::vector<double> integralflux(my::numscal_);

  // loop over all scalars
  for(int k=0;k<my::numscal_;++k)
  {
    integralflux[k] = 0.0;

    // compute scalar values at nodes
    for (int inode=0; inode< my::nen_;++inode)
    {
      phinod(inode) = ephinp[inode*my::numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {

      const double fac = my::EvalShapeFuncAndIntFac(intpoints,iquad,&(this->normal_));

      const double porosity = ComputePorosity(ele);

      // get velocity at integration point
      my::velint_.Multiply(evelnp,my::funct_);

      // normal velocity (note: normal_ is already a unit(!) normal)
      const double normvel = my::velint_.Dot(my::normal_);

      // scalar at integration point
      const double phi = my::funct_.Dot(phinod);

      const double val = porosity*phi*normvel*fac;
      integralflux[k] += val;
      // add contribution to provided vector (distribute over nodes using shape fct.)
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = vi*my::numdofpernode_+k;
        erhs[fvi] += val*my::funct_(vi);
      }
    }
  }

  return integralflux;

} //ScaTraEleBoundaryCalcPoro<distype>::ConvectiveFlux


/*----------------------------------------------------------------------*
 |  get the material constants  (protected)                  vuong 10/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::ComputePorosity(
    const DRT::FaceElement* ele       //!< the element we are dealing with
  )
{
  double porosity=0.0;

  if(isnodalporosity_)
  {
    porosity = eporosity_.Dot(my::funct_);
  }
  else
  {
    dserror("porosity calculation not yet implemented for non-node-based porosity!");
  }

  return porosity;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::nurbs9>;
