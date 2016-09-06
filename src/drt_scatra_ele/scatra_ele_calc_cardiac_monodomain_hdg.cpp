/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_cardiac_monodomain_hdg.cpp

 \brief main file containing routines for calculation of HDG cardiac monodomain element

 <pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "scatra_ele_calc_cardiac_monodomain_hdg.H"

#include "scatra_ele_calc_hdg.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "scatra_ele_parameter_timint.H"

#include "../drt_mat/myocard.H"
#include "../drt_mat/matlist.H"

#include "../drt_fem_general/drt_utils_polynomial.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::ScaTraEleCalcHDGCardiacMonodomain(
    const int numdofpernode,
    const int numscal,
    const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::ScaTraEleCalcHDG(numdofpernode,numscal,disname)
{

}

/*----------------------------------------------------------------------*
 | singleton access method                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim> * DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  bool create )
{
  static std::map<std::string,ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>(numdofpernode,numscal,disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( typename std::map<std::string,ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
 | singleton destruction                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", false );
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::Materials(
  const Teuchos::RCP<const MAT::Material>   material, //!< pointer to current material
  const int                                 k,        //!< id of current scalar
  Epetra_SerialDenseMatrix*                 difftensor,
  Epetra_SerialDenseVector*                 ivecn,
  Epetra_SerialDenseVector*                 ivecnp,
  Epetra_SerialDenseMatrix*                 ivecnpderiv,
  double*                                   diff1          //!< main diffusivity
)
{
  if (material->MaterialType() == INPAR::MAT::m_myocard)
    MatMyocard(material,k,difftensor,ivecn,ivecnp,ivecnpderiv,diff1);
  else dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                      hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::MatMyocard(
  const Teuchos::RCP<const MAT::Material>   material,   //!< pointer to current material
  const int                                 k,          //!< id of current scalar
  Epetra_SerialDenseMatrix*                 difftensor,
  Epetra_SerialDenseVector*                 ivecn,
  Epetra_SerialDenseVector*                 ivecnp,
  Epetra_SerialDenseMatrix*                 ivecnpderiv,
  double*                                   diff1          //!< main diffusivity
  )
{
  const Teuchos::RCP<const MAT::Myocard>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Myocard>(material);

  // get constant diffusivity
  LINALG::Matrix<probdim,probdim> diff(true);
  actmat->Diffusivity(diff);

  *diff1 = actmat->Parameter()->diff1;


  for (unsigned int i=0; i<this->nsd_; ++i)
    for (unsigned int j=0; j<this->nsd_; ++j)
      (*difftensor)(i,j) = diff(i,j);

  // coordinate of material gauss points
  LINALG::Matrix<probdim,1> mat_gp_coord(true);
  // values of shape function at material gauss points
  Epetra_SerialDenseVector values_mat_gp(this->shapes_->ndofs_);

  double imatgpnpderiv(0.);
  double imatgpnp(0.);
  double imatgpn(0.);

  // polynomial space to get the value of the shape function at the material gauss points
  DRT::UTILS::PolynomialSpaceParams params(distype,this->shapes_->degree_,this->usescompletepoly_);
  polySpace_ = DRT::UTILS::PolynomialSpaceCache<probdim>::Instance().Create(params);

  Epetra_SerialDenseMatrix ivecnpderiv_gp(this->shapes_->ndofs_,this->shapes_->ndofs_);
  Epetra_SerialDenseVector ivecnp_gp(this->shapes_->ndofs_);
  Epetra_SerialDenseVector ivecn_gp(this->shapes_->ndofs_);

  const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim> intpoints(SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(actmat->Parameter()->num_gp));

    for (int q = 0; q < intpoints.IP().nquad; q++)
    {
      double phinpgp = 0.0;
      double phingp = 0.0;

      // gaussian points coordinates
      for (int idim=0;idim<DRT::UTILS::DisTypeToDim<distype>::dim;idim++)
        mat_gp_coord(idim) = intpoints.IP().qxg[q][idim];

      //gaussian weight
      double gp_mat_alpha = intpoints.IP().qwgt[q];

      polySpace_->Evaluate(mat_gp_coord,values_mat_gp);

      // loop over shape functions
      for (unsigned int i=0; i<this->shapes_->ndofs_; ++i)
      {
        phingp += values_mat_gp(i)*this->interiorPhin_(i);
        phinpgp += values_mat_gp(i)*this->interiorPhinp_(i);
      }

      // Reaction term at material gauss points
      imatgpnpderiv = actmat->ReaCoeffDeriv(phinpgp,this->Dt(),q);
      imatgpn = actmat->ReaCoeff(phingp, this->Dt(),q);
      imatgpnp = actmat->ReaCoeff(phinpgp, this->Dt(),q);

      // loop over shape functions
      for (unsigned int i=0; i<this->shapes_->ndofs_; i++)
      {
        for (unsigned int j=0; j<this->shapes_->ndofs_; j++)
          ivecnpderiv_gp(i,j) += imatgpnpderiv*values_mat_gp(i)*values_mat_gp(j)*this->shapes_->xji.Invert(this->shapes_->xjm) *gp_mat_alpha;
        ivecnp_gp(i) += imatgpnp*values_mat_gp(i)*this->shapes_->xji.Invert(this->shapes_->xjm) *gp_mat_alpha;
        ivecn_gp(i) += imatgpn*values_mat_gp(i)*this->shapes_->xji.Invert(this->shapes_->xjm) *gp_mat_alpha;
      }
    }

  *ivecnpderiv = ivecnpderiv_gp;
  *ivecnp = ivecnp_gp;
  *ivecn = ivecn_gp;



  return;
} // ScaTraEleCalcHDGCardiacMonodomain<distype>::MatMyocard


/*----------------------------------------------------------------------*
 |  Material Time Update                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::TimeUpdateMaterial(
    const DRT::Element*   ele   //!< the element we are dealing with
    )
{

  std::vector<Teuchos::RCP<MAT::Myocard> > updatemat;

  // access the general material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // first, determine the materials which need a time update, i.e. myocard materials
  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<MAT::MatList> actmat
    = Teuchos::rcp_dynamic_cast<MAT::MatList>(material);
    if (actmat->NumMat() < this->numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<this->numscal_;++k)
    {
      const int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == INPAR::MAT::m_myocard)
      {
        // reference to Teuchos::rcp not possible here, since the material
        // is required to be not const for this application
        updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(singlemat));
      }
    }
  }

  if (material->MaterialType() == INPAR::MAT::m_myocard)
  {
    // reference to Teuchos::rcp not possible here, since the material is required to be
    // not const for this application
    updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(material));
  }

  if (updatemat.size()>0) // found at least one material to be updated
  {

    for (unsigned i=0;i<updatemat.size();i++)
      updatemat[i]->Update(Teuchos::null, 0.0);

  }

  return;

}


/*----------------------------------------------------------------------*
 |  Get Material Internal State for Restart              hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::GetMaterialInternalState(
    const DRT::Element*      ele,   //!< the element we are dealing with
    Teuchos::ParameterList&           params,
    DRT::Discretization &    discretization
    )
{
  // NOTE: add integral values only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    // access the general material
    Teuchos::RCP<MAT::Material> material = ele->Material();
    Teuchos::RCP<Epetra_MultiVector> material_internal_state = params.get< Teuchos::RCP<Epetra_MultiVector> >("material_internal_state");

    if( material->MaterialType() == INPAR::MAT::m_myocard)
    {
      Teuchos::RCP<MAT::Myocard> material = Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
      for (int k = 0; k< material->GetNumberOfInternalStateVariables(); ++k)
      {
        double material_state = 0;
        unsigned int nqpoints = this->shapes_->nqpoints_;
        for (unsigned int q = 0; q < nqpoints; ++q)
        {
          material_state += material->GetInternalState(k,q);

        }
        int err = material_internal_state->ReplaceGlobalValue(ele->Id(),k, material_state/nqpoints);
        if(err != 0) dserror("%i",err);
      }
    }

    params.set< Teuchos::RCP<Epetra_MultiVector> >("material_internal_state", material_internal_state);

  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set Material Internal State after Restart            hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::SetMaterialInternalState(
    const DRT::Element*               ele,   //!< the element we are dealing with
    Teuchos::ParameterList&           params,
    DRT::Discretization &             discretization
    )
{
  // NOTE: add integral values only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
//  if (discretization.Comm().MyPID()==0)
  {
    // access the general material
    Teuchos::RCP<MAT::Material> material = ele->Material();
    Teuchos::RCP<Epetra_MultiVector> material_internal_state = params.get< Teuchos::RCP<Epetra_MultiVector> >("material_internal_state");

    if( material->MaterialType() == INPAR::MAT::m_myocard)
    {
      Teuchos::RCP<MAT::Myocard> material = Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
      for (int k = 0; k< material->GetNumberOfInternalStateVariables(); ++k)
      {
        int nqpoints = this->shapes_->nqpoints_;
        for (int q = 0; q < nqpoints; ++q)
        {
          Teuchos::RCP<Epetra_Vector> material_internal_state_component = Teuchos::rcp((*material_internal_state)(k*nqpoints+q),false);
          material->SetInternalState(k,(*material_internal_state_component)[ele->Id()],q);
        }
      }
    }
  }


  return;
}


// template classes
// 1D elements
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line2,1>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line2,2>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line2,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad4,2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad4,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::quad9,2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::nurbs9,2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::hex8,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::hex27,3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tet4,3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::tet10,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::pyramid5,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<DRT::Element::nurbs27>;
