/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_NP.cpp

\brief evaluation of ScaTra boundary elements for Nernst-Planck formulation

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "../drt_mat/ion.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"

#include "scatra_ele.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_boundary_calc_elch_NP.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    bool create
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchNP<distype>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchNP<distype>(numdofpernode,numscal,disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcElchNP<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",false);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::ScaTraEleBoundaryCalcElchNP(const int numdofpernode,const int numscal, const std::string& disname)
  : // constructor of base class
    myelch::ScaTraEleBoundaryCalcElch(numdofpernode,numscal,disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                       fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateNeumann(
    DRT::FaceElement*                   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    DRT::Condition&                     condition,
    DRT::Element::LocationArray&        la,
    Epetra_SerialDenseVector&           elevec1,
    const double                        scalar
    )
{
  // call base class routine
  my::EvaluateNeumann(ele,params,discretization,condition,la,elevec1,scalar);

  // add boundary flux contributions to potential equation
  switch(myelch::ElchParams()->EquPot())
  {
  case INPAR::ELCH::equpot_enc_pde:
  case INPAR::ELCH::equpot_enc_pde_elim:
  {
    for(int k=0; k<my::numscal_; ++k)
    {
      // get valence
      const double valence_k = GetValence(ele->ParentElement()->Material(),k);

      for (int vi=0; vi<my::nen_; ++vi)
        elevec1[vi*my::numdofpernode_+my::numscal_] += valence_k*elevec1[vi*my::numdofpernode_+k];
    }

    break;
  }

  case INPAR::ELCH::equpot_enc:
  case INPAR::ELCH::equpot_poisson:
  case INPAR::ELCH::equpot_laplace:
    // do nothing in these cases
    break;

  default:
  {
    dserror("Closing equation for electric potential not recognized!");
    break;
  }
  } // switch(myelch::ElchParams()->EquPot())

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate electrode kinetics boundary condition            fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateElchBoundaryKinetics(
    const DRT::Element*                 ele,        ///< the actual boundary element
    Epetra_SerialDenseMatrix&           emat,       ///< element-matrix
    Epetra_SerialDenseVector&           erhs,       ///< element-rhs
    const std::vector<double>&          ephinp,     ///< actual conc. and pot. values
    const std::vector<double>&          ehist,      ///< element history vector
    double                              timefac,    ///< time factor
    Teuchos::RCP<const MAT::Material>   material,   ///< the material
    Teuchos::RCP<DRT::Condition>        cond,       ///< the condition
    const int                           nume,       ///< number of transferred electrons
    const std::vector<int>              stoich,     ///< stoichiometry of the reaction
    const int                           kinetics,   ///< desired electrode kinetics model
    const double                        pot0,       ///< actual electrode potential on metal side
    const double                        frt,        ///< factor F/RT
    const double                        scalar      ///< scaling factor for element matrix and residual contributions
)
{
  // call base class routine
  myelch::EvaluateElchBoundaryKinetics(ele,emat,erhs,ephinp,ehist,timefac,material,cond,nume,stoich,kinetics,pot0,frt,scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch(myelch::ElchParams()->EquPot())
  {
  case INPAR::ELCH::equpot_enc:
  {
    // do nothing, since no boundary integral present
    break;
  }

  case INPAR::ELCH::equpot_enc_pde:
  case INPAR::ELCH::equpot_enc_pde_elim:
  {
    for(int k=0; k<my::numscal_; ++k)
    {
      for(int vi=0; vi<my::nen_; ++vi)
      {
        for(int ui=0; ui<my::nen_; ++ui)
        {
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k) += nume*emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += nume*emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
        }

        erhs[vi*my::numdofpernode_+my::numscal_] += nume*erhs[vi*my::numdofpernode_+k];
      }
    }

    break;
  }

  // need special treatment for Laplace equation due to missing scaling with inverse of Faraday constant
  case INPAR::ELCH::equpot_laplace:
  {
    for(int k=0; k<my::numscal_; ++k)
    {
      for(int vi=0; vi<my::nen_; ++vi)
      {
        for(int ui=0; ui<my::nen_; ++ui)
        {
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k) += INPAR::ELCH::faraday_const*nume*emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += INPAR::ELCH::faraday_const*nume*emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
        }

        erhs[vi*my::numdofpernode_+my::numscal_] += INPAR::ELCH::faraday_const*nume*erhs[vi*my::numdofpernode_+k];
      }
    }

    break;
  }

  case INPAR::ELCH::equpot_poisson:
  {
    dserror("Poisson equation combined with electrode boundary conditions not implemented!");
    break;
  }

  default:
  {
    dserror("Unknown closing equation for electric potential!");
    break;
  }
  } // switch(myelch::ElchParams()->EquPot())

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateElchBoundaryKinetics


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>&   material,   // element material
    const int                                  k           // species number
    ) const
{
  double valence(0.);

  if(material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList> matlist = Teuchos::rcp_static_cast<const MAT::MatList>(material);

    const Teuchos::RCP<const MAT::Material> species = matlist->MaterialById(matlist->MatID(k));

    if(species->MaterialType() == INPAR::MAT::m_ion)
    {
      valence = Teuchos::rcp_static_cast<const MAT::Ion>(species)->Valence();
      if(abs(valence) < 1.e-14)
        dserror("Received zero valence!");
    }
    else
      dserror("Material species is not an ion!");
  }

  else
    dserror("Unknown material!");

  return valence;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<DRT::Element::nurbs9>;
