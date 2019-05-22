/*----------------------------------------------------------------------*/
/*!

\brief evaluation of ScaTra boundary elements for Nernst-Planck formulation

\level 2

\maintainer Christoph Schmidt
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleBoundaryCalcElchNP* delete_me)
{
  static std::map<std::string, ScaTraEleBoundaryCalcElchNP<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleBoundaryCalcElchNP<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleBoundaryCalcElchNP<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
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
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::ScaTraEleBoundaryCalcElchNP(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      myelch::ScaTraEleBoundaryCalcElch(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateAction(
    DRT::FaceElement* ele,                     //!< boundary element
    Teuchos::ParameterList& params,            //!< parameter list
    DRT::Discretization& discretization,       //!< discretization
    SCATRA::BoundaryAction action,             //!< action
    DRT::Element::LocationArray& la,           //!< location array
    Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::bd_calc_elch_boundary_kinetics:
    {
      myelch::CalcElchBoundaryKinetics(
          ele, params, discretization, la, elemat1_epetra, elevec1_epetra, 1.);

      break;
    }

    default:
    {
      myelch::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch action

  return 0;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateAction


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                       fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateNeumann(DRT::FaceElement* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, DRT::Condition& condition,
    DRT::Element::LocationArray& la, Epetra_SerialDenseVector& elevec1, const double scalar)
{
  // call base class routine
  my::EvaluateNeumann(ele, params, discretization, condition, la, elevec1, scalar);

  // add boundary flux contributions to potential equation
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        // get valence
        const double valence_k = GetValence(ele->ParentElement()->Material(), k);

        for (int vi = 0; vi < my::nen_; ++vi)
          elevec1[vi * my::numdofpernode_ + my::numscal_] +=
              valence_k * elevec1[vi * my::numdofpernode_ + k];
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
  }  // switch(myelch::elchparams_->EquPot())

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition         fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateElchBoundaryKinetics(
    const DRT::Element* ele,         ///< current element
    Epetra_SerialDenseMatrix& emat,  ///< element matrix
    Epetra_SerialDenseVector& erhs,  ///< element right-hand side vector
    const std::vector<LINALG::Matrix<my::nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<LINALG::Matrix<my::nen_, 1>>& ehist,  ///< nodal history vector
    double timefac,                                         ///< time factor
    Teuchos::RCP<const MAT::Material> material,             ///< material
    Teuchos::RCP<DRT::Condition> cond,  ///< electrode kinetics boundary condition
    const int nume,                     ///< number of transferred electrons
    const std::vector<int> stoich,      ///< stoichiometry of the reaction
    const int kinetics,                 ///< desired electrode kinetics model
    const double pot0,                  ///< electrode potential on metal side
    const double frt,                   ///< factor F/RT
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // call base class routine
  myelch::EvaluateElchBoundaryKinetics(ele, emat, erhs, ephinp, ehist, timefac, material, cond,
      nume, stoich, kinetics, pot0, frt, scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (int vi = 0; vi < my::nen_; ++vi)
        {
          for (int ui = 0; ui < my::nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] += nume * erhs[vi * my::numdofpernode_ + k];
        }
      }

      break;
    }

    // need special treatment for Laplace equation due to missing scaling with inverse of Faraday
    // constant
    case INPAR::ELCH::equpot_laplace:
    {
      const double faraday = myelch::elchparams_->Faraday();
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (int vi = 0; vi < my::nen_; ++vi)
        {
          for (int ui = 0; ui < my::nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                faraday * nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                faraday * nume *
                emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] +=
              faraday * nume * erhs[vi * my::numdofpernode_ + k];
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
  }  // switch(myelch::elchparams_->EquPot())

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::EvaluateElchBoundaryKinetics


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>& material,  // element material
    const int k                                         // species number
    ) const
{
  double valence(0.);

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList> matlist =
        Teuchos::rcp_static_cast<const MAT::MatList>(material);

    const Teuchos::RCP<const MAT::Material> species = matlist->MaterialById(matlist->MatID(k));

    if (species->MaterialType() == INPAR::MAT::m_ion)
    {
      valence = Teuchos::rcp_static_cast<const MAT::Ion>(species)->Valence();
      if (abs(valence) < 1.e-14) dserror("Received zero valence!");
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
