/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation

\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_boundary_calc_elch_diffcond.hpp"

#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_ion.hpp"
#include "4C_mat_newman.hpp"
#include "4C_scatra_ele_calc_elch_diffcond.hpp"
#include "4C_scatra_ele_parameter_boundary.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>*
Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>>(
            new ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype,
    probdim>::ScaTraEleBoundaryCalcElchDiffCond(const int numdofpernode, const int numscal,
    const std::string& disname)
    :  // constructor of base class
      myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode, numscal, disname),
      // initialization of diffusion manager
      dmedc_(Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCond(my::numscal_)))
{
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 08/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::evaluate_action(
    Core::Elements::FaceElement* ele,                 //!< boundary element
    Teuchos::ParameterList& params,                   //!< parameter list
    Discret::Discretization& discretization,          //!< discretization
    ScaTra::BoundaryAction action,                    //!< action
    Core::Elements::Element::LocationArray& la,       //!< location array
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    Core::LinAlg::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    Core::LinAlg::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    Core::LinAlg::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case ScaTra::BoundaryAction::calc_elch_boundary_kinetics:
    {
      // access material of parent element
      Teuchos::RCP<Core::Mat::Material> material = ele->parent_element()->Material();

      // extract porosity from material and store in diffusion manager
      if (material->MaterialType() == Core::Materials::m_elchmat)
      {
        const auto* elchmat = static_cast<const Mat::ElchMat*>(material.get());

        for (int iphase = 0; iphase < elchmat->NumPhase(); ++iphase)
        {
          Teuchos::RCP<const Core::Mat::Material> phase =
              elchmat->PhaseById(elchmat->PhaseID(iphase));

          if (phase->MaterialType() == Core::Materials::m_elchphase)
          {
            dmedc_->SetPhasePoro(
                (static_cast<const Mat::ElchPhase*>(phase.get()))->Epsilon(), iphase);
          }
          else
            FOUR_C_THROW("Invalid material!");
        }
      }

      else
        FOUR_C_THROW("Invalid material!");

      // process electrode kinetics boundary condition
      myelch::calc_elch_boundary_kinetics(
          ele, params, discretization, la, elemat1_epetra, elevec1_epetra, dmedc_->GetPhasePoro(0));

      break;
    }

    default:
    {
      myelch::evaluate_action(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch action

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                       fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::evaluate_neumann(
    Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
    const double scalar)
{
  // get material of parent element
  Teuchos::RCP<Core::Mat::Material> mat = ele->parent_element()->Material();

  if (mat->MaterialType() == Core::Materials::m_elchmat)
  {
    const auto* actmat = static_cast<const Mat::ElchMat*>(mat.get());

    for (int iphase = 0; iphase < actmat->NumPhase(); ++iphase)
    {
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const Core::Mat::Material> singlemat = actmat->PhaseById(phaseid);

      if (singlemat->MaterialType() == Core::Materials::m_elchphase)
      {
        const auto* actsinglemat = static_cast<const Mat::ElchPhase*>(singlemat.get());

        dmedc_->SetPhasePoro(actsinglemat->Epsilon(), iphase);
      }
      else
        FOUR_C_THROW("Invalid material!");
    }
  }
  else
    FOUR_C_THROW("Invalid material!");

  // call base class routine
  my::evaluate_neumann(
      ele, params, discretization, condition, la, elevec1, dmedc_->GetPhasePoro(0));

  // add boundary flux contributions to potential equation
  if (myelch::elchparams_->boundary_flux_coupling())
  {
    switch (myelch::elchparams_->EquPot())
    {
      case Inpar::ElCh::equpot_divi:
      {
        for (int k = 0; k < my::numscal_; ++k)
        {
          // get valence
          const double valence_k = get_valence(mat, k);

          for (int vi = 0; vi < nen_; ++vi)
            elevec1[vi * my::numdofpernode_ + my::numscal_] +=
                valence_k * elevec1[vi * my::numdofpernode_ + k];
        }  // loop over scalars

        break;
      }

      default:
      {
        FOUR_C_THROW("Closing equation for electric potential not recognized!");
        break;
      }
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition         fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype,
    probdim>::evaluate_elch_boundary_kinetics(const Core::Elements::Element*
                                                  ele,  ///< current element
    Core::LinAlg::SerialDenseMatrix& emat,              ///< element matrix
    Core::LinAlg::SerialDenseVector& erhs,              ///< element right-hand side vector
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
    double timefac,                                           ///< time factor
    Teuchos::RCP<const Core::Mat::Material> material,         ///< material
    Teuchos::RCP<Core::Conditions::Condition> cond,  ///< electrode kinetics boundary condition
    const int nume,                                  ///< number of transferred electrons
    const std::vector<int> stoich,                   ///< stoichiometry of the reaction
    const int kinetics,                              ///< desired electrode kinetics model
    const double pot0,                               ///< electrode potential on metal side
    const double frt,                                ///< factor F/RT
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // call base class routine
  myelch::evaluate_elch_boundary_kinetics(ele, emat, erhs, ephinp, ehist, timefac, material, cond,
      nume, stoich, kinetics, pot0, frt, scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch (myelch::elchparams_->EquPot())
  {
    case Inpar::ElCh::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case Inpar::ElCh::equpot_divi:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (int vi = 0; vi < nen_; ++vi)
        {
          for (int ui = 0; ui < nen_; ++ui)
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

    default:
    {
      FOUR_C_THROW("Unknown closing equation for electric potential!");
      break;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition (electrochemistry)   fang 12/14 |
 *-------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::evaluate_s2_i_coupling(
    const Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& eslavematrix, Core::LinAlg::SerialDenseMatrix& emastermatrix,
    Core::LinAlg::SerialDenseVector& eslaveresidual)
{
  switch (my::scatraparamsboundary_->KineticModel())
  {
    case Inpar::S2I::kinetics_nointerfaceflux:
      break;
    case Inpar::S2I::kinetics_constantinterfaceresistance:
    {
      myelectrode::evaluate_s2_i_coupling(
          ele, params, discretization, la, eslavematrix, emastermatrix, eslaveresidual);
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Evaluation of scatra-scatra interface kinetics for electrochemistry problems with "
          "conforming interface discretization must have an electrode on the slave side and the "
          "electrolyte on the master side.");
      break;
    }
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype,
    probdim>::evaluate_s2_i_coupling_od(const Core::Elements::FaceElement* ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization,
    Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& eslavematrix)
{
  switch (my::scatraparamsboundary_->KineticModel())
  {
    case Inpar::S2I::kinetics_nointerfaceflux:
      break;
    case Inpar::S2I::kinetics_constantinterfaceresistance:
    {
      myelectrode::evaluate_s2_i_coupling_od(ele, params, discretization, la, eslavematrix);
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Evaluation of scatra-scatra interface kinetics for electrochemistry problems with "
          "conforming interface discretization must have an electrode on the slave side and the "
          "electrolyte on the master side.");
      break;
    }
  }
}


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 12/14 |
 *-------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
double Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::get_valence(
    const Teuchos::RCP<const Core::Mat::Material>& material,  // element material
    const int k                                               // species number
) const
{
  double valence(0.);

  if (material->MaterialType() == Core::Materials::m_elchmat)
  {
    const Teuchos::RCP<const Mat::ElchMat> elchmat =
        Teuchos::rcp_dynamic_cast<const Mat::ElchMat>(material);

    // safety check
    if (elchmat->NumPhase() != 1) FOUR_C_THROW("Only one material phase is allowed at the moment!");

    // loop over phases
    for (int iphase = 0; iphase < elchmat->NumPhase(); ++iphase)
    {
      const Teuchos::RCP<const Mat::ElchPhase> phase =
          Teuchos::rcp_dynamic_cast<const Mat::ElchPhase>(
              elchmat->PhaseById(elchmat->PhaseID(iphase)));

      // loop over species within phase
      for (int imat = 0; imat < phase->NumMat(); ++imat)
      {
        const Teuchos::RCP<const Core::Mat::Material> species = phase->MatById(phase->MatID(imat));

        if (species->MaterialType() == Core::Materials::m_newman)
        {
          valence = Teuchos::rcp_static_cast<const Mat::Newman>(species)->Valence();
          if (abs(valence) < 1.e-14) FOUR_C_THROW("Received zero valence!");
        }
        else if (species->MaterialType() == Core::Materials::m_ion)
        {
          valence = Teuchos::rcp_static_cast<const Mat::Ion>(species)->Valence();
          if (abs(valence) < 1.e-14) FOUR_C_THROW("Received zero valence!");
        }
        else
          FOUR_C_THROW("Unknown material species!");
      }
    }
  }

  else
    FOUR_C_THROW("Unknown material!");

  return valence;
}


// template classes
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::quad4, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::quad8, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::quad9, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::tri6, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::line3, 2>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::nurbs3, 2>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<Core::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
