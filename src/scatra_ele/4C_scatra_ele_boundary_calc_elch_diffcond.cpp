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
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>>(
            new ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype,
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
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele,                            //!< boundary element
    Teuchos::ParameterList& params,                   //!< parameter list
    DRT::Discretization& discretization,              //!< discretization
    SCATRA::BoundaryAction action,                    //!< action
    DRT::Element::LocationArray& la,                  //!< location array
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_elch_boundary_kinetics:
    {
      // access material of parent element
      Teuchos::RCP<CORE::MAT::Material> material = ele->ParentElement()->Material();

      // extract porosity from material and store in diffusion manager
      if (material->MaterialType() == CORE::Materials::m_elchmat)
      {
        const auto* elchmat = static_cast<const MAT::ElchMat*>(material.get());

        for (int iphase = 0; iphase < elchmat->NumPhase(); ++iphase)
        {
          Teuchos::RCP<const CORE::MAT::Material> phase =
              elchmat->PhaseById(elchmat->PhaseID(iphase));

          if (phase->MaterialType() == CORE::Materials::m_elchphase)
          {
            dmedc_->SetPhasePoro(
                (static_cast<const MAT::ElchPhase*>(phase.get()))->Epsilon(), iphase);
          }
          else
            FOUR_C_THROW("Invalid material!");
        }
      }

      else
        FOUR_C_THROW("Invalid material!");

      // process electrode kinetics boundary condition
      myelch::CalcElchBoundaryKinetics(
          ele, params, discretization, la, elemat1_epetra, elevec1_epetra, dmedc_->GetPhasePoro(0));

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
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                       fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::EvaluateNeumann(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Condition& condition, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& elevec1, const double scalar)
{
  // get material of parent element
  Teuchos::RCP<CORE::MAT::Material> mat = ele->ParentElement()->Material();

  if (mat->MaterialType() == CORE::Materials::m_elchmat)
  {
    const auto* actmat = static_cast<const MAT::ElchMat*>(mat.get());

    for (int iphase = 0; iphase < actmat->NumPhase(); ++iphase)
    {
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const CORE::MAT::Material> singlemat = actmat->PhaseById(phaseid);

      if (singlemat->MaterialType() == CORE::Materials::m_elchphase)
      {
        const auto* actsinglemat = static_cast<const MAT::ElchPhase*>(singlemat.get());

        dmedc_->SetPhasePoro(actsinglemat->Epsilon(), iphase);
      }
      else
        FOUR_C_THROW("Invalid material!");
    }
  }
  else
    FOUR_C_THROW("Invalid material!");

  // call base class routine
  my::EvaluateNeumann(ele, params, discretization, condition, la, elevec1, dmedc_->GetPhasePoro(0));

  // add boundary flux contributions to potential equation
  if (myelch::elchparams_->BoundaryFluxCoupling())
  {
    switch (myelch::elchparams_->EquPot())
    {
      case INPAR::ELCH::equpot_divi:
      {
        for (int k = 0; k < my::numscal_; ++k)
        {
          // get valence
          const double valence_k = GetValence(mat, k);

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
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype,
    probdim>::EvaluateElchBoundaryKinetics(const DRT::Element* ele,  ///< current element
    CORE::LINALG::SerialDenseMatrix& emat,                           ///< element matrix
    CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
    double timefac,                                           ///< time factor
    Teuchos::RCP<const CORE::MAT::Material> material,         ///< material
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

    case INPAR::ELCH::equpot_divi:
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
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::EvaluateS2ICoupling(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& eslavematrix, CORE::LINALG::SerialDenseMatrix& emastermatrix,
    CORE::LINALG::SerialDenseVector& eslaveresidual)
{
  switch (my::scatraparamsboundary_->KineticModel())
  {
    case INPAR::S2I::kinetics_nointerfaceflux:
      break;
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    {
      myelectrode::EvaluateS2ICoupling(
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
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::EvaluateS2ICouplingOD(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& eslavematrix)
{
  switch (my::scatraparamsboundary_->KineticModel())
  {
    case INPAR::S2I::kinetics_nointerfaceflux:
      break;
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    {
      myelectrode::EvaluateS2ICouplingOD(ele, params, discretization, la, eslavematrix);
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
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::GetValence(
    const Teuchos::RCP<const CORE::MAT::Material>& material,  // element material
    const int k                                               // species number
) const
{
  double valence(0.);

  if (material->MaterialType() == CORE::Materials::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat> elchmat =
        Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

    // safety check
    if (elchmat->NumPhase() != 1) FOUR_C_THROW("Only one material phase is allowed at the moment!");

    // loop over phases
    for (int iphase = 0; iphase < elchmat->NumPhase(); ++iphase)
    {
      const Teuchos::RCP<const MAT::ElchPhase> phase =
          Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(
              elchmat->PhaseById(elchmat->PhaseID(iphase)));

      // loop over species within phase
      for (int imat = 0; imat < phase->NumMat(); ++imat)
      {
        const Teuchos::RCP<const CORE::MAT::Material> species = phase->MatById(phase->MatID(imat));

        if (species->MaterialType() == CORE::Materials::m_newman)
        {
          valence = Teuchos::rcp_static_cast<const MAT::Newman>(species)->Valence();
          if (abs(valence) < 1.e-14) FOUR_C_THROW("Received zero valence!");
        }
        else if (species->MaterialType() == CORE::Materials::m_ion)
        {
          valence = Teuchos::rcp_static_cast<const MAT::Ion>(species)->Valence();
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
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<CORE::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
