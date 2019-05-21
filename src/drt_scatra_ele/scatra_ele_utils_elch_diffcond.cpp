/*----------------------------------------------------------------------*/
/*!

\brief utility class supporting element evaluation for concentrated electrolytes

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_utils_elch_diffcond.H"
#include "scatra_ele_calc_elch_diffcond.H"
#include "scatra_ele_calc_elch_diffcond_multiscale.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/newman_multiscale.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>*
DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::Instance(
    const int numdofpernode,                     ///< number of degrees of freedom per node
    const int numscal,                           ///< number of transported scalars per node
    const std::string& disname,                  ///< name of discretization
    const ScaTraEleUtilsElchDiffCond* delete_me  ///< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleUtilsElchDiffCond<distype>*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleUtilsElchDiffCond<distype>(numdofpernode, numscal, disname);
  }

  // destruct instance
  else
  {
    for (typename std::map<std::string, ScaTraEleUtilsElchDiffCond<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  // return existing or newly created instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::ScaTraEleUtilsElchDiffCond(
    const int numdofpernode,    ///< number of degrees of freedom per node
    const int numscal,          ///< number of transported scalars per node
    const std::string& disname  ///< name of discretization
    )
    : myelectrode::ScaTraEleUtilsElchElectrode(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate electrolyte material                             fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatElchMat(
    const Teuchos::RCP<const MAT::Material>& material,  //!< electrolyte material
    const std::vector<double>& concentrations,          //!< local concentration values
    const INPAR::ELCH::EquPot& equpot,  //!< type of closing equation for electric potential
    const double& ffrt,                 //!< factor F²/RT
    const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& diffmanager,  //!< diffusion manager
    INPAR::ELCH::DiffCondMat& diffcondmat                               //!< ion type
)
{
  // cast material to electrolyte material
  const Teuchos::RCP<const MAT::ElchMat> elchmat =
      Teuchos::rcp_static_cast<const MAT::ElchMat>(material);

  // safety check
  if (elchmat->NumPhase() != 1) dserror("Can only have a single electrolyte phase at the moment!");

  // extract electrolyte phase
  const Teuchos::RCP<const MAT::Material> elchphase = elchmat->PhaseById(elchmat->PhaseID(0));

  if (elchphase->MaterialType() == INPAR::MAT::m_elchphase)
    // evaluate electrolyte phase
    MatElchPhase(elchphase, concentrations, equpot, ffrt, diffmanager, diffcondmat);
  else
    dserror("Invalid material type!");

  return;
}

/*----------------------------------------------------------------------*
 | evaluate electrolyte phase                                fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatElchPhase(
    const Teuchos::RCP<const MAT::Material>& material,  //!< electrolyte phase
    const std::vector<double>& concentrations,          //!< local concentration values
    const INPAR::ELCH::EquPot& equpot,  //!< type of closing equation for electric potential
    const double& ffrt,                 //!< factor F²/RT
    const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& diffmanager,  //!< diffusion manager
    INPAR::ELCH::DiffCondMat& diffcondmat                               //!< ion type
)
{
  // cast material to electrolyte phase
  const Teuchos::RCP<const MAT::ElchPhase> matelchphase =
      Teuchos::rcp_static_cast<const MAT::ElchPhase>(material);

  // set porosity
  diffmanager->SetPhasePoro(matelchphase->Epsilon(), 0);

  // set tortuosity
  diffmanager->SetPhaseTort(matelchphase->Tortuosity(), 0);

  // loop over materials within electrolyte phase
  for (int imat = 0; imat < matelchphase->NumMat(); ++imat)
  {
    const Teuchos::RCP<const MAT::Material> material =
        matelchphase->MatById(matelchphase->MatID(imat));

    switch (material->MaterialType())
    {
      case INPAR::MAT::m_newman:
      case INPAR::MAT::m_newman_multiscale:
      {
        // safety check
        if (matelchphase->NumMat() != 1)
          dserror("Newman material must be the only transported species!");

        // set ion type
        diffcondmat = INPAR::ELCH::diffcondmat_newman;

        // evaluate standard Newman material
        if (material->MaterialType() == INPAR::MAT::m_newman)
          MatNewman(material, concentrations[0], diffmanager);

        // evaluate multi-scale Newman material
        else
          MatNewmanMultiScale(material, concentrations[0], diffmanager);

        break;
      }

      case INPAR::MAT::m_ion:
      {
        // set ion type
        diffcondmat = INPAR::ELCH::diffcondmat_ion;

        myelch::MatIon(material, imat, equpot, diffmanager);

        // calculation of conductivity and transference number based on diffusion coefficient and
        // valence
        if (imat == matelchphase->NumMat() - 1)
        {
          diffmanager->CalcConductivity(matelchphase->NumMat(), ffrt, concentrations);
          diffmanager->CalcTransNum(matelchphase->NumMat(), concentrations);
        }

        break;
      }

      default:
      {
        dserror("Invalid material type!");
        break;
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate standard Newman material                         fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatNewman(
    const Teuchos::RCP<const MAT::Material>& material,  //!< Newman material
    const double& concentration,                        //!< local concentration value
    const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& diffmanager  //!< diffusion manager
)
{
  // cast material to Newman material
  const Teuchos::RCP<const MAT::Newman> matnewman =
      Teuchos::rcp_static_cast<const MAT::Newman>(material);

  // valence of ionic species
  diffmanager->SetValence(matnewman->Valence(), 0);

  // concentration depending diffusion coefficient
  diffmanager->SetIsotropicDiff(matnewman->ComputeDiffusionCoefficient(concentration), 0);
  // derivation of concentration depending diffusion coefficient wrt all ionic species
  diffmanager->SetDerivIsoDiffCoef(matnewman->ComputeFirstDerivDiffCoeff(concentration), 0, 0);

  // concentration depending transference number
  diffmanager->SetTransNum(matnewman->ComputeTransferenceNumber(concentration), 0);
  // derivation of concentration depending transference number wrt all ionic species
  diffmanager->SetDerivTransNum(matnewman->ComputeFirstDerivTrans(concentration), 0, 0);

  // thermodynamic factor of electrolyte solution
  diffmanager->SetThermFac(matnewman->ComputeThermFac(concentration));
  // derivative of conductivity with respect to concentrations
  diffmanager->SetDerivThermFac(matnewman->ComputeFirstDerivThermFac(concentration), 0);

  // conductivity and first derivative can maximally depend on one concentration
  // since time curve is used as input routine
  // conductivity of electrolyte solution
  diffmanager->SetCond(matnewman->ComputeConductivity(concentration));
  // derivative of conductivity with respect to concentrations
  diffmanager->SetDerivCond(matnewman->ComputeFirstDerivCond(concentration), 0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate multi-scale Newman material                      fang 07/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatNewmanMultiScale(
    const Teuchos::RCP<const MAT::Material>& material,  //!< Newman material
    const double& concentration,                        //!< local concentration value
    const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& diffmanager  //!< diffusion manager
)
{
  // evaluate standard Newman material
  MatNewman(material, concentration, diffmanager);

  // cast material and diffusion manager
  const Teuchos::RCP<const MAT::NewmanMultiScale> newmanmultiscale =
      Teuchos::rcp_dynamic_cast<const MAT::NewmanMultiScale>(material);
  if (newmanmultiscale == Teuchos::null) dserror("Invalid material!");
  const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCondMultiScale> diffmanagermultiscale =
      Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElchDiffCondMultiScale>(diffmanager);
  if (diffmanagermultiscale == Teuchos::null) dserror("Invalid diffusion manager!");

  // set electronic conductivity
  diffmanagermultiscale->SetSigma(newmanmultiscale->Sigma());

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<DRT::Element::nurbs27>;
