/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_utils_elch_electrode.cpp

\brief utility class supporting element evaluation for electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_utils_elch_electrode.H"
#include "scatra_ele_calc_elch_electrode.H"

#include "../drt_mat/electrode.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>*
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(
    const int numdofpernode,                      ///< number of degrees of freedom per node
    const int numscal,                            ///< number of transported scalars per node
    const std::string& disname,                   ///< name of discretization
    const ScaTraEleUtilsElchElectrode* delete_me  ///< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleUtilsElchElectrode<distype>*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleUtilsElchElectrode<distype>(numdofpernode, numscal, disname);
  }

  // destruct instance
  else
  {
    for (typename std::map<std::string, ScaTraEleUtilsElchElectrode<distype>*>::iterator i =
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
void DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::ScaTraEleUtilsElchElectrode(
    const int numdofpernode,    ///< number of degrees of freedom per node
    const int numscal,          ///< number of transported scalars per node
    const std::string& disname  ///< name of discretization
    )
    : myelch::ScaTraEleUtilsElch(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate electrode material                               fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::MatElectrode(
    const Teuchos::RCP<const MAT::Material> material,                   //!< electrode material
    const double concentration,                                         //!< concentration
    const Teuchos::RCP<ScaTraEleDiffManagerElchElectrode>& diffmanager  //!< diffusion manager
)
{
  const MAT::Electrode* matelectrode = static_cast<const MAT::Electrode*>(material.get());

  // diffusion coefficient
  diffmanager->SetIsotropicDiff(matelectrode->ComputeDiffusionCoefficient(concentration), 0);

  // derivative of diffusion coefficient with respect to concentration
  diffmanager->SetDerivIsoDiffCoef(matelectrode->ComputeFirstDerivDiffCoeff(concentration), 0, 0);

  // electronic conductivity
  diffmanager->SetCond(matelectrode->ComputeConductivity(concentration));

  // derivative of electronic conductivity with respect to concentration
  diffmanager->SetDerivCond(matelectrode->ComputeFirstDerivCond(concentration), 0);

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::nurbs27>;
