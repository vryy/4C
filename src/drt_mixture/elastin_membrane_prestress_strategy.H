/*----------------------------------------------------------------------*/
/*! \file
\brief Prestress strategy for the membrane part of elastin

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_ELASTIN_MEMBRANE_PRESTRESS_STRATEGY_H
#define BACI_MIXTURE_ELASTIN_MEMBRANE_PRESTRESS_STRATEGY_H

#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations
namespace Teuchos
{
  class ParameterList;
}
namespace MAT
{
  class CylinderCoordinateSystemProvider;
}

namespace MIXTURE
{
  // Forward declaration of mixture classes
  class MixtureConstituent;
  class MixtureRule;
  class ElastinMembraneEvaluation;

  /*!
   * \brief Prestress strategy for the elastin membrane material used for growth and remodeling
   * simulations with homogenized constrained mixtures
   */
  class ElastinMembranePrestressStrategy
  {
   public:
    /*!
     * \brief Evaluates the fraction of the elastin membrane contribution to the total elastin
     * response to ensure equilibrium in the reference configuration
     *
     * \param mixtureRule Mixture rule
     * \param anisotropy Cylinder coordinate system
     * \param constituent Constituent that needs to be prestressed
     * \param membraneEvaluation Evaluator of the membrane material
     * \param params Container for additional information
     * \param gp Gauss point
     * \param eleGID Global element id
     * \return double Fraction of the elastin membrane contribution to the total elastin response
     */
    virtual double EvaluateMueFrac(MixtureRule& mixtureRule,
        const MAT::CylinderCoordinateSystemProvider& anisotropy,
        MIXTURE::MixtureConstituent& constituent, ElastinMembraneEvaluation& membraneEvaluation,
        Teuchos::ParameterList& params, int gp, int eleGID) const = 0;
  };

  /*!
   * \brief Evaluation of the membrane material of elastin
   */
  class ElastinMembraneEvaluation
  {
   public:
    /*!
     * \brief Evaluates the 2D membrane stress
     *
     * \param S (out)     : Membrane 2. Piola-Kirchhoff stress tensor in stress like Voigt notation
     * \param params (in) : Container for additional information
     * \param gp (in)     : Gauss point
     * \param eleGID (in) : Global element id
     */
    virtual void EvaluateMembraneStress(
        LINALG::Matrix<6, 1>& S, Teuchos::ParameterList& params, int gp, int eleGID) = 0;
  };
}  // namespace MIXTURE
#endif  // BACI_MIXTURE_ELASTIN_MEMBRANE_PRESTRESS_STRATEGY_H
