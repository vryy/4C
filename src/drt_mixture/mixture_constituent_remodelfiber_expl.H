/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a remodel constituent with explicit update rule

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_H
#define BACI_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_H

#include <Teuchos_RCPDecl.hpp>
#include "../drt_mat/matpar_material.H"
#include "mixture_constituent_remodelfiber.H"

namespace MIXTURE
{
  class MixtureConstituent;

  namespace PAR
  {
    class MixtureConstituent_RemodelFiberExpl : public MIXTURE::PAR::MixtureConstituent_RemodelFiber
    {
     public:
      explicit MixtureConstituent_RemodelFiberExpl(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      ~MixtureConstituent_RemodelFiberExpl() override = default;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::MixtureConstituent> CreateConstituent() override;
    };
  }  // namespace PAR

  /*!
   * \brief Remodel fiber constituent with an explicit update rule
   */
  class MixtureConstituent_RemodelFiberExpl : public MIXTURE::MixtureConstituent_RemodelFiber
  {
   public:
    explicit MixtureConstituent_RemodelFiberExpl(
        MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl* params);

    /// Returns the material type enum
    INPAR::MAT::MaterialType MaterialType() const override;

    /*!
     * \brief Update rule to evaluate and update the growth evolution
     *
     * \param F Deformation gradient of the finished time step
     * \param iFext inverse of the inelastic external deformation gradient
     * \param params Container for additional information
     * \param dt Time step
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void UpdateElasticPart(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFext,
        Teuchos::ParameterList& params, double dt, int gp, int eleGID) override;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_H
