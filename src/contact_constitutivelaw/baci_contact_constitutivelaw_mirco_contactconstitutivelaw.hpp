/*----------------------------------------------------------------------*/
/*! \file
\brief Implements a default mirco function as contact constitutive law

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_CONTACT_CONSTITUTIVELAW_MIRCO_CONTACTCONSTITUTIVELAW_HPP
#define BACI_CONTACT_CONSTITUTIVELAW_MIRCO_CONTACTCONSTITUTIVELAW_HPP

#include "baci_config.hpp"

#include "baci_contact_constitutivelaw_contactconstitutivelaw.hpp"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "baci_linalg_serialdensematrix.hpp"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    /*----------------------------------------------------------------------*/
    /** \brief constitutive law parameters for a mirco contact law to the contact pressure
     *
     */
    class MircoConstitutiveLawParams : public Parameter
    {
     public:
      /** \brief standard constructor
       * \param[in] container containing the law parameter from the input file
       */
      MircoConstitutiveLawParams(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container);

      /// create constitutive law instance of matching type with my parameters
      Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> CreateConstitutiveLaw() override;

      /// @name get-functions for the Constitutive Law parameters of a mirco function
      //@{

      int GetFirstMatID() const { return firstmatid_; };
      int GetSecondMatID() const { return secondmatid_; };
      double GetLateralLength() const { return lateralLength_; };
      double GetTolerance() const { return tolerance_; };
      double GetMaxIteration() const { return maxIteration_; };
      bool GetWarmStartingFlag() const { return warmStartingFlag_; };
      Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> GetTopology() const { return topology_; };
      double GetCompositeYoungs() const { return compositeYoungs_; };
      double GetGridSize() const { return gridSize_; };
      double GetComplianceCorrection() const { return elasticComplianceCorrection_; };
      double GetFiniteDifferenceFraction() const { return finiteDifferenceFraction_; };
      double GetActiveGapTolerance() const { return activeGapTolerance_; };
      Teuchos::Ptr<std::vector<double>> GetMeshGrid() const { return meshgrid_; };

      void SetParameters();

     private:
      /// @name Constitutive Law parameters of a mirco function
      //@{

      int firstmatid_;
      int secondmatid_;
      double lateralLength_;
      int resolution_;
      bool randomTopologyFlag_;
      bool randomSeedFlag_;
      int randomGeneratorSeed_;
      double tolerance_;
      int maxIteration_;
      bool warmStartingFlag_;
      Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> topology_;
      double compositeYoungs_;
      double gridSize_;
      double elasticComplianceCorrection_;
      Teuchos::Ptr<std::vector<double>> meshgrid_;
      double finiteDifferenceFraction_;
      double activeGapTolerance_;
      std::string topologyFilePath_;
      //@}
    };  // class

    /*----------------------------------------------------------------------*/
    /** \brief implements a mirco contact constitutive law relating the gap to the
     * contact pressure
     */
    class MircoConstitutiveLaw : public ConstitutiveLaw
    {
     public:
      /// construct the constitutive law object given a set of parameters
      explicit MircoConstitutiveLaw(CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams* params);

      //! @name Access methods

      /// contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const override
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_mirco;
      }

      /// Return quick accessible contact constitutive law parameter data
      CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const override { return params_; }

      //! @name Evaluation methods
      //@{
      /** \brief evaluate the constitutive law
       *
       * The pressure response for a gap is calucated using MIRCO, which uses BEM for solving
       * contact between a rigid rough surface and a linear elastic half plane.
       *
       * \param gap contact gap at the mortar node
       * \return The pressure response from MIRCO
       */
      double Evaluate(double gap, CONTACT::RoughNode* cnode) override;

      /** \brief Evaluate derivative of the constitutive law
       *
       * The derivative of the pressure response is approximated using a finite difference approach
       * by calling MIRCO twice at two different gap values and doing a backward difference
       * approximation for the linearization.
       *
       * \param gap contact gap at the mortar node
       * \return Derivative of the pressure responses from MIRCO
       */
      double EvaluateDeriv(double gap, CONTACT::RoughNode* cnode) override;
      //@}

     private:
      /// my constitutive law parameters
      CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams* params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // CONTACT_CONSTITUTIVELAW_MIRCO_CONTACTCONSTITUTIVELAW_H
