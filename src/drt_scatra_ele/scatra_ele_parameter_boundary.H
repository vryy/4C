/*---------------------------------------------------------------------*/
/*! \file
\brief singleton class holding all interface parameters required for boundary element evaluation

\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef SCATRA_ELE_PARAMETER_BOUNDARY_H
#define SCATRA_ELE_PARAMETER_BOUNDARY_H

#include "scatra_ele_parameter_base.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_inpar/inpar_s2i.H"
#include "vector"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    class ScaTraEleParameterBoundary : public ScaTraEleParameterBase
    {
     public:
      //! destructor
      virtual ~ScaTraEleParameterBoundary() { return; };

      //! singleton access method
      static ScaTraEleParameterBoundary* Instance(
          const std::string& disname, const ScaTraEleParameterBoundary* delete_me = nullptr);

      //! singleton destruction
      virtual void Done();

      //! set parameters for specific kintetic model
      void SetParameters(Teuchos::ParameterList& parameters) override;

      /// Set actual parameters
      ///@{
      void SetAlpha(Teuchos::ParameterList& parameters);
      void SetChargeTransferConstant(Teuchos::ParameterList& parameters);
      void SetConvTolIterNum(Teuchos::ParameterList& parameters);
      void SetDensityMolarMass(Teuchos::ParameterList& parameters);
      void SetNumElectrons(Teuchos::ParameterList& parameters);
      void SetNumScal(Teuchos::ParameterList& parameters);
      void SetPeltier(Teuchos::ParameterList& parameters);
      void SetPermeabilities(Teuchos::ParameterList& parameters);
      void SetRegularization(Teuchos::ParameterList& parameters);
      void SetResistance(Teuchos::ParameterList& parameters);
      void SetResistivity(Teuchos::ParameterList& parameters);
      void SetStoichiometries(Teuchos::ParameterList& parameters);
      ///@}

      /// Return parameters
      ///@{
      double AlphaA() const { return alphaa_; };
      double AlphaC() const { return alphac_; };
      double ChargeTransferConstant() const { return kr_; };
      DRT::Condition::ConditionType ConditionType() const { return conditiontype_; };
      double ConvtolimplicitBV() const { return convtolimplicitBV_; };
      double Density() const { return density_; };
      double ItemaximplicitBV() const { return itemaxmimplicitBV_; };
      int KineticModel() const { return kineticmodel_; };
      double MolarMass() const { return molarmass_; };
      int NumElectrons() const { return numelectrons_; };
      int NumScal() const { return numscal_; };
      double Peltier() const { return peltier_; };
      std::vector<double>* Permeabilities() const { return permeabilities_; };
      double RegularizationParameter() const { return regularizationparameter_; };
      INPAR::S2I::RegularizationType RegularizationType() const { return regularizationtype_; };
      double Resistance() const { return resistance_; };
      double Resistivity() const { return resistivity_; };
      std::vector<int>* Stoichiometries() const { return stoichiometries_; };
      ///@}

     private:
      //! private constructor for singletons
      ScaTraEleParameterBoundary(const std::string& disname);

      /// symmetry coefficient of anodic intercalation reaction
      double alphaa_;

      /// symmetry coefficient of cathodic intercalation reaction
      double alphac_;

      /// condition type of stored condition
      DRT::Condition::ConditionType conditiontype_;

      /// convergence tolerance of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      double convtolimplicitBV_;

      /// density of the interface layer material associated to S2ICouplingGrowth condition
      double density_;

      /// maximum number of iterations of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      int itemaxmimplicitBV_;

      /// kinetic model of s2i condition
      int kineticmodel_;

      /// charge transfer constant
      double kr_;

      /// molar mass of the interface layer material associated to S2ICouplingGrowth condition
      double molarmass_;

      /// number of electrons involved in charge transfer
      int numelectrons_;

      /// total number of transported scalars
      int numscal_;

      /// peltier coefficient
      double peltier_;

      /// vector of constant permeabilities
      std::vector<double>* permeabilities_;

      /// regularization factor for S2ICouplingGrowth condition
      double regularizationparameter_;

      /// type of regularization for S2ICouplingGrowth condition
      INPAR::S2I::RegularizationType regularizationtype_;

      /// interface resistance associated with S2ICoupling condition
      double resistance_;

      /// resistivity of the interface layer material associated to S2ICouplingGrowth condition
      double resistivity_;

      /// vector of stoichiometric coefficients for scatra-scatra interface
      std::vector<int>* stoichiometries_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SCATRA_ELE_PARAMETER_BOUNDARY_H */
