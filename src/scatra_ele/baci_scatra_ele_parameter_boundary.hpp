/*---------------------------------------------------------------------*/
/*! \file
\brief singleton class holding all interface parameters required for boundary element evaluation

\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_PARAMETER_BOUNDARY_HPP
#define BACI_SCATRA_ELE_PARAMETER_BOUNDARY_HPP

#include "baci_config.hpp"

#include "baci_inpar_s2i.hpp"
#include "baci_lib_condition.hpp"
#include "baci_scatra_ele_parameter_base.hpp"

#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    class ScaTraEleParameterBoundary : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterBoundary* Instance(const std::string& disname);

      //! set parameters for specific kintetic model
      void SetParameters(Teuchos::ParameterList& parameters) override;

      /// Set actual parameters
      ///@{
      void SetAlpha(Teuchos::ParameterList& parameters);
      void SetCapacitance(Teuchos::ParameterList& parameters);
      void SetChargeTransferConstant(Teuchos::ParameterList& parameters);
      void SetConvTolIterNum(Teuchos::ParameterList& parameters);
      void SetDensityMolarMass(Teuchos::ParameterList& parameters);
      void SetEnergySubstanceRatio(Teuchos::ParameterList& parameters);
      void SetIsPseudoContact(Teuchos::ParameterList& parameters);
      void SetNumElectrons(Teuchos::ParameterList& parameters);
      void SetNumScal(Teuchos::ParameterList& parameters);
      void SetOnOff(Teuchos::ParameterList& parameters);
      void SetPeltier(Teuchos::ParameterList& parameters);
      void SetPermeabilities(Teuchos::ParameterList& parameters);
      void SetRegularization(Teuchos::ParameterList& parameters);
      void SetResistance(Teuchos::ParameterList& parameters);
      void SetResistivity(Teuchos::ParameterList& parameters);
      void SetStoichiometries(Teuchos::ParameterList& parameters);
      void SetThermoPerm(Teuchos::ParameterList& parameters);
      ///@}

      /// Return parameters
      ///@{
      double AlphaA() const { return alphaa_; }
      double AlphaC() const { return alphac_; }
      double Capacitance() const { return capacitance_; }
      double ChargeTransferConstant() const { return kr_; }
      DRT::Condition::ConditionType ConditionType() const { return conditiontype_; }
      double ConvtolimplicitBV() const { return convtolimplicitBV_; }
      double Density() const { return density_; }
      bool IsPseudoContact() const { return is_pseudo_contact_; }
      int ItemaximplicitBV() const { return itemaxmimplicitBV_; }
      int KineticModel() const { return kineticmodel_; }
      double MolarHeatCapacity() const { return molar_heat_capacity_; }
      double MolarMass() const { return molarmass_; }
      int NumElectrons() const { return numelectrons_; }
      int NumScal() const { return numscal_; }
      const std::vector<int>* OnOff() const { return onoff_; }
      double Peltier() const { return peltier_; }
      const std::vector<double>* Permeabilities() const { return permeabilities_; }
      double RegularizationParameter() const { return regularizationparameter_; }
      INPAR::S2I::RegularizationType RegularizationType() const { return regularizationtype_; }
      double Resistance() const { return resistance_; }
      double Resistivity() const { return resistivity_; }
      const std::vector<int>* Stoichiometries() const { return stoichiometries_; }
      double ThermoPerm() const { return thermoperm_; }

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

      /// density of the interface layer material associated to S2IKineticsGrowth condition
      double density_;

      /// ratio between flux of species (amount) and flux of energy
      double molar_heat_capacity_;

      /// indicating if pseudo contact is considered at the s2i interface, i.e. no flux if interface
      /// is under tensile stresses
      bool is_pseudo_contact_;

      /// maximum number of iterations of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      int itemaxmimplicitBV_;

      /// kinetic model of s2i condition
      int kineticmodel_;

      /// charge transfer constant
      double kr_;

      /// molar mass of the interface layer material associated to S2IKineticsGrowth condition
      double molarmass_;

      /// number of electrons involved in charge transfer
      int numelectrons_;

      /// total number of transported scalars
      int numscal_;

      /// which dofs are constrained by this condition
      const std::vector<int>* onoff_;

      /// peltier coefficient
      double peltier_;

      /// vector of constant permeabilities
      const std::vector<double>* permeabilities_;

      /// regularization factor for S2IKineticsGrowth condition
      double regularizationparameter_;

      /// type of regularization for S2IKineticsGrowth condition
      INPAR::S2I::RegularizationType regularizationtype_;

      /// interface resistance associated with S2ICoupling condition
      double resistance_;

      /// resistivity of the interface layer material associated to S2IKineticsGrowth condition
      double resistivity_;

      /// surface capacitance at interface associated with S2ICoupling condition
      double capacitance_;

      /// vector of stoichiometric coefficients for scatra-scatra interface
      const std::vector<int>* stoichiometries_;

      /// constant permeability for heat at interface
      double thermoperm_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // SCATRA_ELE_PARAMETER_BOUNDARY_H
