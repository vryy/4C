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
          const std::string& disname,  //!< name of discretization
          const ScaTraEleParameterBoundary* delete_me =
              nullptr  //!< creation/destruction indication
      );

      //! singleton destruction
      virtual void Done();

      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      /// return symmetry coefficient of anodic intercalation reaction
      double AlphaA() const { return alphaa_; };

      /// return symmetry coefficient of cathodic intercalation reaction
      double AlphaC() const { return alphac_; };

      /// return the condition type
      DRT::Condition::ConditionType ConditionType() const { return conditiontype_; };

      /// return the density associated to the material of the S2ICouplingGrowth condition
      double Density() const { return density_; };

      /// return the kinetic model
      int KineticModel() const { return kineticmodel_; };

      /// return charge transfer constant
      double Kr() const { return kr_; };

      /// return the molar mass associated to the material of the S2ICouplingGrowth condition
      double MolarMass() const { return molarmass_; };

      /// return number of electrons involved in charge transfer
      int NumElectrons() const { return numelectrons_; };

      /// return total number of existing scalars
      int NumScal() const { return numscal_; };

      /// return peltier
      double Peltier() const { return peltier_; };

      /// return vector of constant permeabilities
      std::vector<double>* Permeabilities() const { return permeabilities_; };

      /// return the regularization parameter associated to the S2ICouplingGrowth condition
      double RegularizationParameter() const { return regularizationparameter_; };

      /// return the type of regularization associated to the S2ICouplingGrowth condition
      INPAR::S2I::RegularizationType RegularizationType() const { return regularizationtype_; };

      /// return the resistivity associated to the material of the S2ICouplingGrowth condition
      double Resistivity() const { return resistivity_; };

      /// return vector of stoichiometric coefficients for scatra-scatra interface
      std::vector<int>* Stoichiometries() const { return stoichiometries_; };

      /// return interface resistance associated with S2ICoupling condition
      double Resistance() const { return resistance_; };

      /// return convergence tolerance of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      double ConvtolimplicitBV() const { return convtolimplicitBV_; };

      /// maximum number of iterations of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      double ItemaximplicitBV() const { return itemaxmimplicitBV_; };

     private:
      //! private constructor for singletons
      ScaTraEleParameterBoundary(const std::string& disname  //!< name of discretization
      );

      /// symmetry coefficient of anodic intercalation reaction
      double alphaa_;

      /// symmetry coefficient of cathodic intercalation reaction
      double alphac_;

      /// condition type of stored condition
      DRT::Condition::ConditionType conditiontype_;

      /// density of the interface layer material associated to S2ICouplingGrowth condition
      double density_;

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

      /// resistivity of the interface layer material associated to S2ICouplingGrowth condition
      double resistivity_;

      /// vector of stoichiometric coefficients for scatra-scatra interface
      std::vector<int>* stoichiometries_;

      /// interface resistance associated with S2ICoupling condition
      double resistance_;

      /// convergence tolerance of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      double convtolimplicitBV_;

      /// maximum number of iterations of local Newton-Raphson iteration for modified Butler-Volmer
      /// equation
      double itemaxmimplicitBV_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SCATRA_ELE_PARAMETER_BOUNDARY_H */
