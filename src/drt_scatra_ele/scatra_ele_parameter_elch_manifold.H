/*----------------------------------------------------------------------*/
/*! \file

\brief parameter class for elch on manifold kinetics

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_ELCH_MANIFOLD_H
#define SCATRA_ELE_PARAMETER_ELCH_MANIFOLD_H

#include "scatra_ele_parameter_base.H"
#include "../drt_inpar/inpar_ssi.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterElchManifold : public ScaTraEleParameterBase
    {
     public:
      //!  singleton access method
      //! \param disname       name of discretization
      //! \param delete_me     creation/destruction indication
      static ScaTraEleParameterElchManifold* Instance(const std::string& disname);

      void Done() override;

      void SetParameters(Teuchos::ParameterList& parameters) override;

      //! return read-only parameters
      //@{
      int KineticModel() const { return kinetic_model_; }
      int NumElectrons() const { return num_electrons_; }
      double Resistance() const { return resistance_; }
      bool EvaluateMasterSide() const { return evaluate_master_side_; }
      //@}

     private:
      //! private constructor for singletons
      ScaTraEleParameterElchManifold(const std::string& disname);

      //! parameters
      //@{
      bool evaluate_master_side_;
      int kinetic_model_;
      int num_electrons_;
      double resistance_;
      //@}

      //! map of all instance
      static std::map<std::string, std::unique_ptr<ScaTraEleParameterElchManifold>> instances_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
