/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all predictors.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_GENERIC_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_GENERIC_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

#include <Teuchos_RCP.hpp>

// forward declaration ...
class Epetra_Vector;
namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  class Dbc;
  namespace IMPLICIT
  {
    class Generic;
  }  // namespace IMPLICIT
  namespace TimeInt
  {
    class BaseDataGlobalState;
    class BaseDataIO;
  }  // namespace TimeInt
  namespace Predict
  {
    class Generic
    {
     public:
      //! constructor
      Generic();

      //! destructor
      virtual ~Generic() = default;

      //! initialize the base class variables
      virtual void Init(const enum Inpar::STR::PredEnum& type,
          const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr,
          const Teuchos::RCP<STR::Dbc>& dbc_ptr,
          const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TimeInt::BaseDataIO>& iodata_ptr,
          const Teuchos::RCP<Teuchos::ParameterList>& noxparams_ptr);

      //! setup of the specific predictor
      virtual void Setup() = 0;

      //! Get the predictor type enum
      const Inpar::STR::PredEnum& get_type() const { return type_; };

      //! returns the name of the used predictor
      virtual std::string Name() const;

      //! Preprocess the predictor step
      virtual void pre_predict(::NOX::Abstract::Group& grp);

      //! Pre-/Postprocess the specific predictor step
      void Predict(::NOX::Abstract::Group& grp);

      //! Calculate the specific predictor step
      virtual void Compute(::NOX::Abstract::Group& grp) = 0;

      //! Postprocess the predictor step
      virtual void post_predict(::NOX::Abstract::Group& grp);

      //! return a constant reference to the global state object (read only)
      const STR::TimeInt::BaseDataGlobalState& global_state() const;

      //! print the result of the predictor step
      void Print() const;

      //! Run before the external force are computed and assembled
      virtual bool pre_apply_force_external(Epetra_Vector& fextnp) const;

     protected:
      //! returns init state
      const bool& is_init() const { return isinit_; };

      //! returns setup state
      const bool& is_setup() const { return issetup_; };

      void check_init() const;

      void check_init_setup() const;

      Teuchos::RCP<STR::IMPLICIT::Generic>& impl_int_ptr();
      STR::IMPLICIT::Generic& impl_int();

      Teuchos::RCP<STR::Dbc>& dbc_ptr();
      STR::Dbc& dbc();

      Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& global_state_ptr();
      STR::TimeInt::BaseDataGlobalState& global_state();

      Teuchos::RCP<STR::TimeInt::BaseDataIO>& io_data_ptr();
      STR::TimeInt::BaseDataIO& io_data();

      Teuchos::RCP<Teuchos::ParameterList>& nox_params_ptr();
      Teuchos::ParameterList& nox_params();

     protected:
      //! indicates if the Init() function has been called
      bool isinit_;

      //! indicates if the Setup() function has been called
      bool issetup_;

     private:
      //! predictor type
      enum Inpar::STR::PredEnum type_;

      //! pointer to the implicit integrator
      Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr_;

      //! pointer to the dirichlet boundary condition object
      Teuchos::RCP<STR::Dbc> dbc_ptr_;

      //! global state pointer
      Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> gstate_ptr_;

      //! input/output data pointer
      Teuchos::RCP<STR::TimeInt::BaseDataIO> iodata_ptr_;

      Teuchos::RCP<Teuchos::ParameterList> noxparams_ptr_;
    };  // class  Generic
  }     // namespace Predict
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
