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

namespace Solid
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
      virtual void init(const enum Inpar::Solid::PredEnum& type,
          const Teuchos::RCP<Solid::IMPLICIT::Generic>& implint_ptr,
          const Teuchos::RCP<Solid::Dbc>& dbc_ptr,
          const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<Solid::TimeInt::BaseDataIO>& iodata_ptr,
          const Teuchos::RCP<Teuchos::ParameterList>& noxparams_ptr);

      //! setup of the specific predictor
      virtual void setup() = 0;

      //! Get the predictor type enum
      const Inpar::Solid::PredEnum& get_type() const { return type_; };

      //! returns the name of the used predictor
      virtual std::string name() const;

      //! Preprocess the predictor step
      virtual void pre_predict(::NOX::Abstract::Group& grp);

      //! Pre-/Postprocess the specific predictor step
      void predict(::NOX::Abstract::Group& grp);

      //! Calculate the specific predictor step
      virtual void compute(::NOX::Abstract::Group& grp) = 0;

      //! Postprocess the predictor step
      virtual void post_predict(::NOX::Abstract::Group& grp);

      //! return a constant reference to the global state object (read only)
      const Solid::TimeInt::BaseDataGlobalState& global_state() const;

      //! print the result of the predictor step
      void print() const;

      //! Run before the external force are computed and assembled
      virtual bool pre_apply_force_external(Epetra_Vector& fextnp) const;

     protected:
      //! returns init state
      const bool& is_init() const { return isinit_; };

      //! returns setup state
      const bool& is_setup() const { return issetup_; };

      void check_init() const;

      void check_init_setup() const;

      Teuchos::RCP<Solid::IMPLICIT::Generic>& impl_int_ptr();
      Solid::IMPLICIT::Generic& impl_int();

      Teuchos::RCP<Solid::Dbc>& dbc_ptr();
      Solid::Dbc& dbc();

      Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& global_state_ptr();
      Solid::TimeInt::BaseDataGlobalState& global_state();

      Teuchos::RCP<Solid::TimeInt::BaseDataIO>& io_data_ptr();
      Solid::TimeInt::BaseDataIO& io_data();

      Teuchos::RCP<Teuchos::ParameterList>& nox_params_ptr();
      Teuchos::ParameterList& nox_params();

     protected:
      //! indicates if the init() function has been called
      bool isinit_;

      //! indicates if the setup() function has been called
      bool issetup_;

     private:
      //! predictor type
      enum Inpar::Solid::PredEnum type_;

      //! pointer to the implicit integrator
      Teuchos::RCP<Solid::IMPLICIT::Generic> implint_ptr_;

      //! pointer to the dirichlet boundary condition object
      Teuchos::RCP<Solid::Dbc> dbc_ptr_;

      //! global state pointer
      Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState> gstate_ptr_;

      //! input/output data pointer
      Teuchos::RCP<Solid::TimeInt::BaseDataIO> iodata_ptr_;

      Teuchos::RCP<Teuchos::ParameterList> noxparams_ptr_;
    };  // class  Generic
  }     // namespace Predict
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
