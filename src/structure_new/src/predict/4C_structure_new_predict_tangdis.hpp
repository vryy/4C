/*-----------------------------------------------------------*/
/*! \file

\brief Tangential displacement predictor.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_TANGDIS_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_TANGDIS_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_abstract_prepostoperator.hpp"
#include "4C_structure_new_predict_generic.hpp"

// forward declaration
class Epetra_Vector;

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    class Group;
  }  // namespace Nln
}  // namespace NOX
namespace STR
{
  namespace Predict
  {
    class TangDis : public Generic
    {
     public:
      TangDis();
      //! setup class specific stuff
      void setup() override;

      //! do the class specific predictor step
      void Compute(::NOX::Abstract::Group& grp) override;

      //! return the dbc increment
      const Epetra_Vector& get_dbc_incr() const;

      //! return the switch for the pre/post operator
      const bool& is_apply_linear_reaction_forces() const;

      //! derived
      bool pre_apply_force_external(Epetra_Vector& fextnp) const override;

     private:
      Teuchos::RCP<Epetra_Vector> dbc_incr_ptr_;

      bool apply_linear_reaction_forces_;
    };  // class TangDis
  }     // namespace Predict
}  // namespace STR

namespace NOX
{
  namespace Nln
  {
    namespace GROUP
    {
      namespace PrePostOp
      {
        /*! \brief Tangential Displacement helper class
         *
         *  This class is an implementation of the NOX::Nln::Abstract::PrePostOperator
         *  and is used to modify the computeF() routines of the given NOX::Nln::Group
         *  (see STR::Predict::TangDis). It's called by the wrapper class
         *  NOX::Nln::GROUP::PrePostOperator.
         *
         *  \author Michael Hiermeier */
        class TangDis : public NOX::Nln::Abstract::PrePostOperator
        {
         public:
          //! constructor
          TangDis(const Teuchos::RCP<const STR::Predict::TangDis>& tang_predict_ptr);

          //! add the linear reaction forces
          void runPostComputeF(Epetra_Vector& F, const NOX::Nln::Group& grp) override;

         private:
          //! pointer to the tangdis object (read-only)
          Teuchos::RCP<const STR::Predict::TangDis> tang_predict_ptr_;
        };  // class TangDis
      }     // namespace PrePostOp
    }       // namespace GROUP
  }         // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
