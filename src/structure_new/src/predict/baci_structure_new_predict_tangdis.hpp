/*-----------------------------------------------------------*/
/*! \file

\brief Tangential displacement predictor.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_PREDICT_TANGDIS_HPP
#define BACI_STRUCTURE_NEW_PREDICT_TANGDIS_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_abstract_prepostoperator.hpp"
#include "baci_structure_new_predict_generic.hpp"

// forward declaration
class Epetra_Vector;

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    class Group;
  }  // namespace NLN
}  // namespace NOX
namespace STR
{
  namespace PREDICT
  {
    class TangDis : public Generic
    {
     public:
      TangDis();
      //! setup class specific stuff
      void Setup() override;

      //! do the class specific predictor step
      void Compute(::NOX::Abstract::Group& grp) override;

      //! return the dbc increment
      const Epetra_Vector& GetDbcIncr() const;

      //! return the switch for the pre/post operator
      const bool& IsApplyLinearReactionForces() const;

      //! derived
      bool PreApplyForceExternal(Epetra_Vector& fextnp) const override;

     private:
      Teuchos::RCP<Epetra_Vector> dbc_incr_ptr_;

      bool applyLinearReactionForces_;
    };  // class TangDis
  }     // namespace PREDICT
}  // namespace STR

namespace NOX
{
  namespace NLN
  {
    namespace GROUP
    {
      namespace PrePostOp
      {
        /*! \brief Tangential Displacement helper class
         *
         *  This class is an implementation of the NOX::NLN::Abstract::PrePostOperator
         *  and is used to modify the computeF() routines of the given NOX::NLN::Group
         *  (see STR::PREDICT::TangDis). It's called by the wrapper class
         *  NOX::NLN::GROUP::PrePostOperator.
         *
         *  \author Michael Hiermeier */
        class TangDis : public NOX::NLN::Abstract::PrePostOperator
        {
         public:
          //! constructor
          TangDis(const Teuchos::RCP<const STR::PREDICT::TangDis>& tang_predict_ptr);


          //! add the linear reaction forces
          void runPostComputeF(Epetra_Vector& F, const NOX::NLN::Group& grp) override;

         private:
          //! pointer to the tangdis object (read-only)
          Teuchos::RCP<const STR::PREDICT::TangDis> tang_predict_ptr_;
        };  // class TangDis
      }     // namespace PrePostOp
    }       // namespace GROUP
  }         // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_PREDICT_TANGDIS_H
