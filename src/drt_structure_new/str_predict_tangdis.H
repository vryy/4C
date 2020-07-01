/*-----------------------------------------------------------*/
/*! \file

\brief Tangential displacement predictor.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_PREDICT_TANGDIS_H_
#define STR_PREDICT_TANGDIS_H_

#include "str_predict_generic.H"
#include "../solver_nonlin_nox/nox_nln_abstract_prepostoperator.H"

// forward declaration
class Epetra_Vector;
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

      virtual ~TangDis(){};

      //! setup class specific stuff
      virtual void Setup();

      //! do the class specific predictor step
      virtual void Compute(NOX::Abstract::Group& grp);

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
          TangDis(const Teuchos::RCP<const ::STR::PREDICT::TangDis>& tang_predict_ptr);

          //! destructor
          virtual ~TangDis(){};

          //! add the linear reaction forces
          virtual void runPostComputeF(Epetra_Vector& F, const NOX::NLN::Group& grp);

         private:
          //! pointer to the tangdis object (read-only)
          Teuchos::RCP<const ::STR::PREDICT::TangDis> tang_predict_ptr_;
        };  // class TangDis
      }     // namespace PrePostOp
    }       // namespace GROUP
  }         // namespace NLN
}  // namespace NOX

#endif /* STR_PREDICT_TANGDIS_H_ */
