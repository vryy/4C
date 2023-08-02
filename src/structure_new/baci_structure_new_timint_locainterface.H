/*-----------------------------------------------------------*/
/*! \file

\brief  interface to Nox and LOCA for structural dynamics


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_LOCAINTERFACE_H
#define BACI_STRUCTURE_NEW_TIMINT_LOCAINTERFACE_H

#include "baci_loca_continuation_nln_interface_required.H"
#include "baci_structure_new_timint_noxinterface.H"

namespace STR
{
  namespace TIMINT
  {
    class LocaInterface : public virtual STR::TIMINT::NoxInterface,
                          public virtual LOCA::NLN::Interface::Required
    {
     public:
      //! constructor
      LocaInterface();

      //! destructor
      ~LocaInterface() override = default;

      //! setup function
      void Setup() override;

      /*! compute the right hand side entries
       *  (derived from NOX::Epetra::Interface::Required) */
      bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override;

      /*! compute jacobian
       *  ( derived from NOX::Epetra::Inteface::Jacobian) */
      bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;

      /*! compute right hand side and jacobian
       *  (derived from NOX::NLN::Interface::Jacobian) */
      bool computeFandJacobian(
          const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac) override;

      /*! Set the continuation parameters for the compute routines
       *  (derived from LOCA::Epetra::Interface::Required) */
      void setParameters(const LOCA::ParameterVector& p) override;

     protected:
      double displ_fac_;
    };
  }  // namespace TIMINT
}  // namespace STR


#endif  // STRUCTURE_NEW_TIMINT_LOCAINTERFACE_H
