/*-----------------------------------------------------------*/
/*! \file

\brief  interface to Nox and LOCA for structural dynamics

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_LOCAINTERFACE_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_LOCAINTERFACE_H_

#include "str_timint_noxinterface.H"
#include "../loca_continuation/loca_nln_interface_required.H"

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
      virtual ~LocaInterface(){};

      //! setup function
      virtual void Setup();

      /*! compute the right hand side entries
       *  (derived from NOX::Epetra::Interface::Required) */
      virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

      /*! compute jacobian
       *  ( derived from NOX::Epetra::Inteface::Jacobian) */
      virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

      /*! compute right hand side and jacobian
       *  (derived from NOX::NLN::Interface::Jacobian) */
      virtual bool computeFandJacobian(
          const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac);

      /*! Set the continuation parameters for the compute routines
       *  (derived from LOCA::Epetra::Interface::Required) */
      virtual void setParameters(const LOCA::ParameterVector& p);

     protected:
      double displ_fac_;
    };
  }  // namespace TIMINT
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_LOCAINTERFACE_H_ */
