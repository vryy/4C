/*-----------------------------------------------------------*/
/*! \file

\brief LOCA nonlinear interface required defintion


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_LOCA_CONTINUATION_NLN_INTERFACE_REQUIRED_H
#define BACI_LOCA_CONTINUATION_NLN_INTERFACE_REQUIRED_H

#include <LOCA_Epetra_Interface_Required.H>
#include "baci_solver_nonlin_nox_interface_required.H"

namespace LOCA
{
  namespace NLN
  {
    namespace Interface
    {
      class Required : public virtual LOCA::Epetra::Interface::Required,
                       public virtual NOX::NLN::Interface::Required
      {
       public:
        //! constructor
        Required(){};

        //! destructor
        virtual ~Required(){};
      };  // class Required
    }     // namespace Interface
  }       // namespace NLN
}  // namespace LOCA


#endif  // LOCA_CONTINUATION_NLN_INTERFACE_REQUIRED_H
