/*
 * loca_nln_interface_required.H
 *
 *  Created on: Nov 19, 2015
 *      Author: hiermeier
 */

#ifndef SRC_LOCA_CONTINUATION_LOCA_NLN_INTERFACE_REQUIRED_H_
#define SRC_LOCA_CONTINUATION_LOCA_NLN_INTERFACE_REQUIRED_H_

#include <LOCA_Epetra_Interface_Required.H>
#include "../solver_nonlin_nox/nox_nln_interface_required.H"

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


#endif /* SRC_LOCA_CONTINUATION_LOCA_NLN_INTERFACE_REQUIRED_H_ */
