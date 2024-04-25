/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for the Lagrange solving strategy of the augmented
       framework.

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_LAGRANGE_INTERFACE_HPP
#define FOUR_C_CONTACT_AUG_LAGRANGE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
  {
    namespace LAGRANGE
    {
      class Interface : public CONTACT::AUG::Interface
      {
       public:
        /** \brief Alternative constructor
         *
         *  A prerequisite for this constructor is, that the passed
         *  shared interface data object has been filled/initialized already.
         *
         *  \param interfaceData_ptr (in) : filled shared augmented contact interface
         *                          data container object
         *
         *  \author hiermeier \date 03/17 */
        Interface(const Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer>& idata_ptr);

        /// constructor
        Interface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr,
            const int id, const Epetra_Comm& comm, const int dim,
            const Teuchos::ParameterList& icontact, const bool selfcontact);

       protected:
        /** \brief Assemble the linearization matrix contributions for the
         *         augmentation term [derived]
         *
         *  This term vanishes in the standard Lagrange formulation.
         *
         *  \author hiermeier \date 03/17 */
        void AssembleDGGLinMatrix(
            CORE::LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const override{};

      };  // class Interface
    }     // namespace LAGRANGE
  }       // namespace AUG
}  // namespace CONTACT



FOUR_C_NAMESPACE_CLOSE

#endif
