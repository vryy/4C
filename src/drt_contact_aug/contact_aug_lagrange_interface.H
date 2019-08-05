/*----------------------------------------------------------------------*/
/*!
\brief Interface class for the Lagrange solving strategy of the augmented
       framework.

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_INTERFACE_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_INTERFACE_H_

#include "contact_augmented_interface.H"

namespace CONTACT
{
  namespace AUG
  {
    namespace LAGRANGE
    {
      class Interface : public ::CONTACT::AUG::Interface
      {
       public:
        /** \brief Alternative constructor
         *
         *  A prerequisite for this constructor is, that the passed
         *  shared interface data object has been filled/initialized already.
         *
         *  \param idata_ptr (in) : filled shared augmented contact interface
         *                          data container object
         *
         *  \author hiermeier \date 03/17 */
        Interface(const Teuchos::RCP<CONTACT::AUG::IDataContainer>& idata_ptr);

        /// constructor
        Interface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr,
            const int id, const Epetra_Comm& comm, const int dim,
            const Teuchos::ParameterList& icontact, const bool selfcontact,
            INPAR::MORTAR::RedundantStorage redundant);

        /// destructor
        virtual ~Interface(){};

       protected:
        /** \brief Assemble the linearization matrix contributions for the
         *         augmentation term [derived]
         *
         *  This term vanishes in the standard Lagrange formulation.
         *
         *  \author hiermeier \date 03/17 */
        virtual void AssembleDGGLinMatrix(
            LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const;

      };  // class Interface
    }     // namespace LAGRANGE
  }       // namespace AUG
}  // namespace CONTACT



#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_LAGRANGE_INTERFACE_H_ */
