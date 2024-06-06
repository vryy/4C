/*---------------------------------------------------------------------*/
/*! \file
\brief Steepest ascent interface based on the augmented contact
       formulation.

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_STEEPEST_ASCENT_INTERFACE_HPP
#define FOUR_C_CONTACT_AUG_STEEPEST_ASCENT_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    namespace SteepestAscent
    {
      /** \brief Augmented steepest ascent contact interface class
       *
       *  \author hiermeier \date 03/17 */
      class Interface : public CONTACT::Aug::Interface
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
        Interface(const Teuchos::RCP<CONTACT::Aug::InterfaceDataContainer>& interfaceData_ptr);

        /// constructor
        Interface(const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData_ptr, int id,
            const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact,
            bool selfcontact);

       protected:
        Teuchos::RCP<CONTACT::Aug::INTERFACE::AssembleStrategy>
        create_node_based_assemble_strategy() override;

      };  // class Interface

      namespace INTERFACE
      {
        /*--------------------------------------------------------------------------*/
        /** \brief Node based assemble strategy
         *
         *  Assembly of nodal stored quantities. This a derived class which overloads
         *  some few methods of the base class in order to establish the desired
         *  condensed system contributions.
         *
         *  \author hiermeier \date 07/17 */
        template <typename assemble_policy>
        class NodeBasedAssembleStrategy
            : public CONTACT::Aug::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>
        {
         public:
          /// constructor
          explicit NodeBasedAssembleStrategy(Interface* inter)
              : CONTACT::Aug::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>(inter)
          { /* empty */
          }

          /// derived
          void Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const override;

          /// derived
          void assemble_dgg_lin_matrix(
              Core::LinAlg::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const override;
        };

      }  // namespace INTERFACE
    }    // namespace SteepestAscent
  }      // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
