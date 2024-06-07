/*----------------------------------------------------------------------*/
/*! \file

\brief Type definitions for elements

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_ELEMENTTYPE_HPP
#define FOUR_C_FEM_GENERAL_ELEMENTTYPE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseOperator;
}

namespace Discret
{
  class Discretization;
}
/// Subclass of ParObjectType that adds element type specific methods
/*!
  Element types need to be initialized. Furthermore, there is a pre_evaluate
  method and the ability to read elements from dat files. And finally the
  element specific setup of null spaces for multi grid preconditioning is
  here, too.

  \note There are boundary elements that do not need all of this
  functionality.

  \author u.kue
  \date 06/10
 */

namespace Core::Elements
{
  class ElementType : public Core::Communication::ParObjectType
  {
   protected:
    // only derived classes might create an instance
    ElementType();

   public:
    /// setup the dat file input line definitions for this type of element
    virtual void setup_element_definition(
        std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
    {
    }

    /// create an element from a dat file specifier
    virtual Teuchos::RCP<Element> Create(
        const std::string eletype, const std::string eledistype, const int id, const int owner)
    {
      return Teuchos::null;
    }

    /// create an empty element
    virtual Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) = 0;

    /// initialize the element type
    virtual int Initialize(Discret::Discretization& dis);

    /// preevaluation
    virtual inline void pre_evaluate(Discret::Discretization& dis, Teuchos::ParameterList& p,
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3)
    {
      return;
    }

    /*!
    \brief Get nodal block information to create a null space description

    \note All elements will fill \c nv, but \c np is only filled by some elements.

    @param[in] dwele Element pointer
    @param[out] numdf Number of degrees of freedom per node
    @param[out] dimns Nullspace dimension
    @param[out] nv Number of degrees of freedom for balance of linear momentum (e.g. solid
                displacement, fluid velocity)
    @param[out] np Number of degrees of freedom for local constraints (e.g. fluid pressure)
    */
    virtual void nodal_block_information(
        Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) = 0;

    /// do the null space computation
    virtual Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
        Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) = 0;
  };

}  // namespace Core::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
