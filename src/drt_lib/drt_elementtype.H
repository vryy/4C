/*----------------------------------------------------------------------*/
/*! \file

\brief Type definitions for elements

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_ELEMENTTPYE_H
#define DRT_ELEMENTTPYE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "drt_parobjectfactory.H"
#include "drt_element.H"

// forward declarations
namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{
  class Element;
  class Discretization;

  /// Subclass of ParObjectType that adds element type specific methods
  /*!
    Element types need to be initialized. Furthermore, there is a PreEvaluate
    method and the ability to read elements from dat files. And finally the
    element specific setup of null spaces for multi grid preconditioning is
    here, too.

    \note There are boundary elements that do not need all of this
    functionality.

    \author u.kue
    \date 06/10
   */
  class ElementType : public ParObjectType
  {
   protected:
    // only derived classes might create an instance
    ElementType();

   public:
    /// setup the dat file input line definitions for this type of element
    virtual void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
    {
    }

    /// create an element from a dat file specifier
    virtual Teuchos::RCP<Element> Create(
        const std::string eletype, const std::string eledistype, const int id, const int owner)
    {
      return Teuchos::null;
    }

    /// create an empty element
    virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner) = 0;

    /// initialize the element type
    virtual int Initialize(DRT::Discretization& dis);

    /// preevaluation
    virtual inline void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
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
    virtual void NodalBlockInformation(
        DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) = 0;

    /// do the null space computation
    virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
        const double* x0, int numdf, int dimns) = 0;
  };

}  // namespace DRT

#endif
