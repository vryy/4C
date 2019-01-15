/*!
\file geometry_pair_element_types.H

\brief Element types that can be part of a pair. This types can be used as a template argument.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_TYPES_H_
#define GEOMETRY_PAIR_ELEMENT_TYPES_H_


#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Base class for the element type. The template parameters are stored in static const
   * member, so they can be accessed from derived classes.
   * @param discretization Type of shape function.
   * @param values_per_node Number of nodal values per node (standard elements have 1, Hermitian
   * shape functions have 2)
   */
  template <DRT::Element::DiscretizationType discretization, unsigned int values_per_node>
  class ElementDiscretizationBase
  {
   public:
    //! Type of shape function that will be used when evaluating the shape functions.
    static const DRT::Element::DiscretizationType discretization_ = discretization;

    //! Dimension of element (curve=1, surface=2, volume=3).
    static const unsigned int dim_ = DRT::UTILS::DisTypeToDim<discretization_>::dim;

    //! Number of values per node.
    static const unsigned int n_val_ = values_per_node;

    //! Number of nodes for this element.
    static const unsigned int n_nodes_ =
        DRT::UTILS::DisTypeToNumNodePerEle<discretization_>::numNodePerElement;

    //! Number of unknowns for this element.
    static const unsigned int n_dof_ = 3 * n_val_ * n_nodes_;

   public:
    /**
     * \brief Standard Constructor.
     */
    ElementDiscretizationBase(){};

    /**
     * \brief Destructor.
     */
    virtual ~ElementDiscretizationBase(){};
  };


  /**
   * \brief Diescretization based on standard (in most cases Lagrangian) shape functions.
   * @param discretization Type of shape function.
   */
  template <DRT::Element::DiscretizationType discretization>
  class ElementDiscretizationStandard : public ElementDiscretizationBase<discretization, 1>
  {
   public:
    /**
     * \brief Standard Constructor.
     */
    ElementDiscretizationStandard(){};

    /**
     * \brief Destructor.
     */
    virtual ~ElementDiscretizationStandard(){};
  };


  class ElementDiscretizationHermite
      : public ElementDiscretizationBase<DRT::Element::DiscretizationType::line2, 2>
  {
   public:
    /**
     * \brief Standard Constructor.
     */
    ElementDiscretizationHermite(){};

    /**
     * \brief Destructor.
     */
    virtual ~ElementDiscretizationHermite(){};
  };

}  // namespace GEOMETRYPAIR

#endif
