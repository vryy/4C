/*!
\file geometry_pair_element_types.H

\brief Element types that can be part of a pair. This types can be used as a template argument.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_TYPES_H_
#define GEOMETRY_PAIR_ELEMENT_TYPES_H_


#include "geometry_pair_utility_classes.H"

#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief This structure "converts" the type of a volume element. This is done in this way, so
   * the type conversion can be done at compile time.
   */
  template <DRT::Element::DiscretizationType discretization>
  struct ElementDiscretizationToVolumeType
  {
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex8>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex20>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex27>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::tet4>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::tet10>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  };


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
  };


  /**
   * \brief Discretization based on standard (in most cases Lagrangian) shape functions.
   * @param discretization Type of shape function.
   */
  template <DRT::Element::DiscretizationType discretization>
  class ElementDiscretizationStandard : public ElementDiscretizationBase<discretization, 1>
  {
   public:
    //! Volume type of the element. This member will be defined for hex8/20/27 tet4/10. If it is
    //! accessed for a different element type, there will be a compiler error.
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        ElementDiscretizationToVolumeType<discretization>::volume_type_;
  };


  /**
   * \brief Discretization based on hermitian shape functions.
   * @param discretization Type of shape function.
   */
  class ElementDiscretizationHermite
      : public ElementDiscretizationBase<DRT::Element::DiscretizationType::line2, 2>
  {
  };


  /**
   * Shortcuts to element types are created here, so the explicit template initialisations are
   * better readable.
   */
  typedef ElementDiscretizationHermite t_hermite;
  typedef ElementDiscretizationStandard<DRT::Element::hex8> t_hex8;
  typedef ElementDiscretizationStandard<DRT::Element::hex20> t_hex20;
  typedef ElementDiscretizationStandard<DRT::Element::hex27> t_hex27;
  typedef ElementDiscretizationStandard<DRT::Element::tet4> t_tet4;
  typedef ElementDiscretizationStandard<DRT::Element::tet10> t_tet10;

}  // namespace GEOMETRYPAIR

#endif
