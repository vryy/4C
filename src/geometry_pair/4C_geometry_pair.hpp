/*----------------------------------------------------------------------*/
/*! \file

\brief base class for two or more elements that interact geometrically, e.g. contact or meshtying
problems.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_HPP
#define FOUR_C_GEOMETRY_PAIR_HPP


#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace CORE::Elements
{
  class Element;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Base class for all geometry pairs.
   */
  class GeometryPair
  {
   public:
    /**
     * \brief Constructor.
     * @param element1 Pointer to the first element.
     * @param element2 Pointer to the second element.
     */
    GeometryPair(const CORE::Elements::Element* element1, const CORE::Elements::Element* element2)
        : element1_(element1), element2_(element2){};

    /**
     * \brief Destructor.
     */
    virtual ~GeometryPair() = default;

    /**
     * \brief Pointer to the first element.
     */
    inline const CORE::Elements::Element* Element1() const { return element1_; };

    /**
     * \brief Pointer to the second element.
     */
    inline const CORE::Elements::Element* Element2() const { return element2_; };

    /**
     * \brief Set the pointer to the second element.
     *
     * For surface elements the pairs need the face element representing the surface, not the volume
     * element.
     */
    virtual inline void SetElement2(const CORE::Elements::Element* element2)
    {
      element2_ = element2;
    };

   private:
    //! Link to the element containing the first geometry.
    const CORE::Elements::Element* element1_;

    //! Link to the element containing the second geometry.
    const CORE::Elements::Element* element2_;
  };
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
