/*! \file
\brief Interface for every material that can evaluate solid material laws

\level 3

*/

#ifndef MAT_TRAIT_SOLID_H
#define MAT_TRAIT_SOLID_H
namespace MAT
{
  namespace TRAIT
  {
    /*!
     * This interface should define all core functionality for solid materials. Currently this is
     * done in So3Material.
     */
    class Solid
    {
     public:
      virtual ~Solid() = default;
    };
  }  // namespace TRAIT
}  // namespace MAT

#endif
