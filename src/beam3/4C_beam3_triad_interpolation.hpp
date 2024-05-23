/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief base class for a triad interpolation scheme


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAM3_TRIAD_INTERPOLATION_HPP
#define FOUR_C_BEAM3_TRIAD_INTERPOLATION_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace LARGEROTATIONS
{
  /**
   * \brief abstract base class for a triad interpolation scheme
   */
  template <typename T>
  class TriadInterpolation
  {
   public:
    //! @name Constructors and destructors and related methods

    /** \brief Standard Constructor
     *
     *  \author grill
     *  \date 01/17 */
    TriadInterpolation();

    /** \brief Destructor
     *
     *  \author grill
     *  \date 01/17 */
    virtual ~TriadInterpolation() = default;

    /** \brief return appropriate derived (templated) class (acts as a simple factory)
     *
     *  \author grill
     *  \date 01/17 */
    static Teuchos::RCP<TriadInterpolation<T>> Create(unsigned int numnodes);
    //@}


    //! @name Public evaluation methods

    /** \brief reset interpolation scheme with nodal quaternions
     *
     *  \author grill
     *  \date 01/17 */
    virtual void Reset(std::vector<CORE::LINALG::Matrix<4, 1, T>> const& nodal_quaternions) = 0;

    /** \brief reset interpolation scheme with nodal triads
     *
     *  \author grill
     *  \date 01/17 */
    virtual void Reset(std::vector<CORE::LINALG::Matrix<3, 3, T>> const& nodal_triads) = 0;

    /** \brief compute the interpolated triad at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void get_interpolated_triad_at_xi(
        CORE::LINALG::Matrix<3, 3, T>& triad, const double xi) const = 0;

    /** \brief compute the interpolated quaternion at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void get_interpolated_quaternion_at_xi(
        CORE::LINALG::Matrix<4, 1, T>& quaternion, const double xi) const = 0;
    //@}

   private:
    //! @name Private evaluation methods

    //@}

   private:
    //! @name member variables

    //@}
  };

}  // namespace LARGEROTATIONS

FOUR_C_NAMESPACE_CLOSE

#endif
