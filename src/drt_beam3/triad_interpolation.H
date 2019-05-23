/*-----------------------------------------------------------------------------------------------*/
/*!

\brief base class for a triad interpolation scheme

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef TRIAD_INTERPOLATION
#define TRIAD_INTERPOLATION

#include "../linalg/linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>

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
    virtual ~TriadInterpolation(){};

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
    virtual void Reset(std::vector<LINALG::TMatrix<T, 4, 1>> const& nodal_quaternions) = 0;

    /** \brief reset interpolation scheme with nodal triads
     *
     *  \author grill
     *  \date 01/17 */
    virtual void Reset(std::vector<LINALG::TMatrix<T, 3, 3>> const& nodal_triads) = 0;

    /** \brief compute the interpolated triad at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void GetInterpolatedTriadAtXi(
        LINALG::TMatrix<T, 3, 3>& triad, const double xi) const = 0;

    /** \brief compute the interpolated quaternion at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void GetInterpolatedQuaternionAtXi(
        LINALG::TMatrix<T, 4, 1>& quaternion, const double xi) const = 0;
    //@}

   private:
    //! @name Private evaluation methods

    //@}

   private:
    //! @name member variables

    //@}
  };

}  // namespace LARGEROTATIONS

#endif
