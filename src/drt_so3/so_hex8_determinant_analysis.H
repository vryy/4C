/*-----------------------------------------------------------*/
/*!
\file so_hex8_determinant_analysis.H

\maintainer Michael Hiermeier

\date Sep 5, 2018

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_SO3_SO_HEX8_DETERMINANT_ANALYSIS_H_
#define SRC_DRT_SO3_SO_HEX8_DETERMINANT_ANALYSIS_H_

#include "so_hex8.H"
#include <functional>

//#define DEBUG_SO_HEX8_DET_ANALYSIS
namespace DRT
{
  namespace ELEMENTS
  {
    /** \brief Analyse the Jacobian determiant of Hex8 element
     *
     *  Based on the work
     *  [1] Johnen, et al, "Robust and efficient validation of the linear
     *  hexahedral element", 2017
     *
     *  \author hiermeier \date 09/18 */
    class So_Hex8_Determinant_Analysis
    {
      struct BezierCube;

     public:
      /** \brief Create a new instance of this class
       *
       *  If it is called the very first time, the related static members are
       *  set-up as well.
       *
       *  \return RCP to the created class object. */
      static Teuchos::RCP<So_Hex8_Determinant_Analysis> create();

      /// destructor
      ~So_Hex8_Determinant_Analysis() = default;

      /** \brief Test the validity of the current element
       *
       *  \param x_curr  Current nodal coordinates
       *  \param rc      recursion counter pointer (optional)
       *
       *  \return TRUE if the element is valid, otherwise return FALSE. */
      bool isValid(
          const LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr, unsigned* rc = NULL) const;

     private:
      /** \brief constructor (private)
       *
       *  Use the create function instead.
       *
       *  \author hiermeier \date 09/18 */
      So_Hex8_Determinant_Analysis() = default;

      /** \brief Fill the static matrix \c map_q_ with the values given in
       *  reference [1]
       *
       *  The matrix is later used to map the 20 Jacobian determinant values given
       *  by the 20 computed TET4 volumes (multiplied by 6) to the 27 Bézier
       *  coefficients.
       *
       *  \author hiermeier \date 09/18 */
      static void buildMapLagrange20ToBezier27();

      /** \brief Build the full map from the 2nd order Lagrange coefficients
       *  to the 2nd-order Bézier coefficients
       *
       *  This is a one-time cost and is only done once during the very first
       *  create call. The computation makes it necessary to invert a 27x27
       *  matrix. This is achieved by the suitable LAPACK routine. The result is
       *  identical to the matrix T in reference [1].
       *
       *  \author hiermeier \date 09/18 */
      static void buildMapLagrangeToBezier();

      /** \brief Build the map from 27 Bézier coefficients to 27 Lagrange
       *  coefficients (both for 2nd order polynomials)
       *
       *  This can be performed cheaply, since only the evaluation of the
       *  2nd order Bezier basis functions at the different points is necessary.
       *
       *  \author hiermeier \date 09/18 */
      static void buildMapBezierToLagrange(LINALG::Matrix<27, 27>& map_b2l);

      /** \brief Alternative call which allows an additional shift and scale of
       *  the parametric coordinates
       *
       *  This routine is important for the potentially necessary subdivision.
       *
       *  \param map_b2l  Map from Bézier to Lagrange coefficients
       *  \param scale    Array with three scaling values for xi, eta and zeta,
       *                  respectively.
       *  \param shift    Array with three shifting values for xi, eta and zeta,
       *                  respectively.
       *
       *  \author hiermeier \date 09/18 */
      static void buildMapBezierToLagrange(
          LINALG::Matrix<27, 27>& map_b2l, const double* scale, const double* shift);

      /** \brief Build Sub-map matrix from Bézier to Lagrange coefficients
       *
       *  \param left   Array containing the three left border points of the
       *                sub-cube
       *  \param right  Array containing the three right border points of the
       *                sub-cube
       *  \param sub_map_b2l  resulting sub-map matrix
       *
       *  \author hiermeier \date 09/18 */
      void buildSubMapBezierToLagrange(
          const double* left, const double* right, LINALG::Matrix<27, 27>& sub_map_b2l) const;

      /** \brief Compute the Bézier coefficients of the sub-domain
       *
       *  \param bcoeffs      Bézier coefficients of the original domain
       *  \param subcube      Current sub-domain
       *  \param sub_bcoeffs  Bézier coefficients of the current sub-domain
       *
       *  \author hiermeier \date 09/18 */
      void getBezierCoeffsOfSubdomain(const LINALG::Matrix<27, 1>& bcoeffs, BezierCube& subcube,
          LINALG::Matrix<27, 1>& sub_bcoeffs) const;


      /** \brief Compute the sub-cube border points based on the parent cube
       *   borders
       *
       *   \param l   Array containing the three left border points of the parent
       *              cube
       *   \param r   Array containing the three right border points of the parent
       *              cube
       *   \param subcubes  List containing all eight sub-cubes. This method fills
       *              the left_ and right_ variables of these sub-cubes
       *              corresponding to arrays containing the 3 left and right
       *              border coordinates, respectively.
       *
       *   \author hiermeier \date 09/18 */
      void getSubCubeBorders(
          const double* l, const double* r, std::list<BezierCube>& subcubes) const;

      /** \brief Returns TRUE if the tested entries contain one which is invalid
       *
       *  Invalid means in this context that the value is zero or negative.
       *
       *  \param entries  ptr to the first entry of the array
       *  \param length   number of entries in the array
       *
       *  \author hiermeier \date 09/18 */
      bool hasInvalidEntry(const double* entries, const unsigned length) const;

      /** \brief Perform a recursive subdivision of the cubes and refine
       *  the estimates for the Jacobian determinant
       *
       *  \param bcoeffs  Bézier coefficients of the initial Hex8 element
       *  \param left     Array containing the three left border coordinates
       *                  of the parent (sub-)cube
       *  \param right    Array containing the three left border coordinates
       *                  of the parent (sub-)cube
       *
       *  \author hiermeier \date 09/18 */
      bool recursiveSubdivision(const LINALG::Matrix<27, 1>& bcoeffs, const double* left,
          const double* right, unsigned& rcount) const;


      /** \brief Compute the 20 independent scaled TET4 volumes
       *
       *  The computed volumes are multiplied by a factor of 6 such that the
       *  first 8 (related to the corners of the Hex8 element) correspond to
       *  Jacobian determinant values. The subsequent 12 TET4 volumes correspond
       *  to the center points of the edges. See reference [1] for more
       *  information.
       *
       *  \param tet4_volumes  result of the calculation.
       *  \param x_curr        matrix containing the current coordinates of the
       *                       HEX8 corners
       *
       *  \author hiermeier \date 09/18 */
      void compute20Tet4Volumes(LINALG::Matrix<20, 1>& tet4_volumes,
          const LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr) const;

      /** \brief Compute scaled TET4 volume at the corners of the HEX8
       *
       *  \param tet4_volumes  Contains the results
       *  \param x_curr        Current nodal coordinates of the HEX8
       *  \param f_index0      Function ptr to compute the node id for the
       *                       first corner point of the TET4
       *  \param f_index1      Function ptr to compute the node id for the
       *                       second corner point of the TET4
       *  \param f_index2      Function ptr to compute the node id for the
       *                       third corner point of the TET4
       *  \param f_index3      Function ptr to compute the node id for the
       *                       fourth corner point of the TET4
       *  \param offset        offset for the index of the result array
       *
       *  \return new offset value (increased by the four computed TET4 volumes)
       *
       *  \author hiermeier \date 09/18 */
      unsigned computeTet4VolAtCorners(LINALG::Matrix<20, 1>& tet4_volumes,
          const LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr,
          std::function<unsigned(unsigned i)> f_index0,
          std::function<unsigned(unsigned i)> f_index1,
          std::function<unsigned(unsigned i)> f_index2,
          std::function<unsigned(unsigned i)> f_index3, unsigned offset) const;

      /** \brief Compute scaled TET4 volume at the edges of the HEX8
       *
       *  \param tet4_volumes  Contains the results
       *  \param x_curr        Current nodal coordinates of the HEX8
       *  \param f_index0      Function ptr to compute the node id for the
       *                       first corner point of the TET4
       *  \param f_index1      Function ptr to compute the node id for the
       *                       second corner point of the TET4
       *  \param f_index20     Function ptr to compute the node id for the
       *                       left corner point of the third TET4 corner
       *  \param f_index21     Function ptr to compute the node id for the
       *                       right corner point of the third TET4 corner
       *  \param f_index30     Function ptr to compute the node id for the
       *                       left corner of the fourth TET4 corner
       *  \param f_index31     Function ptr to compute the node id for the
       *                       right corner of the fourth TET4 corner
       *  \param offset        offset for the index of the result array
       *
       *  \return new offset value (increased by the four computed TET4 volumes)
       *
       *  \author hiermeier \date 09/18 */
      unsigned computeTet4VolAtEdges(LINALG::Matrix<20, 1>& tet4_volumes,
          const LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr,
          std::function<unsigned(unsigned i)> f_index0,
          std::function<unsigned(unsigned i)> f_index1,
          std::function<unsigned(unsigned i)> f_index20,
          std::function<unsigned(unsigned i)> f_index21,
          std::function<unsigned(unsigned i)> f_index30,
          std::function<unsigned(unsigned i)> f_index31, unsigned offset) const;

      double computeTet4Volume(const LINALG::Matrix<NUMDIM_SOH8, 4>& tet4_ncoords) const;

      /** \brief Evaluate the second order bezier function
       *
       *  \param t  parameteric coordinate \f$ t \in (0,1) \f$
       *  \param n  desired bezier polynomial. n must be between 0 and 2.
       *
       *  \author hiermeier \date 09/18 */
      static double bezierFunc2(const double t, unsigned n)
      {
        switch (n)
        {
          case 0:
            return (1.0 - t) * (1.0 - t);
          case 1:
            return 2.0 * t * (1.0 - t);
          case 2:
            return t * t;
          default:
          {
            dserror("The desired bezier function #%d is not defined.", n);
            exit(EXIT_FAILURE);
          }
        }
      };

      /** \brief struct containing important necessary information of a Bézier
       *  (sub-)cube */
      struct BezierCube
      {
        /// initialize members
        void init()
        {
          std::fill(left_, left_ + 3, 0.0);
          std::fill(right_, right_ + 3, 1.0);
        }

        /// print content to output stream
        void print(std::ostream& os) const;

        /// array containing the left border coordinates
        double left_[3];

        /// array containing the right border coordinates
        double right_[3];
      };

     private:
      /** \brief Bézier point coordinates following the scheme in reference [1]
       *
       *  The first 8 correspond to the eight corners.
       *  The following 12 to the center points of the 12 edges of the Hex8.
       *  The next 6 to the center points of the faces.
       *  The last one corresponds to the spatial center of the Hex8.
       *
       *  The parameter space goes from 0.0 to 1.0, since this is numerically
       *  more stable.
       *
       *  \author hiermeier \date 09/18 */
      static const double bezier_points_[27][3];

      /// Bézier function indices corresponding to the Bézier point coordinates
      static const unsigned bezier_indices_[27][3];

      /// TRUE if all static members have been successfully set-up.
      static bool issetup_;

      /// mapping from 20 Jacobian entries to 27 Bézier coefficients
      static LINALG::Matrix<27, 20> map_q_;

      // mapping from the Lagrange to the Bézier space
      static LINALG::Matrix<27, 27> map_l2b_;
    };

    /** \brief Fast modulus routine
     *
     *  This routine is supposed to be remarkably faster than the default
     *  modulus operator. According to https://youtu.be/nXaxk27zwlk?t=56m34s
     *
     *  \author hiermeier \date 09/18 */
    inline unsigned fast_mod(const unsigned input, const unsigned ceil)
    {
      /* apply the modulo operator only when needed
       * (i.e. when the input is greater than the ceiling) */
      return input >= ceil ? input % ceil : input;
    }
  }  // namespace ELEMENTS
}  // namespace DRT



#endif /* SRC_DRT_SO3_SO_HEX8_DETERMINANT_ANALYSIS_H_ */
