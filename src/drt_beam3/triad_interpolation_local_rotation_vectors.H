/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief a triad interpolation scheme based on local rotation vectors


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef TRIAD_INTERPOLATION_LOCAL_ROTVECS
#define TRIAD_INTERPOLATION_LOCAL_ROTVECS

#include "triad_interpolation.H"

#include "../drt_lib/drt_element.H"

#include "../linalg/linalg_fixedsizematrix.H"

namespace LARGEROTATIONS
{
  /**
   * \brief a triad interpolation scheme based on local rotation vectors
   *        see Shoemake (1985) and Crisfield/Jelenic (1999) for formulae and details
   */

  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors : public LARGEROTATIONS::TriadInterpolation<T>
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}


    //! @name Constructors and destructors and related methods

    /** \brief Standard Constructor
     *
     *  \author grill
     *  \date 01/17 */
    TriadInterpolationLocalRotationVectors();

    /** \brief Destructor
     *
     *  \author grill
     *  \date 01/17 */
    virtual ~TriadInterpolationLocalRotationVectors(){};
    //@}

    //! @name Accessors

    /** \brief get node I which is part of the definition of the reference triad
     *
     *  \author grill
     *  \date 01/17 */
    inline unsigned int NodeI() const { return nodeI_; }

    /** \brief get node J which is part of the definition of the reference triad
     *
     *  \author grill
     *  \date 01/17 */
    inline unsigned int NodeJ() const { return nodeJ_; }

    //@}


    //! @name Derived methods

    /** \brief reset interpolation scheme with nodal quaternions
     *
     *  \author grill
     *  \date 01/17 */
    virtual void Reset(std::vector<LINALG::Matrix<4, 1, T>> const& nodal_quaternions);

    /** \brief reset interpolation scheme with nodal triads
     *
     *  \author grill
     *  \date 01/17 */
    virtual void Reset(std::vector<LINALG::Matrix<3, 3, T>> const& nodal_triads);

    /** \brief compute the interpolated triad at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void GetInterpolatedTriadAtXi(LINALG::Matrix<3, 3, T>& triad, const double xi) const;

    /** \brief compute the interpolated quaternion at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void GetInterpolatedQuaternionAtXi(
        LINALG::Matrix<4, 1, T>& quaternion, const double xi) const;

    //@}

    //! @name specific methods of this triad interpolation scheme (based on local rotation vectors)

    /** \brief compute the interpolated triad based on given local rotation vector
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedTriad(
        LINALG::Matrix<3, 3, T>& triad, const LINALG::Matrix<3, 1, T>& Psi_l) const;

    /** \brief compute the interpolated quaternion based on given local rotation vector
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedQuaternion(
        LINALG::Matrix<4, 1, T>& quaternion, const LINALG::Matrix<3, 1, T>& Psi_l) const;

    /** \brief compute the local rotation vector at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVectorAtXi(
        LINALG::Matrix<3, 1, T>& Psi_l, const double xi) const;

    /** \brief compute the local rotation vector based on given shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVector(
        LINALG::Matrix<3, 1, T>& Psi_l, const LINALG::Matrix<1, numnodes, double>& I_i) const;


    /** \brief compute the arc-length derivative of the local rotation vector at any point
     *         \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVectorDerivativeAtXi(
        LINALG::Matrix<3, 1, T>& Psi_l_s, const double jacobifac, const double xi) const;

    /** \brief compute the arc-length derivative of the local rotation vector based on given
     *         shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVectorDerivative(LINALG::Matrix<3, 1, T>& Psi_l_s,
        const LINALG::Matrix<1, numnodes, double>& I_i_xi, const double jacobifac) const;


    /** \brief compute the generalized rotational interpolation matrices for all nodes at
     *         any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatricesAtXi(
        std::vector<LINALG::Matrix<3, 3, T>>& Itilde, const double xi) const;

    /** \brief compute the generalized rotational interpolation matrices for all nodes
     *         based on given local rotation vector and shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatrices(
        std::vector<LINALG::Matrix<3, 3, T>>& Itilde, const LINALG::Matrix<3, 1, T>& Psi_l,
        const LINALG::Matrix<1, numnodes, double>& I_i) const;


    /** \brief compute the arc-length derivative of generalized rotational interpolation
     *         matrices for all nodes based on given local rotation vector and shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatricesDerivative(
        std::vector<LINALG::Matrix<3, 3, T>>& Itilde_prime, const LINALG::Matrix<3, 1, T>& Psi_l,
        const LINALG::Matrix<3, 1, T>& Psi_l_s, const LINALG::Matrix<1, numnodes, double>& I_i,
        const LINALG::Matrix<1, numnodes, double>& I_i_xi, const double jacobifac) const;

    /** \brief compute the arc-length derivative of generalized rotational interpolation
     *         matrices for all nodes based on given local rotation vector and shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatricesDerivative(
        std::vector<LINALG::Matrix<3, 3, T>>& Itilde_prime, const LINALG::Matrix<3, 1, T>& Psi_l,
        const LINALG::Matrix<3, 1, T>& Psi_l_s, const LINALG::Matrix<1, numnodes, double>& I_i,
        const LINALG::Matrix<1, numnodes, double>& I_i_s) const;

    //@}

   private:
    //! @name Private methods

    /** \brief set the two nodes I and J that are used to define the reference triad later on
     *
     *  \author grill
     *  \date 01/2017 */
    void SetNodeIandJ();

    //! get the interpolation scheme from the given number of nodes
    DRT::Element::DiscretizationType GetDisType() const;

    //! compute quaternion corresponding to reference triad Lambda_r according to (3.9), Jelenic
    //! 1999
    void CalcRefQuaternion(const LINALG::Matrix<4, 1, T>& Q_nodeI,
        const LINALG::Matrix<4, 1, T>& Q_nodeJ, LINALG::Matrix<4, 1, T>& Q_r) const;

    //! compute angle of relative rotation between node I and J according to (3.10), Jelenic 1999
    void CalcPhi_IJ(const LINALG::Matrix<4, 1, T>& Q_nodeI, const LINALG::Matrix<4, 1, T>& Q_nodeJ,
        LINALG::Matrix<3, 1, T>& Phi_IJ) const;

    //! compute nodal local rotations according to (3.8), Jelenic 1999
    void CalcPsi_li(const LINALG::Matrix<4, 1, T>& Q_i, const LINALG::Matrix<4, 1, T>& Q_r,
        LINALG::Matrix<3, 1, T>& Psi_li) const;

    //! compute interpolated local relative rotation \Psi^l according to (3.11), Jelenic 1999
    void Calc_Psi_l(const std::vector<LINALG::Matrix<3, 1, T>>& Psi_li,
        const LINALG::Matrix<1, numnodes, double>& func, LINALG::Matrix<3, 1, T>& Psi_l) const;

    //! compute derivative of interpolated local relative rotation \Psi^l with respect to reference
    //! arc-length parameter s according to (3.11), Jelenic 1999
    void Calc_Psi_l_s(const std::vector<LINALG::Matrix<3, 1, T>>& Psi_li,
        const LINALG::Matrix<1, numnodes, double>& deriv_xi, const double& jacobi,
        LINALG::Matrix<3, 1, T>& Psi_l_s) const;

    //! compute local triad \Lambda from Crisfield 1999, eq. (4.7)
    void Calc_Lambda(const LINALG::Matrix<3, 1, T>& Psi_l, const LINALG::Matrix<4, 1, T>& Q_r,
        LINALG::Matrix<3, 3, T>& Lambda) const;

    //! compute quaternion equivalent to local triad \Lambda from Crisfield 1999, eq. (4.7)
    void Calc_Qgauss(const LINALG::Matrix<3, 1, T>& Psi_l, const LINALG::Matrix<4, 1, T>& Q_r,
        LINALG::Matrix<4, 1, T>& Qgauss) const;

    //! compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certain Gauss
    //! point
    void computeItilde(const LINALG::Matrix<3, 1, T>& Psil,
        std::vector<LINALG::Matrix<3, 3, T>>& Itilde, const LINALG::Matrix<3, 1, T>& phiIJ,
        const LINALG::Matrix<3, 3, T>& Lambdar, const std::vector<LINALG::Matrix<3, 1, T>>& Psili,
        const LINALG::Matrix<1, numnodes, double>& funct) const;

    //! compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999 for all nodes i at a certain Gauss
    //! point
    void computeItildeprime(const LINALG::Matrix<3, 1, T>& Psil,
        const LINALG::Matrix<3, 1, T>& Psilprime, std::vector<LINALG::Matrix<3, 3, T>>& Itildeprime,
        const LINALG::Matrix<3, 1, T>& phiIJ, const LINALG::Matrix<3, 3, T>& Lambdar,
        const std::vector<LINALG::Matrix<3, 1, T>>& Psili,
        const LINALG::Matrix<1, numnodes, double>& funct,
        const LINALG::Matrix<1, numnodes, double>& deriv_s) const;

    //! compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999
    void Calc_vI(LINALG::Matrix<3, 3, T>& vI, const LINALG::Matrix<3, 1, T>& phiIJ) const;

    //! compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999
    void Calc_vJ(LINALG::Matrix<3, 3, T>& vJ, const LINALG::Matrix<3, 1, T>& phiIJ) const;

    //@}


   private:
    //! @name member variables

    //! node I for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
    unsigned int nodeI_;

    //! node J for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
    unsigned int nodeJ_;

    //! this determines the kind of shape functions which are to be applied
    DRT::Element::DiscretizationType distype_;

    //! nodal triads stored as quaternions
    std::vector<LINALG::Matrix<4, 1, T>> Qnode_;

    //! reference quaternion Q_r corresponding to reference triad Lambda_r
    LINALG::Matrix<4, 1, T> Q_r_;

    //! local rotation angles at nodes: angles between nodal triads and reference triad
    std::vector<LINALG::Matrix<3, 1, T>> Psi_li_;

    //@}
  };

}  // namespace LARGEROTATIONS

#endif
