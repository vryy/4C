/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief a triad interpolation scheme based on local rotation vectors

\maintainer Maximilian Grill

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
    virtual void Reset(std::vector<LINALG::TMatrix<T, 4, 1>> const& nodal_quaternions);

    /** \brief reset interpolation scheme with nodal triads
     *
     *  \author grill
     *  \date 01/17 */
    virtual void Reset(std::vector<LINALG::TMatrix<T, 3, 3>> const& nodal_triads);

    /** \brief compute the interpolated triad at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void GetInterpolatedTriadAtXi(LINALG::TMatrix<T, 3, 3>& triad, const double xi) const;

    /** \brief compute the interpolated quaternion at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    virtual void GetInterpolatedQuaternionAtXi(
        LINALG::TMatrix<T, 4, 1>& quaternion, const double xi) const;

    //@}

    //! @name specific methods of this triad interpolation scheme (based on local rotation vectors)

    /** \brief compute the interpolated triad based on given local rotation vector
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedTriad(
        LINALG::TMatrix<T, 3, 3>& triad, const LINALG::TMatrix<T, 3, 1>& Psi_l) const;

    /** \brief compute the interpolated quaternion based on given local rotation vector
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedQuaternion(
        LINALG::TMatrix<T, 4, 1>& quaternion, const LINALG::TMatrix<T, 3, 1>& Psi_l) const;

    /** \brief compute the local rotation vector at any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVectorAtXi(
        LINALG::TMatrix<T, 3, 1>& Psi_l, const double xi) const;

    /** \brief compute the local rotation vector based on given shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVector(
        LINALG::TMatrix<T, 3, 1>& Psi_l, const LINALG::TMatrix<double, 1, numnodes>& I_i) const;


    /** \brief compute the arc-length derivative of the local rotation vector at any point
     *         \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVectorDerivativeAtXi(
        LINALG::TMatrix<T, 3, 1>& Psi_l_s, const double jacobifac, const double xi) const;

    /** \brief compute the arc-length derivative of the local rotation vector based on given
     *         shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetInterpolatedLocalRotationVectorDerivative(LINALG::TMatrix<T, 3, 1>& Psi_l_s,
        const LINALG::TMatrix<double, 1, numnodes>& I_i_xi, const double jacobifac) const;


    /** \brief compute the generalized rotational interpolation matrices for all nodes at
     *         any point \xi \in [-1,1] in parameter space
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatricesAtXi(
        std::vector<LINALG::TMatrix<T, 3, 3>>& Itilde, const double xi) const;

    /** \brief compute the generalized rotational interpolation matrices for all nodes
     *         based on given local rotation vector and shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatrices(
        std::vector<LINALG::TMatrix<T, 3, 3>>& Itilde, const LINALG::TMatrix<T, 3, 1>& Psi_l,
        const LINALG::TMatrix<double, 1, numnodes>& I_i) const;


    /** \brief compute the arc-length derivative of generalized rotational interpolation
     *         matrices for all nodes based on given local rotation vector and shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatricesDerivative(
        std::vector<LINALG::TMatrix<T, 3, 3>>& Itilde_prime, const LINALG::TMatrix<T, 3, 1>& Psi_l,
        const LINALG::TMatrix<T, 3, 1>& Psi_l_s, const LINALG::TMatrix<double, 1, numnodes>& I_i,
        const LINALG::TMatrix<double, 1, numnodes>& I_i_xi, const double jacobifac) const;

    /** \brief compute the arc-length derivative of generalized rotational interpolation
     *         matrices for all nodes based on given local rotation vector and shape function values
     *
     *  \author grill
     *  \date 01/17 */
    void GetNodalGeneralizedRotationInterpolationMatricesDerivative(
        std::vector<LINALG::TMatrix<T, 3, 3>>& Itilde_prime, const LINALG::TMatrix<T, 3, 1>& Psi_l,
        const LINALG::TMatrix<T, 3, 1>& Psi_l_s, const LINALG::TMatrix<double, 1, numnodes>& I_i,
        const LINALG::TMatrix<double, 1, numnodes>& I_i_s) const;

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
    void CalcRefQuaternion(const LINALG::TMatrix<T, 4, 1>& Q_nodeI,
        const LINALG::TMatrix<T, 4, 1>& Q_nodeJ, LINALG::TMatrix<T, 4, 1>& Q_r) const;

    //! compute angle of relative rotation between node I and J according to (3.10), Jelenic 1999
    void CalcPhi_IJ(const LINALG::TMatrix<T, 4, 1>& Q_nodeI,
        const LINALG::TMatrix<T, 4, 1>& Q_nodeJ, LINALG::TMatrix<T, 3, 1>& Phi_IJ) const;

    //! compute nodal local rotations according to (3.8), Jelenic 1999
    void CalcPsi_li(const LINALG::TMatrix<T, 4, 1>& Q_i, const LINALG::TMatrix<T, 4, 1>& Q_r,
        LINALG::TMatrix<T, 3, 1>& Psi_li) const;

    //! compute interpolated local relative rotation \Psi^l according to (3.11), Jelenic 1999
    void Calc_Psi_l(const std::vector<LINALG::TMatrix<T, 3, 1>>& Psi_li,
        const LINALG::TMatrix<double, 1, numnodes>& func, LINALG::TMatrix<T, 3, 1>& Psi_l) const;

    //! compute derivative of interpolated local relative rotation \Psi^l with respect to reference
    //! arc-length parameter s according to (3.11), Jelenic 1999
    void Calc_Psi_l_s(const std::vector<LINALG::TMatrix<T, 3, 1>>& Psi_li,
        const LINALG::TMatrix<double, 1, numnodes>& deriv_xi, const double& jacobi,
        LINALG::TMatrix<T, 3, 1>& Psi_l_s) const;

    //! compute local triad \Lambda from Crisfield 1999, eq. (4.7)
    void Calc_Lambda(const LINALG::TMatrix<T, 3, 1>& Psi_l, const LINALG::TMatrix<T, 4, 1>& Q_r,
        LINALG::TMatrix<T, 3, 3>& Lambda) const;

    //! compute quaternion equivalent to local triad \Lambda from Crisfield 1999, eq. (4.7)
    void Calc_Qgauss(const LINALG::TMatrix<T, 3, 1>& Psi_l, const LINALG::TMatrix<T, 4, 1>& Q_r,
        LINALG::TMatrix<T, 4, 1>& Qgauss) const;

    //! compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certain Gauss
    //! point
    void computeItilde(const LINALG::TMatrix<T, 3, 1>& Psil,
        std::vector<LINALG::TMatrix<T, 3, 3>>& Itilde, const LINALG::TMatrix<T, 3, 1>& phiIJ,
        const LINALG::TMatrix<T, 3, 3>& Lambdar, const std::vector<LINALG::TMatrix<T, 3, 1>>& Psili,
        const LINALG::TMatrix<double, 1, numnodes>& funct) const;

    //! compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999 for all nodes i at a certain Gauss
    //! point
    void computeItildeprime(const LINALG::TMatrix<T, 3, 1>& Psil,
        const LINALG::TMatrix<T, 3, 1>& Psilprime,
        std::vector<LINALG::TMatrix<T, 3, 3>>& Itildeprime, const LINALG::TMatrix<T, 3, 1>& phiIJ,
        const LINALG::TMatrix<T, 3, 3>& Lambdar, const std::vector<LINALG::TMatrix<T, 3, 1>>& Psili,
        const LINALG::TMatrix<double, 1, numnodes>& funct,
        const LINALG::TMatrix<double, 1, numnodes>& deriv_s) const;

    //! compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999
    void Calc_vI(LINALG::TMatrix<T, 3, 3>& vI, const LINALG::TMatrix<T, 3, 1>& phiIJ) const;

    //! compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999
    void Calc_vJ(LINALG::TMatrix<T, 3, 3>& vJ, const LINALG::TMatrix<T, 3, 1>& phiIJ) const;

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
    std::vector<LINALG::TMatrix<T, 4, 1>> Qnode_;

    //! reference quaternion Q_r corresponding to reference triad Lambda_r
    LINALG::TMatrix<T, 4, 1> Q_r_;

    //! local rotation angles at nodes: angles between nodal triads and reference triad
    std::vector<LINALG::TMatrix<T, 3, 1>> Psi_li_;

    //@}
  };

}  // namespace LARGEROTATIONS

#endif
