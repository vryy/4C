/*----------------------------------------------------------------------*/
/*!
\file condif2_impl.H

\brief Internal implementation of Condif2 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID2
#ifdef CCADISCRET

#ifndef CONDIF2_IMPL_H
#define CONDIF2_IMPL_H


#include "condif2.H"

namespace DRT
{
namespace ELEMENTS
{

  /// Internal Condif2 implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Condif2 element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Condif2)
    from the mathematical contents (this class). Of course there are
    different implementations of the Condif3 element, this is just one
    such implementation.

    The Condif2 will allocate exactly one object of this class for all
    condif elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>History</h3>

    The implementation here is the standard convection-diffusion element
    capable of dealing with systems of transported scalars.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version.

    \author vg
    \date 08/08
  */
  class Condif2Impl
  {
  public:
    /// Constructor with number of nodes and number of degree of freedom
    explicit Condif2Impl(int iel, int numdofpernode, int numscal);

    /// Internal implementation class for specific scalar transport element
    static Condif2Impl* Impl(DRT::ELEMENTS::Condif2* c2);

    //! Calculate matrix and rhs. Here the whole thing is hidden.
    void Sysmat(
        const DRT::ELEMENTS::Condif2*   ele, ///< the element those matrix is calculated
        const vector<double>&           ehist, ///< rhs from beginning of time step
        const vector<double>&           edens, ///< density
        Epetra_SerialDenseMatrix*       sys_mat,///< element matrix to calculate
        Epetra_SerialDenseVector*       residual, ///< element rhs to calculate
        Epetra_SerialDenseVector&       subgrdiff, ///< subgrid-diff.-scaling vector
        const struct _MATERIAL*         material, ///< material pointer
        const double                    time, ///< current simulation time
        const double                    timefac, ///< time discretization factor
        const Epetra_SerialDenseVector& evel, ///< nodal velocities at n+1
        bool                            temperature, ///< temperature flag
        string                          fssgd, ///< subgrid-diff. flag
        const bool                      is_stationary ///< flag indicating stationary formulation
    );

    //! calculate mass matrix and rhs for initializing One-Step-Theta time integration
    void InitializeOST(
        const DRT::ELEMENTS::Condif2*   ele,
        const vector<double>&           ephi0,
        const vector<double>&           edens,
        Epetra_SerialDenseMatrix&       massmat,
        Epetra_SerialDenseVector&       rhs,
        Epetra_SerialDenseVector&       subgrdiff,
        const struct _MATERIAL*         material,
        const double                    time,
        const double                    timefac,
        const Epetra_SerialDenseVector& evel,
        bool                            temperature,
        string                          fssgd
    );

    //! calculate normalized subgrid-diffusivity matrix
    void CalcSubgridDiffMatrix(
        const DRT::ELEMENTS::Condif2*   ele,
        Epetra_SerialDenseMatrix&       sys_mat_sd,
        const double                    timefac,
        const bool                      is_stationary
    );


  private:

    //! get the body force
    void BodyForce(
        const DRT::ELEMENTS::Condif2* ele, ///< the element we are dealing with
        const double time                  ///< current simulation time
    );

    //! calculate stabilization parameter
    void CalTau(
        const DRT::ELEMENTS::Condif2*&          ele,
        Epetra_SerialDenseVector&               subgrdiff,
        const Epetra_SerialDenseVector&         evel,
        const DRT::Element::DiscretizationType& distype,
        const double&                           timefac,
        string                                  fssgd,
        const bool&                             is_stationary,
        const bool                              initial
    );

    //! evaluate shape functions and their derivatives at current integration point
    void EvalShapeFuncAndDerivsAtIntPoint(
        const DRT::UTILS::IntegrationPoints2D& intpoints, ///< integration points
        const int& iquad,                                 ///< id of current Gauss point
        const DRT::Element::DiscretizationType& distype,  ///< distinguish between DiscretizationType
        const bool& higher_order_ele,                     ///< are second derivatives needed?
        const DRT::ELEMENTS::Condif2*   ele               ///< the element
    );

    //! calculate 2nd derivatives of shape functions (local and global)
    void CalSecondDeriv(
        const double&                            e1,      ///< first coordinate of GP
        const double&                            e2,      ///< second coordinate of GP
        const DRT::Element::DiscretizationType&  distype  ///< distinguish between DiscretizationType
    );

    //! calculate element matrix and rhs vector (One-Step-Theta, BDF2)
    void CalMat(
        Epetra_SerialDenseMatrix& estif,
        Epetra_SerialDenseVector& eforce,
        const bool                higher_order_ele,
        const double&             timefac,
        const int&                dofindex
    );

    //! calculate element matrix and rhs vector (Stationary problem formulation)
    void CalMatStationary(
        Epetra_SerialDenseMatrix& estif,
        Epetra_SerialDenseVector& eforce,
        const bool                higher_order_ele,
        const int&                dofindex
    );

    //! number of nodes
    const int iel_;
    //! number of dof per node
    const int numdofpernode_;
    //! number of transported scalars (numscal_ <= numdofpernode_)
    const int numscal_;
    //! node coordinates
    Epetra_SerialDenseMatrix xyze_;
    //! bodyforce in element nodes
    Epetra_SerialDenseVector bodyforce_;
    //! diffusivity / diffusivities (in case of systems) / thermal conductivity
    vector<double> diffus_;
    //! charge numbers of all ion species (ELCH specific)
    vector<double> valence_;
    //! specific heat capacity at constant pressure (in case of temperature eq.)
    double shcacp_;
    //! array for shape functions
    Epetra_SerialDenseVector funct_;
    //! array for density-weighted shape functions
    Epetra_SerialDenseVector densfunct_;
    //! array for shape function derivatives w.r.t r,s,t
    Epetra_SerialDenseMatrix deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    Epetra_SerialDenseMatrix deriv2_;
    //! transposed jacobian "dx/ds"
    Epetra_SerialDenseMatrix xjm_;
    //! inverse of transposed jacobian "ds/dx"
    Epetra_SerialDenseMatrix xij_;
    //! global derivatives of shape functions w.r.t x,y,z
    Epetra_SerialDenseMatrix derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    Epetra_SerialDenseMatrix derxy2_;
    //! rhs force in Gauss point
    vector<double> rhs_;
    //! vector containing all values from previous timelevel n in gausspoint
    vector<double> hist_;
    //! velocity vector in gausspoint
    vector<double> velint_;
    //! the stabilisation parameters (one per transported scalar)
    vector<double> tau_;
    //! artificial diffusivities (one per transported scalar)
    vector<double> kart_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    Epetra_SerialDenseMatrix xder2_;
    //! integration factor for current GP: fac = GaussWeight * det(J)
    double fac_;
    //! convective part used in stabilization terms
    vector<double> conv_;
    //! diffusive part used in stabilization terms
    vector<double> diff_;
  };

}
}

#endif

#endif
#endif
