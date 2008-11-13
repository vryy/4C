/*----------------------------------------------------------------------*/
/*!
\file condif3_impl.H

\brief Internal implementation of Condif3 element

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef CONDIF3_IMPL_H
#define CONDIF3_IMPL_H


#include "condif3.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

namespace DRT
{
namespace ELEMENTS
{

/// Interface base class for Condif3Impl
/*!
  This class exists to provide a common interface for all template
  versions of Condif3Impl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of Condif3Impl.
 */
class Condif3ImplInterface
{

public:
  /// Empty constructor
  Condif3ImplInterface() {}
  /// Empty destructor
  virtual ~Condif3ImplInterface() {}
  /// Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in Condif3Impl.
   */
  virtual int Evaluate(Condif3*                   ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra,
                       RefCountPtr<MAT::Material> mat,
                       _MATERIAL*                 actmat
                       ) = 0;


  /// Internal implementation class for condif element
  static Condif3ImplInterface* Impl(DRT::ELEMENTS::Condif3* c3);

};

  /// Internal Condif3 implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Condif3 element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Condif3)
    from the mathematical contents (this class). Of course there are
    different implementations of the Condif3 element, this is just one
    such implementation.

    The Condif3 will allocate exactly one object of this class for all
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

    \author gjb
    \date 08/08
  */
template<DRT::Element::DiscretizationType distype>
class Condif3Impl: public Condif3ImplInterface
{
public:
  /// Constructor
  Condif3Impl(int numdofpernode, int numscal);

  //! number of nodes
  static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! Evaluate
  virtual int Evaluate(Condif3*                   ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra,
                       RefCountPtr<MAT::Material> mat,
                       _MATERIAL*                 actmat);


  private:

    //! Calculate matrix and rhs. Here the whole thing is hidden.
    virtual void Sysmat(
        const DRT::ELEMENTS::Condif3*                             ele, ///< the element those matrix is calculated
        const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ephinp,///< current scalar field
        const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ehist, ///< rhs from beginning of time step
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>&          edensnp, ///< current density field
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>&          epotnp, ///< el. potential at element nodes
        Epetra_SerialDenseMatrix&                                 sys_mat,///< element matrix to calculate
        Epetra_SerialDenseVector&                                 residual, ///< element rhs to calculate
        Epetra_SerialDenseVector&                                 subgrdiff, ///< subgrid-diff.-scaling vector
        const struct _MATERIAL*                                   material, ///< material pointer
        const double                                              time, ///< current simulation time
        const double                                              dt, ///< current time-step length
        const double                                              timefac, ///< time discretization factor
        const LINALG::FixedSizeSerialDenseMatrix<3,iel>&          evelnp,///< nodal velocities at t_{n+1}
        const bool                                                temperature, ///< temperature flag
        const enum Condif3::TauType                               whichtau, ///< flag for stabilization parameter definition
        const string                                              fssgd, ///< subgrid-diff. flag
        const bool                                                is_stationary, ///< flag indicating stationary formulation
        const double                                              frt ///< factor F/RT needed for ELCH calculations
    );

    //! calculate mass matrix and rhs for initializing One-Step-Theta time integration
    virtual void InitializeOST(
        const DRT::ELEMENTS::Condif3*                             ele,  ///< the element
        const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ephi0, ///< initial scalar field
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>&          edens0, ///< initial density field
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>&          epot0,  ///< el. potential at element nodes
        Epetra_SerialDenseMatrix&                                 massmat, ///< mass matrix to calculate
        Epetra_SerialDenseVector&                                 rhs, ///< element rhs to calculate
        Epetra_SerialDenseVector&                                 subgrdiff, ///< subgrid-diff.-scaling vector
        const struct _MATERIAL*                                   material, ///< material pointer
        const double                                              time, ///< current simulation time
        const double                                              dt, ///< current time-step length
        const double                                              timefac, ///< time discretization factor
        const LINALG::FixedSizeSerialDenseMatrix<3,iel>&          evel0, ///< initial nodal velocities
        const bool                                                temperature, ///< temperature flag
        const enum Condif3::TauType                               whichtau, ///< flag for stabilization parameter definition
        const string                                              fssgd, ///< subgrid-diff. flag
        const double                                              frt ///< factor F/RT needed for ELCH calculations
    );

    //! calculate normalized subgrid-diffusivity matrix
    virtual void CalcSubgridDiffMatrix(
        const DRT::ELEMENTS::Condif3* ele,
        Epetra_SerialDenseMatrix&     sys_mat_sd,
        const double                  timefac,
        const bool                    is_stationary
    );

    //! get the body force
    virtual void BodyForce(
        const DRT::ELEMENTS::Condif3* ele, ///< the element we are dealing with
        const double                  time ///< current simulation time
    );

    //! get the material parameters
    virtual void GetMaterialParams(
        const struct _MATERIAL* material,    ///< material pointer
        const bool&             temperature  ///< temperature flag
    );

    //! calculate stabilization parameter
    virtual void CalTau(
        const DRT::ELEMENTS::Condif3*&                   ele, ///< the element
        Epetra_SerialDenseVector&                        subgrdiff, ///< subgrid-diff.-scaling vector
        const LINALG::FixedSizeSerialDenseMatrix<3,iel>& evel,///< nodal velocities
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>& edens, ///< density
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>& epot,  ///< el. potential at element nodes
        const double                                     dt, ///< current time-step length
        const double&                                    timefac, ///< time discretization factor
        const enum Condif3::TauType                      whichtau, ///< flag for stabilization parameter definition
        string                                           fssgd, ///< subgrid-diff. flag
        const bool&                                      is_stationary, ///< flag indicating stationary formulation
        const bool                                       initial, ///< flag for intitial calculation
        const double&                                    frt ///< factor F/RT needed for ELCH calculations
    );

    //! evaluate shape functions and their derivatives at current integration point
    virtual void EvalShapeFuncAndDerivsAtIntPoint(
        const DRT::UTILS::IntegrationPoints3D& intpoints, ///< integration points
        const int&                             iquad,     ///< id of current Gauss point
        const bool&                            higher_order_ele,///< are second derivatives needed?
        const DRT::ELEMENTS::Condif3*          ele        ///< the element
    );

    //! calculate 2nd derivatives of shape functions (local and global)
    virtual void CalSecondDeriv(
        const double& e1, ///< first coordinate of GP
        const double& e2, ///< second coordinate of GP
        const double& e3  ///< third coordinate of GP
    );

    //! calculate element matrix and rhs vector (One-Step-Theta, BDF2)
    virtual void CalMat(
        Epetra_SerialDenseMatrix& estif,
        Epetra_SerialDenseVector& eforce,
        const bool                higher_order_ele,
        const double&             timefac,
        const int&                dofindex
    );

    //! calculate element matrix and rhs vector (Stationary problem formulation)
    virtual void CalMatStationary(
        Epetra_SerialDenseMatrix& estif,
        Epetra_SerialDenseVector& eforce,
        const bool                higher_order_ele,
        const int&                dofindex
    );

    //! calculate element matrix and rhs vector (Electrochemistry problems)
    virtual void CalMatElch(
        Epetra_SerialDenseMatrix&                                 emat,   ///< element matrix
        Epetra_SerialDenseVector&                                 erhs,   ///< element rhs
        const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ephinp,///< current scalar field
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>&          epotnp, ///< el. potential at element nodes
        const bool&                                               higher_order_ele, ///< are second derivatives needed?
        const double&                                             frt,    ///< factor F/RT
        const bool&                                               is_stationary, ///< flag indicating stationary formulation
        const double&                                             timefac ///< time discretization factor
    );

    //! calculate error compared to analytical solutions
    virtual void CalErrorComparedToAnalytSolution(
        const DRT::ELEMENTS::Condif3*                    ele, ///< the element
        ParameterList&                                   params, ///< parameter list
        const LINALG::FixedSizeSerialDenseMatrix<iel,2>& ephinp, ///< actual numerical sol. of concentrations
        const LINALG::FixedSizeSerialDenseMatrix<iel,1>& epotnp, ///< actual numerical solution of el. potential
        Epetra_SerialDenseVector&                        errors, ///< vector of error values
        struct _MATERIAL*                                material ///< material pointer
    );

    //! number of dof per node
    const int numdofpernode_;
    //! number of transported scalars (numscal_ <= numdofpernode_)
    const int numscal_;
    //! node coordinates
    LINALG::FixedSizeSerialDenseMatrix<3,iel> xyze_;
    //! bodyforce in element nodes
    vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> > bodyforce_;
    //! diffusivity / diffusivities (in case of systems) / thermal conductivity
    vector<double> diffus_;
    //! charge numbers of all ion species (ELCH specific)
    vector<double> valence_;
    //! specific heat capacity at constant pressure (in case of temperature eq.)
    double shcacp_;
    //! array for shape functions
    LINALG::FixedSizeSerialDenseMatrix<iel,1> funct_;
    //! array for density-weighted shape functions
    LINALG::FixedSizeSerialDenseMatrix<iel,1> densfunct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::FixedSizeSerialDenseMatrix<3,iel> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::FixedSizeSerialDenseMatrix<6,iel> deriv2_;
    //! transposed jacobian "dx/ds"
    LINALG::FixedSizeSerialDenseMatrix<3,3> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::FixedSizeSerialDenseMatrix<3,3> xij_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::FixedSizeSerialDenseMatrix<3,iel> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::FixedSizeSerialDenseMatrix<6,iel> derxy2_;
    //! rhs force in Gauss point
    vector<double> rhs_;
    //! vector containing all values from previous timelevel n in gausspoint
    vector<double> hist_;
    //! velocity vector in gausspoint
    LINALG::FixedSizeSerialDenseMatrix<3,1> velint_;
    //! migration velocity vector (divided by D_k*z_k) in gausspoint (ELCH)
    LINALG::FixedSizeSerialDenseMatrix<3,1> migvelint_;
    //! the stabilisation parameters (one per transported scalar)
    vector<double> tau_;
    //! artificial diffusivities (one per transported scalar)
    vector<double> kart_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::FixedSizeSerialDenseMatrix<6,3> xder2_;
    //! integration factor for current GP: fac = GaussWeight * det(J)
    double fac_;
    //! convective part used in stabilization terms
    LINALG::FixedSizeSerialDenseMatrix<iel,1> conv_;
    //! diffusive part used in stabilization terms
    LINALG::FixedSizeSerialDenseMatrix<iel,1> diff_;
    //! migration part (divided by D_k*z_k) used in stabilization terms (ELCH)
    LINALG::FixedSizeSerialDenseMatrix<iel,1> mig_;
    //! (global) gradient of el. potential at integration point (ELCH)
    LINALG::FixedSizeSerialDenseMatrix<3,1> gradpot_;
    //! values of all transported scalars in gausspoint (ELCH)
    vector<double> conint_;
    //! (global) gradient of transported scalar at integration point
    LINALG::FixedSizeSerialDenseMatrix<3,1> gradphi_;

  };

}
}

#endif

#endif // CCADISCRET
#endif // D_FLUID3

