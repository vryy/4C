/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_impl.H

\brief Internal implementation of scalar transport elements

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#ifndef SCATRA_ELE_IMPL_H
#define SCATRA_ELE_IMPL_H


#include "scatra_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


// forward declarations
struct _MATERIAL;


namespace DRT
{
namespace ELEMENTS
{

/// Interface base class for ScaTraImpl
/*!
  This class exists to provide a common interface for all template
  versions of ScaTraImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of ScaTraImpl.
 */
class ScaTraImplInterface
{

public:
  /// Empty constructor
  ScaTraImplInterface() {}
  /// Empty destructor
  virtual ~ScaTraImplInterface() {}
  /// Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in ScaTraImpl.
   */
  virtual int Evaluate(DRT::Element*              ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  //! calculate mass flux (variable number of nodes)
  virtual void CalculateFluxSerialDense(
      LINALG::SerialDenseMatrix&      flux,
      DRT::Element*&            ele,
      vector<double>&           ephinp,
      struct _MATERIAL*         material,
      bool&                     temperature,
      double&                   frt,
      Epetra_SerialDenseVector& evel,
      std::string&              fluxtypestring,
      int&                      dofindex
  ) = 0;

  /// Internal implementation class for scalar transport elements
  static ScaTraImplInterface* Impl(DRT::Element* ele);

};

  /// Internal Scalar transport element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the transport element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Transport)
    from the mathematical contents (this class). Of course there are
    different implementations of the Transport element, this is just one
    such implementation.

    The Transport element will allocate exactly one object of this class 
    for all transport elements with the same number of nodes in the mesh. 
    This allows us to use exactly matching working arrays (and keep them
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
class ScaTraImpl: public ScaTraImplInterface
{
public:
  /// Constructor
  ScaTraImpl(int numdofpernode, int numscal);

  //! number of nodes
  static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of components necessary to store second derivatives
  // 1 component  for nsd=1:  (N,xx)
  // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
  // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
  static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

  //! Evaluate
  virtual int Evaluate(DRT::Element*              ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       );

  //! calculate mass flux (variable number of nodes)
  virtual void CalculateFluxSerialDense(
      LINALG::SerialDenseMatrix&      flux,
      DRT::Element*&            ele,
      vector<double>&           ephinp,
      struct _MATERIAL*         material,
      bool&                     temperature,
      double&                   frt,
      Epetra_SerialDenseVector& evel,
      std::string&              fluxtypestring,
      int&                      dofindex
  );

  private:

    //! Calculate matrix and rhs. Here the whole thing is hidden.
    virtual void Sysmat(
        DRT::Element*                         ele, ///< the element those matrix is calculated
        const vector<LINALG::Matrix<iel,1> >& ephinp,///< current scalar field
        const vector<LINALG::Matrix<iel,1> >& ehist, ///< rhs from beginning of time step
        const LINALG::Matrix<iel,1>&          edensnp, ///< current density field
        const LINALG::Matrix<iel,1>&          epotnp, ///< el. potential at element nodes
        Epetra_SerialDenseMatrix&             sys_mat,///< element matrix to calculate
        Epetra_SerialDenseVector&             residual, ///< element rhs to calculate
        Epetra_SerialDenseVector&             subgrdiff, ///< subgrid-diff.-scaling vector
        struct _MATERIAL*                     material, ///< material pointer
        const double                          time, ///< current simulation time
        const double                          dt, ///< current time-step length
        const double                          timefac, ///< time discretization factor
        const double                          alphaF, ///< factor for generalized-alpha time integration
        const LINALG::Matrix<nsd_,iel>&       evelnp,///< nodal velocities at t_{n+1}
        const bool                            temperature, ///< temperature flag
        const bool                            conservative, ///< flag for conservative form
        const enum SCATRA::TauType            whichtau, ///< flag for stabilization parameter definition
        const string                          fssgd, ///< subgrid-diff. flag
        const bool                            is_stationary, ///< flag indicating stationary formulation
        const bool                            is_genalpha, ///< generalized-alpha flag
        const bool                            islinear, ///< flag for linear/nonlinear problem
        const double                          frt, ///< factor F/RT needed for ELCH calculations
        const double                          thermpressdt ///< time deriv. of therm. press. (for temp. eq.)
    );

    //! calculate mass matrix and rhs for determining initial time derivative
    virtual void InitialTimeDerivative(
        DRT::Element*                         ele,  ///< the element
        const vector<LINALG::Matrix<iel,1> >& ephi0, ///< initial scalar field
        const LINALG::Matrix<iel,1>&          edens0, ///< initial density field
        const LINALG::Matrix<iel,1>&          epot0,  ///< el. potential at element nodes
        Epetra_SerialDenseMatrix&             massmat, ///< mass matrix to calculate
        Epetra_SerialDenseVector&             rhs, ///< element rhs to calculate
        Epetra_SerialDenseVector&             subgrdiff, ///< subgrid-diff.-scaling vector
        struct _MATERIAL*                     material, ///< material pointer
        const double                          time, ///< current simulation time
        const double                          dt, ///< current time-step length
        const double                          timefac, ///< time discretization factor
        const LINALG::Matrix<nsd_,iel>&       evel0, ///< initial nodal velocities
        const bool                            temperature, ///< temperature flag
        const bool                            conservative, ///< flag for conservative form
        const enum SCATRA::TauType            whichtau, ///< flag for stabilization parameter definition
        const string                          fssgd, ///< subgrid-diff. flag
        const double                          frt, ///< factor F/RT needed for ELCH calculations
        const double                          thermpressdt ///< time deriv. of therm. press. (for temp. eq.)
    );

    //! calculate normalized subgrid-diffusivity matrix
    virtual void CalcSubgridDiffMatrix(
        const DRT::Element*           ele,
        Epetra_SerialDenseMatrix&     sys_mat_sd,
        const double                  timefac,
        const bool                    is_stationary
    );

    //! get the body force
    virtual void BodyForce(
        const DRT::Element*    ele, ///< the element we are dealing with
        const double           time ///< current simulation time
    );

    //! get the material parameters
    virtual void GetMaterialParams(
        struct _MATERIAL*       material,    ///< material pointer
        const bool&             temperature  ///< temperature flag
    );

    //! calculate stabilization parameter
    virtual void CalTau(
        DRT::Element*                    ele, ///< the element
        Epetra_SerialDenseVector&        subgrdiff, ///< subgrid-diff.-scaling vector
        const LINALG::Matrix<nsd_,iel>&  evel,///< nodal velocities
        const LINALG::Matrix<iel,1>&     edens, ///< density
        const LINALG::Matrix<iel,1>&     epot,  ///< el. potential at element nodes
        const double                     dt, ///< current time-step length
        const double&                    timefac, ///< time discretization factor
        const enum SCATRA::TauType       whichtau, ///< flag for stabilization parameter definition
        const string                     fssgd, ///< subgrid-diff. flag
        const bool&                      is_stationary, ///< flag indicating stationary formulation
        const bool                       initial, ///< flag for intitial calculation
        const double&                    frt ///< factor F/RT needed for ELCH calculations
    );

    //! evaluate shape functions and their derivatives at current integration point
    virtual void EvalShapeFuncAndDerivsAtIntPoint(
        const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, ///< integration points
        const int&                                   iquad,     ///< id of current Gauss point
        const bool&                                  higher_order_ele,///< are second derivatives needed?
        const int&                                   eleid       ///< the element id
    );

    //! calculate element matrix and rhs vector (One-Step-Theta, BDF2, Gen.-alpha)
    virtual void CalMat(
        Epetra_SerialDenseMatrix&            estif,
        Epetra_SerialDenseVector&            eforce,
        const vector<LINALG::Matrix<iel,1> >& ephinp,
        const bool                           higher_order_ele,
        const bool                           conservative,
        const bool                           is_genalpha,
        const double                         timefac,
        const double                         alphaF,
        const int                            dofindex
    );

    //! calculate element matrix and rhs vector (Stationary problem formulation)
    virtual void CalMatStationary(
        Epetra_SerialDenseMatrix& estif,
        Epetra_SerialDenseVector& eforce,
        const bool                higher_order_ele,
        const bool                conservative,
        const int                 dofindex
    );

    //! calculate matrix and rhs vector (incremental formulation)
    virtual void CalMatInc(
        Epetra_SerialDenseMatrix&             emat,
        Epetra_SerialDenseVector&             erhs,
        const vector<LINALG::Matrix<iel,1> >& ephinp,
        const bool                            higher_order_ele,
        const bool                            is_stationary,
        const double                          timefac
    );

    //! calculate element matrix and rhs vector (Electrochemistry problems)
    virtual void CalMatElch(
        Epetra_SerialDenseMatrix&             emat,   ///< element matrix
        Epetra_SerialDenseVector&             erhs,   ///< element rhs
        const vector<LINALG::Matrix<iel,1> >& ephinp, ///< current scalar field
        const LINALG::Matrix<iel,1>&          epotnp, ///< el. potential at element nodes
        const bool                            higher_order_ele, ///< are second derivatives needed?
        const double                          frt,    ///< factor F/RT
        const bool                            is_stationary, ///< flag indicating stationary formulation
        const double                          timefac ///< time discretization factor
    );

    //! calculate error compared to analytical solutions
    virtual void CalErrorComparedToAnalytSolution(
        const DRT::Element*          ele, ///< the element
        ParameterList&               params, ///< parameter list
        const LINALG::Matrix<iel,2>& ephinp, ///< actual numerical sol. of concentrations
        const LINALG::Matrix<iel,1>& epotnp, ///< actual numerical solution of el. potential
        Epetra_SerialDenseVector&    errors, ///< vector of error values
        struct _MATERIAL*            material ///< material pointer
    );

    //! calculate flux vector for actual solution
    virtual void CalculateFlux(
        LINALG::Matrix<3,iel>&          flux,
        const DRT::Element*             ele,
        const vector<double>&           ephinp,
        struct _MATERIAL*               material,
        const bool                      temperature,
        const double                    frt,
        const Epetra_SerialDenseVector& evel,
        const SCATRA::FluxType          fluxtype,
        const int                       dofindex
    );

    //! calculate integrals of temperature/concentrations, density and domain
   virtual void CalculateTempAndDens(
        const DRT::Element*       ele,     ///< pointer to the element
        const vector<double>&     ephinp,  ///< current scalar field
        const vector<double>&     edensnp, ///< current density field
        Epetra_SerialDenseVector& scalars  ///< vector for the computed integrals
        );

    //! calculate domain and bodyforce integral
   virtual void CalculateDomainAndBodyforce(
        Epetra_SerialDenseVector&  scalars, ///< vector for the computed integrals
        const DRT::Element*        ele,     ///< pointer to the element
        const double               time     ///< current simulation time
        );

    //! calculate the Laplacian for all shape functions(strong form)
    inline void GetLaplacianStrongForm
    (
      LINALG::Matrix<iel,1>& diff,
      const LINALG::Matrix<numderiv2_,iel>& deriv2
     )
    {
      diff.Clear();
      // compute N,xx  +  N,yy +  N,zz for each shape function at integration point
      for (int i=0; i<iel; ++i)
      {
        for (int j = 0; j<nsd_; ++j)
        {
          diff(i) += deriv2(j,i);
        }
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    inline void GetLaplacianWeakForm(
        double& val,
        const LINALG::Matrix<nsd_,iel>& derxy,
        const int vi,
        const int ui)
    {
      val = 0.0;
      for (int j = 0; j<nsd_; j++)
      {
        val += derxy(j, vi)*derxy(j, ui);
      }
      return;
    };

    //! calculate momentum divergence
    void GetMomentumDivergence(double&                          mdiv,
                               const LINALG::Matrix<nsd_,nsd_>& mderxy)
    {
      mdiv = 0.0;
      // compute (rho*vel),x  + (rho*vel),y +  (rho*vel),z at integration point
      for (int j = 0; j<nsd_; ++j)
      {
        mdiv += mderxy(j,j);
      }
      return;
    };


    //! number of dof per node
    const int numdofpernode_;
    //! number of transported scalars (numscal_ <= numdofpernode_)
    const int numscal_;
    //! node coordinates
    LINALG::Matrix<nsd_,iel> xyze_;
    //! bodyforce in element nodes
    vector<LINALG::Matrix<iel,1> > bodyforce_;
    //! diffusivity / diffusivities (in case of systems) / thermal conductivity
    vector<double> diffus_;
    //! charge numbers of all ion species (ELCH specific)
    vector<double> valence_;
    //! specific heat capacity at constant pressure (in case of temperature eq.)
    double shcacp_;
    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_,1> xsi_;
    //! array for shape functions
    LINALG::Matrix<iel,1> funct_;
    //! array for density-weighted shape functions
    LINALG::Matrix<iel,1> densfunct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,iel> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<numderiv2_,iel> deriv2_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xij_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,iel> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<numderiv2_,iel> derxy2_;
    //! rhs force in Gauss point
    vector<double> rhs_;
    //! vector containing all values from previous timelevel n in gausspoint
    vector<double> hist_;
    //! velocity vector in gausspoint
    LINALG::Matrix<nsd_,1> velint_;
    //! migration velocity vector (divided by D_k*z_k) in gausspoint (ELCH)
    LINALG::Matrix<nsd_,1> migvelint_;
    //! momentum (density-weighted velocity) derivatives
    LINALG::Matrix<nsd_,nsd_> mderxy_;
    //! momentum divergence
    double mdiv_;
    //! the stabilisation parameters (one per transported scalar)
    vector<double> tau_;
    //! artificial diffusivities (one per transported scalar)
    vector<double> kart_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::Matrix<numderiv2_,nsd_> xder2_;
    //! integration factor for current GP: fac = GaussWeight * det(J)
    double fac_;
    //! convective part used in stabilization terms
    LINALG::Matrix<iel,1> conv_;
    //! diffusive part used in stabilization terms
    LINALG::Matrix<iel,1> diff_;
    //! convective part of migration term (divided by D_k*z_k) used in stabilization terms (ELCH)
    LINALG::Matrix<iel,1> migconv_;
    // reactive part of migration term used in stabilization terms
    LINALG::Matrix<iel,1> migrea_;
    //! (global) gradient of el. potential at integration point (ELCH)
    LINALG::Matrix<nsd_,1> gradpot_;
    //! values of all transported scalars in gausspoint (ELCH)
    vector<double> conint_;
    //! (global) gradient of transported scalar at integration point
    LINALG::Matrix<nsd_,1> gradphi_;
    //! (global) 2nd gradient (Laplacian) of transported scalar at integration point
    LINALG::Matrix<numderiv2_,1> lapphi_;
    //! (global) Laplacian of shape functions at integration point
    LINALG::Matrix<iel,1> laplace_;

  };

}
}

#endif

#endif // CCADISCRET
#endif // D_FLUID3 or D_FLUID2

