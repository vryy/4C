/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_impl.H

\brief Internal implementation of scalar transport boundary elements

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

#ifndef SCATRA_ELE_BOUNDARY_IMPL_H
#define SCATRA_ELE_BOUNDARY_IMPL_H


#include "scatra_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "scatra_element.H"


// forward declarations
//struct _MATERIAL;


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
class ScaTraBoundaryImplInterface
{

public:
  /// Empty constructor
  ScaTraBoundaryImplInterface() {}
  /// Empty destructor
  virtual ~ScaTraBoundaryImplInterface() {}
  /// Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in ScaTraBoundaryImpl.
   */
  virtual int Evaluate(DRT::ELEMENTS::TransportBoundary* ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  /// Evaluate a Neumann boundary condition
  /*!
    This class does not provide a definition for this function, it
    must be defined in ScaTraBoundaryImpl.
   */
  virtual int EvaluateNeumann(DRT::Element*              ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1_epetra
                              ) = 0;

  /// Internal implementation class for scalar transport elements
  static ScaTraBoundaryImplInterface* Impl(DRT::Element* ele);

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
class ScaTraBoundaryImpl: public ScaTraBoundaryImplInterface
{
public:
  /// Constructor
  ScaTraBoundaryImpl(int numdofpernode, int numscal);

  //! number of nodes
  static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! Evaluate
  virtual int Evaluate(DRT::ELEMENTS::TransportBoundary* ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       );

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(DRT::Element*              ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1_epetra
                              );

  private:

    //! evaluate shape functions and derivatives at int. point
    void EvalShapeFuncAndIntFac(
        const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
        const int&                                   iquad,      ///< id of current Gauss point
        const int&                                   eleid       ///< the element id
    );

    //! evaluate an electrode kinetics boundary condition
    void EvaluateElectrodeKinetics(
       const DRT::Element*        ele,   ///< the actual boundary element
       Epetra_SerialDenseMatrix& emat,  ///< element-matrix
       Epetra_SerialDenseVector& erhs,  ///< element-rhs
       const vector<double>&   ephinp,  ///< actual conc. and pot. values
       Teuchos::RCP<const MAT::Material> material,  ///< the material
       const int           reactantid,  ///< MATERIAL id of reactive species
       const std::string     kinetics,  ///< desired electrode kinetics model
       const double              pot0,  ///< actual electrode potential on metal side
       const double            alphaa,  ///< anodic factor in Butler-Volmer model
       const double            alphac,  ///< cathodic factor in Butler-Volmer model
       const double                i0,  ///< exchange current density
       const double               frt,  ///< factor F/RT
       const double             gamma,  ///< exponent for concentration dependency
       const bool              iselch   ///< problemtype elch or scatra
       );

    //! provide electrode kinetics status information
    void ElectrodeStatus(
        const DRT::Element*        ele,  ///< the actual boundary element
        ParameterList&          params,  ///< the parameter list
        const vector<double>&   ephinp,  ///< actual conc. and pot. values
        const std::string     kinetics,  ///< desired electrode kinetics model
        const double              pot0,  ///< actual electrode potential on metal side
        const double            alphaa,  ///< anodic factor in Butler-Volmer model
        const double            alphac,  ///< cathodic factor in Butler-Volmer model
        const double                i0,  ///< exchange current density
        const double               frt,  ///< factor F/RT
        const double             gamma,  ///< exponent for concentration dependency
        const bool              iselch   ///< problemtype elch or scatra
    );

     //! integral of shape functions over boundary surface
     void IntegrateShapeFunctions(
         const DRT::Element*       ele,    ///< the actual boundary element
         ParameterList&            params, ///< the parameter list
         Epetra_SerialDenseVector& elevec1,///< result vector (to be assembled)
         const bool                addarea ///< flag for area calculation
     );

    //! integral of normal diffusive flux and velocity div. over boundary surface
    void DifffluxAndDivuIntegral(
        const DRT::Element*             ele,
        ParameterList&                  params,
        const vector<double>&           ediffflux,
        const vector<double>&           edivu);

    //! Compute a constant normal vector for a boundary element
    void GetConstNormal(
        LINALG::Matrix<nsd_+1,1>&            normal, ///< the constant normal vector
        const LINALG::Matrix<nsd_+1,iel>&    xyze    ///< element node coordinates
        );

    //! calculate potential Neumann inflow terms
    void NeumannInflow(
        const DRT::Element*                   ele,
        const vector<LINALG::Matrix<iel,1> >& ephinp,
        const vector<LINALG::Matrix<iel,1> >& edensnp,
        const LINALG::Matrix<nsd_+1,iel>&     evelnp,
        Epetra_SerialDenseMatrix&             emat,
        Epetra_SerialDenseVector&             erhs,
        const bool                            is_stationary,
        const bool                            is_genalpha,
        const bool                            is_incremental,
        const double                          timefac,
        const double                          alphaF);

    //! number of dof per node
    const int numdofpernode_;
    //! number of transported scalars (numscal_ <= numdofpernode_)
    const int numscal_;
    //! flag for using ALE formulation
    bool isale_;
    //! node coordinates
    LINALG::Matrix<nsd_+1,iel> xyze_;
    //! nodal displacement values for ALE
    LINALG::Matrix<nsd_+1,iel> edispnp_;
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
    //! array for density-weighted shape functions at n+1/n+alpha_F
    LINALG::Matrix<iel,1> densfunct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,iel> deriv_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,iel> derxy_;
    //! unit normal vector at integration point
    LINALG::Matrix<nsd_+1,1> normal_;
    //! velocity vector in gausspoint
    LINALG::Matrix<nsd_+1,1> velint_;
    //! metric tensor at integration point
    LINALG::Matrix<nsd_,nsd_>  metrictensor_;
    //! integration factor for current GP: fac = GaussWeight * drs
    double fac_;

  };

}
}

#endif

#endif // CCADISCRET
#endif // D_FLUID3 or D_FLUID2

