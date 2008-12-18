/*----------------------------------------------------------------------*/
/*!
\file fluid3_stationary.H

\brief Internal implementation of Fluid3 element -- stationary formulation

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

#ifndef FLUID3_STATIONARY_H
#define FLUID3_STATIONARY_H

#include "fluid3.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
//#include <blitz/array.h>


namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid3StationaryImpl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid3StationaryImpl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid3StationaryImpl.
   */
  class Fluid3StationaryImplInterface
  {
  public:
    /// Empty constructor
    Fluid3StationaryImplInterface() {}
    /// Empty destructor
    virtual ~Fluid3StationaryImplInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3StationaryImpl.
     */
    virtual int Evaluate(Fluid3*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat,
                         _MATERIAL* actmat) = 0;

    /// Internal implementation class for fluid element -- stationary formulation
    static Fluid3StationaryImplInterface* Impl(DRT::ELEMENTS::Fluid3* f3);

  };

  /// Internal Fluid3 implementation -- stationary formulation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3 element in stationary formulation.
    Additionally the method Sysmat() provides a clean and fast
    element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid3)
    from the mathematical contents (this class). Of course there are
    different implementations of the Fluid3 element, this is just one
    such implementation.

    The Fluid3 will allocate exactly one object of this class for all
    fluid elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>History</h3>
    Right now we do not read any stabilization parameters from the
    input file but have a fixed version with precompiled values.

    \author gjb
    \date 11/07
  */
  template <DRT::Element::DiscretizationType distype>
  class Fluid3StationaryImpl : public Fluid3StationaryImplInterface
  {
  public:
    /// Constructor
    explicit Fluid3StationaryImpl();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /// Evaluate
    /*!
      The evaluate function for the stationary fluid case.
      Why doesn't someone who actually knows what this function does
      write a useful comment?
     */
    virtual int Evaluate(Fluid3*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat,
                         _MATERIAL*                 actmat);

    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs for stationary problem formulation
      Here the whole thing is hidden.

      \param ele     		  (i) the element those matrix is calculated
      \param evelnp  		  (i) nodal velocities at n+1
      \param fsevelnp         (i) fine-scale nodal velocities at n+1
      \param eprenp  		  (i) nodal pressure
      \param edensnp          (i) nodal density
      \param estif   		  (o) element matrix to calculate
      \param eforce  		  (o) element rhs to calculate
      \param material 		  (i) fluid material
      \param pseudotime       (i) current simulation pseudotime
      \param newton     	  (i) full Newton or fixed-point-like
      \param loma             (i) boolean flag for potential low-Mach-number solver
      \param higher_order_ele (i) keep or drop second derivatives
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
      \param pspg     		  (i) boolean flag for stabilisation
      \param supg     		  (i) boolean flag for stabilisation
      \param vstab     		  (i) boolean flag for stabilisation
      \param cstab     		  (i) boolean flag for stabilisation
      \param cross     		  (i) boolean flag for stabilisation
      \param reynolds     	  (i) boolean flag for stabilisation
      \param Cs               (i) Smagorinsky model parameter
     */
    void Sysmat(Fluid3*                                          ele,
                const LINALG::Matrix<3,iel>& evelnp,
                const LINALG::Matrix<3,iel>& fsevelnp,
                const LINALG::Matrix<iel,1>& eprenp,
                const LINALG::Matrix<iel,1>& edensnp,
                LINALG::Matrix<4*iel,4*iel>& estif,
                LINALG::Matrix<4*iel,1>&     eforce,
                struct _MATERIAL*                                material,
                double                                           pseudotime,
                const bool                                       newton,
                const bool                                       loma,
                const bool                                       conservative,
                const bool                                       higher_order_ele,
                const enum Fluid3::StabilisationAction           fssgv,
                const enum Fluid3::StabilisationAction           pspg,
                const enum Fluid3::StabilisationAction           supg,
                const enum Fluid3::StabilisationAction           vstab,
                const enum Fluid3::StabilisationAction           cstab,
                const enum Fluid3::StabilisationAction           cross,
                const enum Fluid3::StabilisationAction           reynolds,
                const double                                     Cs);

  private:

    /*!
      \brief calculate stabilization parameter for stationary formulation
    */
    void CalTauStationary(
      Fluid3*                                          ele,
      const LINALG::Matrix<3,iel>& evelnp,
      const LINALG::Matrix<3,iel>& fsevelnp,
      const LINALG::Matrix<iel,1>& edensnp,
      const double                                     visc,
      const enum Fluid3::StabilisationAction           fssgv,
      const double                                     Cs);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid3*      ele,
                   const double pseudotime);


  private:
    double vart_;                      /// fine-scale subgrid viscosity
    LINALG::Matrix<3,iel> xyze_;      /// nodal position array
    LINALG::Matrix<3,iel> edeadng_;   /// nodal dead load
    LINALG::Matrix<iel,1> funct_;     /// nodal shape functions
    LINALG::Matrix<iel,1> densfunct_; /// density-weighted nodal shape functions
    LINALG::Matrix<iel,1> functdens_; /// nodal shape functions divided by density
    LINALG::Matrix<3,iel> deriv_;     /// shape derivatives with respect to \xi
    LINALG::Matrix<6,iel> deriv2_;    /// second shape derivatives with respect to \xi
    LINALG::Matrix<3,3> xjm_;       /// transposed of the jacobian matrix d x / d \xi
    LINALG::Matrix<3,3> xji_;       /// inverse of the jacobian matrix
    LINALG::Matrix<3,3> vderxy_;    /// velocity vector derivative
    LINALG::Matrix<3,3> mderxy_;    /// momentum vector derivative
    LINALG::Matrix<3,3> fsvderxy_;  /// fine-scale velocity vector derivative
    LINALG::Matrix<3,iel> derxy_;     /// shape derivative with respect to x
    LINALG::Matrix<3,iel> densderxy_; /// density-weighted shape derivative w.r.t. x
    LINALG::Matrix<6,iel> derxy2_;    /// second shape derivative with respect to x
    LINALG::Matrix<3,1> bodyforce_; /// volumetric load at gausspoint
    LINALG::Matrix<3,1> velino_;
    LINALG::Matrix<3,1> velint_;    /// density-weighted velocity
    LINALG::Matrix<3,1> ndwvelint_; /// non-density-weighted velocity
    LINALG::Matrix<3,1> fsvelint_;  /// fine-scale velocity
    LINALG::Matrix<3,1> gradp_;     /// pressure gradient
    LINALG::Matrix<3,1> tau_;       /// tau_Mu, tau_Mp, tau_C
    LINALG::Matrix<3*3,iel> viscs2_;    /// viscous term including 2nd derivatives
    LINALG::Matrix<iel,1> conv_c_;    /// linearisation of convect, convective part
    LINALG::Matrix<iel,1> ndwconv_c_; /// non-density-weighted linearisation of conv.
    double mdiv_;                      /// momentum divergence
    double vdiv_;                      /// velocity divergence
    LINALG::Matrix<3,1> rhsmom_;    /// momentum equation right hand side
    LINALG::Matrix<3,1> conv_old_;
    LINALG::Matrix<3,1> visc_old_;
    LINALG::Matrix<3,1> res_old_;
    LINALG::Matrix<iel,1> conv_resM_;
    LINALG::Matrix<6,3> xder2_;
  };

}
}

#endif

#endif
#endif
