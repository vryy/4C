/*----------------------------------------------------------------------*/
/*!
\file fluid2_stationary.H

\brief Internal implementation of Fluid2 element (stationary formulation)

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

#ifndef FLUID2_STATIONARY_H
#define FLUID2_STATIONARY_H

#include "fluid2.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
//#include <blitz/array.h>


namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid2Stationary
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid2Stationary. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid2Stationary.
   */
  class Fluid2StationaryInterface
  {
  public:
    /// Empty constructor
    Fluid2StationaryInterface() {}
    /// Empty destructor
    virtual ~Fluid2StationaryInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int Evaluate(Fluid2*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat) = 0;

    /// Internal implementation class for fluid element
    static Fluid2StationaryInterface* Impl(DRT::ELEMENTS::Fluid2* f2);

  };

  /// Internal Fluid2 implementation -- stationary formulation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid2 element in stationary formulation.
    Additionally the method Sysmat() provides a clean and fast
    element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid2)
    from the mathematical contents (this class). Of course there are
    different implementations of the Fluid2 element, this is just one
    such implementation.

    The Fluid2 will allocate exactly one object of this class for all
    fluid elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    For matrix calculation the blitz library is used.

    \author vg
    \date 08/08
  */
  template <DRT::Element::DiscretizationType distype>
  class Fluid2Stationary: public Fluid2StationaryInterface
  {
  public:
    /// Constructor
    explicit Fluid2Stationary();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /// Evaluate
    /*!
      The evaluate function for the general fluid case.
      Why doesn't someone who actually knows what this function does
      write a useful comment?
     */
    virtual int Evaluate(Fluid2*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat);
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
    void Sysmat(Fluid2*                                  ele,
                const LINALG::Matrix<2,iel>&             evelnp,
                const LINALG::Matrix<2,iel>&             fsevelnp,
                const LINALG::Matrix<iel,1>&             eprenp,
                LINALG::Matrix<3*iel,3*iel>&             estif,
                LINALG::Matrix<3*iel,    1>&             eforce,
                Teuchos::RCP<const MAT::Material>        material,
                double                                   pseudotime,
                const bool                               newton,
                const bool                               conservative,
                const bool                               higher_order_ele,
                const enum Fluid2::FineSubgridVisc       fssgv,
                const enum Fluid2::StabilisationAction   pspg,
                const enum Fluid2::StabilisationAction   supg,
                const enum Fluid2::StabilisationAction   vstab,
                const enum Fluid2::StabilisationAction   cstab,
                const enum Fluid2::StabilisationAction   cross,
                const enum Fluid2::StabilisationAction   reynolds,
                const double                             Cs);

    /// Internal implementation class for specific fluid element
    static Fluid2Stationary* StationaryImpl(DRT::ELEMENTS::Fluid2* f2);

  private:

    /*!
      \brief calculate stabilization parameter for stationary formulation
    */
    void CalTauStationary(Fluid2*                             ele,
                          const LINALG::Matrix<2,iel>&        evelnp,
                          const LINALG::Matrix<2,iel>&        fsevelnp,
                          const double                        visc,
                          const enum Fluid2::FineSubgridVisc  fssgv,
                          const double                        Cs);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid2*      ele,
                   const double pseudotime);

  private:
    double vart_;                                         /// fine-scale subgrid viscosity
    LINALG::Matrix<2,iel> xyze_;      /// nodal position array
    LINALG::Matrix<2,iel> edeadng_;   /// nodal dead load
    LINALG::Matrix<iel,1> funct_;     /// nodal shape functions
    LINALG::Matrix<2,iel> deriv_;     /// shape derivatives with respect to \xi
    LINALG::Matrix<3,iel> deriv2_;    /// second shape derivatives with respect to \xi
    LINALG::Matrix<2,2>   xjm_;       /// transposed of the jacobian matrix d x / d \xi
    LINALG::Matrix<2,2>   xji_;       /// inverse of the jacobian matrix
    LINALG::Matrix<2,2>   vderxy_;    /// velocity vector derivative
    LINALG::Matrix<2,2>   fsvderxy_;  /// fine-scale velocity vector derivative
    LINALG::Matrix<2,iel> derxy_;     /// shape derivative with respect to x
    LINALG::Matrix<3,iel> derxy2_;    /// second shape derivative with respect to x
    LINALG::Matrix<2,1>   bodyforce_; /// volumetric load at gausspoint
    //LINALG::Matrix<,1> velino_;
    LINALG::Matrix<2,1>   velint_;    /// density-weighted velocity
    LINALG::Matrix<2,1>   fsvelint_;  /// fine-scale velocity
    LINALG::Matrix<2,1>   gradp_;     /// pressure gradient
    LINALG::Matrix<3,1>   tau_;       /// tau_Mu, tau_Mp, tau_C
    LINALG::Matrix<4,iel> viscs2_;    /// viscous term including 2nd derivatives
    LINALG::Matrix<iel,1> conv_c_;    /// linearisation of convect, convective part
    double vdiv_;                      /// velocity divergence
    LINALG::Matrix<2,1>   rhsmom_;    /// momentum equation right hand side
    LINALG::Matrix<2,1>   conv_old_;
    LINALG::Matrix<2,1>   visc_old_;
    LINALG::Matrix<2,1>   res_old_;
    LINALG::Matrix<iel,1> conv_resM_;
    LINALG::Matrix<3,2>   xder2_;
  };

}
}

#endif

#endif
#endif
