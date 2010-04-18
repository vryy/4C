/*----------------------------------------------------------------------*/
/*!
\file fluid3_lin_impl.H

\brief Internal implementation of linearised Fluid3 element

<pre>
Maintainer: Christiane Foerster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_LIN_IMPL_H
#define FLUID3_LIN_IMPL_H

#include "fluid3.H"
#include "fluid3_ele_impl_utils.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
//#include <blitz/array.h>


namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid3lin_Impl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid3lin_Impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid3lin_Impl.
   */
  class Fluid3lin_ImplInterface
  {
  public:
    /// Empty constructor
    Fluid3lin_ImplInterface() {}
    /// Empty destructor
    virtual ~Fluid3lin_ImplInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3lin_Impl.
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
                         RefCountPtr<MAT::Material> mat) = 0;

    /// Internal implementation class for linear fluid element
    static Fluid3lin_ImplInterface* Impl(DRT::ELEMENTS::Fluid3* f3);
  };

  /// Internal Fluid3 implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the linearised Fluid3 element. Additionally the method
    Sysmat() provides a clean and fast element implementation.

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

    In particular the present implementation is based on the information
    obtained from Miguel Fernandez (personal communication 02/08). Thus
    the fluid element is linearised by using a previous velocity within
    the convective term. This also naturally linearises all stabilisation
    terms. An additional stabilisation term is required to ensure stability
    of the linearised matrix.

    <h3>Properties</h3>

    The linearisation reduces the element generally to first order accuracy
    in time.
    The linearisation DOES NOT work for ALE. Thus the present element
    is for pure fluid calculations only!!!

    <h3>History</h3>

    The implementation here is based on the default fluid3 element. Further
    historic comment can therefore be found in fluid3_impl.H.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version with precompiled values.

    \author chfoe
    \date 02/08
  */
  template <DRT::Element::DiscretizationType distype>
  class Fluid3lin_Impl : public Fluid3lin_ImplInterface
  {
  public:
    /// Constructor
    explicit Fluid3lin_Impl(int numdofpernode);

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of space dimensions
    static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

    /// Evaluate
    /*!
      The evaluate function for the linear fluid case.
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
                         RefCountPtr<MAT::Material> mat);

    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs. Here the whole thing is hidden.

      \note Stationary calculations are not (yet) supported.

      \param ele     (i) the element those matrix is calculated
      \param evelnp  (i) nodal velocities at n+1
      \param eprenp  (i) nodal pressure
      \param edensnp (i) nodal density
      \param emhist  (i) rhs for momentum equation from beginning of time step
      \param echist  (i) rhs for continuity equation from beginning of time step
      \param edispnp (i) nodal displacements (on moving mesh)
      \param egridv  (i) grid velocity (on moving mesh)
      \param estif   (o) element matrix to calculate
      \param esv     (o) fine-scale subgrid-viscosity element matrix
      \param eforce  (o) element rhs to calculate
      \param material (i) fluid material
      \param time     (i) current simulation time
      \param timefac  (i) time discretization factor
      \param fssgv    (i) flag for type of fine-scale subgrid viscosity
      \param Cs_fs    (i) Smagorinsky model param. for fine-scale subgrid viscosity
      */
    void Sysmat(Fluid3* ele,
                const LINALG::Matrix<nsd_,iel>&     evelnp,
                const LINALG::Matrix<iel,1>&     eprenp,
                const LINALG::Matrix<iel,1>&     edensnp,
                const LINALG::Matrix<nsd_,iel>&     emhist,
                const LINALG::Matrix<iel,1>&     echist,
                LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&     estif,
                LINALG::Matrix<(nsd_+1)*iel,1>&         eforce,
                Teuchos::RCP<const MAT::Material>                    material,
                double                                               time,
                double                                               timefac,
                INPAR::FLUID::PhysicalType							 physicaltype);

  private:

    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_,1> xsi_;

    //! number of components necessary to store second derivatives
    // 1 component  for nsd=1:  (N,xx)
    // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
    // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
    static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;
    /*!
      \brief calculate stabilization parameter
    */
    void Caltau(
      Fluid3*                                          ele,
      const LINALG::Matrix<nsd_,iel>& evelnp,
      const LINALG::Matrix<iel,1>& edensnp,
      const double                                     visc,
      const double                                     timefac);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid3* ele, const double time);


  private:
    //number of dof's
    const int numdofpernode_;

    LINALG::Matrix<nsd_,iel> xyze_;
    LINALG::Matrix<nsd_,iel> edeadng_;
    LINALG::Matrix<iel,1> funct_;
    LINALG::Matrix<iel,1> densfunct_;
    LINALG::Matrix<iel,1> functdens_;
    LINALG::Matrix<nsd_,iel> deriv_;
    LINALG::Matrix<numderiv2_,iel> deriv2_;
    LINALG::Matrix<nsd_,nsd_> xjm_;
    LINALG::Matrix<nsd_,nsd_> xji_;
    LINALG::Matrix<nsd_,nsd_> vderxy_;
    LINALG::Matrix<nsd_,nsd_> mderxy_;
    LINALG::Matrix<nsd_,1> pderxy_;
    LINALG::Matrix<nsd_,numderiv2_> vderxy2_;
    LINALG::Matrix<nsd_,iel> derxy_;
    LINALG::Matrix<nsd_,iel> densderxy_;
    LINALG::Matrix<numderiv2_,iel> derxy2_;
    LINALG::Matrix<nsd_,1> bodyforce_;
    LINALG::Matrix<nsd_,1> histmom_;
    double histcon_;
    LINALG::Matrix<nsd_,1> velino_;
    LINALG::Matrix<nsd_,1> velint_; /* newest, i.e. previous time leve velocity*/
    LINALG::Matrix<nsd_,1> gradp_;
    LINALG::Matrix<3,1> tau_;
    // this array once had three dimensions, I combined the first two
    // to one.
    LINALG::Matrix<nsd_*nsd_,iel> viscs2_; /* viscous term including 2nd derivatives */
    LINALG::Matrix<iel,1> conv_;   /* linearisation of convect, convective part */
    double mdiv_;
    LINALG::Matrix<nsd_,1> rhsmom_;
    double rhscon_;
    LINALG::Matrix<numderiv2_,nsd_> xder2_;
    // this was a blitz::Array<int,1>
    std::vector<int> numepn_;
  };

}
}

#endif

#endif
#endif
