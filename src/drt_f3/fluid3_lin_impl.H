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
#include <blitz/array.h>


namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid3lin_Impl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid3StationaryImpl. The only function
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
                         RefCountPtr<MAT::Material> mat,
                         _MATERIAL* actmat) = 0;

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
  template <int iel>
  class Fluid3lin_Impl : public Fluid3lin_ImplInterface
  {
  public:
    /// Constructor
    explicit Fluid3lin_Impl();

    /// Evaluate
    /*!
      The evaluate function for the linear fluid case.
      Why doesn't someone who actually knows what this function does
      write a usefull comment?
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
                         _MATERIAL* actmat);

    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs. Here the whole thing is hidden.

      \note Stationary calculations are not (yet) supported.

      \param ele     (i) the element those matrix is calculated
      \param evelnp  (i) nodal velocities at n+1
      \param eprenp  (i) nodal pressure
      \param evhist  (i) rhs from beginning of time step
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
                const LINALG::FixedSizeSerialDenseMatrix<3,iel>&     evelnp,
                const LINALG::FixedSizeSerialDenseMatrix<iel,1>&     eprenp,
                const LINALG::FixedSizeSerialDenseMatrix<3,iel>&     evhist,
                LINALG::FixedSizeSerialDenseMatrix<4*iel,4*iel>&           estif,
                LINALG::FixedSizeSerialDenseMatrix<4*iel,1>&           eforce,
                struct _MATERIAL*       material,
                double                  time,
                double                  timefac);

  private:

    /*!
      \brief calculate stabilization parameter
    */
    void Caltau(
      Fluid3* ele,
      const LINALG::FixedSizeSerialDenseMatrix<3,iel>&           evelnp,
      const DRT::Element::DiscretizationType  distype,
      const double                            visc,
      const double                            timefac);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid3* ele, const double time, struct _MATERIAL* material);

    /*!
     * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
     */
    void gder2(Fluid3* ele);

  private:
    LINALG::FixedSizeSerialDenseMatrix<3,iel> xyze_;
    LINALG::FixedSizeSerialDenseMatrix<3,iel> edeadng_;
    LINALG::FixedSizeSerialDenseMatrix<iel,1> funct_;
    LINALG::FixedSizeSerialDenseMatrix<3,iel> deriv_;
    LINALG::FixedSizeSerialDenseMatrix<6,iel> deriv2_;
    LINALG::FixedSizeSerialDenseMatrix<3,3> xjm_;
    LINALG::FixedSizeSerialDenseMatrix<3,3> xji_;
    LINALG::FixedSizeSerialDenseMatrix<3,3> vderxy_;
    LINALG::FixedSizeSerialDenseMatrix<3,1> pderxy_;
    LINALG::FixedSizeSerialDenseMatrix<3,6> vderxy2_;
    LINALG::FixedSizeSerialDenseMatrix<3,iel> derxy_;
    LINALG::FixedSizeSerialDenseMatrix<6,iel> derxy2_;
    LINALG::FixedSizeSerialDenseMatrix<3,1> bodyforce_;
    LINALG::FixedSizeSerialDenseMatrix<3,1> histvec_;
    LINALG::FixedSizeSerialDenseMatrix<3,1> velino_;
    LINALG::FixedSizeSerialDenseMatrix<3,1> velint_; /* newest, i.e. previous time leve velocity*/
    LINALG::FixedSizeSerialDenseMatrix<3,1> gradp_;
    LINALG::FixedSizeSerialDenseMatrix<3,1> tau_;
    // this array once had three dimensions, I combined the first two
    // to one.
    LINALG::FixedSizeSerialDenseMatrix<3*3,iel> viscs2_; /* viscous term including 2nd derivatives */
    LINALG::FixedSizeSerialDenseMatrix<iel,1> conv_;   /* linearisation of convect, convective part */
    LINALG::FixedSizeSerialDenseMatrix<3,1> rhsint_; /* total right hand side terms at int.-point */
    LINALG::FixedSizeSerialDenseMatrix<6,3> xder2_;
    // this was a blitz::Array<int,1>
    std::vector<int> numepn_;
  };

}
}

#endif

#endif
#endif
