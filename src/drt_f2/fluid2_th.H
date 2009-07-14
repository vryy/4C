/*
 * fluid2_th.H
 *
 *  Created on: Jun 3, 2009
 *      Author: wiesner
 */

#ifdef D_FLUID2
#ifdef CCADISCRET

#ifndef FLUID2_TH_H_
#define FLUID2_TH_H_


#include "fluid2.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

namespace DRT
{
namespace ELEMENTS
{
/// Interface class for Fluid2TH
class Fluid2THInterface
{
public:

    Fluid2THInterface() {}
    virtual ~Fluid2THInterface() {}

    virtual int CalcGradPAndMassMatrix(Fluid2*            ele,
            ParameterList&        params,
            DRT::Discretization&    discretization,
            vector<int>&        lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseMatrix& elemat2_epetra,
            Epetra_SerialDenseVector&   elevec1_epetra) = 0;
    virtual int CalcSysmatAndResidual(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseVector& elevec1_epetra,
            RefCountPtr< MAT::Material >   mat) = 0;

    virtual int CalcImpulseEqnImplicit(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseVector& elevec1_epetra,
            RefCountPtr< MAT::Material >   mat) = 0;
    virtual int CalcImpulseEqnSemiImplicit(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseVector& elevec1_epetra,
            RefCountPtr< MAT::Material >   mat) = 0;

    virtual int CalcResidual(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseVector& elevec_epetra,
            RefCountPtr<MAT::Material> mat) = 0;

    static Fluid2THInterface* Impl(DRT::ELEMENTS::Fluid2* f2);
};

template<DRT::Element::DiscretizationType DISTYPE>
struct DisTypeToDOFNum {};	//< struct with number of DOFs per element (template based)
template<> struct DisTypeToDOFNum<DRT::Element::quad9>   {static const int numDOFsPerElement = 22;};
template<> struct DisTypeToDOFNum<DRT::Element::quad8>   {static const int numDOFsPerElement = 20;};
template<> struct DisTypeToDOFNum<DRT::Element::quad4>   {static const int numDOFsPerElement = 12;};

/// Internal implementation of Taylor-Hood (Q2Q1) fluid element
/*!
This internal class keeps all the working arrays needed to
calculate the Taylor-Hood fluid element.

The idea is to separate the element maintenance (class Fluid2)
from the mathematical contents (this class). Of course there are
different implementations of the Fluid2 element, this is just one
such implementation.

This code supports pressure correction methods as fluid solver
as well as implicit standard fluid solvers.

\author tw
\date 06/09
*/
template <DRT::Element::DiscretizationType distype>
class Fluid2TH: public Fluid2THInterface
{
public:
    explicit Fluid2TH();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of dofs per element
    static const int idofs = DisTypeToDOFNum<distype>::numDOFsPerElement;

    /// for pressure correction method
    /*!
    \brief calculate mass matrix, inverse lumped mass matrix and pressure gradient operator
    */
    virtual int CalcGradPAndMassMatrix(Fluid2*            ele,
            ParameterList&        params,
            DRT::Discretization&    discretization,
            vector<int>&        lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseMatrix& elemat2_epetra,
            Epetra_SerialDenseVector&   elevec1_epetra);

    /// for pressure correction method
    /*!
    \brief solve impulse equation implicitly for intermediate velocity only
    */
    virtual int CalcImpulseEqnImplicit(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseVector& elevec1_epetra,
            RefCountPtr< MAT::Material >   material);

    /// for pressure correction method
    /*!
    \brief solve impulse equation semi-implicitly for intermediate velocity only
    */
    virtual int CalcImpulseEqnSemiImplicit(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseVector& elevec1_epetra,
            RefCountPtr< MAT::Material >   material);

    /// for implicit standard solver
    /*!
    \brief calculate system matrix and RHS
    */
    virtual int CalcSysmatAndResidual(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseMatrix& elemat1_epetra,
            Epetra_SerialDenseVector& elevec1_epetra,
            RefCountPtr< MAT::Material >   material);

    /// necessary for partitioned FSI
    /*!
    calculate RHS only
    */
    virtual int CalcResidual(Fluid2* ele,
            ParameterList& params,
            DRT::Discretization& discretization,
            vector<int> lm,
            Epetra_SerialDenseVector& elevec_epetra,
            RefCountPtr<MAT::Material> mat);

    static Fluid2TH* THImpl(DRT::ELEMENTS::Fluid2* f2);

private:

    void Caltau(Fluid2*                                 ele,
            const LINALG::Matrix<2,iel>&           evelnp,
            const LINALG::Matrix<2,iel>&           fsevelnp,
            const LINALG::Matrix<iel,1>&           edensnp,
            const enum Fluid2::TauType              whichtau,
            Teuchos::RCP<const MAT::Material>       material,
            double&                                 visc,
            const double                            timefac,
            const double                            dt,
            const enum Fluid2::TurbModelAction      turb_mod_action,
            double&                                 Cs,
            double&                                 visceff,
            const enum Fluid2::FineSubgridVisc      fssgv);

    void CalVisc(Teuchos::RCP<const MAT::Material> material,
            double&                 visc);

    void BodyForce(Fluid2*      ele,
            const double time);
private:
    //! artificial viscosity
    double vart_;
    //! node coordinates
    LINALG::Matrix<2,iel> xyze_;
    //! dead load in nodes
    LINALG::Matrix<2,iel> edeadng_;
    //! array for shape functions
    LINALG::Matrix<iel,1> funct_;
    //! array for shape functions (pressure)
    LINALG::Matrix<iel,1>   pfunct_;
    //! array for density-weighted shape functions
    LINALG::Matrix<iel,1> densfunct_;
    //! array for shape functions divided by density
    LINALG::Matrix<iel,1> functdens_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<2,iel> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<3,iel> deriv2_;
    LINALG::Matrix<2,iel>   pderiv_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<2,2>   xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<2,2>   xji_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<2,2>   vderxy_;
    //! global momentum derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<2,2>   mderxy_;
    //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<2,2>   fsvderxy_;

    LINALG::Matrix<2,3>   vderxy2_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<2,iel> derxy_;

    LINALG::Matrix<2,iel> pderxy_;
    //! density-weighted global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<2,iel> densderxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<3,iel> derxy2_;
    //! bodyforce in gausspoint
    LINALG::Matrix<2,1>   bodyforce_;
    //! vector containing all values from previous timelevel n for momentum equation
    LINALG::Matrix<2,1> histmom_;
    //! value from previous timelevel n for continuity equation
    double histcon_;
    //LINALG::Matrix<,1> velino_;
    //! density-weighted velocity vector in gausspoint
    LINALG::Matrix<2,1>   velint_;
    LINALG::Matrix<2,1>   velints_;
    //! fine-scale velocity vector in gausspoint
    LINALG::Matrix<2,1>   fsvelint_;
    //! pressure gradient in gausspoint
    LINALG::Matrix<2,1>   gradp_;
    LINALG::Matrix<3,1>   tau_;
    LINALG::Matrix<4,iel> viscs2_;
    //! linearisation of convection, convective part
    LINALG::Matrix<iel,1> conv_c_;
    LINALG::Matrix<2,1> convvelint_;
    double mdiv_;
    double vdiv_;
    LINALG::Matrix<2,1>   rhsmom_;
    double rhscon_;
    LINALG::Matrix<2,1>   conv_old_;
    LINALG::Matrix<2,1>   visc_old_;
    LINALG::Matrix<2,1>   res_old_;
    LINALG::Matrix<iel,1> conv_resM_;
    LINALG::Matrix<3,2>   xder2_;


};

}
}


#endif /* FLUID2_TH_H_ */

#endif /* CCADISCRET */
#endif /* D_FLUID2 */
