
/*----------------------------------------------------------------------*/
/*!
\file artery_lin_exp.H

\brief Internal implementation of artery_lin_exp element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_ARTNET
#ifdef CCADISCRET

#ifndef ARTERY_LIN_EXP_H
#define ARTERY_LIN_EXP_H

#include "artery.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_fixedsizematrix.H"


namespace DRT
{
namespace ELEMENTS
{
  /// Interface base class for artery_lin_exp
  /*!
    This class exists to provide a common interface for all template
    versions of artery_lin_exp. The only function
    this class actually defines is Expl, which returns a pointer to
    the appropriate version of artery_lin_exp.
   */
  class ArteryExpInterface
  {
  public:
    /// Empty constructor
    ArteryExpInterface() {}
    /// Empty destructor
    virtual ~ArteryExpInterface() {}
    /// Evaluate the element
    virtual int Evaluate(Artery*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat) = 0;

    virtual void Initial(Artery*                                  ele,
                         ParameterList&                           params,
                         DRT::Discretization&                     discretization,
                         vector<int>&                             lm,
                         Teuchos::RCP<const MAT::Material>        material)=0;

    virtual void EvaluateTerminalBC(Artery*                      ele,
                                    ParameterList&               params,
                                    DRT::Discretization&         discretization,
                                    vector<int>&                 lm,
                                    Epetra_SerialDenseMatrix&    elemat1,
                                    Epetra_SerialDenseVector&    elevec1,
                                    RefCountPtr<MAT::Material>   mat) = 0;
    /// Internal implementation class for artery element
    static ArteryExpInterface* Expl(DRT::ELEMENTS::Artery* art);

  };


  /// Internal artery implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the artery element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    \author ismail
    \date 01/09
  */

  template <DRT::Element::DiscretizationType distype>
  class ArteryLinExp: public ArteryExpInterface
  {
  public:
    /// Constructor
    explicit ArteryLinExp();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /// Evaluate
    /*!
      The evaluate function for the general artery case.
     */
    virtual int Evaluate(Artery*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat);

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param eqnp             (i) nodal volumetric flow rate at n+1
      \param evelnp           (i) nodal velocity at n+1
      \param eareanp          (i) nodal cross-sectional area at n+1
      \param eprenp           (i) nodal pressure at n+1
      \param estif            (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) artery material/dimesion
      \param time             (i) current simulation time
      \param dt               (i) timestep
      */
    void Sysmat(Artery*                                  ele,
                const LINALG::Matrix<iel,1>&             eqnp,
                const LINALG::Matrix<iel,1>&             eareanp,
                LINALG::Matrix<2*iel,2*iel>&             estif,
                LINALG::Matrix<2*iel,    1>&             eforce,
                Teuchos::RCP<const MAT::Material>        material,
                double                                   dt);

    bool SolveRiemann(Artery*                             ele,
                      ParameterList&                      params,
                      const LINALG::Matrix<iel,1>&        eqnp,
                      const LINALG::Matrix<iel,1>&        eareanp,
                      Teuchos::RCP<const MAT::Material>   material,
                      double                              dt);

    virtual void EvaluateTerminalBC(Artery*                      ele,
                                    ParameterList&               params,
                                    DRT::Discretization&         disctretization,
                                    vector<int>&                 lm,
                                    Epetra_SerialDenseMatrix&    elemat1,
                                    Epetra_SerialDenseVector&    elevec1,
                                    RefCountPtr<MAT::Material>   mat);

    /*!
      \brief get the initial values of the degrees of freedome at the node

      \param ele              (i) the element those matrix is calculated
      \param eqnp             (i) nodal volumetric flow rate at n+1
      \param evelnp           (i) nodal velocity at n+1
      \param eareanp          (i) nodal cross-sectional area at n+1
      \param eprenp           (i) nodal pressure at n+1
      \param estif            (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) artery material/dimesion
      \param time             (i) current simulation time
      \param dt               (i) timestep
      */
    virtual void Initial(Artery*                                  ele,
                         ParameterList&                           params,
                         DRT::Discretization&                     discretization,
                         vector<int>&                             lm,
                         Teuchos::RCP<const MAT::Material>        material);

    /*!
     \Essential functions to compute the results of essentail matrices 
    */
#if 0
    LINALG::Matrix<2,1> F    (double xi);
    LINALG::Matrix<2,1> B    (double xi);
    LINALG::Matrix<2,2> Bu   (double xi);
    LINALG::Matrix<2,2> H    (double xi);
    LINALG::Matrix<2,1> dFdxi(double xi);
    LINALG::Matrix<2,1> FLW  (double xi, double dt);
    LINALG::Matrix<2,1> BLW  (double xi, double dt);
#endif
    
  private:
    //! node coordinates
    //    LINALG::Matrix<3,iel> xyze_;
    //! node 1D coordinates on the "s" axis
    LINALG::Matrix<iel,1> se_;
    //! array for shape functions
    LINALG::Matrix<iel,1> funct_;
    //! array for shape function derivatives w.r.t s
    LINALG::Matrix<1,iel> deriv_;
    //! transposed array for shape function derivatives w.r.t s
    LINALG::Matrix<iel,1> tderiv_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<1,1> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<1,1> xji_;
    //! nodal volumetric flow rate at time step "n"
    LINALG::Matrix<iel,1> qn_;
    //! nodal cross-sectional area at time step "n"
    LINALG::Matrix<iel,1> an_;
    //! global derivatives of shape functions w.r.t s
    LINALG::Matrix<2,iel> derxy_;
    //! vector containing the initial cross-sectional area at the element nodes
    LINALG::Matrix<iel,1> area0_;
    //! vector containing the initial thickness at the element nodes
    LINALG::Matrix<iel,1> th_;
    //! vector containing the initial Youngs modulus at the element nodes
    LINALG::Matrix<iel,1> young_;
    //! vector containing the fixed external pressure
    LINALG::Matrix<iel,1> pext_;
    //! Element Length
    //    double L_;
    //! Initial Cross-sectional Area of a the artery at element node 1
    double area1_0_;
    //! Initial Cross-sectional Area of a the artery at element node 2
    double area2_0_;
    //! Artery Thickness at element node 1
    double th1_;
    //! Artery Thickness at element node 2
    double th2_;
    //! Artery's Youngs Modulus of Elasticity at element node 1
    double young1_;
    //! Artery's Youngs Modulus of Elasticity at element node 2
    double young2_;
    //! Blood viscosity
    double visc_;
    //! Blood density
    double dens_;
    //! Poisson's ration of arterial wall
    double nue_;
    //! Viscosity parameter
    double Kr_;
    //! Time indipendent external artery wall pressure at terminal 1
    double pext1_;
    //! Time indipendent external artery wall pressure at terminal 2
    double pext2_;
    //! Forward Characteristic wave at time step n+1 at terminal 1
    double Wf1np_;
    //! Backward Characteristic wave at time step n+1 at terminal 1
    double Wb1np_;
    //! Forward Characteristic wave at time step n+1 at terminal 2
    double Wf2np_;
    //! Backward Characteristic wave at time step n+1 at terminal 2
    double Wb2np_;
  };

}
}

#endif

#endif
#endif 
