/*--------------------------------------------------------------------------*/
/*!
\file acou_visc_ele_calc.H

\brief Routines for acoustic elements

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef ACOU_VISC_ELE_CALC_H
#define ACOU_VISC_ELE_CALC_H

#include "acou_ele_interface.H"
#include "acou_visc_ele.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_inpar/inpar_acou.H"

namespace DRT
{
namespace ELEMENTS
{

/// Acoustic element implementation
/*!

  \author schoeder
  \date 02/14
*/

template<DRT::Element::DiscretizationType distype>
class AcouViscEleCalc : public AcouEleInterface
{
public:
  //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
  static const unsigned int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of scalar dofs per cell (quad/hex assumption right now)
  static const unsigned int ndofs_ = DRT::UTILS::FixedPower<DRT::ELEMENTS::AcouVisc::degree+1,nsd_>::value;

  //! number of scalar dofs per face (quad/hex assumption right now)
  static const unsigned int nfdofs_ = DRT::UTILS::FixedPower<DRT::ELEMENTS::AcouVisc::degree+1,nsd_-1>::value;

  ///! number of faces on element
  static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

  virtual int IntegrateShapeFunction(
      DRT::ELEMENTS::Acou*     ele,
      DRT::Discretization&      discretization,
      const std::vector<int>&   lm            ,
      Epetra_SerialDenseVector& elevec1)
  { dserror("Not implemented"); return 1; };

  void ComputeError(
      DRT::ELEMENTS::Acou*          ele,
      Teuchos::ParameterList&       params,
      Teuchos::RCP<MAT::Material>&  mat,
      DRT::Discretization&          discretization,
      const std::vector<int>&       lm,
      Epetra_SerialDenseVector&     elevec);

  /// projection of function field
  virtual int ProjectField(DRT::ELEMENTS::Acou*                 ele,
                           Teuchos::ParameterList&              params,
                           Teuchos::RCP<MAT::Material>&         mat,
                           DRT::Discretization&                 discretization,
                           const std::vector<int>&              lm,
                           Epetra_SerialDenseVector&            elevec1,
                           Epetra_SerialDenseVector&            elevec2);

  /// projection of optical field
  virtual int ProjectOpticalField(DRT::ELEMENTS::Acou*                 ele,
                                  Teuchos::ParameterList&              params,
                                  Teuchos::RCP<MAT::Material>&         mat,
                                  DRT::Discretization&                 discretization,
                                  const std::vector<int>&              lm,
                                  Epetra_SerialDenseVector&            elevec1,
                                  Epetra_SerialDenseVector&            elevec2);

  /// map solution values to nodes
  void NodeBasedValues(DRT::ELEMENTS::Acou*                 ele,
                       DRT::Discretization&                 discretization,
                       const std::vector<int>&              lm,
                       Epetra_SerialDenseVector&            elevec1);

  /// map solution values to nodes
  void NodeBasedPsi(const Teuchos::RCP<MAT::Material>    &mat,
                    DRT::ELEMENTS::Acou*                 ele,
                    DRT::Discretization&                 discretization,
                    const std::vector<int>&              lm,
                    Epetra_SerialDenseVector&            elevec1,
                    double                               dt);

  /// Evaluate the element
  /*!
    Generic virtual interface function. Called via base pointer.
   */
  virtual int Evaluate(DRT::ELEMENTS::Acou*          ele,
                       DRT::Discretization &         discretization,
                       const std::vector<int> &      lm,
                       Teuchos::ParameterList&       params,
                       Teuchos::RCP<MAT::Material> & mat,
                       Epetra_SerialDenseMatrix&     elemat1_epetra,
                       Epetra_SerialDenseMatrix&     elemat2_epetra,
                       Epetra_SerialDenseVector&     elevec1_epetra,
                       Epetra_SerialDenseVector&     elevec2_epetra,
                       Epetra_SerialDenseVector&     elevec3_epetra,
                       bool                          offdiag = false);

  /// Evaluate the element at specified gauss points
  virtual int Evaluate(DRT::ELEMENTS::Acou*                 ele,
                       DRT::Discretization &                discretization,
                       const std::vector<int> &             lm,
                       Teuchos::ParameterList&              params,
                       Teuchos::RCP<MAT::Material> &        mat,
                       Epetra_SerialDenseMatrix&            elemat1_epetra,
                       Epetra_SerialDenseMatrix&            elemat2_epetra,
                       Epetra_SerialDenseVector&            elevec1_epetra,
                       Epetra_SerialDenseVector&            elevec2_epetra,
                       Epetra_SerialDenseVector&            elevec3_epetra,
                       const DRT::UTILS::GaussIntegration & intpoints,
                       bool                                 offdiag = false);


  virtual ~AcouViscEleCalc() {}

  /// Singleton access method
  static AcouViscEleCalc<distype> * Instance( bool create=true );

  /// called upon destruction
  virtual void Done();


private:
  /// private Constructor since we are a Singleton.
  AcouViscEleCalc();

  /// holds shape values on the present element and all the faces
  struct ShapeValues
  {
    static const unsigned int nen_    = AcouViscEleCalc<distype>::nen_;
    static const unsigned int nfn_    = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;
    static const unsigned int nsd_    = AcouViscEleCalc<distype>::nsd_;
    static const unsigned int ndofs_  = AcouViscEleCalc<distype>::ndofs_;
    static const unsigned int nfdofs_ = AcouViscEleCalc<distype>::nfdofs_;
    static const unsigned int nfaces_ = AcouViscEleCalc<distype>::nfaces_;

    ShapeValues ();

    void Evaluate            (const DRT::Element &ele);
    void EvaluateFace        (const DRT::Element &ele,
                              const unsigned int  face);

    /// underlying polynomial space for element interior, created in constructor
    const DRT::UTILS::LagrangeBasis<nsd_> polySpace_;

    /// underlying polynomial space for faces
    const DRT::UTILS::LagrangeBasis<nsd_-1> polySpaceFace_;

    LINALG::Matrix<nen_,ndofs_>   funct;              /// values of mapping shape functions on all quadrature points
    LINALG::Matrix<nfn_,nfdofs_>  functF;             /// values of mapping shape functions on all face quadrature points
    LINALG::Matrix<nsd_,nen_>     deriv;              /// gradients of mapping shape functions
    LINALG::Matrix<nsd_-1,nfn_>   derivF;             /// gradients of mapping shape functions on face
    LINALG::Matrix<nsd_-1,nsd_-1> metricTensor;       /// metric tensor on face
    LINALG::Matrix<nsd_,1>        normal;             /// normal vector
    LINALG::Matrix<nsd_,ndofs_>   xyzreal;            /// coordinates of quadrature points in real space
    LINALG::Matrix<nsd_,nfdofs_>  xyzFreal;           /// coordinates of face quadrature points in real space

    Epetra_SerialDenseMatrix      shfunct;            /// evaluated HDG shape functions on all quadrature points
    Epetra_SerialDenseVector      shfunctAvg;         /// average of shfunctF on cell
    Epetra_SerialDenseMatrix      shderiv;            /// evaluated HDG shape function gradients in unit coordinates
    Epetra_SerialDenseMatrix      shderxy;            /// evaluated HDG shape function gradients in real coordinates
    Epetra_SerialDenseMatrix      shfunctF;           /// evaluated shape functions for HDG face polynomials, permuted to account for face orientation
    Epetra_SerialDenseMatrix      shfunctFNoPermute;  /// evaluated shape functions for HDG face polynomials in natural ordering
    std::vector<Epetra_SerialDenseMatrix> shfunctI;   /// evaluated shape functions on face for interior HDG polynomials
    Epetra_SerialDenseMatrix      normals;            /// normal vectors on a single face for all quadrature points

    LINALG::Matrix<nsd_,1>        xsi;                /// quadrature points
    LINALG::Matrix<nsd_-1,1>      xsiF;               /// face quadrature points
    LINALG::Matrix<nsd_,nsd_>     xjm;                /// Jacobi matrix of transformation
    LINALG::Matrix<nsd_,nsd_>     xji;                /// inverse of Jacobi matrix of transformation
    LINALG::Matrix<nsd_,nen_>     xyze;               /// element nodes
    LINALG::Matrix<nsd_,nfn_>     xyzeF;              /// face nodes
    Epetra_SerialDenseVector      jfac;               /// Jacobian determinant times quadrature weight
    Epetra_SerialDenseVector      jfacF;              /// face Jacobian determinant times quadrature weight

    std::vector<std::vector<int> > faceNodeOrder;     /// numbering of nodes belonging to faces

    Teuchos::RCP<DRT::UTILS::GaussPoints> quadrature_;
    Teuchos::RCP<DRT::UTILS::GaussPoints> fquadrature_;
  };

  /// local solver that inverts local problem on an element and can solve with various vectors
  struct LocalSolver
  {
    static const unsigned int nsd_    = AcouViscEleCalc<distype>::nsd_;
    static const unsigned int ndofs_  = AcouViscEleCalc<distype>::ndofs_;
    static const unsigned int nfdofs_ = AcouViscEleCalc<distype>::nfdofs_;
    static const unsigned int nfaces_ = AcouViscEleCalc<distype>::nfaces_;

    LocalSolver (const ShapeValues &shapeValues);

    const ShapeValues        &shapes_;  /// evaluated shape values

    /// add terms corresponding to the absorbing boundary condition
    void ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                            Teuchos::ParameterList&     params,
                            Teuchos::RCP<MAT::Material> & mat,
                            int                         face,
                            Epetra_SerialDenseMatrix    &elemat,
                            Epetra_SerialDenseVector    &elevec);

    /// add terms corresponding to the absorbing boundary condition
    void ComputeSourcePressureMonitor(DRT::ELEMENTS::Acou*        ele,
                            Teuchos::ParameterList&     params,
                            Teuchos::RCP<MAT::Material> & mat,
                            int                         face,
                            Epetra_SerialDenseMatrix    &elemat,
                            Epetra_SerialDenseVector    &elevec);

    /// add terms corresponding to the absorbing boundary condition
    void ComputeSourcePressureMonitorLine3D(DRT::ELEMENTS::Acou*        ele,
                              Teuchos::ParameterList&     params,
                              Teuchos::RCP<MAT::Material> & mat,
                              int                         face,
                              Epetra_SerialDenseMatrix    &elemat,
                              Epetra_SerialDenseVector    &elevec);

    /// set up interior matrices
    void ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material> &mat,
                                 double                            dt,
                                 INPAR::ACOU::DynamicType          dyna);
    /// set up face matrices
    void ComputeFaceMatrices(const int                          face,
                             const Teuchos::RCP<MAT::Material>& mat,
                             double                             dt);
    // condense the local matrix (involving cell velocity gradients, velocities and pressure) into
    // the element matrix for the trace and similarly for the residuals
    void CondenseLocalPart(Epetra_SerialDenseMatrix &elemat,
                           INPAR::ACOU::DynamicType dyna);

    /// compute the residual
    void ComputeResidual(Teuchos::ParameterList&     params,
                         DRT::ELEMENTS::Acou &                ele,
                         Epetra_SerialDenseVector          & eleVec,
                         Epetra_SerialDenseVector          & interiorGradVeln,
                         Epetra_SerialDenseVector          & interiorVeln,
                         Epetra_SerialDenseVector          & interiorPressn,
                         Epetra_SerialDenseVector          & interiorDensn,
                         std::vector<double>               traceVal,
                         INPAR::ACOU::DynamicType          dyna);


    /// evaluates the source term for the adjoint integration!
    void EvaluateSourceAdjoint(Teuchos::ParameterList&     params,
                               DRT::ELEMENTS::Acou &                ele,
                               Epetra_SerialDenseVector          & sourcevec);

    /// evaluates source term for given coordinates
    void EvaluateFaceAdjoint(double fnodexyz[][nsd_],
                             double values[],
                             int numfnode,
                             const double (&xyz)[nsd_],
                             double &val) const;

    Epetra_SerialDenseMatrix amat;
    Epetra_SerialDenseMatrix invamat;
    Epetra_SerialDenseMatrix bmat;
    Epetra_SerialDenseMatrix cmat;
    Epetra_SerialDenseMatrix dmat;
    Epetra_SerialDenseMatrix emat;
    Epetra_SerialDenseMatrix ehatmat;
    Epetra_SerialDenseMatrix fmat;
    Epetra_SerialDenseMatrix gmat;
    Epetra_SerialDenseMatrix hmat;
    Epetra_SerialDenseMatrix imat;
    Epetra_SerialDenseMatrix invimat;
    Epetra_SerialDenseMatrix jmat;
    Epetra_SerialDenseMatrix kmat;
    Epetra_SerialDenseMatrix lmat;
    Epetra_SerialDenseMatrix lhatmat;
    Epetra_SerialDenseMatrix mmat;
    Epetra_SerialDenseMatrix nmat;
    Epetra_SerialDenseMatrix omat;
    Epetra_SerialDenseMatrix pmat;

  };
  /// updates interior variables and calculates residual
  void UpdateInteriorVariablesAndComputeResidual(DRT::Discretization &             discretization,
                                                 Teuchos::ParameterList&           params,
                                                 DRT::ELEMENTS::Acou &             ele,
                                                 const Teuchos::RCP<MAT::Material> &mat,
                                                 Epetra_SerialDenseVector          & elevec,
                                                 double                            dt,
                                                 bool                              errormaps,
                                                 bool                              updateonly);

  /// calls local solver to compute matrices: internal and face
  void ComputeMatrices(const Teuchos::RCP<MAT::Material> &mat,
                       DRT::ELEMENTS::Acou &             ele,
                       double                            dt,
                       INPAR::ACOU::DynamicType          dyna,
                       bool                              adjoint);


  /// add terms corresponding to the absorbing boundary condition
  void ComputePMonNodeVals(DRT::ELEMENTS::Acou*        ele,
                          Teuchos::ParameterList&     params,
                          Teuchos::RCP<MAT::Material> & mat,
                          int                         face,
                          Epetra_SerialDenseMatrix    &elemat,
                          Epetra_SerialDenseVector    &elevec);

  /// reads from global vectors
  void ReadGlobalVectors(const DRT::Element     & ele,
                         DRT::Discretization    & discretization,
                         const std::vector<int> & lm);

  void EvaluateAll(const int start_func,
                   const double (&xyz)[nsd_],
                   double  &p) const;

  void EvaluateLight(double lightxyz[][nsd_],
                     double values[],
                     int    numnode,
                     const double (&xyz)[nsd_],
                     double  &p,
                     double absorptioncoeff) const;

  /// local data object
  ShapeValues shapes_;

  /// local solver object
  LocalSolver localSolver_;

  std::vector<double> traceVal_;          /// extracted values from trace solution vector
  std::vector<double> traceValm_;         /// extracted values from old trace

  /// local values from interior solution vector (velocity gradients, velocities, pressure) at n+1
  std::vector<double> interiorValnp_;
  /// local values from interior solution vector (velocity gradients, velocities, pressure) at n
  std::vector<double> interiorValn_;
  /// local values from interior solution vector (velocity gradients, velocities, pressure) at n-1
  std::vector<double> interiorValnm_;
  /// local values from interior solution vector (velocity gradients, velocities, pressure) at n-2
  std::vector<double> interiorValnmm_;
  /// local values from interior solution vector (velocity gradients, velocities, pressure) at n-2
  std::vector<double> interiorValnmmm_;

  /// local values from interior solution vector (gradient, velocity) at n+1
  Epetra_SerialDenseVector interiorGradVelnp_;
  Epetra_SerialDenseVector interiorVelnp_;
  Epetra_SerialDenseVector interiorPressnp_;
  Epetra_SerialDenseVector interiorDensnp_;

  /// local values from interior solution vector (gradient, velocity) at n
  Epetra_SerialDenseVector interiorGradVeln_;
  Epetra_SerialDenseVector interiorVeln_;
  Epetra_SerialDenseVector interiorPressn_;
  Epetra_SerialDenseVector interiorDensn_;

  /// local values from interior solution vector (gradient, velocity) at n-1
  Epetra_SerialDenseVector interiorGradVelnm_;
  Epetra_SerialDenseVector interiorVelnm_;
  Epetra_SerialDenseVector interiorPressnm_;
  Epetra_SerialDenseVector interiorDensnm_;

  /// local values from interior solution vector (gradient, velocity) at n-2
  Epetra_SerialDenseVector interiorGradVelnmm_;
  Epetra_SerialDenseVector interiorVelnmm_;
  Epetra_SerialDenseVector interiorPressnmm_;
  Epetra_SerialDenseVector interiorDensnmm_;

  /// local values from interior solution vector (gradient, velocity) at n-3
  Epetra_SerialDenseVector interiorGradVelnmmm_;
  Epetra_SerialDenseVector interiorVelnmmm_;
  Epetra_SerialDenseVector interiorPressnmmm_;
  Epetra_SerialDenseVector interiorDensnmmm_;

  /// time integration scheme
  INPAR::ACOU::DynamicType dyna_;

};


} // namespace ELEMENTS
} // namespace DRT

#endif /* ACOU_VISC_ELE_CALC_H */
