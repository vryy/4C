/*--------------------------------------------------------------------------*/
/*!
\file acou_sol_ele_calc.H

\brief Routines for acoustic elements

<pre>
\level 2

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef ACOU_SOL_ELE_CALC_H
#define ACOU_SOL_ELE_CALC_H

#include "acou_ele_interface.H"
#include "acou_sol_ele.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_inpar/inpar_acou.H"
#include "../drt_fem_general/drt_utils_shapevalues_hdg.H"

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
class AcouSolEleCalc : public AcouEleInterface
{
public:
  //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
  static const unsigned int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

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
      Epetra_SerialDenseVector&     elevec);

  /// projection of function field
  virtual int ProjectField(DRT::ELEMENTS::Acou*                 ele,
                           Teuchos::ParameterList&              params,
                           Epetra_SerialDenseVector&            elevec1,
                           Epetra_SerialDenseVector&            elevec2);

  /// projection of Dirichlet function field
  virtual int ProjectDirichField(DRT::ELEMENTS::Acou*                 ele,
                                 Teuchos::ParameterList&              params,
                                 Epetra_SerialDenseVector&            elevec1);

  /// projection of optical field
  virtual int ProjectOpticalField(DRT::ELEMENTS::Acou*                 ele,
                                  Teuchos::ParameterList&              params,
                                  Epetra_SerialDenseVector&            elevec2);

  /// projection of optical field
  virtual void ElementInit(DRT::ELEMENTS::Acou*                 ele,
                           Teuchos::ParameterList&              params);

  /// map solution values to nodes
  void NodeBasedValues(const Teuchos::RCP<MAT::Material>    &mat,
                       DRT::ELEMENTS::Acou*                 ele,
                       Epetra_SerialDenseVector&            elevec1,
                       bool                                 writestress);

  /// map solution values to nodes
  void NodeBasedPsi(const Teuchos::RCP<MAT::Material>    &mat,
                    DRT::ELEMENTS::Acou*                 ele,
                    Epetra_SerialDenseVector&            elevec1,
                    double                               dt);


  ///Initialize the shape functions and solver to the given element (degree is runtime parameter)
  void InitializeShapes(const DRT::ELEMENTS::Acou* ele);

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


  virtual ~AcouSolEleCalc() {}

  /// Singleton access method
  static AcouSolEleCalc<distype> * Instance( bool create=true );

  /// called upon destruction
  virtual void Done();


private:
  /// private Constructor since we are a Singleton.
  AcouSolEleCalc();

  /// local solver that inverts local problem on an element and can solve with various vectors
  struct LocalSolver
  {
    static const unsigned int nsd_    = AcouSolEleCalc<distype>::nsd_;
    static const unsigned int nfaces_ = AcouSolEleCalc<distype>::nfaces_;

    LocalSolver (const DRT::ELEMENTS::Acou* ele,
                 const DRT::UTILS::ShapeValues<distype> &shapeValues,
                 DRT::UTILS::ShapeValuesFace<distype> &shapeValuesFace,
                 bool completepoly);

    /// add terms corresponding to the absorbing boundary condition
    void ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                            Teuchos::ParameterList&     params,
                            Teuchos::RCP<MAT::Material> & mat,
                            int                         face,
                            Epetra_SerialDenseMatrix    &elemat,
                            Epetra_SerialDenseVector    &elevec);

    void ComputeSource(Teuchos::ParameterList&           params,
                       Epetra_SerialDenseVector          & interiorSourcen,
                       Epetra_SerialDenseVector          & interiorSourcenp);

    /// set up interior matrices
    void ComputeInteriorMatrices();

    /// set up face matrices
    void ComputeFaceMatrices(const int face);

    // condense the local matrix (involving cell velocity gradients, velocities and pressure) into
    // the element matrix for the trace and similarly for the residuals
    void CondenseLocalPart(Epetra_SerialDenseMatrix &elemat,
                           INPAR::ACOU::DynamicType dyna);

    /// compute the residual
    void ComputeResidual(Teuchos::ParameterList&           params,
                         Epetra_SerialDenseVector          & eleVec,
                         Epetra_SerialDenseVector          & interiorGradVeln,
                         Epetra_SerialDenseVector          & interiorVeln,
                         Epetra_SerialDenseVector          & interiorPressn,
                         std::vector<double>               traceVal,
                         INPAR::ACOU::DynamicType          dyna);

    const unsigned int ndofs_;

    const DRT::UTILS::ShapeValues<distype> &shapes_; /// evaluated shape values
    DRT::UTILS::ShapeValuesFace<distype> &shapesface_; /// evaluated shape values

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
    Epetra_SerialDenseMatrix jmat;
    Epetra_SerialDenseMatrix kmat;
    Epetra_SerialDenseMatrix lmat;
    Epetra_SerialDenseMatrix mmat;
    Epetra_SerialDenseMatrix nmat;
  };
  /// updates interior variables and calculates residual
  void UpdateInteriorVariablesAndComputeResidual(Teuchos::ParameterList&           params,
                                                 DRT::ELEMENTS::Acou *             ele,
                                                 Epetra_SerialDenseVector          & elevec,
                                                 double                            dt,
                                                 bool                              errormaps,
                                                 bool                              updateonly);

  /// projects field to new polynomial space
  void ProjectPadapField(DRT::ELEMENTS::Acou *             ele,
                         int                               newdeg);

  /// calls local solver to compute matrices: internal and face
  void ComputeMatrices(const Teuchos::RCP<MAT::Material> &mat,
                       DRT::ELEMENTS::Acou &             ele,
                       double                            dt,
                       INPAR::ACOU::DynamicType          dyna);

  /// reads from global vectors
  void ReadGlobalVectors(DRT::Element     * ele,
                         DRT::Discretization    & discretization,
                         const std::vector<int> & lm);


  /// writes internal field to global vectors
  void FillRestartVectors(DRT::Element           * ele,
                          DRT::Discretization    & discretization);

  /// writes internal field to global vectors
  void ElementInitFromRestart(DRT::Element           * ele,
                              DRT::Discretization    & discretization);

  /// calculate error maps with local postprocessing
  double EstimateError(DRT::ELEMENTS::AcouSol & ele, Epetra_SerialDenseVector & L);

  void EvaluateAll(const int start_func,
                   const double (&xyz)[nsd_],
                   double  &p,
                   double  (&v)[nsd_],
                   double  (&gv)[nsd_*nsd_],
                   double t) const;

  /// local data object
  Teuchos::RCP<DRT::UTILS::ShapeValues<distype> > shapes_;
  Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype> > shapesface_;

  /// local solver object
  Teuchos::RCP<LocalSolver> localSolver_;

  std::vector<double> traceVal_;          /// extracted values from trace solution vector

  /// local values from interior solution vector (gradient, velocity) at n+1
  Epetra_SerialDenseVector interiorGradVelnp_;
  Epetra_SerialDenseVector interiorVelnp_;
  Epetra_SerialDenseVector interiorPressnp_;

  /// time integration scheme
  INPAR::ACOU::DynamicType dyna_;

  bool usescompletepoly_;
};


} // namespace ELEMENTS
} // namespace DRT

#endif /* ACOU_SOL_ELE_CALC_H */
