/*--------------------------------------------------------------------------*/
/*!
\file acou_ele_calc.H

\brief Routines for acoustic elements

<pre>
\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef ACOU_ELE_CALC_H
#define ACOU_ELE_CALC_H

#include "acou_ele_interface.H"
#include "acou_ele.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_inpar/inpar_acou.H"
#include "../drt_fem_general/drt_utils_shapevalues_hdg.H"

namespace DRT
{

namespace ELEMENTS
{

/// Acoustic element implementation

  template<DRT::Element::DiscretizationType distype>
  class AcouEleCalc : public AcouEleInterface
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

    /// projection of optical field
    virtual void ElementInit(DRT::ELEMENTS::Acou*                 ele,
                             Teuchos::ParameterList&              params);

    /// projection of optical field
    virtual void ElementAcouOptInit(DRT::ELEMENTS::Acou*                 ele,
                             Teuchos::ParameterList&              params);

    /// map solution values to nodes
    void NodeBasedValues(DRT::ELEMENTS::Acou*                 ele,
                         Epetra_SerialDenseVector&            elevec1,
                         bool                                 padaptivity);

    /// map solution values to nodes
    void ComputePressureAverage(Epetra_SerialDenseVector& elevec);

    /// map solution values to nodes
    void NodeBasedPsi(Epetra_SerialDenseVector&            elevec1,
                      double                               dt,
                      double                               fac);


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


    virtual ~AcouEleCalc() {}

    /// Singleton access method
    static AcouEleCalc<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();


  private:
    /// private Constructor since we are a Singleton.
    AcouEleCalc();

    /// local solver that inverts local problem on an element and can solve with various vectors
    struct LocalSolver
    {
      static const unsigned int nsd_    = AcouEleCalc<distype>::nsd_;
      static const unsigned int nfaces_ = AcouEleCalc<distype>::nfaces_;

      LocalSolver (DRT::UTILS::ShapeValues<distype> &shapeValues, INPAR::ACOU::DynamicType &dyna);

      void FaceSpecificConstruction(const DRT::ELEMENTS::Acou* ele,
                                    bool completepoly);

      /// compute the residual
      void ComputeResidual(Teuchos::ParameterList&             params,
                           Epetra_SerialDenseVector          & eleVec,
                           Epetra_SerialDenseVector          & interiorVeln,
                           Epetra_SerialDenseVector          & interiorPressn);

      void ComputeSource(Teuchos::ParameterList&           params,
                         Epetra_SerialDenseVector          & interiorSourcen,
                         Epetra_SerialDenseVector          & interiorSourcenp);

      /// add terms corresponding to the absorbing boundary condition
      void ComputeAbsorbingBC(DRT::Discretization &       discretization,
                              DRT::ELEMENTS::Acou*        ele,
                              Teuchos::ParameterList&     params,
                              Teuchos::RCP<MAT::Material> & mat,
                              int                         face,
                              Epetra_SerialDenseMatrix    &elemat,
                              int                          indexstart);

      /// add terms corresponding to the absorbing boundary condition
      void ComputeBoundaryIntegral(DRT::ELEMENTS::Acou*   ele,
                              Teuchos::ParameterList&     params,
                              int                         face);

      /// add terms corresponding to the absorbing boundary condition
      void ComputeSourcePressureMonitor(DRT::ELEMENTS::Acou*        ele,
                              Teuchos::ParameterList&     params,
                              Teuchos::RCP<MAT::Material> & mat,
                              int                         face,
                              Epetra_SerialDenseMatrix    &elemat,
                              Epetra_SerialDenseVector    &elevec,
                              int                          indexstart);

      /// add terms corresponding to the absorbing boundary condition
      void ComputeSourcePressureMonitorLine3D(DRT::ELEMENTS::Acou*        ele,
                                Teuchos::ParameterList&     params,
                                Teuchos::RCP<MAT::Material> & mat,
                                int                         face,
                                Epetra_SerialDenseMatrix    &elemat,
                                Epetra_SerialDenseVector    &elevec,
                                int                          indexstart);

      /// calls local solver to compute matrices: internal and face
      void ComputeMatrices(DRT::Discretization &             discretization,
                           const Teuchos::RCP<MAT::Material> &mat,
                           DRT::ELEMENTS::Acou &             ele,
                           double                            dt,
                           INPAR::ACOU::DynamicType          dyna,
                           bool                              adjoint);

      /// set up interior matrices
      void ComputeInteriorMatrices(double                            dt,
                                   double                            rho,
                                   double                            c);

      /// set up face matrices
      void ComputeFaceMatrices(const int                          face,
                               double                             dt,
                               int                                indexstart,
                               double                             c);

      // condense the local matrix (involving cell velocity gradients, velocities and pressure) into
      // the element matrix for the trace and similarly for the residuals
      void CondenseLocalPart(Epetra_SerialDenseMatrix &elemat);

      void EvaluateFaceAdjoint(double fnodexyz[][nsd_],
                               double values[],
                               int numfnode,
                               const double (&xyz)[nsd_],
                               double &val) const;

      /// projection of function field
      int ProjectField(DRT::ELEMENTS::Acou*                 ele,
                       Teuchos::ParameterList&              params,
                       Epetra_SerialDenseVector&            elevec1,
                       Epetra_SerialDenseVector&            elevec2);

      /// projection of Dirichlet function field
      int ProjectDirichField(DRT::ELEMENTS::Acou*                 ele,
                             Teuchos::ParameterList&              params,
                             Epetra_SerialDenseVector&            elevec1);

      /// projection of optical field
      int ProjectOpticalField(DRT::ELEMENTS::Acou*                 ele,
                                      Teuchos::ParameterList&              params,
                                      Epetra_SerialDenseVector&            elevec2);

      void EvaluateAll(const int start_func,
                       const double (&xyz)[nsd_],
                       double  &p,
                       double  (&v)[nsd_]) const;

      void EvaluateLight(double lightxyz[][nsd_],
                         double values[],
                         int    numnode,
                         const double (&xyz)[nsd_],
                         double  &p,
                         double absorptioncoeff) const;

      // convention: we sort the entries in the matrices the following way:
      // first come the velocity gradients, then the velocities, and finally the pressure
      // we also build the matrix in a block-fashion, keeping the dofs for individual components
      // closest to each other. I.e. the blocks are in 2D for g_00, g_01, g_10, g_11, v_0, v_1, p
      // and similarly for 3D

      const unsigned int ndofs_;

      DRT::UTILS::ShapeValues<distype> &shapes_; /// evaluated shape values
      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype> > shapesface_; /// evaluated shape values

      Epetra_SerialDenseMatrix  Amat;     /// velocity velocity
      Epetra_SerialDenseMatrix  invAmat;  /// inverse of Amat
      Epetra_SerialDenseMatrix  Bmat;     /// pressure gradient
      Epetra_SerialDenseMatrix  Mmat;     /// pressure pressure
      Epetra_SerialDenseMatrix  Dmat;     /// pressure pressure

      Epetra_SerialDenseMatrix  Cmat;     /// velocity trace
      Epetra_SerialDenseMatrix  Emat;     /// pressure trace
      Epetra_SerialDenseMatrix  Gmat;     /// trace trace

      Epetra_SerialDenseMatrix  Hmat;     /// velocity pressure
      Epetra_SerialDenseMatrix  Imat;     /// velocity trace
      Epetra_SerialDenseMatrix  Jmat;     /// pressure trace

      // auxiliary stuff
      Epetra_SerialDenseMatrix  massPart;
      Epetra_SerialDenseMatrix  massPartW;

      INPAR::ACOU::DynamicType &dyna_;
    };

    /// updates interior variables and calculates residual
    void UpdateInteriorVariablesAndComputeResidual(Teuchos::ParameterList&           params,
                                                   DRT::ELEMENTS::Acou &             ele,
                                                   const Teuchos::RCP<MAT::Material> &mat,
                                                   Epetra_SerialDenseVector          & elevec,
                                                   double                            dt,
                                                   bool                              errormaps,
                                                   bool                              updateonly);

    void ComputeGradientContribution(DRT::Discretization &             discretization,
                                     DRT::ELEMENTS::Acou &             ele,
                                     const Teuchos::RCP<MAT::Material> &mat,
                                     int                               offset);

    /// projects field to new polynomial space
    void ProjectPadapField(DRT::ELEMENTS::Acou &             ele,
                           int                               newdeg);

    /// add terms corresponding to the absorbing boundary condition
    void ComputePMonNodeVals(DRT::ELEMENTS::Acou*        ele,
                            Teuchos::ParameterList&     params,
                            Teuchos::RCP<MAT::Material> & mat,
                            int                         face,
                            Epetra_SerialDenseMatrix    &elemat,
                            Epetra_SerialDenseVector    &elevec,
                            int                          indexstart);

    /// reads from global vectors
    void ReadGlobalVectors(DRT::Element     * ele,
                           DRT::Discretization    & discretization,
                           const std::vector<int> & lm,
                           const bool               padaptivity,
                           const int                acouopt);

    /// writes internal field to global vectors
    void FillRestartVectors(DRT::Element           * ele,
                            DRT::Discretization    & discretization);

    /// writes internal field to global vectors
    void ElementInitFromRestart(DRT::Element           * ele,
                                DRT::Discretization    & discretization);

    /// writes internal field to global vectors
    void ElementInitFromAcouOptRestart(DRT::Element    * ele,
                                DRT::Discretization    & discretization,
                                const std::vector<int> & lm);

    /// calculate error maps with local postprocessing
    double EstimateError(DRT::ELEMENTS::Acou & ele, Epetra_SerialDenseVector & p);

    /// local data object
    Teuchos::RCP<DRT::UTILS::ShapeValues<distype> > shapes_;

    /// local solver object
    Teuchos::RCP<LocalSolver> localSolver_;

    std::vector<double> traceVal_;          /// extracted values from trace solution vector

    /// local values from interior solution vector (gradient, velocity) at n+1
    Epetra_SerialDenseVector interiorVelnp_;
    Epetra_SerialDenseVector interiorPressnp_;

    /// time integration scheme
    INPAR::ACOU::DynamicType dyna_;

    bool usescompletepoly_;
  };
} // namespace ELEMENTS
} // namespace DRT


#endif /* ACOU_ELE_CALC_H */
