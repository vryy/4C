/*--------------------------------------------------------------------------*/
/*!
\file acou_ele_calc.H

\brief Routines for acoustic elements

<pre>
Maintainer: Svenja Schoeder
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

    /// projection of optical field
    virtual void ElementInit(DRT::ELEMENTS::Acou*                 ele,
                             Teuchos::ParameterList&              params);

    /// map solution values to nodes
    void NodeBasedValues(DRT::ELEMENTS::Acou*                 ele,
                         DRT::Discretization&                 discretization,
                         const std::vector<int>&              lm,
                         Epetra_SerialDenseVector&            elevec1,
                         bool                                 padaptivity);

    /// map solution values to nodes
    void NodeBasedPsi(const Teuchos::RCP<MAT::Material>    &mat,
                      DRT::ELEMENTS::Acou*                 ele,
                      DRT::Discretization&                 discretization,
                      const std::vector<int>&              lm,
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

      LocalSolver (const DRT::ELEMENTS::Acou* ele,
                   const DRT::UTILS::ShapeValues<distype> &shapeValues,
                   DRT::UTILS::ShapeValuesFace<distype> &shapeValuesFace,
                   bool completepoly);

      /// compute the residual
      void ComputeResidual(Epetra_SerialDenseVector          & eleVec,
                           Epetra_SerialDenseVector          & interiorVeln,
                           Epetra_SerialDenseVector          & interiorPressn,
                           std::vector<double>               traceVal,
                           INPAR::ACOU::DynamicType          dyna);

      /// add terms corresponding to the absorbing boundary condition
      void ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                              Teuchos::ParameterList&     params,
                              Teuchos::RCP<MAT::Material> & mat,
                              int                         face,
                              Epetra_SerialDenseMatrix    &elemat,
                              Epetra_SerialDenseVector    &elevec,
                              int                          indexstart);

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

      /// set up interior matrices
      void ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material> &mat,
                                   double                            dt,
                                   INPAR::ACOU::DynamicType          dyna);

      /// set up face matrices
      void ComputeFaceMatrices(const int                          face,
                               const Teuchos::RCP<MAT::Material>& mat,
                               double                             dt,
                               int                                indexstart);

      // condense the local matrix (involving cell velocity gradients, velocities and pressure) into
      // the element matrix for the trace and similarly for the residuals
      void CondenseLocalPart(Epetra_SerialDenseMatrix &elemat,
                             INPAR::ACOU::DynamicType dyna);

      void EvaluateFaceAdjoint(double fnodexyz[][nsd_],
                               double values[],
                               int numfnode,
                               const double (&xyz)[nsd_],
                               double &val) const;

      // convention: we sort the entries in the matrices the following way:
      // first come the velocity gradients, then the velocities, and finally the pressure
      // we also build the matrix in a block-fashion, keeping the dofs for individual components
      // closest to each other. I.e. the blocks are in 2D for g_00, g_01, g_10, g_11, v_0, v_1, p
      // and similarly for 3D

      const unsigned int ndofs_;

      const DRT::UTILS::ShapeValues<distype> &shapes_; /// evaluated shape values
      DRT::UTILS::ShapeValuesFace<distype> &shapesface_; /// evaluated shape values

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

    /// projects field to new polynomial space
    void ProjectPadapField(DRT::ELEMENTS::Acou &             ele,
                           int                               newdeg);

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
                            Epetra_SerialDenseVector    &elevec,
                            int                          indexstart);

    /// reads from global vectors
    void ReadGlobalVectors(DRT::Element     * ele,
                           DRT::Discretization    & discretization,
                           const std::vector<int> & lm,
                           const bool               padaptivity);

    /// writes internal field to global vectors
    void FillRestartVectors(DRT::Element           * ele,
                            DRT::Discretization    & discretization);

    /// writes internal field to global vectors
    void ElementInitFromRestart(DRT::Element           * ele,
                                DRT::Discretization    & discretization);

    /// reads from global vectors
    void VectorHandling(DRT::ELEMENTS::Acou*     ele,
                        Teuchos::ParameterList&       params,
                        double& dt);

    void EvaluateAll(const int start_func,
                     const double (&xyz)[nsd_],
                     double  &p,
                     double rho) const;

    void EvaluateLight(double lightxyz[][nsd_],
                       double values[],
                       int    numnode,
                       const double (&xyz)[nsd_],
                       double  &p,
                       double rho,
                       double absorptioncoeff) const;

    /// calculate error maps
    double CalculateError(DRT::ELEMENTS::Acou & ele,Epetra_SerialDenseMatrix & h, Epetra_SerialDenseVector & rhs, Epetra_SerialDenseVector & p);

    /// local data object
    Teuchos::RCP<DRT::UTILS::ShapeValues<distype> > shapes_;
    Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype> > shapesface_;

    /// local solver object
    Teuchos::RCP<LocalSolver> localSolver_;

    std::vector<double> traceVal_;          /// extracted values from trace solution vector
    std::vector<double> traceValm_;         /// extracted values from old trace

    /// local values from interior solution vector (gradient, velocity) at n+1
    Epetra_SerialDenseVector interiorVelnp_;
    Epetra_SerialDenseVector interiorPressnp_;

    /// local values from interior solution vector (gradient, velocity) at n
    Epetra_SerialDenseVector interiorVeln_;
    Epetra_SerialDenseVector interiorPressn_;

    /// local values from interior solution vector (gradient, velocity) at n-1
    Epetra_SerialDenseVector interiorVelnm_;
    Epetra_SerialDenseVector interiorPressnm_;

    /// local values from interior solution vector (gradient, velocity) at n-2
    Epetra_SerialDenseVector interiorVelnmm_;
    Epetra_SerialDenseVector interiorPressnmm_;

    /// local values from interior solution vector (gradient, velocity) at n-3
    Epetra_SerialDenseVector interiorVelnmmm_;
    Epetra_SerialDenseVector interiorPressnmmm_;

    /// time integration scheme
    INPAR::ACOU::DynamicType dyna_;

    bool usescompletepoly_;
  };
} // namespace ELEMENTS
} // namespace DRT


#endif /* ACOU_ELE_CALC_H */
