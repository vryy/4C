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

    //! number of scalar dofs per cell (quad/hex assumption right now)
    static const unsigned int ndofs_ = DRT::UTILS::FixedPower<DRT::ELEMENTS::Acou::degree+1,nsd_>::value;

    //! number of scalar dofs per face (quad/hex assumption right now)
    static const unsigned int nfdofs_ = DRT::UTILS::FixedPower<DRT::ELEMENTS::Acou::degree+1,nsd_-1>::value;

    ///! number of faces on element
    static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

    virtual int IntegrateShapeFunction(
        DRT::ELEMENTS::Acou*     ele,
        DRT::Discretization&      discretization,
        const std::vector<int>&   lm            ,
        Epetra_SerialDenseVector& elevec1)
    { dserror("Not implemented"); return 1; };

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

    /// holds shape values on the present element and all the faces
    struct ShapeValues
    {
      static const unsigned int nen_    = AcouEleCalc<distype>::nen_;
      static const unsigned int nfn_    = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;
      static const unsigned int nsd_    = AcouEleCalc<distype>::nsd_;
      static const unsigned int ndofs_  = AcouEleCalc<distype>::ndofs_;
      static const unsigned int nfdofs_ = AcouEleCalc<distype>::nfdofs_;
      static const unsigned int nfaces_ = AcouEleCalc<distype>::nfaces_;

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
      static const unsigned int nsd_    = AcouEleCalc<distype>::nsd_;
      static const unsigned int ndofs_  = AcouEleCalc<distype>::ndofs_;
      static const unsigned int nfdofs_ = AcouEleCalc<distype>::nfdofs_;
      static const unsigned int nfaces_ = AcouEleCalc<distype>::nfaces_;

      LocalSolver (const ShapeValues &shapeValues);

      /// compute the residual
      void ComputeResidual(Epetra_SerialDenseVector          & eleVec,
                           const Teuchos::RCP<MAT::Material> &mat,
                           Epetra_SerialDenseVector          & interiorGradnp,
                           Epetra_SerialDenseVector          & interiorVelnp,
                           std::vector<double>               traceVal,
                           INPAR::ACOU::DynamicType          dyna);
      /// compute the residual for adjoint time integration
      void ComputeResidualAdjoint(Epetra_SerialDenseVector          & eleVec,
                                  const Teuchos::RCP<MAT::Material> &mat,
                                  Epetra_SerialDenseVector          & interiorGradnp,
                                  Epetra_SerialDenseVector          & interiorVelnp,
                                  std::vector<double>               traceVal,
                                  INPAR::ACOU::DynamicType dyna);

      /// add terms corresponding to the absorbing boundary condition
      void ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                              Teuchos::ParameterList&     params,
                              Teuchos::RCP<MAT::Material> & mat,
                              int                         face,
                              Epetra_SerialDenseMatrix    &elemat,
                              Epetra_SerialDenseVector    &elevec);

      /// integrate the objective function over the given boundary
      void ComputeObjfIntegral(DRT::ELEMENTS::Acou*        ele,
                              Teuchos::ParameterList&     params,
                              Teuchos::RCP<MAT::Material> & mat,
                              int                         face,
                              double                      dt);

      /// set up interior matrices
      void ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material> &mat, double dt);

      /// set up face matrices
      void ComputeFaceMatrices(const int                          face,
                               const Teuchos::RCP<MAT::Material>& mat,
                               double                             dt);

      // condense the local matrix (involving cell velocity gradients, velocities and pressure) into
      // the element matrix for the trace and similarly for the residuals
      void CondenseLocalPart(Epetra_SerialDenseMatrix &elemat,
                             const Teuchos::RCP<MAT::Material>& mat,
                             INPAR::ACOU::DynamicType dyna);
      void CondenseLocalPartAdjoint(Epetra_SerialDenseMatrix &elemat,
                                    const Teuchos::RCP<MAT::Material>& mat,
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

      const ShapeValues        &shapes_;  /// evaluated shape values

      Epetra_SerialDenseMatrix  Amat;     /// gradient gradient
      Epetra_SerialDenseMatrix  invAmat;  /// inverse of Amat
      Epetra_SerialDenseMatrix  Bmat;     /// velocity gradient
      Epetra_SerialDenseMatrix  Dmat;     /// velocity velocity
      Epetra_SerialDenseMatrix  Mmat;     /// velocity velocity
      Epetra_SerialDenseMatrix  DmMmat;   /// velocity velocity

      Epetra_SerialDenseMatrix  Cmat;     /// gradient trace
      Epetra_SerialDenseMatrix  Emat;     /// velocity trace
      Epetra_SerialDenseMatrix  Gmat;     /// trace trace

    };

    /// updates interior variables and calculates residual
    void UpdateInteriorVariablesAndComputeResidual(DRT::Discretization &             discretization,
                                                   Teuchos::ParameterList&           params,
                                                   DRT::ELEMENTS::Acou &             ele,
                                                   Epetra_SerialDenseVector          & elevec,
                                                   const Teuchos::RCP<MAT::Material> &mat,
                                                   double                            dt,
                                                   bool                              errormaps,
                                                   bool                              updateonly);
    void UpdateInteriorVariablesAndComputeResidualAdjoint(DRT::Discretization &          discretization,
                                                       DRT::ELEMENTS::Acou &             ele,
                                                       Epetra_SerialDenseVector          & elevec,
                                                       const Teuchos::RCP<MAT::Material> &mat,
                                                       bool                              updateonly);

    /// reads from global vectors
    void ReadGlobalVectors(const DRT::Element     & ele,
                           DRT::Discretization    & discretization,
                           const std::vector<int> & lm);

    void EvaluateAll(const int start_func,
                     const double (&xyz)[nsd_],
                     double (&u)[nsd_],
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
    ShapeValues shapes_;

    /// local solver object
    LocalSolver localSolver_;

    LINALG::Matrix<nsd_,nen_> ebofoaf_;     /// body force (see fluid_ele_calc.cpp)
    LINALG::Matrix<nsd_,nen_> eprescpgaf_;  /// pressure gradient body force
    LINALG::Matrix<nen_,1>    escabofoaf_;  /// scalar body force for loma

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
    Epetra_SerialDenseVector interiorGradnp_;
    Epetra_SerialDenseVector interiorVelnp_;

    /// local values from interior solution vector (gradient, velocity) at n
    Epetra_SerialDenseVector interiorGradn_;
    Epetra_SerialDenseVector interiorVeln_;

    /// local values from interior solution vector (gradient, velocity) at n-1
    Epetra_SerialDenseVector interiorGradnm_;
    Epetra_SerialDenseVector interiorVelnm_;

    /// local values from interior solution vector (gradient, velocity) at n-2
    Epetra_SerialDenseVector interiorGradnmm_;
    Epetra_SerialDenseVector interiorVelnmm_;

    /// local values from interior solution vector (gradient, velocity) at n-3
    Epetra_SerialDenseVector interiorGradnmmm_;
    Epetra_SerialDenseVector interiorVelnmmm_;

    /// time integration scheme
    INPAR::ACOU::DynamicType dyna_;
  };
} // namespace ELEMENTS
} // namespace DRT


#endif /* ACOU_ELE_CALC_H */
