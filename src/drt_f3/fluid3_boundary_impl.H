/*----------------------------------------------------------------------*/
/*!
\file fluid3_boundary_impl.H

\brief Internal implementation of fluid3 boundary elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#ifndef FLUID3_BOUNDARY_IMPL_H
#define FLUID3_BOUNDARY_IMPL_H

#include "fluid3.H"
#include "fluid3_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_scatra.H"



namespace DRT
{
namespace ELEMENTS
{

/// Interface base class for Fluid3BoundaryImpl
/*!
  This class exists to provide a common interface for all template
  versions of Fluid3BoundaryImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of Fluid3BoundaryImpl.
 */
class Fluid3BoundaryImplInterface
{

public:
  /// Empty constructor
  Fluid3BoundaryImplInterface() {}
  /// Empty destructor
  virtual ~Fluid3BoundaryImplInterface() {}
  /// Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in Fluid3BoundaryImpl.
   */
  virtual int Evaluate(DRT::ELEMENTS::Fluid3Boundary* ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  /// Evaluate a Neumann boundary condition
  /*!
    This class does not provide a definition for this function, it
    must be defined in Fluid3BoundaryImpl.
   */
  virtual int EvaluateNeumann(DRT::ELEMENTS::Fluid3Boundary*      ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1_epetra,
                              Epetra_SerialDenseMatrix*  elemat1
                              ) = 0;

  /*!
    \brief Evaluate the surface elements weak Dirichlet boundary
           conditions

    This class does not provide a definition for this function, it
    is defined in the implementation.

    \author gammi (gamnitzer@lnm.mw.tum.de), 11/09

    \param surfele             (in): surface element in geometry
                                     associated with weak dbc
    \param params              (in): parameters (time integration
                                     parent lms, condition)
    \param discretization      (in): the discretization
    \param lm                  (in): the location vector of surfele
    \param elemat_epetra   (in/out): element matrix (due to wdbc)
    \param elevec_epetra   (in/out): element force vector (due to wdbc)

  */
  /*
  virtual int EvaluateWeakDBC(
    Fluid3Boundary*             surfele       ,
    ParameterList&             params        ,
    DRT::Discretization&       discretization,
    vector<int>&               lm            ,
    Epetra_SerialDenseMatrix&  elemat_epetra ,
    Epetra_SerialDenseVector&  elevec_epetra ) = 0;
*/
  /// Internal implementation class for Fluid3Boundary elements
  static Fluid3BoundaryImplInterface* Impl(const DRT::Element* ele);

};

  /// Internal Fluid3Boundary element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3Boundary element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid3Boundary)
    from the mathematical contents (this class). Of course there are
    different implementations of the Fluid3Boundary element, this is just one
    such implementation.

    The Fluid3Boundary element will allocate exactly one object of this class
    for all Fluid3Boundary elements with the same number of nodes in the mesh.
    This allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>History</h3>

    The implementation here is the standard convection-diffusion element
    capable of dealing with systems of transported scalars.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version.

    \author gjb
    \date 08/08
    \author ehrl
    \date 03/10
  */
template<DRT::Element::DiscretizationType distype>
class Fluid3BoundaryImpl: public Fluid3BoundaryImplInterface
{
public:
  /// Constructor
  Fluid3BoundaryImpl(int numdofpernode);

  //! number of element nodes
  static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions of the Fluid3Boundary element
  static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  static const int nsd_ = bdrynsd_+1;

  //! Evaluate
  virtual int Evaluate(DRT::ELEMENTS::Fluid3Boundary* ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       );

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(DRT::ELEMENTS::Fluid3Boundary*    ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1_epetra,
                              Epetra_SerialDenseMatrix*  elemat1
                              );

private:

    //! number of dof's per node
    const int numdofpernode_;

    //! number of space dimensions of the parent element

    void CreateLinesTri(const int& nline,
                        const int& nnode);

    void CreateLinesQuad(const int& nline,
                         const int& nnode);

    /*!
    \brief integrate shapefunction over surface element

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,
    */

    virtual void IntegrateShapeFunction(
                                    DRT::ELEMENTS::Fluid3Boundary*    ele,
                                    ParameterList&                    params,
                                    DRT::Discretization&              discretization,
                                    vector<int>&                      lm,
                                    Epetra_SerialDenseVector&         elevec1,
                                    const std::vector<double>&        edispnp);

    /*!
      \brief Calculate node normals

      Acc. to Wall (7.13), surface integral.

      \param elevec1 (out)      : Weighted element normal for assembly
      \param edispnp (in)       : Displacement-vector
    */
    void ElementNodeNormal(
                         DRT::ELEMENTS::Fluid3Boundary*   ele,
                         ParameterList&                 params,
                         DRT::Discretization&           discretization,
                         vector<int>&                   lm,
                         Epetra_SerialDenseVector&      elevec1,
                         const std::vector<double>&     edispnp);
    /*!
      \brief Calculate mean curvature H. Interpolate the results to achieve better
      results in the surface tension algorithm (c0 field -> c1 field).

      \param elevec1  (out)     : Nodal values of mean curvature
      \param edispnp  (in)      : Displacement-vector
      \param enormals (in)      : Node normals
    */
    void ElementMeanCurvature(
                              DRT::ELEMENTS::Fluid3Boundary*    ele,
                              ParameterList&                    params,
                              DRT::Discretization&              discretization,
                              vector<int>&                      lm,
                              Epetra_SerialDenseVector&         elevec1,
                              const std::vector<double>&        edispnp,
                              std::vector<double>&              enormals);


    /*!
      \brief Integrate surface tension

      \param elevec1   (out)    : RHS contribution of surface tension effect
      \param edispnp   (in)     : Displacement-vector
      \param enormals  (in)     : Node normals
      \param curvature (in)     : Interpolated curvature
    */
    void ElementSurfaceTension(
                               DRT::ELEMENTS::Fluid3Boundary*   ele,
                               ParameterList&                   params,
                               DRT::Discretization&             discretization,
                               vector<int>&                     lm,
                               Epetra_SerialDenseVector&        elevec1,
                               const std::vector<double>&       edispnp,
                               std::vector<double>&             enormals,
                               std::vector<double>&             ecurvature);

    /*!
    brief integrate elemental areas over a surface

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    */
    void AreaCaculation(DRT::ELEMENTS::Fluid3Boundary*    ele,
                        ParameterList&                    params,
                        DRT::Discretization&              discretization,
                        vector<int>&                      lm);

   /*!
    brief integrate elemental flow rates over a surface

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    */
    void FlowRateParameterCalculation(
                    DRT::ELEMENTS::Fluid3Boundary*    ele,
                    ParameterList&                    params,
                    DRT::Discretization&              discretization,
                    vector<int>&                      lm);

    /*!
     brief integrate elemental flow rates over a surface

     \param params (in/out)    : ParameterList for communication between control routine
                                 and elements
     \param discretization (in): A reference to the underlying discretization
     \param lm (in)            : location vector of this element
     */
     void ComputeFlowRate(
                     DRT::ELEMENTS::Fluid3Boundary*    ele,
                     ParameterList&                    params,
                     DRT::Discretization&              discretization,
                     vector<int>&                      lm,
                     Epetra_SerialDenseVector&         elevec1);

    /*!
    \brief determine elemental flow rate and its derivatives wrt dofs

    \param params (in/out)    : parameterList for communication between control routine
                                and elements
    \param discretization (in): reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : 2nd derivative of flowrates wrt
                                velocities and displacements
    \param elemat2 (out)      : 2nd derivative of flowrates wrt displacements
    \param elevec1 (out)      : derivative of flowrates wrt velocities
    \param elevec2 (out)      : derivative of flowrate wrt displacements
    \param elevec3 (out)      : flowrate

    */
    void FlowRateDeriv(
                       DRT::ELEMENTS::Fluid3Boundary*   ele,
                       ParameterList&                   params,
                       DRT::Discretization&             discretization,
                       vector<int>&                     lm,
                       Epetra_SerialDenseMatrix&        elemat1,
                       Epetra_SerialDenseMatrix&        elemat2,
                       Epetra_SerialDenseVector&        elevec1,
                       Epetra_SerialDenseVector&        elevec2,
                       Epetra_SerialDenseVector&        elevec3);


    /*!
    \brief apply impedance boundary condition (outlet pressures)

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,
    */
    void ImpedanceIntegration(
            DRT::ELEMENTS::Fluid3Boundary*    ele,
            ParameterList&                    params,
            DRT::Discretization&              discretization,
            vector<int>&                      lm,
            Epetra_SerialDenseVector&         elevec1);

    /*!
    \brief apply outflow boundary condition which is necessary for the conservative
    element formulation (since the convective term was partially integrated)

    \param params (in)        : ParameterList for communication between control routine
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : matrix to be filled by element.
    \param elevec1 (out)      : vector to be filled by element.

    */
    void ConservativeOutflowConsistency(
      DRT::ELEMENTS::Fluid3Boundary* ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseVector&  elevec1);

    /*!
    \brief compute potential Neumann inflow

    \param params (in)        : ParameterList for communication between control routine
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                                the controling method does not epxect the element to fill
                                this matrix.                          and elements
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,

    */
    void NeumannInflow(
      DRT::ELEMENTS::Fluid3Boundary*  ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseVector&  elevec1);

    enum ActionType
    {
      none,
      integrate_Shapefunction,
      areacalc,
      flowratecalc,
      calc_flowrate,
      flowratederiv,
      Outletimpedance,
      calc_node_normal,
      calc_node_curvature,
      calc_surface_tension,
      enforce_weak_dbc,
      conservative_outflow_bc,
      calc_Neumann_inflow
    };

  //! evaluate shape functions and derivatives at int. point
  double EvalShapeFuncAndIntFac(
      const DRT::UTILS::IntPointsAndWeights<bdrynsd_>&  intpoints,  ///< integration points
      const int                                         iquad,
      const LINALG::Matrix<nsd_,iel>&                   xyze,
      const std::vector<Epetra_SerialDenseVector>*      myknots,
      const Epetra_SerialDenseVector*                   weights,
      LINALG::Matrix<bdrynsd_,1>&                       xsi,
      LINALG::Matrix<iel,1>&                            funct,
      LINALG::Matrix<bdrynsd_,iel>&                     deriv,
      LINALG::Matrix<nsd_,1>*                           normal = NULL,
      double*                                           drs = NULL
  );

  //! evaluate Nurbs specific stuff
  bool GetKnotVectorAndWeightsForNurbs(
      DRT::ELEMENTS::Fluid3Boundary*              ele,
      DRT::Discretization&                        discretization,
      std::vector<Epetra_SerialDenseVector>&      mypknots,
      std::vector<Epetra_SerialDenseVector>&      myknots,
      Epetra_SerialDenseVector&                   weights,
      double&                                     normalfac);

  };  // end class fluid3boundaryImpl

}
}

#endif

#endif // CCADISCRET
#endif // D_FLUID3 or D_FLUID2

