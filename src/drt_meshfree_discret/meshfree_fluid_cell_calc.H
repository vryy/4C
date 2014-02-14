/*!---------------------------------------------------------------------------

\file meshfree_fluid_cell_calc.H

\brief main file containing routines for calculation of meshfree fluid cell

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#ifndef MESHFREE_FLUID_CELL_CALC_H
#define MESHFREE_FLUID_CELL_CALC_H

#include "meshfree_fluid_cell_interface.H"
#include "drt_meshfree_cell_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // DRT::UTILS::DisTypeToDim
#include "../drt_lib/drt_element.H"                                   // Element::DiscretizationType

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace DRT
{
  class Discretization;

namespace MESHFREE
{
  class MeshfreeDiscretization;
  class Cell;
}

namespace ELEMENTS
{
  class MeshfreeFluid;
  class FluidEleParameter;
  class FluidEleParameterTimInt;

/*--------------------------------------------------------------------------*/
/*!
 * \brief Internal meshfree fluid cell implementation for
 *        backgroundmesh integration
 *
 * This internal class is based on Uli Kuettler's implementations of standard
 * scatra or fluid elements. It keeps all the working arrays needed to
 * calculate the transport element. Additionally the method Sysmat()
 * provides a clean and fast element implementation.
 *
 * <h3>Purpose</h3>
 *
 * The idea is to separate the element maintenance (class Transport)
 * from the mathematical contents (this class). Of course there are
 * different implementations of the Transport element, this is just one
 * such implementation.
 *
 * The Transport element will allocate exactly one object of this class
 * for all transport elements with the same number of nodes in the mesh.
 * This allows us to use exactly matching working arrays (and keep them
 * around.)
 *
 * The code is meant to be as clean as possible. This is the only way
 * to keep it fast. The number of working arrays has to be reduced to
 * a minimum so that the element fits into the cache. (There might be
 * room for improvements.)
 *
 * <h3>History</h3>
 *
 * The implementation here is the standard convection-diffusion element
 * capable of dealing with systems of transported scalars.
 *
 * Right now we do not read any stabilization parameters from the
 * input file but have a fixed version.
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2013
 */
/*--------------------------------------------------------------------------*/

template<DRT::Element::DiscretizationType distype>
class MeshfreeFluidCellCalc : public MeshfreeFluidCellInterface
{
  /*========================================================================*/
  //! @name friends
  /*========================================================================*/

  friend class FluidEleParameter;
  friend class FluidEleParameterTimInt;

public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  virtual ~MeshfreeFluidCellCalc() {}

  /*========================================================================*/
  //! @name Evaluate method
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Generic evaluation of the element. Called via base pointer.
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int Evaluate(
    DRT::ELEMENTS::MeshfreeFluid* cell,
    DRT::Discretization &         discretization,
    const std::vector<int> &      lm,
    Teuchos::ParameterList&       params,
    Teuchos::RCP<MAT::Material> & mat,
    Epetra_SerialDenseMatrix&     elemat1_epetra,
    Epetra_SerialDenseMatrix&     elemat2_epetra,
    Epetra_SerialDenseVector&     elevec1_epetra,
    Epetra_SerialDenseVector&     elevec2_epetra,
    Epetra_SerialDenseVector&     elevec3_epetra,
    bool                          offdiag = false
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate supporting methods of the element. Called via base
   * pointer.
   *
   * Interface function for supporting methods of the element
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int EvaluateService(
    DRT::ELEMENTS::MeshfreeFluid* ele,
    Teuchos::ParameterList&       params,
    Teuchos::RCP<MAT::Material>&  mat,
    DRT::Discretization&          discretization,
    std::vector<int>&             lm,
    Epetra_SerialDenseMatrix&     elemat1,
    Epetra_SerialDenseMatrix&     elemat2,
    Epetra_SerialDenseVector&     elevec1,
    Epetra_SerialDenseVector&     elevec2,
    Epetra_SerialDenseVector&     elevec3
    );


  /*========================================================================*/
  //! @name Query methods (inherited form base class)
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief general function to compute the error (analytical solution) for
   *        particular problem type at specific Gauss points
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int ComputeError(
    DRT::ELEMENTS::MeshfreeFluid* cell,
    Teuchos::ParameterList&       params,
    Teuchos::RCP<MAT::Material>&  mat,
    DRT::Discretization&          discretization,
    std::vector<int>&             lm,
    Epetra_SerialDenseVector&     elevec,
    const DRT::UTILS::GaussIntegration& intpoints
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief general function to compute the error (analytical solution) for
   *        particular problem type
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int ComputeError(
    DRT::ELEMENTS::MeshfreeFluid* cell,
    Teuchos::ParameterList&       params,
    Teuchos::RCP<MAT::Material>&  mat,
    DRT::Discretization&          discretization,
    std::vector<int>&             lm,
    Epetra_SerialDenseVector&     elevec
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief general function to compute the error (analytical solution) for
   *        particular problem type
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  int IntegrateShapeFunction(
    DRT::ELEMENTS::MeshfreeFluid* cell,
    DRT::Discretization&          discretization,
    const std::vector<int>&       lm,
    Epetra_SerialDenseVector&     elevec1_epetra
    );

//  /*------------------------------------------------------------------------*/
//  /*!
//   * \brief Finite difference check for debugging
//   *///                                                   (public) nis Jan13
//  /*------------------------------------------------------------------------*/
//  virtual void FDcheck(const LINALG::Matrix<nsd_,nen_>&                    evelaf,
//                       const LINALG::Matrix<nsd_,nen_>&                    eveln,
//                       const LINALG::Matrix<nsd_,nen_>&                    fsevelaf,
//                       const LINALG::Matrix<nen_,1>&                       epreaf,
//                       const LINALG::Matrix<nsd_,nen_>&                    eaccam,
//                       const LINALG::Matrix<nen_,1>&                       escaaf,
//                       const LINALG::Matrix<nen_,1>&                       escaam,
//                       const LINALG::Matrix<nen_,1>&                       escadtam,
//                       const LINALG::Matrix<nsd_,nen_>&                    emhist,
//                       const LINALG::Matrix<nsd_,nen_>&                    edispnp,
//                       const LINALG::Matrix<nsd_,nen_>&                    egridv,
//                       const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
//                       const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
//                       const LINALG::Matrix<(nsd_+1)*nen_,    1>&          eforce,
//                       const Teuchos::RCP<const MAT::Material>             material,
//                       const double                                        timefac,
//                       const double&                                       Cs,
//                       const double&                                       Cs_delta_sq,
//                       const double&                                       l_tau);

protected:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Protected constructor since we are a Singleton
   *///                                                (protected) nis Jan13
  /*------------------------------------------------------------------------*/
  MeshfreeFluidCellCalc();

  /*========================================================================*/
  //! @name Static member variables
  /*========================================================================*/

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of points
  static const int nek_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of Gauss points
  static const int ngp_ = DRT::MESHFREE::DisTypeToNumGaussPoints<distype>::ngp;

  //! number of degrees of freedomd per node (nsd_ for velocity, 1 for pressure)
  static const int numdofpernode_ = nsd_+1;

  /*========================================================================*/
  //! @name Auxiliary evaluate method
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief evaluate function for meshfree fluid cell
   *
   * Specific evaluate function without any knowledge about DRT objects. This
   * way the element evaluation is independent of the specific mesh storage.
   *
   *///                                                (protected) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int Evaluate(
    const DRT::ELEMENTS::MeshfreeFluid* cell,     //< (i) pointer to meshfree cell
    const Teuchos::ParameterList&       params,   //< (i) paramater list
    LINALG::SerialDenseMatrix &         elemat1,  //< (o) element matrix to calculate
    LINALG::SerialDenseMatrix &         elemat2,  //< (o) linearization wrt mesh motion
    LINALG::SerialDenseVector &         elevec1,  //< (o) element rhs to calculate
    const LINALG::SerialDenseMatrix &   evelaf,   //< (i) nodal velocities at n+alpha_F/n+1
    const LINALG::SerialDenseVector &   epreaf,   //< (i) nodal pressure at n+alpha_F/n+1
    const LINALG::SerialDenseVector &   eprenp,   //< (i) nodal pressure at n+1
    const LINALG::SerialDenseMatrix &   evelnp,   //< (i) nodal velocities at n+1 (np_genalpha)
    const LINALG::SerialDenseMatrix &   emhist,   //< (i) time rhs for momentum equation
    const LINALG::SerialDenseMatrix &   eaccam,   //< (i) nodal accelerations at n+alpha_M
    const LINALG::SerialDenseMatrix &   edispnp,  //< (i) nodal displacements (on moving mesh)
    const LINALG::SerialDenseMatrix &   egridv,   //< (i) grid velocity (on moving mesh)
    Teuchos::RCP<const MAT::Material>   mat,      //< (i) fluid material
    const bool                          isale,    //< (i) ALE flag
    const bool                          offdiag   //< (i) flag for off-diagonal matrix evaluation
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate element matrix and rhs
   *///                                                (protected) nis Jan13
  /*------------------------------------------------------------------------*/
  //! calculate element matrix and rhs
  virtual void Sysmat(
    const DRT::ELEMENTS::MeshfreeFluid* cell,           //< (i) pointer to meshfree cell
    const LINALG::SerialDenseMatrix&    evelaf,         //< (i) nodal velocities at n+alpha_F/n+1
    const LINALG::SerialDenseMatrix&    evelnp,         //< (i) nodal velocities at n+1 (np_genalpha)
    const LINALG::SerialDenseVector&    epreaf,         //< (i) nodal pressure at n+alpha_F/n+1
    const LINALG::SerialDenseVector&    eprenp,         //< (i) nodal pressure at n+1
    const LINALG::SerialDenseMatrix&    eaccam,         //< (i) nodal accelerations at n+alpha_M
    const LINALG::SerialDenseMatrix&    emhist,         //< (i) time rhs for momentum equation
    const LINALG::SerialDenseMatrix&    edispnp,        //< (i) nodal displacements (on moving mesh)
    const LINALG::SerialDenseMatrix&    egridv,         //< (i) grid velocity (on moving mesh)
    LINALG::SerialDenseMatrix&          estif,          //< (o) element matrix to calculate
    LINALG::SerialDenseMatrix&          emesh,          //< (o) linearization wrt mesh motion
    LINALG::SerialDenseVector&          eforce,         //< (o) element rhs to calculate
    Teuchos::RCP<const MAT::Material>   mat,            //< (i) fluid material
    const bool                          isale           //< (i) ALE flag
    );

  /*========================================================================*/
  //! @name Auxiliary method
  /*========================================================================*/

  //! calculate body force from nodal conditions
  void BodyForce(
    const DRT::ELEMENTS::MeshfreeFluid * cell,  //< pointer to meshfree cell
    const double *                       cgxyz  //< global position of current Gauss point
    );

  //! get material parameters
  void GetMaterialParams(
    Teuchos::RCP<const MAT::Material> mat,      //< reference pointer to material
    const LINALG::SerialDenseMatrix&  evelaf    //< velocity at time n+alpha_f / n+1
    );

  //! Provide linearization of Garlerkin momentum residual with respect to the velocities
  void LinGalMomResU(
    LINALG::SerialDenseMatrix & lin_resM_Du,    //< linearisation of the Garlerkin momentum residual
    const double                timefacfac      //< = timefac x fac
    );

  //! Compute element matrix and rhs entries: inertia, convective and
  //! reactive terms of the Galerkin part
  void InertiaConvectionReactionGalPart(
    LINALG::SerialDenseMatrix       & estif_u,        //< block (weighting function v x u)
    LINALG::SerialDenseMatrix       & velforce,       //< rhs forces velocity
    const LINALG::SerialDenseMatrix & lin_resM_Du,    //< linearisation of the Garlerkin momentum residual
    LINALG::Matrix<nsd_,1>          & resM_Du,        //< linearisation of the Garlerkin momentum residual
    const double                      rhsfac          //< right-hand-side factor
    );

  //! Compute element matrix entries: for the viscous terms of the Galerkin part
  void ViscousGalPart(
    LINALG::SerialDenseMatrix & estif_u,        //< block (weighting function v x u)
    LINALG::SerialDenseMatrix & velforce,       //< rhs forces velocity
    LINALG::Matrix<nsd_,nsd_> & viscstress,     //< viscous stresses
    const double                timefacfac,     //< = timefac x fac
    const double                rhsfac          //< right-hand-side factor
    );

  //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
  void PressureGalPart(
    LINALG::SerialDenseMatrix & estif_p_v,      //< block (weighting function v x pressure)
    LINALG::SerialDenseMatrix & velforce,       //< rhs forces velocity
    const double                timefacfac,     //< = timefac x fac
    const double                timefacfacpre,
    const double                rhsfac,         //< right-hand-side factor
    const double                press           //< pressure at integration point
    );

  //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
  void ContinuityGalPart(
    LINALG::SerialDenseMatrix & estif_q_u,      //< block (weighting function q x u)
    LINALG::SerialDenseVector & preforce,       //< rhs forces pressure
    const double                timefacfac,     //< = timefac x fac
    const double                timefacfacpre,
    const double                rhsfac          //< right-hand-side factor
    );

  //! Compute element matrix entries: pressure projection terms
  void PressureProjection(
    LINALG::SerialDenseMatrix & ppmat
    );

  //! Finalize pressure projection terms
  void PressureProjectionFinalize(
    LINALG::SerialDenseMatrix & estif_p_v,
    LINALG::SerialDenseMatrix & estif_q_u,
    LINALG::SerialDenseMatrix & ppmat,
    LINALG::SerialDenseMatrix & velforce,
    LINALG::SerialDenseVector & preforce,
    const LINALG::SerialDenseMatrix & evel,
    const LINALG::SerialDenseVector & epre
    );

  //! Compute element matrix entries: body force terms on rhs
  void BodyForceRhsTerm(
    LINALG::SerialDenseMatrix & velforce,       //< rhs forces velocity
    const double                rhsfac          //< right-hand-side factor for residuals
    );

  //! Compute element matrix entries: conservative formulation
  void ConservativeFormulation(
    LINALG::SerialDenseMatrix & estif_u,        //< block (weighting function v x u)
    LINALG::SerialDenseMatrix & velforce,       //< rhs forces velocity
    const double                timefacfac,     //< = timefac x fac
    const double                rhsfac          //< right-hand-side factor
    );

  //! fill cell matrix and vectors with the global values
  void ExtractValuesFromGlobalVector(
    const DRT::Discretization & discretization, //< discretization
    const std::vector<int>    & lm,             //<
    LINALG::SerialDenseMatrix * matrixtofill,   //< vector field
    LINALG::SerialDenseVector * vectortofill,   //< scalar field
    const std::string           state           //< state of the global vector
    );

  //! identify elements of inflow section
  void InflowElement(DRT::Element* ele);

  /*========================================================================*/
  //! @name Dofs and nodes
  /*========================================================================*/

  //! number of element nodes (not static for meshfree)
  int nen_;

  //! reference to meshfree discretisation
  DRT::MESHFREE::MeshfreeDiscretization* discret_;

  /*========================================================================*/
  //! @name other - te be kept clean!!!
  /*========================================================================*/

  //! element id
  int eid_;
  //! pointer to parameter list
  DRT::ELEMENTS::FluidEleParameter* fldpara_;
  //! pointer to parameter list for time integration
  DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
  //! identify elements of inflow section
  //! required for turbulence modeling
  bool is_inflow_ele_;

  //! node coordinates
  LINALG::SerialDenseMatrix nxyz_; // size: <nsd_,nen_>
  //! point coordinates
  LINALG::SerialDenseMatrix kxyz_; // size: <nsd_,nek_>
  //! Gauss point coordinates for boundary cell
  LINALG::SerialDenseMatrix gxyz_;
  //! Gauss point weights for boundary cell
  LINALG::SerialDenseVector gw_;
  //! array for solution basis functions: (nen_)
  Teuchos::RCP<LINALG::SerialDenseVector> sfunct_;
  //! global derivatives of solution basis functions w.r.t x,y,z: (nsd_ x nen_)
  Teuchos::RCP<LINALG::SerialDenseMatrix> sderiv_;
  //! array for weighting basis functions: (nen_)
  Teuchos::RCP<LINALG::SerialDenseVector> wfunct_;
  //! global derivatives of weighting basis functions w.r.t x,y,z: (nsd_ x nen_)
  Teuchos::RCP<LINALG::SerialDenseMatrix> wderiv_;
  //! global velocity derivatives in gausspoint w.r.t x,y,z
  LINALG::Matrix<nsd_,nsd_> vderxy_;
  //! bodyforce in gausspoint
  LINALG::Matrix<nsd_,1> bodyforce_;
  //! vector containing all values from previous timelevel n for momentum equation
  LINALG::Matrix<nsd_,1> histmom_;
  //! velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> velint_;
  //! grid velocity u_G at integration point
  LINALG::Matrix<nsd_,1> gridvelint_;
  //! ale convective velocity c=u-u_G at integration point
  LINALG::Matrix<nsd_,1> convvelint_;
  //! acceleration vector in gausspoint (for gen_alpha)
  LINALG::Matrix<nsd_,1> accint_;
  //! pressure gradient in gausspoint
  LINALG::Matrix<nsd_,1> gradp_;
  //! linearisation of convection, convective part
  LINALG::SerialDenseVector conv_c_; // size: <nen_,1>
  //! velocity divergence at t_(n+alpha_F) or t_(n+1)
  double vdiv_;
  //! total right hand side terms at int.-point for momentum equation
  LINALG::Matrix<nsd_,1> rhsmom_;
  //! (u_old*nabla)u_old
  LINALG::Matrix<nsd_,1> conv_old_;
  //! Jacobian determinant
  double det_;
  //! integration factor
  double fac_;
  //! physical viscosity
  double visc_;
  //! reaction coefficient
  double reacoeff_;

  /*========================================================================*/
  //! @name Polynomial pressure projection variables assuming constant mode a=1
  /*========================================================================*/

  //! "mass matrix" of projection modes - here element volume
  double D_;
  //! prolongator(?) of projection with solution basis functions
  LINALG::SerialDenseVector Eu_;
  //! integrated discrete velocity divergence with solution basis functions
  LINALG::SerialDenseMatrix Fu_;

  //! prolongator(?) of projection with weighting basis functions
  LINALG::SerialDenseVector Ev_;
  //! integrated discrete velocity divergence with weighting basis functions
  LINALG::SerialDenseMatrix Fv_;

  /*========================================================================*/
  //! @name LOMA-specific variables (for preparation)
  /*========================================================================*/

  //! density at t_(n+alpha_F) or t_(n+1)
  double densaf_;
  //! density at t_(n+alpha_M)
  double densam_;
  //! density at t_(n)
  double densn_;

};
} // end ELEMENTS
} // end DRT

#endif // end MESHFREE_FLUID_CELL_CALC_H
