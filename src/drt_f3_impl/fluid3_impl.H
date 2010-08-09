/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl.H

\brief Internal implementation of Fluid3 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_IMPL_H
#define FLUID3_IMPL_H

#include "fluid3_interface.H"

#include "../drt_f3/fluid3_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_fluid/fluid_rotsym_periodicbc.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
namespace ELEMENTS
{

  class Fluid3ImplParameter;

  /// Internal Fluid3 implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3 element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

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

    <h3>History</h3>

    The implementation here is based on the work by Christiane
    Foerster. The core element code has been linearized using MuPAD.
    The reaction terms have been removed from the lhs.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version with precompiled values.

    \author u.kue
    \date 07/07
  */


  template<DRT::Element::DiscretizationType distype>
  class Fluid3Impl: public Fluid3ImplInterface
  {
  public:

    friend class Fluid3ImplParameter;
    /// Constructor
    Fluid3Impl();

    //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
    static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of space dimensions
    static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

    static const int numdofpernode_ = nsd_+1;

    /*!
      \brief evaluate function for Fluid3 element
     */
    virtual int Evaluate(DRT::ELEMENTS::Fluid3*     ele,
                         ParameterList&             params,
                         DRT::Discretization&       discretization,
                         vector<int>&               lm,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra,
                         RefCountPtr<MAT::Material> mat);

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1
      \param eveln            (i) nodal velocities at n
      \param fsevelaf         (i) fine-scale nodal velocities at n+alpha_F/n+1
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1
      \param            (i) nodal accelerations at n+alpha_M
      \param escaaf           (i) nodal scalar at n+alpha_F/n+1
      \param escaam           (i) nodal scalar at n+alpha_M/n
      \param escadtam         (i) nodal scalar derivatives at n+alpha_M/n+1
      \param emhist           (i) time rhs for momentum equation
      \param edispnp          (i) nodal displacements (on moving mesh)
      \param egridv           (i) grid velocity (on moving mesh)
      \param estif            (o) element matrix to calculate
      \param emesh            (o) linearization wrt mesh motion
      \param eforce           (o) element rhs to calculate
      \param material         (i) fluid material
      \param time             (i) current simulation time
      \param timefac          (i) time discretization factor
      \param newton           (i) boolean flag for linearisation
      \param loma             (i) boolean flag for potential low-Mach-number solver
      \param conservative     (i) boolean flag for conservative form
      \param is_genalpha      (i) boolean flag for generalized-alpha time integration
      \param higher_order_ele (i) keep or drop second derivatives
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
      \param pspg             (i) boolean flag for stabilisation
      \param supg             (i) boolean flag for stabilisation
      \param vstab            (i) boolean flag for stabilisation
      \param cstab            (i) boolean flag for stabilisation
      \param cross            (i) boolean flag for stabilisation
      \param reynolds         (i) boolean flag for stabilisation
      \param whichtau         (i) selecting available tau definition
      \param turb_mod_action  (i) selecting turbulence model (none, Smagorisky,
                                  dynamic Smagorinsky, Smagorinsky with van Driest
                                  damping for channel flows)
      \param Cs               (i) Smagorinsky model parameter
      \param Cs_delta_sq      (i) Model parameter computed by dynamic Smagorinsky
                                  approach (Cs*h*h)
      \param l_tau            (i) viscous length scale, required for van driest
                                  damping function and defined on input
      */
    // 3D sysmat
#if 0
    void Sysmat(Fluid3*                                       ele,
                const LINALG::Matrix<nsd_,nen_>&               evelaf,
                const LINALG::Matrix<nsd_,nen_>&               eveln,
                const LINALG::Matrix<nsd_,nen_>&               fsevelaf,
                const LINALG::Matrix<nen_,1>&                  epreaf,
                const LINALG::Matrix<nsd_,nen_>&               eaccam,
                const LINALG::Matrix<nen_,1>&                  escaaf,
                const LINALG::Matrix<nen_,1>&                  escaam,
                const LINALG::Matrix<nen_,1>&                  escadtam,
                const LINALG::Matrix<nsd_,nen_>&               emhist,
                const LINALG::Matrix<nsd_,nen_>&               edispnp,
                const LINALG::Matrix<nsd_,nen_>&               egridv,
                LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    estif,
                LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    emesh,
                LINALG::Matrix<(nsd_+1)*nen_,    1>&           eforce,
                const double                                   thermpressaf,
                const double                                   thermpressam,
                const double                                   thermpressdtam,
                Teuchos::RCP<const MAT::Material>             material,
                const double                                  timefac,
                double&                                       Cs,
                double&                                       Cs_delta_sq,
                double&                                       l_tau);

    // 2D sysmat
    void Sysmat2D(Fluid3*                                  ele,
                const LINALG::Matrix<nsd_,nen_>&             evelaf,
                const LINALG::Matrix<nsd_,nen_>&             eveln,
                const LINALG::Matrix<nsd_,nen_>&             fsevelaf,
                const LINALG::Matrix<nen_,1>&             epreaf,
                const LINALG::Matrix<nsd_,nen_>&             ,
                const LINALG::Matrix<nen_,1>&             escaaf,
                const LINALG::Matrix<nen_,1>&             escaam,
                const LINALG::Matrix<nen_,1>&             escadtam,
                const LINALG::Matrix<nsd_,nen_>&             emhist,
                const LINALG::Matrix<nsd_,nen_>&             edispnp,
                const LINALG::Matrix<nsd_,nen_>&             egridv,
                LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&             estif,
                LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&             emesh,
                LINALG::Matrix<(nsd_+1)*nen_,    1>&             eforce,
                const double                                   thermpressaf,
                const double                                   thermpressam,
                const double                                   thermpressdtam,
                Teuchos::RCP<const MAT::Material>        material,
                const double                             timefac,
                double&                                  Cs);
#endif
    // 2D/3D mixed sysmat
    void Sysmat2D3D(Fluid3*                                    ele,
                const LINALG::Matrix<nsd_,nen_>&               evelaf,
                const LINALG::Matrix<nsd_,nen_>&               eveln,
                const LINALG::Matrix<nsd_,nen_>&               fsevelaf,
                const LINALG::Matrix<nen_,1>&                  epreaf,
                const LINALG::Matrix<nsd_,nen_>&               eaccam,
                const LINALG::Matrix<nen_,1>&                  escaaf,
                const LINALG::Matrix<nen_,1>&                  escaam,
                const LINALG::Matrix<nen_,1>&                  escadtam,
                const LINALG::Matrix<nsd_,nen_>&               emhist,
                const LINALG::Matrix<nsd_,nen_>&               edispnp,
                const LINALG::Matrix<nsd_,nen_>&               egridv,
                LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&   estif,
                LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&   emesh,
                LINALG::Matrix<(nsd_+1)*nen_,    1>&           eforce,
                const double                                   thermpressaf,
                const double                                   thermpressam,
                const double                                   thermpressdtam,
                Teuchos::RCP<const MAT::Material>              material,
                double&                                        Cs_delta_sq);

    // FD-check for debugging
    void FDcheck(
      Fluid3*                                               ele,
      const LINALG::Matrix<nsd_,nen_>&                      evelaf,
      const LINALG::Matrix<nsd_,nen_>&                      eveln,
      const LINALG::Matrix<nsd_,nen_>&                      fsevelaf,
      const LINALG::Matrix<nen_,1>&                         epreaf,
      const LINALG::Matrix<nsd_,nen_>&                      eaccam,
      const LINALG::Matrix<nen_,1>&                         escaaf,
      const LINALG::Matrix<nen_,1>&                         escaam,
      const LINALG::Matrix<nen_,1>&                         escadtam,
      const LINALG::Matrix<nsd_,nen_>&                      emhist,
      const LINALG::Matrix<nsd_,nen_>&                      edispnp,
      const LINALG::Matrix<nsd_,nen_>&                      egridv,
      const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    estif,
      const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    emesh,
      const LINALG::Matrix<(nsd_+1)*nen_,    1>&            eforce,
      const double                                          thermpressaf,
      const double                                          thermpressam,
      const double                                          thermpressdtam,
      const Teuchos::RCP<const MAT::Material>               material,
      const double                                          timefac,
      const double&                                         Cs,
      const double&                                         Cs_delta_sq,
      const double&                                         l_tau);

  private:

    /*!
     * \brief calculate body force from node condition
     */
    void BodyForce(Fluid3*      ele);

    //! number of components necessary to store second derivatives
    // 1 component  for nsd=1:  (N,xx)
    // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
    // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
    static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

    /*!
     * \brief evaluate shape functions and their derivatives at element center
     */
    void EvalShapeFuncAndDerivsAtEleCenter(const int  eleid);

    /*!
     * \brief evaluate shape functions and their derivatives at integration point
     */
    void EvalShapeFuncAndDerivsAtIntPoint(
        const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,
        const int                              iquad,
        const int                              eleid);

    /*!
       \brief set solution flags and parameters
     */
    //const double SetSolutionParameter(Fluid3*                    ele,
    //                          ParameterList&             params);

    /*!
     * \brief fill elment matrix and vectors with the global values
     */

    // TODO: think about location
    void ExtractValuesFromGlobalVector( const DRT::Discretization&          discretization, ///< discretization
                                        const vector<int>&                  lm,             ///<
                                        LINALG::Matrix<nsd_,nen_> *         matrixtofill,   ///< vector field
                                        LINALG::Matrix<nen_,1> *            vectortofill,   ///< scalar field
                                        const string                        state)          ///< state of the global vector
    {
      // get state of the global vector
      RefCountPtr<const Epetra_Vector> matrix_state = discretization.GetState(state);
      if(matrix_state == null)
        dserror("Cannot get state vector %s", state.c_str());

      // extract local values of the global vectors
      vector<double> mymatrix(lm.size());
      DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

      // rotate the vector field in the case of rotationally symmetric boundary conditions
      if(matrixtofill != NULL)
        rotsymmpbc_->RotateMyValuesIfNecessary(mymatrix);

      for (int inode=0; inode<nen_; ++inode)  // number of nodes
      {
        // fill a vector field via a pointer
        if (matrixtofill != NULL)
        {
          for(int idim=0; idim<nsd_; ++idim) // number of dimensions
            {
              (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
            }  // end for(idim)
        }
        // fill a scalar field via a pointer
        if (vectortofill != NULL)
          (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*numdofpernode_)];
      } // end for(inodes)

      return;
   }  // End FillElementVector


    /*!
       \brief get material parameters
     */
    void GetMaterialParams(Teuchos::RCP<const MAT::Material>  material,
                           const LINALG::Matrix<nsd_,nen_>&        evelaf,
                           const LINALG::Matrix<nen_,1>&           escaaf,
                           const LINALG::Matrix<nen_,1>&           escaam,
                           const double                            thermpressaf,
                           const double                            thermpressam,
                           const double                            thermpressdtam);

    /*!
       \brief get material parameters
     */
    void GetTurbulenceParams(Fluid3*                    ele,
                             ParameterList&             turbmodelparams,
                             //double&                    Cs,
                             double&                    Cs_delta_sq,
                             //double&                    l_tau,
                             int&                       nlayer);

    /*!
      \brief calculate (all-scale) subgrid viscosity
    */
    void CalcSubgrVisc(const LINALG::Matrix<nsd_,nen_>&     evelaf,
                       const double                        vol,
                       double&                             Cs,
                       double&                             Cs_delta_sq,
                       double&                             l_tau);

    /*!
      \brief calculate fine-scale subgrid viscosity
    */
    void CalcFineScaleSubgrVisc(const LINALG::Matrix<nsd_,nen_>&     evelaf,
                                const LINALG::Matrix<nsd_,nen_>&     fsevelaf,
                                const double                        vol,
                                double&                             Cs);

    /*!
      \brief calculate stabilization parameter

      Options for tau:

      o Franca-Barrenechea-Valentin-Wall:
           tauMp: Barrenechea, Valentin
           tauMu: Franca, Valentin
           tauC : Wall

      o Bazilevs: bazilevs
           tauM : Bazilevs, Whiting
           tauC : Bazilevs

      o Codina: codina
           tauMp: Barrenechea, Valentin
           tauMu: Franca, Valentin
           tauC : Codina
    */

    void CalcStabParameter(
        const double                            timefac,
        const double                            vol);


    void CalcDivEps(
        const LINALG::Matrix<nsd_,nen_>&        evelaf,
        const double &                          timefac);


    void GetResidualMomentumEq(
        const LINALG::Matrix<nsd_,nen_>&              eaccam,
        const double &                                timefac,
        double &                                      rhsfac);

    void UpdateSubscaleVelocity(
        Fluid3*       ele,
        double &      fac1,
        double &      fac2,
        double &      fac3,
        double &      facMtau,
        int &         iquad);


    void GetResidualContinuityEq(
        const LINALG::Matrix<nsd_,nen_>&          eveln,
        const LINALG::Matrix<nen_,1>&             escaaf,
        const LINALG::Matrix<nen_,1>&             escaam,
        const LINALG::Matrix<nen_,1>&             escadtam,
        const double & timefac);


    void LinGalMomResU(
        LINALG::Matrix<nsd_*nsd_,nen_> &        lin_resM_Du,
        const double &                          timefacfac);

    void LinGalMomResU_subscales(
        LINALG::Matrix<nen_*nsd_,nen_>          estif_p_v,
        LINALG::Matrix<nsd_*nsd_,nen_> &        lin_resM_Du,
        LINALG::Matrix<nsd_,1> &                resM_Du,
        const double &                          timefacfac,
        const double &                          facMtau);

    void InertiaAndConvectionGalPart(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
        LINALG::Matrix<nsd_,nen_> &             velforce,
        LINALG::Matrix<nsd_*nsd_,nen_> &        lin_resM_Du,
        LINALG::Matrix<nsd_,1> &                resM_Du,
        const double &                          rhsfac);


    void ViscousGalPart(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
        LINALG::Matrix<nsd_,nsd_> &             viscstress,
        const double &                          timefacfac,
        const double &                          rhsfac);


    void ContStab_and_ViscousTermRhs(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
        LINALG::Matrix<nsd_,nen_> &               velforce,
        LINALG::Matrix<nsd_,nsd_> &               viscstress,
        const double &                            timefac,
        const double &                            timefacfac,
        const double &                            tau_C,
        const double &                            rhsfac);


    void PressureGalPart(
        LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
        LINALG::Matrix<nsd_,nen_> &               velforce,
        const double &                            timefacfac,
        const double &                            rhsfac,
        const double &                            press);


    void ContinuityGalPart(
        LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,
        LINALG::Matrix<nen_,1> &                  preforce,
        const double &                            timefacfac,
        const double &                            rhsfac);


    void BodyForceRhsTerm(
        LINALG::Matrix<nsd_,nen_> &               velforce);


    void ConservativeFormulation(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
        LINALG::Matrix<nsd_,nen_> &               velforce,
        const double&                             timefacfac,
        const double &                            rhsfac);


    void LomaGalPart(
        LINALG::Matrix<nen_, nen_*nsd_> &       estif_q_u,
        LINALG::Matrix<nen_,1> &                preforce,
        const double &                          timefacfac,
        const double &                          rhsfac);


    void StabLinGalMomResU(
        LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
        const double&                             timefacfac);


    void PSPG(
        LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,
        LINALG::Matrix<nen_,nen_> &               ppmat,
        LINALG::Matrix<nen_,1> &                  preforce,
        LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
        const double &                            fac3,
        const double &                            timefacfac,
        const double &                            tau_Mp,
        const double &                            rhsfac);

    void SUPG(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
        LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
        LINALG::Matrix<nsd_,nen_> &               velforce,
        LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
        const double &                            fac3,
        const double &                            timefacfac,
        const double &                            timetauM);


    void ViscStab(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
        LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
        LINALG::Matrix<nsd_,nen_> &               velforce,
        LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
        const double &                            timefac,
        const double &                            tau_Mp,
        const double &                            vstabfac,
        const double &                            fac3);


    void CrossStressStab(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
        LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
        LINALG::Matrix<nsd_,nen_> &               velforce,
        LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
        const double &                            timefac,
        const double &                            timefacfac,
        const double &                            tau_Mp,
        const double &                            fac3);

    void ReynoldsStressStab(
        LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
        LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
        LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
        const double &                            timefacfac,
        const double &                            fac3);

    void FineScaleSubGridViscosityTerm(
        LINALG::Matrix<nsd_,nen_> &               velforce,
        const double &                            fssgviscfac);

    void LinMeshMotion_2D(
        LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
        const LINALG::Matrix<nsd_,nen_>&              evelaf,
        const double &                                press,
        const double &                                timefac,
        const double &                                timefacfac);

    void LinMeshMotion_3D(
        LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
        const LINALG::Matrix<nsd_,nen_>&              evelaf,
        const double &                                press,
        const double &                                timefac,
        const double &                                timefacfac);

    /*!
      \brief calculate rate of strain of (fine-scale) velocity

      \param evel       (i) nodal velocity values
      \param derxy      (i) shape function derivatives
      \param velderxy   (o) velocity derivatives

      \return computed rate of strain
     */
    double GetStrainRate(const LINALG::Matrix<nsd_,nen_>& evel,
                         const LINALG::Matrix<nsd_,nen_>& derxy,
                         LINALG::Matrix<nsd_,nsd_>        & velderxy)
    {
      double rateofstrain=0;

      // TODO: velderxy is calculated two times, probably ;)
      // get velocity derivatives at integration point
      //
      //              +-----  dN (x)
      //   dvel (x)    \        k
      //   -------- =   +     ------ * vel
      //      dx       /        dx        k
      //        j     +-----      j
      //              node k
      //
      // j : direction of derivative x/y/z
      //
      velderxy.MultiplyNT(evel,derxy);

      // compute (resolved) rate of strain
      //
      //          +-                                 -+ 1
      //          |          /   \           /   \    | -
      //          | 2 * eps | vel |   * eps | vel |   | 2
      //          |          \   / ij        \   / ij |
      //          +-                                 -+
      //
      LINALG::Matrix<nsd_,nsd_> two_epsilon;
      for(int rr=0;rr<nsd_;++rr)
      {
        for(int mm=0;mm<nsd_;++mm)
        {
          two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
        }
      }

      for(int rr=0;rr<nsd_;++rr)
      {
        for(int mm=0;mm<nsd_;++mm)
        {
          rateofstrain += two_epsilon(rr,mm)*two_epsilon(mm,rr);
        }
      }

      // sqrt(two_epsilon(rr,mm)*two_epsilon(mm,rr)/4.0*2.0)

      return(sqrt(rateofstrain/2.0));
    }


  private:
    //! DOF's per node
    //const int numdofpernode_;
    //! node coordinates
    LINALG::Matrix<nsd_,nen_> xyze_;
    //! dead load in nodes at n+alpha_F/n+1
    LINALG::Matrix<nsd_,nen_> edeadaf_;
    //! array for shape functions
    LINALG::Matrix<nen_,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,nen_> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<numderiv2_,nen_> deriv2_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xji_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<nsd_,nsd_> vderxy_;
    //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<nsd_,nsd_> fsvderxy_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,nen_> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<numderiv2_,nen_> derxy2_;
    //! bodyforce in gausspoint
    LINALG::Matrix<nsd_,1> bodyforce_;
    //! vector containing all values from previous timelevel n for momentum equation
    LINALG::Matrix<nsd_,1> histmom_;
    //! direction of flow (normed velocity vector)
    LINALG::Matrix<nsd_,1> velino_;
    //! velocity vector in gausspoint
    LINALG::Matrix<nsd_,1> velint_;
    //! fine-scale velocity vector in gausspoint
    LINALG::Matrix<nsd_,1> fsvelint_;
    //! subgrid-scale velocity vector in gausspoint
    LINALG::Matrix<nsd_,1> sgvelint_;
    //! density-weighted ale convective velocity in gausspoint, c=rho*u-rho*u_G
    LINALG::Matrix<nsd_,1> convvelint_;
    //! acceleration vector in gausspoint
    LINALG::Matrix<nsd_,1> accint_;
    //! pressure gradient in gausspoint
    LINALG::Matrix<nsd_,1> gradp_;
    //! the stabilisation parameters -> it is a (3,1) vector for 2D and 3D
    LINALG::Matrix<3,1> tau_;
    //! viscous term including 2nd derivatives
    //! (This array once had three dimensions, now the first two are combined to one.)
    LINALG::Matrix<nsd_*nsd_,nen_> viscs2_;
    //! linearisation of convection, convective part
    LINALG::Matrix<nen_,1> conv_c_;
    //! linearisation of subgrid-scale convection, convective part
    LINALG::Matrix<nen_,1> sgconv_c_;
    //! velocity divergence
    double vdiv_;
    //! total right hand side terms at int.-point for momentum equation
    LINALG::Matrix<nsd_,1> rhsmom_;
    //! (u_old*nabla)u_old
    LINALG::Matrix<nsd_,1> conv_old_;
    //! div epsilon(u_old)
    LINALG::Matrix<nsd_,1> visc_old_;
    //! old residual of momentum equation
    LINALG::Matrix<nsd_,1> momres_old_;
    //! old residual of continuity equation
    double conres_old_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::Matrix<numderiv2_,nsd_> xder2_;

    LINALG::Matrix<nsd_,nsd_> vderiv_;


    //! residual of momentum equation, velocity part which is not integrated by parts (precomputed)
    LINALG::Matrix<nsd_,1> resM_Du_;
    //! linearisation of residual of momentum equation wrt velocities (precomputed)
    LINALG::Matrix<nsd_*nsd_,nen_> lin_resM_Du_;

    LINALG::Matrix<nen_*nsd_,nen_*nsd_>     estif_;
    LINALG::Matrix<nen_*nsd_,nen_>          estif_p_v_;
    LINALG::Matrix<nen_, nen_*nsd_>         estif_q_u_;
    LINALG::Matrix<nen_,nen_>               ppmat;

    LINALG::Matrix<nen_,1>      preforce;
    LINALG::Matrix<nsd_,nen_>   velforce;


    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_,1> xsi_;

    //! flag for material evaluation at Gaussian integration points
    //bool mat_gp_;
    //! flag for stabilization parameter evaluation at Gaussian integration points
    //bool tau_gp_;
    //! Jacobian determinant
    double det_;
    //! integration factor
    double fac_;
    //! physical viscosity
    double visc_;
    //! (all-scale) subgrid viscosity
    double sgvisc_;
    //! effective viscosity = physical viscosity + (all-scale) subgrid viscosity
    double visceff_;
    //! fine-scale subgrid viscosity
    double fssgvisc_;

#if 0
    //! time-step length
    double dt_;
    //! factor for right-hand side due to time-integration scheme
    double omtheta_;
    //! generalised-alpha parameter (connecting velocity and acceleration)
    double gamma_;
    //! generalised-alpha parameter (velocity)
    double alphaF_;
    //! generalised-alpha parameter (acceleration)
    double alphaM_;
    //! generalised-alpha parameter, alphaF_*gamma_*dt_
    double afgdt_;
#endif

    //! LOMA-specific variables:
    //! right-hand-side term at int.-point for continuity equation
    double rhscon_;
    //! density at t_(n+alpha_F) or t_(n+1)
    double densaf_;
    //! density at t_(n+alpha_M)
    double densam_;
    //! density at t_(n)
    double densn_;
    //! delta density for Boussinesq Approximation
    double deltadens_;
    //! factor for scalar time derivative
    double scadtfac_;
    //! factor for convective scalar term at t_(n+alpha_F) or t_(n+1)
    double scaconvfacaf_;
    //! factor for convective scalar term at t_(n)
    double scaconvfacn_;
    //! addition to continuity equation due to thermodynamic pressure
    double thermpressadd_;
    //! velocity vector in gausspoint at t_(n)
    LINALG::Matrix<nsd_,1> velintn_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z at t_(n)
    LINALG::Matrix<nsd_,nsd_> vderxyn_;
    //! scalar gradient at t_(n+alpha_F) or t_(n+1)
    LINALG::Matrix<nsd_,1> grad_scaaf_;
    //! scalar gradient at t_(n)
    LINALG::Matrix<nsd_,1> grad_scan_;
    //! convective scalar term at t_(n+alpha_F) or t_(n+1)
    double conv_scaaf_;
    //! convective scalar term at t_(n)
    double conv_scan_;
    //! thermodynamic pressure at t_(n+alpha_F) or t_(n+1)
    //double thermpressaf_;
    //! thermodynamic pressure at t_(n+alpha_M) or t_(n)
    //double thermpressam_;
    //! time derivative of thermodynamic pressure at t_(n+alpha_M) or t_(n+1)
    //double thermpressdtam_;
    //! for the handling of rotationally symmetric periodic boundary conditions
    FLD::RotationallySymmetricPeriodicBC<distype>* rotsymmpbc_;

    // flags
    //! Flag to (de)activate higher order elements
    //! elements with only mixed second order derivatives are not counted as higher order elements
    //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
    bool is_higher_order_ele_;
#if 0
    //! actual time to evaluate the body BC
    double time_;
    //! flag to (de)activate generalized-alpha time-integration scheme
    bool is_genalpha_;
    //! flag to (de)activate conservative formulation
    bool is_conservative_;
    //! flag to (de)activate stationary formulation
    bool is_stationary_;
    //! flag to (de)activate Newton linearization
    bool is_newton_;
    //! Flag to (de)activate higher order elements
    //! elements with only mixed second order derivatives are not counted as higher order elements
    //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
    bool is_higher_order_ele_;
    //! Flag for physical type of the fluid flow (incompressible, loma, varying_density, Boussinesq)
    INPAR::FLUID::PhysicalType physicaltype_;

    //! Flag to (de)activate time-dependent subgrid stabilization
    INPAR::FLUID::SubscalesTD tds_;
    //! Flag to (de)activate time-dependent term in large-scale momentum equation
    INPAR::FLUID::Transient transient_;
    //! Flag to (de)activate PSPG stabilization
    INPAR::FLUID::PSPG pspg_;
    //! Flag to (de)activate SUPG stabilization
    INPAR::FLUID::SUPG supg_;
    //! Flag to (de)activate viscous term in residual-based stabilization
    INPAR::FLUID::VStab vstab_;
    //! Flag to (de)activate least-squares stabilization of continuity equation
    INPAR::FLUID::CStab cstab_;
    //! Flag to (de)activate cross-stress term -> residual-based VMM
    INPAR::FLUID::CrossStress cross_;
    //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
    INPAR::FLUID::ReynoldsStress reynolds_;

    //! Flag to define tau
    INPAR::FLUID::TauType whichtau_;
    //! flag to (de)activate fine-scale subgrid viscosity
    Fluid3::FineSubgridVisc fssgv_;
    //! flag to define turbulence model
    Fluid3::TurbModelAction turb_mod_action_;


    double timefac;

    double Cs;
    double l_tau;
#endif
/*
    bool boussinesq_;
    bool loma_;
    bool varyingdensity_;
*/

    DRT::ELEMENTS::Fluid3ImplParameter* f3Parameter_;
  };
}
}

#endif

#endif
#endif
