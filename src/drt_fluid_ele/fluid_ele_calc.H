/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of fluid element

\level 1

\maintainer  Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_H
#define FLUID_ELE_CALC_H

#include "fluid_ele_interface.H"
#include "../drt_inpar/inpar_fluid.H"
#include "fluid_ele.H"

// forward declarations
namespace FLD
{
  template <DRT::Element::DiscretizationType distype, int numdofpernode,
      DRT::ELEMENTS::Fluid::EnrichmentType enrtype = DRT::ELEMENTS::Fluid::none>
  class RotationallySymmetricPeriodicBC;

  class TDSEleData;
}  // namespace FLD

namespace MAT
{
  class Material;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameter;
    class FluidEleParameterTimInt;

    /// Fluid element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the fluid element. Additionally, the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class Fluid) from the
      mathematical contents (this class). There are different
      implementations of the fluid element, this is just one such
      implementation.

      The fluid element will allocate exactly one object of this class for all
      fluid elements with the same number of nodes in the mesh. This
      allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>Usability</h3>

      The calculations are done by the Evaluate() method. There are two
      version. The virtual method that is inherited from FluidEleInterface
      (and called from Fluid) and the non-virtual one that does the actual
      work. The non-virtual Evaluate() method must be callable without an actual
      Fluid object.

      \author u.kue
      \date 07/07
    */

    template <DRT::Element::DiscretizationType distype,
        DRT::ELEMENTS::Fluid::EnrichmentType enrtype = DRT::ELEMENTS::Fluid::none>
    class FluidEleCalc : public FluidEleInterface
    {
     public:
      /// dtor
      virtual ~FluidEleCalc() {}

      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const int nen_ = DRT::ELEMENTS::MultipleNumNode<enrtype>::multipleNode *
                              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      static const int numdofpernode_ = nsd_ + 1;

      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1);


      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, const DRT::UTILS::GaussIntegration& intpoints);


      virtual int IntegrateShapeFunctionXFEM(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells)
      {
        dserror("Implemented in derived xfem class!");
        return 1;
      };


      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */
      virtual int EvaluateService(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      /*! \brief Calculate a integrated divergence operator in vector form
       *
       *   The vector valued operator \f$B\f$ is constructed such that
       *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
       *
       *   \author mayr.mt
       *   \date   04/2012
       */
      virtual int CalcDivOp(DRT::ELEMENTS::Fluid* ele,  //< current fluid element
          DRT::Discretization& discretization,          //< fluid discretization
          std::vector<int>& lm,                         //< some DOF management
          Epetra_SerialDenseVector& elevec1             //< reference to element vector to be filled
      );

      /*! \brief Calculate element mass matrix
       *
       *  \author mayr.mt \date 05/2014
       */
      virtual int CalcMassMatrix(DRT::ELEMENTS::Fluid* ele,
          //    Teuchos::ParameterList&              params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra
          //    const DRT::UTILS::GaussIntegration & intpoints
      );

      /*! \brief Interpolate velocity gradient and pressure to given point
       *
       *  \author rauch \date 05/2014
       */
      int InterpolateVelocityGradientAndPressure(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra);

      /*! \brief Interpolate velocity to given point
       *
       *  \author rauch \date 05/2014
       */
      int InterpolateVelocityToNode(Teuchos::ParameterList& params, DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra);

      /*! \brief Interpolate velocity to given point
       *
       *  \author rauch \date 07/2015
       */
      int CorrectImmersedBoundVelocities(Teuchos::ParameterList& params, DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra);

      /*---------------------------------------------------------------------*
       | Action type: interpolate_velocity_to_given_point                    |
       | calculate velocity at given point                       ghamm 12/15 |
       *---------------------------------------------------------------------*/
      int InterpolateVelocityToPoint(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2);

      /*! \brief Interpolate pressure to given point
       *
       *  \author ghamm \date 06/2015
       */
      int InterpolatePressureToPoint(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra);

      /*! \brief Reset debug output of immersed element
       *
       *  \author rauch \date 05/2014
       */
      int ResetImmersedEle(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params);

      /*! \brief Calculate coordinates and velocities and element center
       *
       *  \author bk \date 01/2015
       */
      virtual int CalcVelGradientEleCenter(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2);


      /*! \brief Calculate properties for adaptive time step based on CFL number
       *
       *  \author bk \date 08/2014
       */
      virtual int CalcTimeStep(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1);

      /*! \brief Calculate channel statistics
       *
       *  \author bk \date 05/2014
       */
      virtual int CalcChannelStatistics(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<MAT::Material>& mat);

      /*! \brief Calculate mass flow for periodic hill
       *
       *  \author bk \date 12/2014
       */
      virtual int CalcMassFlowPeriodicHill(DRT::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Teuchos::RCP<MAT::Material>& mat);

      /*! \brief Project velocity gradient to nodal level
       *
       *   \author ghamm
       *   \date   06/2014
       */
      virtual int VelGradientProjection(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2);


      /*! \brief Project pressure gradient to nodal level
       *
       *   \author mwinter
       *   \date   09/2015
       */
      virtual int PresGradientProjection(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2);

      /*! \brief Calculate a divergence of velocity at the element center
       *
       *   \author ehrl
       *   \date   12/2012
       */
      virtual int ComputeDivU(DRT::ELEMENTS::Fluid* ele,  //< current fluid element
          DRT::Discretization& discretization,            //< fluid discretization
          std::vector<int>& lm,                           //< location vector for DOF management
          Epetra_SerialDenseVector& elevec1  //< reference to element vector to be filled
      );

      /// Evaluate element ERROR
      /*!
          general function to compute the error (analytical solution) for particular problem type
       */
      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec);

      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          const DRT::UTILS::GaussIntegration& intpoints2);

      /*!
       \brief Evaluates the analytic solution in the given point
       */
      static void EvaluateAnalyticSolutionPoint(const LINALG::Matrix<nsd_, 1>& xyzint,
          const double t, const INPAR::FLUID::CalcError calcerr, const int calcerrfunctno,
          const Teuchos::RCP<MAT::Material>& mat, LINALG::Matrix<nsd_, 1>& u, double& p,
          LINALG::Matrix<nsd_, nsd_>& dervel, bool isFullImplPressure = false, double deltat = 0.0);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false);

      /// Evaluate the element at specified gauss points
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints, bool offdiag = false);

      virtual int ComputeErrorInterface(DRT::ELEMENTS::Fluid* ele,   ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<MAT::Material>& mat,                          ///< material
          Epetra_SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,                               ///< boundary integration points
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                ///< parameter list
      )
      {
        dserror("Implemented in derived xfem class!");
        return 1;
      };

      /// Evaluate the XFEM cut element
      virtual int EvaluateXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells, bool offdiag = false)
      {
        dserror("Implemented in derived xfem class!");
        return 1;
      }

      /*!
        \brief calculate dissipation of various terms (evaluation of turbulence models)
      */
      virtual int CalcDissipation(Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      /*!
        \brief finite difference check for debugging
      */
      virtual void FDcheck(const LINALG::Matrix<nsd_, nen_>& evelaf,
          const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nsd_, nen_>& fsevelaf,
          const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nsd_, nen_>& eaccam,
          const LINALG::Matrix<nen_, 1>& escaaf, const LINALG::Matrix<nen_, 1>& escaam,
          const LINALG::Matrix<nen_, 1>& escadtam, const LINALG::Matrix<nsd_, nen_>& emhist,
          const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& egridv,
          const LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          const LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          const LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          const Teuchos::RCP<const MAT::Material> material, const double timefac, const double& Cs,
          const double& Cs_delta_sq, const double& l_tau);


      virtual void ElementXfemInterfaceHybridLM(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,  ///< element gauss points
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,                         ///< side coupling matrices
          Teuchos::ParameterList& params,            ///< parameter list
          Teuchos::RCP<MAT::Material>& mat,          ///< material
          Epetra_SerialDenseMatrix& elemat1_epetra,  ///< local system matrix of intersected element
          Epetra_SerialDenseVector&
              elevec1_epetra,               ///< local element vector of intersected element
          Epetra_SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const GEO::CUT::plain_volumecell_set& vcSet  ///< set of plain volume cells
      )
      {
        dserror("Implemented in derived xfem class!");
        return;
      }


      virtual void ElementXfemInterfaceNIT(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,               ///< parameter list
          Teuchos::RCP<MAT::Material>& mat_master,      ///< material for the coupled side
          Teuchos::RCP<MAT::Material>& mat_slave,       ///< material for the coupled side
          Epetra_SerialDenseMatrix& elemat1_epetra,     ///< element matrix
          Epetra_SerialDenseVector& elevec1_epetra,     ///< element vector
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,                ///< side coupling matrices
          Epetra_SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
      )
      {
        dserror("Implemented in derived xfem class!");
        return;
      }

      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra,
          const DRT::UTILS::GaussIntegration& intpoints)
      {
        dserror("Implemented in derived xfem class!");
        return;
      }

      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra)
      {
        dserror("Implemented in derived xfem class!");
        return;
      }

      /// calculate body force from nodal conditions. Static function interface to allow
      /// for its use in FluidEleCalcHDG.
      static void BodyForce(DRT::ELEMENTS::Fluid* ele,    //< pointer to element
          const double time,                              //< current time
          const INPAR::FLUID::PhysicalType physicaltype,  //< physical type
          LINALG::Matrix<nsd_, nen_>& ebofoaf,            //< body force at nodes
          LINALG::Matrix<nsd_, nen_>&
              eprescpgaf,  //< prescribed pressure gradient (required for turbulent channel flow!)
          LINALG::Matrix<nen_, 1>& escabofoaf  //< scatra body force at nodes
      );

      /// calculate correction term at nodes
      static void CorrectionTerm(DRT::ELEMENTS::Fluid* ele,  //< pointer to element
          LINALG::Matrix<1, nen_>& ecorrectionterm           //<correction term at nodes
      );


     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalc();

      /*!
        \brief evaluate function for fluid element

        Specific evaluate function without any knowledge about DRT objects. This
        way the element evaluation is independent of the specific mesh storage.
       */
      virtual int Evaluate(Teuchos::ParameterList& params,
          const LINALG::Matrix<nsd_, nen_>& ebofoaf, const LINALG::Matrix<nsd_, nen_>& eprescpgaf,
          const LINALG::Matrix<nsd_, nen_>& ebofon, const LINALG::Matrix<nsd_, nen_>& eprescpgn,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& elemat1,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& elemat2,
          LINALG::Matrix<(nsd_ + 1) * nen_, 1>& elevec1, const LINALG::Matrix<nsd_, nen_>& evelaf,
          const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nsd_, nen_>& evelam,
          const LINALG::Matrix<nen_, 1>& epream, const LINALG::Matrix<nen_, 1>& eprenp,
          const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nen_, 1>& escaaf,
          const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nsd_, nen_>& eaccam,
          const LINALG::Matrix<nen_, 1>& escadtam, const LINALG::Matrix<nsd_, nen_>& eveldtam,
          const LINALG::Matrix<nen_, 1>& epredtam, const LINALG::Matrix<nen_, 1>& escabofoaf,
          const LINALG::Matrix<nen_, 1>& escabofon, const LINALG::Matrix<nsd_, nen_>& eveln,
          const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& escaam,
          const LINALG::Matrix<nsd_, nen_>& edispnp, const LINALG::Matrix<nsd_, nen_>& egridv,
          const LINALG::Matrix<nsd_, nen_>& egridvn, const LINALG::Matrix<nsd_, nen_>& fsevelaf,
          const LINALG::Matrix<nen_, 1>& fsescaaf, const LINALG::Matrix<nsd_, nen_>& evel_hat,
          const LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const LINALG::Matrix<nen_, 1>& eporo, const LINALG::Matrix<nsd_, 2 * nen_>& egradphi,
          const LINALG::Matrix<nen_, 2 * 1>& ecurvature, Teuchos::RCP<MAT::Material> mat,
          bool isale, bool isowned, double CsDeltaSq, double CiDeltaSq, double* saccn,
          double* sveln, double* svelnp, const DRT::UTILS::GaussIntegration& intpoints,
          bool offdiag);

      /*!
        \brief calculate element matrix and rhs

        \param ebofoaf          (i) body force at n+alpha_F/n+1
        \param eprescpgaf       (i) prescribed pressure gradient at n+alpha_F/n+1 (required for
        turbulent channel flow) \param ebofon           (i) body force at n \param eprescpgaf (i)
        prescribed pressure gradient at n (required for turbulent channel flow) \param evelaf (i)
        nodal velocities at n+alpha_F/n+1 \param evelam           (i) nodal velocities at
        n+alpha_M/n \param eveln            (i) nodal velocities at n \param evelnp           (i)
        nodal velocities at n+1 (np_genalpha) \param fsevelaf         (i) fine-scale nodal
        velocities at n+alpha_F/n+1 \param fsescaaf         (i) fine-scale nodal scalar at
        n+alpha_F/n+1 \param epreaf           (i) nodal pressure at n+alpha_F/n+1 \param epream (i)
        nodal pressure at n+alpha_M/n \param eprenp           (i) nodal pressure at n+1 \param
        eaccam           (i) nodal accelerations at n+alpha_M \param escaaf           (i) nodal
        scalar at n+alpha_F/n+1 \param escaam           (i) nodal scalar at n+alpha_M/n \param
        escadtam         (i) nodal scalar derivatives at n+alpha_M/n+1 \param eveldtam         (i)
        nodal velocity derivatives at n+alpha_M/n+1 \param epredtam         (i) nodal pressure
        derivatives at n+alpha_M/n+1 \param emhist           (i) time rhs for momentum equation
        \param edispnp          (i) nodal displacements (on moving mesh)
        \param egridv           (i) grid velocity (on moving mesh)
        \param estif            (o) element matrix to calculate
        \param emesh            (o) linearization wrt mesh motion
        \param eforce           (o) element rhs to calculate
        \param egradphi         (i) gradient of nodal scalar at nodes
        \param ecurvature       (i) curvature of scalar at nodes
        \param thermpressaf     (i) thermodynamic pressure at n+alpha_F/n+1
        \param thermpressam     (i) thermodynamic pressure at n+alpha_M/n
        \param thermpressdtaf   (i) thermodynamic pressure derivative at n+alpha_F/n+1
        \param thermpressdtam   (i) thermodynamic pressure derivative at n+alpha_M/n+1
        \param material         (i) fluid material
        \param Cs_delta_sq      (i) parameter for dynamic Smagorinsky model (Cs*h*h)
        \param isale            (i) ALE flag
        \param intpoints        (i) Gaussian integration points

        */
      virtual void Sysmat(const LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const LINALG::Matrix<nsd_, nen_>& eprescpgaf, const LINALG::Matrix<nsd_, nen_>& ebofon,
          const LINALG::Matrix<nsd_, nen_>& eprescpgn, const LINALG::Matrix<nsd_, nen_>& evelaf,
          const LINALG::Matrix<nsd_, nen_>& evelam, const LINALG::Matrix<nsd_, nen_>& eveln,
          const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& fsevelaf,
          const LINALG::Matrix<nen_, 1>& fsescaaf, const LINALG::Matrix<nsd_, nen_>& evel_hat,
          const LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& epream,
          const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& eprenp,
          const LINALG::Matrix<nsd_, nen_>& eaccam, const LINALG::Matrix<nen_, 1>& escaaf,
          const LINALG::Matrix<nen_, 1>& escaam, const LINALG::Matrix<nen_, 1>& escadtam,
          const LINALG::Matrix<nsd_, nen_>& eveldtam, const LINALG::Matrix<nen_, 1>& epredtam,
          const LINALG::Matrix<nen_, 1>& escabofoaf, const LINALG::Matrix<nen_, 1>& escabofon,
          const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nsd_, nen_>& edispnp,
          const LINALG::Matrix<nsd_, nen_>& egridv,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const LINALG::Matrix<nen_, 1>& eporo,
          const LINALG::Matrix<nsd_, 2 * nen_>& egradphi,
          const LINALG::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
          double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
          const DRT::UTILS::GaussIntegration& intpoints);


      virtual void SysmatOSTNew(const LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const LINALG::Matrix<nsd_, nen_>& eprescpgaf, const LINALG::Matrix<nsd_, nen_>& ebofon,
          const LINALG::Matrix<nsd_, nen_>& eprescpgn, const LINALG::Matrix<nsd_, nen_>& evelaf,
          const LINALG::Matrix<nsd_, nen_>& evelam, const LINALG::Matrix<nsd_, nen_>& eveln,
          const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& fsevelaf,
          const LINALG::Matrix<nen_, 1>& fsescaaf, const LINALG::Matrix<nsd_, nen_>& evel_hat,
          const LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& epream,
          const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& eprenp,
          const LINALG::Matrix<nsd_, nen_>& eaccam, const LINALG::Matrix<nen_, 1>& escaaf,
          const LINALG::Matrix<nen_, 1>& escaam, const LINALG::Matrix<nen_, 1>& escadtam,
          const LINALG::Matrix<nen_, 1>& escabofoaf, const LINALG::Matrix<nen_, 1>& escabofon,
          const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nsd_, nen_>& edispnp,
          const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const LINALG::Matrix<nen_, 1>& eporo,
          const LINALG::Matrix<nsd_, 2 * nen_>& egradphi,
          const LINALG::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
          double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
          const DRT::UTILS::GaussIntegration& intpoints);



      //! number of components necessary to store second derivatives
      /*!
       1 component  for nsd=1:  (N,xx)

       3 components for nsd=2:  (N,xx ; N,yy ; N,xy)

       6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
      */
      static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

      //! calculate body force from nodal conditions
      void BodyForce(DRT::ELEMENTS::Fluid* ele,  //< pointer to element
          LINALG::Matrix<nsd_, nen_>& ebofoaf,   //< body force at nodes
          LINALG::Matrix<nsd_, nen_>&
              eprescpgaf,  //< prescribed pressure gradient (required for turbulent channel flow!)
          LINALG::Matrix<nen_, 1>& escabofoaf  //< scatra body force at nodes
      );

      //! calculate body force contribution from surface tension
      void AddSurfaceTensionForce(
          const LINALG::Matrix<nen_, 1>& escaaf,           ///< scalar at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& escaam,           ///< scalar at time n+alpha_m / n
          const LINALG::Matrix<nsd_, 2 * nen_>& egradphi,  //<gradient of scalar function at nodes
          const LINALG::Matrix<nen_, 2 * 1>& ecurvature    //<curvature of scalar function at nodes
      );


      //! evaluate shape functions and their derivatives at element center
      virtual void EvalShapeFuncAndDerivsAtEleCenter();

      //! brief evaluate shape functions and their derivatives at integration point
      virtual void EvalShapeFuncAndDerivsAtIntPoint(
          const double* gpcoords,  ///< actual integration point (coords)
          double gpweight          ///< actual integration point (weight)
      );

      //! get ALE grid displacements and grid velocity for element
      void GetGridDispVelALE(DRT::Discretization& discretization, const std::vector<int>& lm,
          LINALG::Matrix<nsd_, nen_>& edispnp, LINALG::Matrix<nsd_, nen_>& egridv);

      //! get ALE grid displacements and grid velocity for element
      void GetGridDispVelALEOSTNew(DRT::Discretization& discretization, const std::vector<int>& lm,
          LINALG::Matrix<nsd_, nen_>& edispnp, LINALG::Matrix<nsd_, nen_>& egridvnp,
          LINALG::Matrix<nsd_, nen_>& egridvn);

      //! get ALE grid displacements for element
      virtual void GetGridDispALE(DRT::Discretization& discretization, const std::vector<int>& lm,
          LINALG::Matrix<nsd_, nen_>& edispnp);

      //! set the (relative) convective velocity at integration point for various physical types
      void SetConvectiveVelint(const bool isale);

      //! set the (relative) convective velocity at integration point for various physical types
      void SetConvectiveVelintN(const bool isale);

      //! set element advective field for Oseen problems
      void SetAdvectiveVelOseen(DRT::ELEMENTS::Fluid* ele);

      //! get material parameters
      void GetMaterialParams(
          Teuchos::RCP<const MAT::Material> material,  ///< reference pointer to material
          const LINALG::Matrix<nsd_, nen_>& evelaf,    ///< velocity at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& epreaf,       ///< pressure at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& epream,       ///< pressure at time n+alpha_m / n
          const LINALG::Matrix<nen_, 1>& escaaf,       ///< scalar at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& escaam,       ///< scalar at time n+alpha_m / n
          const LINALG::Matrix<nen_, 1>&
              escabofoaf,               ///< body force for scalar transport at time n+alpha_f / n+1
          const double thermpressaf,    ///< thermodynamic pressure at time n+alpha_f / n+1
          const double thermpressam,    ///< thermodynamic pressure at time n+alpha_m / n
          const double thermpressdtaf,  ///< time derivative of thermodynamic pressure at time
                                        ///< n+alpha_f / n+1
          const double thermpressdtam,  ///< time derivative of thermodynamic pressure at time
                                        ///< n+alpha_m / n+1
          const double vol              ///< element volume
      );

      //! get material parameters
      void GetMaterialParams(Teuchos::RCP<const MAT::Material> material,
          const LINALG::Matrix<nsd_, nen_>& evelaf, const LINALG::Matrix<nen_, 1>& epreaf,
          const LINALG::Matrix<nen_, 1>& epream, const LINALG::Matrix<nen_, 1>& escaaf,
          const LINALG::Matrix<nen_, 1>& escaam, const LINALG::Matrix<nen_, 1>& escabofoaf,
          const double thermpressaf, const double thermpressam, const double thermpressdtaf,
          const double thermpressdtam, const double vol, double& densam, double& densaf,
          double& densan, double& visc, double& viscn, double& gamma);

      //! return constant mk for stabilization parameters
      virtual double GetMK();

      //! calculate stabilization parameter
      void CalcStabParameter(const double vol);  ///< volume

      //! calculate characteristic element length
      void CalcCharEleLength(const double vol,  ///< volume
          const double vel_norm,                ///< norm of velocity vector
          double& h_u,                          ///< length for tau_Mu
          double& h_p);                         ///< length for tau_Mp/tau_C


      //! calculate div(epsilon(u))
      void CalcDivEps(
          const LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const LINALG::Matrix<nsd_, nen_>& eveln);  ///< velocity at time n+alpha_f / n+1

      //! compute residual of momentum equation and subgrid-scale velocity
      virtual void ComputeSubgridScaleVelocity(
          const LINALG::Matrix<nsd_, nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                              ///< factor for old s.-s. velocities
          double& fac2,                              ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
      );


      //! compute residual of momentum equation and subgrid-scale velocity
      void ComputeSubgridScaleVelocityOSTNew(
          const LINALG::Matrix<nsd_, nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                              ///< factor for old s.-s. velocities
          double& fac2,                              ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      virtual void LinGalMomResU(
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac  ///< = timefac x fac
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      void LinGalMomResUOSTNew(
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac  ///< = timefac x fac
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities in the
      //! case if subscales
      void LinGalMomResU_subscales(LINALG::Matrix<nen_ * nsd_, nen_>&
                                       estif_p_v,  ///< block (weighting function v x pressure)
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,                   ///< linearisation of the Garlerkin momentum residual
          LINALG::Matrix<nsd_, 1>& resM_Du,  ///< residual of the fluid momentum equation
          const double& timefacfac,          ///< (time factor) x (integration factor)
          const double& facMtau              ///< facMtau = modified tau_M (see code)
      );

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      void InertiaConvectionReactionGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,                   ///< linearisation of the Garlerkin momentum residual
          LINALG::Matrix<nsd_, 1>& resM_Du,  ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac,              ///< right-hand-side factor
          const double& rhsfacn              ///< right-hand-side factor time step n
      );

      //! Compute element matrix entries: for the viscous terms of the Galerkin part
      void ViscousGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_, nsd_>& viscstress,             ///< viscous stresses
          const double& timefacfac,                           ///< = timefac x fac
          const double& rhsfac,                               ///< right-hand-side factor
          const double& rhsfacn  ///< right-hand-side factor time step n
      );

      //! Compute element matrix entries: div-grad stabilization and the rhs of the viscous term
      void ContStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefac,                              ///< time factor
          const double& timefacfac,                           ///< = timefac x fac
          const double& timefacfacpre,                        ///< = timefac x fac
          const double& rhsfac,                               ///< right-hand-side factor
          const double& rhsfacn  ///< right-hand-side factor time step n
      );

      //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
      void PressureGalPart(LINALG::Matrix<nen_ * nsd_, nen_>&
                               estif_p_v,        ///< block (weighting function v x pressure)
          LINALG::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,              ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,   ///< right-hand-side factor
          const double& rhsfacn,  ///< right-hand-side factor time step n
          const double& press,    ///< pressure at integration point
          const double& pressn    ///< pressure at integration point
      );

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      void ContinuityGalPart(
          LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                      ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor
          const double& rhsfacn  ///< right-hand-side factor at time step n
      );

      //! Compute element matrix entries: pressure projection terms
      void PressureProjection(LINALG::Matrix<nen_, nen_>& ppmat);

      //! Finalize pressure projection terms
      void PressureProjectionFinalize(LINALG::Matrix<nen_, nen_>& ppmat,
          LINALG::Matrix<nen_, 1>& preforce, const LINALG::Matrix<nen_, 1>& eprenp);

      //! Compute element matrix entries: body force terms on rhs
      void BodyForceRhsTerm(LINALG::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const double rhsfacn   //= 0.0 ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: conservative formulation
      virtual void ConservativeFormulation(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefacfac,                           ///< = timefac x fac
          const double& rhsfac                                ///< right-hand-side factor
      );

      //! Provide linearization of stabilization residual with respect to the velocities
      void StabLinGalMomResU(LINALG::Matrix<nsd_ * nsd_, nen_>&
                                 lin_resM_Du,  ///< linearisation of the stabilization residual
          const double& timefacfac             ///< = timefac x fac
      );

      //! Compute element matrix entries: PSPG
      void PSPG(LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, nen_>& ppmat,                   ///< block (weighting function q x p)
          LINALG::Matrix<nen_, 1>& preforce,                   ///< rhs forces pressure
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const int iquad        ///< index of current integration point
      );

      //! Compute element matrix entries: PSPG
      void PSPGOSTNew(
          LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, nen_>& ppmat,             ///< block (weighting function q x p)
          LINALG::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const int iquad        ///< index of current integration point
      );

      //! Compute element matrix entries: SUPG
      void SUPG(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nen_, 1>& preforce,                  ///< rhs forces pressure
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: SUPG
      void SUPGOSTNew(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nen_, 1>& preforce,                  ///< rhs forces pressure
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: reactive stabilization
      void ReacStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,              ///< linearisation of the stabilization residual
          const double& timefacfac,     ///< = timefac x fac
          const double& timefacfacpre,  ///< = timefacpre x fac
          const double& rhsfac,         ///< right-hand-side factor for residuals
          const double& fac3            ///< factor for residual in current subgrid velocities
      );

      //! Compute element matrix entries: viscous stabilization
      void ViscStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,              ///< linearisation of the stabilization residual
          const double& timefacfac,     ///< = timefac x fac
          const double& timefacfacpre,  ///< = timefac x fac
          const double& rhsfac,         ///< right-hand-side factor for residuals
          const double& fac3            ///< factor for residual in current subgrid velocities
      );

      //! Compute element matrix entries: convective divergence stabilization for XFEM
      void ConvDivStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefacfac,                           ///< = timefac x fac
          const double& rhsfac  ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: cross stress stabilization
      void CrossStressStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const double& fac3     ///< factor for residual in current subgrid velocities
      );

      //! Compute element matrix entries: Reynolds stress stabilization
      void ReynoldsStressStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& timefacfac,  ///< timefac x fac
          const double& timefacfacpre,
          const double& fac3  ///< factor for residual in current subgrid velocities
      );

      //! turbulence related methods
      //! definition in fluid_impl_turbulence_service.cpp

      //! get parameters for multifractal subgrid scales
      void PrepareMultifractalSubgrScales(
          LINALG::Matrix<nsd_, 1>& B_mfs,  ///< coefficient multifractal subgrid scales velocity
          double& D_mfs,  ///< coefficient multifractal subgrid scales scalar (loma only)
          const LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const LINALG::Matrix<nsd_, nen_>&
              fsevelaf,     ///< fine scale velocity at time n+alpha_f / n+1
          const double vol  ///< volume
      );

      //! get turbulence parameter
      void GetTurbulenceParams(
          Teuchos::ParameterList& turbmodelparams,  ///< pointer general turbulence parameter list
          double& Cs_delta_sq,                      ///< parameter CS in dynamic Smagorinsky
          double& Ci_delta_sq,  ///< parameter CI in dynamic Smagorinsky for loma
          int& nlayer,  ///< number of layers for computation of parameter CS in dynamic Smagorinsky
          double CsDeltaSq,  ///< parameter CS in dynamic Smagorinsky computed in DynSmagFilter()
          double CiDeltaSq   ///< parameter CI in dynamic Smagorinsky for loma computed in
                             ///< DynSmagFilter()
      );

      //! calculate (all-scale) subgrid viscosity
      void CalcSubgrVisc(
          const LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double vol,                          ///< volume
          double&
              Cs_delta_sq,     ///< parameter CS in dynamic Smagorinsky // !or Cv in dynamic Vreman!
          double& Ci_delta_sq  ///< parameter CS in dynamic Smagorinsky for loma
      );

      //! calculate fine-scale subgrid viscosity
      void CalcFineScaleSubgrVisc(
          const LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const LINALG::Matrix<nsd_, nen_>&
              fsevelaf,     ///< fine scale velocity at time n+alpha_f / n+1
          const double vol  ///< volume
      );

      //! get coefficient for multifractal subgrid scales (velocity)
      void CalcMultiFracSubgridVelCoef(const double Csgs,  ///< model coefficient
          const double alpha,                              ///< filter width ratio
          const std::vector<double> Nvel,                  ///< number of casacade steps
          LINALG::Matrix<nsd_, 1>& B_mfs                   ///< final coefficient
      );

      //! get coefficient for multifractal subgrid scales (scalar) (loma only)
      void CalcMultiFracSubgridScaCoef(const double Csgs,  ///< model coefficient
          const double alpha,                              ///< model coefficient
          const double Pr,                                 ///< Prandtl number
          const double Pr_limit,  ///< Prandtl number to distinguish between low and high Prandtl
                                  ///< number regime
          const std::vector<double> Nvel,  ///< number of casacade steps (velocity)
          double Nphi,                     ///< number of casacade steps (scalar)
          double& D_mfs                    ///< final coefficient
      );

      //! Compute element matrix entries: fine scale subgrid viscousity rhs term
      void FineScaleSubGridViscosityTerm(
          LINALG::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& fssgviscfac              ///< = (fine scale subgrid viscousity) x timefacfac
      );

      void MultfracSubGridScalesCross(LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,
          LINALG::Matrix<nsd_, nen_>& velforce, const double& timefacfac, const double& rhsfac);

      void MultfracSubGridScalesReynolds(LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,
          LINALG::Matrix<nsd_, nen_>& velforce, const double& timefacfac, const double& rhsfac);

      void MultfracSubGridScalesConsistentResidual();

      //! loma related methods
      //! definition in fluid_impl_loma_service.cpp

      //! update material parameters including subgrid-scale part of scalar
      void UpdateMaterialParams(
          Teuchos::RCP<const MAT::Material> material,  ///< reference pointer to material
          const LINALG::Matrix<nsd_, nen_>& evelaf,    ///< velocity at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& epreaf,       ///< pressure at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& epream,       ///< pressure at time n+alpha_m / n
          const LINALG::Matrix<nen_, 1>& escaaf,       ///< scalar at time n+alpha_f / n+1
          const LINALG::Matrix<nen_, 1>& escaam,       ///< scalar at time n+alpha_m / n
          const double thermpressaf,  ///< thermodynamic pressure at time n+alpha_f / n+1
          const double thermpressam,  ///< thermodynamic pressure at time n+alpha_m / n
          const double scsca          ///< subgrid scalar at integration point
      );

      //! compute additional Galerkin terms on right-hand side of continuity equation
      //! (only required for variable-density flow at low Mach number)
      void ComputeGalRHSContEq(const LINALG::Matrix<nsd_, nen_>& eveln,  ///< velocity at time n
          const LINALG::Matrix<nen_, 1>& escaaf,    ///< scalar at time n+alpha_F/n+1
          const LINALG::Matrix<nen_, 1>& escaam,    ///< scalar at time n+alpha_M/n
          const LINALG::Matrix<nen_, 1>& escadtam,  ///< acceleration at time n+alpha_M/n
          bool isale                                ///< flag for ALE case
      );

      //! compute additional Galerkin terms on right-hand side of continuity equation
      //! (only required for weakly compressibility)
      void ComputeGalRHSContEqWeakComp(
          const LINALG::Matrix<nen_, 1>& epreaf,    ///< pressure at time n+alpha_F/n+1
          const LINALG::Matrix<nen_, 1>& epredtam,  ///< derivative of pressure at time n+alpha_M/n
          bool isale                                ///< flag for ALE case
      );

      //! compute additional Galerkin terms on right-hand side of continuity equation
      //! (only required for artificial compressibility)
      void ComputeGalRHSContEqArtComp(
          const LINALG::Matrix<nen_, 1>& epreaf,   ///< pressure at time n+alpha_F/n+1
          const LINALG::Matrix<nen_, 1>& epren,    ///< pressure at time n
          const LINALG::Matrix<nen_, 1>& escadtam  ///< acceleration at time n+alpha_M/n
      );

      //! compute residual of scalar equation and subgrid-scale part of scalar
      //! (only required for variable-density flow at low Mach number)
      void ComputeSubgridScaleScalar(
          const LINALG::Matrix<nen_, 1>& escaaf,  ///< scalar at time n+alpha_F/n+1
          const LINALG::Matrix<nen_, 1>& escaam   ///< scalar at time n+alpha_M/n
      );

      //! recompute Galerkin terms based on updated material parameters
      //! including s.-s. part of scalar and compute cross-stress term on
      //! right-hand side of continuity equation
      //! (only required for variable-density flow at low Mach number)
      void RecomputeGalAndComputeCrossRHSContEq();

      //! Compute element matrix entries: LOMA
      void LomaGalPart(
          LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                      ///< = timefac x fac
          const double& rhsfac                           ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: artificial compressibility
      void ArtCompPressureInertiaGalPartandContStab(
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          LINALG::Matrix<nen_, nen_>& ppmat              ///< block (weighting function q x p)
      );

      //! Compute element matrix entries: weak compressibility
      void WeakCompPressureInertiaGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          LINALG::Matrix<nen_, nen_>& ppmat              ///< block (weighting function q x p)
      );

      //! ale related methods
      //! definition in fluid_impl_ale_service.cpp

      //! linearisation in the case of mesh motion 2-D
      virtual void LinMeshMotion_2D(
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,  ///< mesh motion
          const LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double& press,                       ///< pressure at integration point
          const double& timefac,                     ///< time factor
          const double& timefacfac                   ///< = timefac x fac
      );

      //! linearisation in the case of mesh motion 3-D
      virtual void LinMeshMotion_3D(
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,  ///< mesh motion
          const LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double& press,                       ///< pressure at integration point
          const double& timefac,                     ///< time factor
          const double& timefacfac                   ///< = timefac x fac
      );

      void GetPorosityAtGP(const LINALG::Matrix<nen_, 1>& eporo);

      /*!
        \brief calculate rate of strain of (fine-scale) velocity

        \param evel       (i) nodal velocity values
        \param derxy      (i) shape function derivatives
        \param velderxy   (o) velocity derivatives

        \return computed rate of strain
      */
      double GetStrainRate(const LINALG::Matrix<nsd_, nen_>& evel)
      {
        double rateofstrain = 0.0;

        // velderxy is computed here since the evaluation of the strain rate can be performed
        // at the element center before the gauss loop

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
        LINALG::Matrix<nsd_, nsd_> velderxy;
        velderxy.MultiplyNT(evel, derxy_);

        // compute (resolved) rate of strain
        //
        //          +-                                 -+ 1
        //          |          /   \           /   \    | -
        //          | 2 * eps | vel |   * eps | vel |   | 2
        //          |          \   / ij        \   / ij |
        //          +-                                 -+
        //
        LINALG::Matrix<nsd_, nsd_> two_epsilon;
        for (int rr = 0; rr < nsd_; ++rr)
        {
          for (int mm = 0; mm < nsd_; ++mm)
          {
            two_epsilon(rr, mm) = velderxy(rr, mm) + velderxy(mm, rr);
          }
        }

        for (int rr = 0; rr < nsd_; ++rr)
        {
          for (int mm = 0; mm < nsd_; ++mm)
          {
            rateofstrain += two_epsilon(rr, mm) * two_epsilon(mm, rr);
          }
        }

        // sqrt(two_epsilon(rr,mm)*two_epsilon(mm,rr)/4.0*2.0)

        return (sqrt(rateofstrain / 2.0));
      }

      //! output values of Cs, visceff and Cs_delta_sq for statistics
      void StoreModelParametersForOutput(const double Cs_delta_sq, const double Ci_delta_sq,
          const int nlayer, const bool isowned, Teuchos::ParameterList& turbmodelparams);

      /*!
       * \brief fill elment matrix and vectors with the global values
       */
      void ExtractValuesFromGlobalVector(
          const DRT::Discretization& discretization,  ///< discretization
          const std::vector<int>& lm,                 ///<
          FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 1, enrtype>& rotsymmpbc,  ///<
          LINALG::Matrix<nsd_, nen_>* matrixtofill,  ///< vector field
          LINALG::Matrix<nen_, 1>* vectortofill,     ///< scalar field
          const std::string state);                  ///< state of the global vector

      //! identify elements of inflow section
      void InflowElement(DRT::Element* ele);

      // FLD::RotationallySymmetricPeriodicBC<distype> & rotsymmpbc, ///<
      //  {
      //    // get state of the global vector
      //    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
      //    if(matrix_state == Teuchos::null)
      //      dserror("Cannot get state vector %s", state.c_str());
      //
      //    // extract local values of the global vectors
      //    std::vector<double> mymatrix(lm.size());
      //    DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);
      //
      //    // rotate the vector field in the case of rotationally symmetric boundary conditions
      //    if(matrixtofill != NULL)
      //      rotsymmpbc.RotateMyValuesIfNecessary(mymatrix);
      //
      //    for (int inode=0; inode<nen_; ++inode)  // number of nodes
      //    {
      //      // fill a vector field via a pointer
      //      if (matrixtofill != NULL)
      //      {
      //        for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      //        {
      //          (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
      //        }  // end for(idim)
      //      }
      //      // fill a scalar field via a pointer
      //      if (vectortofill != NULL)
      //        (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*numdofpernode_)];
      //    }
      //  }



      //==================================================================================
      // OLD FLUID ELE CALC ROUTINES BEFORE OST-HIST MIGRATION.

      //! calculate div(epsilon(u))
      void CalcDivEps(
          const LINALG::Matrix<nsd_, nen_>& evelaf);  ///< velocity at time n+alpha_f / n+1

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      virtual void InertiaConvectionReactionGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,                   ///< linearisation of the Garlerkin momentum residual
          LINALG::Matrix<nsd_, 1>& resM_Du,  ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac               ///< right-hand-side factor
      );

      //! Compute element matrix entries: for the viscous terms of the Galerkin part
      void ViscousGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          LINALG::Matrix<nsd_, nsd_>& viscstress,             ///< viscous stresses
          const double& timefacfac,                           ///< = timefac x fac
          const double& rhsfac                                ///< right-hand-side factor
      );

      //! Compute element matrix entries: div-grad stabilization and the rhs of the viscous term
      void ContStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefac,                              ///< time factor
          const double& timefacfac,                           ///< = timefac x fac
          const double& timefacfacpre,                        ///< = timefac x fac
          const double& rhsfac                                ///< right-hand-side factor
      );

      //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
      void PressureGalPart(LINALG::Matrix<nen_ * nsd_, nen_>&
                               estif_p_v,        ///< block (weighting function v x pressure)
          LINALG::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,              ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor
          const double& press    ///< pressure at integration point
      );

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      virtual void ContinuityGalPart(
          LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                      ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor
      );

      //! Compute element matrix entries: body force terms on rhs
      void BodyForceRhsTerm(LINALG::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& rhsfac);

      //==================================================================================



      //! for the handling of rotationally symmetric periodic boundary conditions
      Teuchos::RCP<FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 1, enrtype>> rotsymmpbc_;
      //! element id
      int eid_;
      //! Flag to (de)activate higher order elements
      //! elements with only mixed second order derivatives are not counted as higher order elements
      //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
      bool is_higher_order_ele_;
      //! pointer to parameter lists
      DRT::ELEMENTS::FluidEleParameter* fldpara_;
      //! pointer to parameter list for time integration
      DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
      //! element type: nurbs
      bool isNurbs_;
      //! weights for nurbs elements
      LINALG::Matrix<nen_, 1> weights_;
      //! knot vector for nurbs elements
      std::vector<Epetra_SerialDenseVector> myknots_;
      //! Gaussian integration points
      DRT::UTILS::GaussIntegration intpoints_;
      //! identify elements of inflow section
      //! required for turbulence modeling
      bool is_inflow_ele_;

      //========================================================

      //! element stiffness block velocity-velocity
      LINALG::Matrix<nen_ * nsd_, nen_ * nsd_> estif_u_;
      //! element stiffness block pressure-velocity
      LINALG::Matrix<nen_ * nsd_, nen_> estif_p_v_;
      //! element stiffness block velocity-pressure
      LINALG::Matrix<nen_, nen_ * nsd_> estif_q_u_;
      //! element stiffness block pressure-pressure
      LINALG::Matrix<nen_, nen_> ppmat_;

      // definition of vectors
      //! element rhs blocks pressure
      LINALG::Matrix<nen_, 1> preforce_;
      //! element rhs blocks velocity
      LINALG::Matrix<nsd_, nen_> velforce_;

      //! definition of velocity-based momentum residual vectors
      LINALG::Matrix<nsd_ * nsd_, nen_> lin_resM_Du_;
      LINALG::Matrix<nsd_, 1> resM_Du_;


      //========================================================

      //! nodal based quantities for an element
      LINALG::Matrix<nsd_, nen_> ebofoaf_;
      LINALG::Matrix<nsd_, nen_> eprescpgaf_;
      LINALG::Matrix<nen_, 1> escabofoaf_;

      LINALG::Matrix<nsd_, nen_> ebofon_;
      LINALG::Matrix<nsd_, nen_> eprescpgn_;
      LINALG::Matrix<nen_, 1> escabofon_;

      LINALG::Matrix<nsd_, nen_> evelaf_;
      LINALG::Matrix<nen_, 1> epreaf_;
      LINALG::Matrix<nsd_, nen_> evelam_;
      LINALG::Matrix<nen_, 1> epream_;
      LINALG::Matrix<nsd_, nen_> evelnp_;
      LINALG::Matrix<nen_, 1> eprenp_;
      LINALG::Matrix<nsd_, nen_> eveln_;
      LINALG::Matrix<nen_, 1> epren_;
      LINALG::Matrix<nsd_, nen_> eaccam_;
      LINALG::Matrix<nen_, 1> escadtam_;
      LINALG::Matrix<nsd_, nen_> eveldtam_;
      LINALG::Matrix<nen_, 1> epredtam_;
      LINALG::Matrix<nen_, 1> escaaf_;
      LINALG::Matrix<nen_, 1> escaam_;
      LINALG::Matrix<nsd_, nen_> emhist_;
      LINALG::Matrix<nen_, 1> eporo_;

      LINALG::Matrix<nsd_, nen_> gradphiele_;
      LINALG::Matrix<nen_, 1> curvatureele_;
      LINALG::Matrix<nsd_, nen_> gradphielen_;
      LINALG::Matrix<nen_, 1> curvatureelen_;

      LINALG::Matrix<nsd_, 2 * nen_> gradphieletot_;
      LINALG::Matrix<nen_, 2> curvatureeletot_;

      LINALG::Matrix<nsd_, nen_> edispnp_;
      LINALG::Matrix<nsd_, nen_> egridv_;
      LINALG::Matrix<nsd_, nen_> egridvn_;

      LINALG::Matrix<nsd_, nen_> fsevelaf_;
      LINALG::Matrix<nen_, 1> fsescaaf_;

      LINALG::Matrix<nsd_, nen_> evel_hat_;
      LINALG::Matrix<nsd_ * nsd_, nen_> ereynoldsstress_hat_;


      //! node coordinates
      LINALG::Matrix<nsd_, nen_> xyze_;
      //! array for shape functions
      LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      LINALG::Matrix<nsd_, nen_> deriv_;
      //! array for second derivatives of shape function w.r.t r,s,t
      LINALG::Matrix<numderiv2_, nen_> deriv2_;
      //! transposed jacobian "dx/ds"
      LINALG::Matrix<nsd_, nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd_, nsd_> xji_;
      //! global velocity derivatives in gausspoint w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> vderxy_;
      //! global derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<nsd_, nen_> derxy_;
      //! global second derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<numderiv2_, nen_> derxy2_;
      //! bodyforce in gausspoint
      LINALG::Matrix<nsd_, 1> bodyforce_;
      // New One Step Theta variables
      //======================================================
      //! denisty multiplied with instationary term for OST implementations
      double dens_theta_;
      //! bodyforce in gausspoint (n)
      LINALG::Matrix<nsd_, 1> bodyforcen_;
      //! (u^n_old*nabla)u^n_old
      LINALG::Matrix<nsd_, 1> conv_oldn_;
      //! div epsilon(u^n_old)
      LINALG::Matrix<nsd_, 1> visc_oldn_;
      //! pressure gradient in gausspoint (n)
      LINALG::Matrix<nsd_, 1> gradpn_;
      //! velocity vector in gausspoint (n)
      LINALG::Matrix<nsd_, 1> velintn_;
      //! physical viscosity (n)
      double viscn_;
      //! old residual of continuity equation (n)
      double conres_oldn_;
      //! prescribed pressure gradient (required for turbulent channel flow!)
      LINALG::Matrix<nsd_, 1> generalbodyforcen_;
      //======================================================
      //! prescribed pressure gradient (required for turbulent channel flow!)
      LINALG::Matrix<nsd_, 1> generalbodyforce_;
      //! vector containing all values from previous timelevel n for momentum equation
      LINALG::Matrix<nsd_, 1> histmom_;
      //! velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> velint_;
      //! subgrid-scale velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> sgvelint_;
      //! grid velocity u_G at integration point
      LINALG::Matrix<nsd_, 1> gridvelint_;
      //! grid velocity u_G at integration point for new ost
      LINALG::Matrix<nsd_, 1> gridvelintn_;
      //! ale convective velocity c=u-u_G at integration point
      LINALG::Matrix<nsd_, 1> convvelint_;
      //! Oseen advective velocity at element nodes
      LINALG::Matrix<nsd_, nen_> eadvvel_;
      //! acceleration vector in gausspoint
      LINALG::Matrix<nsd_, 1> accint_;
      //! pressure gradient in gausspoint
      LINALG::Matrix<nsd_, 1> gradp_;
      //! the stabilisation parameters -> it is a (3,1) vector for 2D and 3D
      LINALG::Matrix<3, 1> tau_;
      //! viscous term including 2nd derivatives
      //! (This array once had three dimensions, now the first two are combined to one.)
      LINALG::Matrix<nsd_ * nsd_, nen_> viscs2_;
      //! linearisation of convection, convective part
      LINALG::Matrix<nen_, 1> conv_c_;
      //! linearisation of subgrid-scale convection, convective part
      LINALG::Matrix<nen_, 1> sgconv_c_;
      //! velocity divergenceat at t_(n+alpha_F) or t_(n+1)
      double vdiv_;
      //! total right hand side terms at int.-point for momentum equation
      LINALG::Matrix<nsd_, 1> rhsmom_;
      //! (u_old*nabla)u_old
      LINALG::Matrix<nsd_, 1> conv_old_;
      //! div epsilon(u_old)
      LINALG::Matrix<nsd_, 1> visc_old_;
      //! old residual of momentum equation
      LINALG::Matrix<nsd_, 1> momres_old_;
      //! old residual of continuity equation
      double conres_old_;
      //! 2nd derivatives of coord.-functions w.r.t r,s,t
      LINALG::Matrix<numderiv2_, nsd_> xder2_;
      //! global velocity second derivatives in gausspoint w.r.t local coordinates
      LINALG::Matrix<nsd_, nsd_> vderiv_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<nsd_, 1> xsi_;
      //! Jacobian determinant
      double det_;
      //! integration factor
      double fac_;
      //! physical viscosity
      double visc_;
      //! effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      double visceff_;
      //! reaction coefficient
      double reacoeff_;

      //! Two-phase specific variables:
      //! surface tension
      double gamma_;

      //! LOMA-specific variables:
      //! physical diffusivity of scalar equation
      double diffus_;
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
      //! convective velocity vector in gausspoint at t_(n)
      LINALG::Matrix<nsd_, 1> convvelintn_;
      //! global velocity derivatives in gausspoint w.r.t x,y,z at t_(n)
      LINALG::Matrix<nsd_, nsd_> vderxyn_;
      //! velocity divergence at at t_(n)
      double vdivn_;
      //! scalar gradient at t_(n+alpha_F) or t_(n+1)
      LINALG::Matrix<nsd_, 1> grad_scaaf_;
      //! scalar gradient at t_(n)
      LINALG::Matrix<nsd_, 1> grad_scan_;
      //! scalar at t_(n+alpha_F) or t_(n+1)
      double scaaf_;
      //! scalar at t_(n)
      double scan_;
      //! time derivative of scalar term (only required for generalized-alpha scheme)
      double tder_sca_;
      //! convective scalar term at t_(n+alpha_F) or t_(n+1)
      double conv_scaaf_;
      //! convective scalar term at t_(n)
      double conv_scan_;
      //! right-hand side of scalar equation
      double scarhs_;
      //! subgrid-scale part of scalar at integration point
      double sgscaint_;

      //! weakly_compressible-specific variables:
      //! pressure at t_(n+alpha_F) or t_(n+1)
      double preaf_;
      //! pressure at t_(n+alpha_M) or t_(n)
      double pream_;
      //! factor for convective pressure term at t_(n+alpha_F) or t_(n+1)
      double preconvfacaf_;
      //! time derivative of pressure
      double tder_pre_;
      //! factor for pressure time derivative
      double predtfac_;
      //! pressure gradient at t_(n+alpha_F) or t_(n+1)
      LINALG::Matrix<nsd_, 1> grad_preaf_;
      //! convective pressure term at t_(n+alpha_F) or t_(n+1)
      double conv_preaf_;
      //! // element correction term
      LINALG::Matrix<1, nen_> ecorrectionterm_;

      //! turbulence-specific variables:
      //! fine-scale velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> fsvelint_;
      //! fine-scale velocity vector in gausspoint for multifractal subgrid-scale modeling
      LINALG::Matrix<nsd_, 1> mffsvelint_;
      //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> fsvderxy_;
      //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z for multifractal
      //! subgrid-scale modeling
      LINALG::Matrix<nsd_, nsd_> mffsvderxy_;
      //! fine scale velocity divergence for multifractal subgrid-scale modeling
      double mffsvdiv_;
      //! (all-scale) subgrid viscosity
      double sgvisc_;
      //! fine-scale subgrid viscosity
      double fssgvisc_;
      //! model parameter for isotropic part of subgrid-stress tensor (dyn Smag for loma)
      double q_sq_;
      //! multifractal subgrid-scale part of scalar at integration point
      double mfssgscaint_;
      //! gradient of multifractal subgrid-scale scalar (for loma)
      LINALG::Matrix<nsd_, 1> grad_fsscaaf_;

      //! norm of velocity at integration point at time t^{n+1}
      double vel_normnp_;
      //! time-dependent subgrid-scales (pointer to element-specific data)
      Teuchos::RCP<FLD::TDSEleData> tds_;

      // polynomial pressure projection matrices
      double D_;
      LINALG::Matrix<nen_, 1> E_;

      LINALG::Matrix<nsd_ * nsd_, nen_> evelafgrad_;
      LINALG::Matrix<nsd_ * nsd_, nen_> evelngrad_;


      // protected:
      // static std::map<int,std::map<int,DRT::ELEMENTS::FluidEleCalc<distype>* >* > instances_;
    };
    // template<DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType
    // enrtype = DRT::ELEMENTS::Fluid::none>
    // std::map<int,std::map<int,DRT::ELEMENTS::FluidEleCalc<distype>* >* >
    // FluidEleCalc<distype>::instances_;
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
