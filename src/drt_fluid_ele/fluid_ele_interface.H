/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_interface.H

\brief Interface of fluid elements

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_INTERFACE_H
#define FLUID_ELE_INTERFACE_H

#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "../drt_cut/cut_utils.H"

namespace MAT
{
  class Material;
}

namespace GEO
{
  namespace CUT
  {
    class BoundaryCell;
    class VolumeCell;
  }
}

namespace DRT
{
  class Discretization;

namespace ELEMENTS
{
  class Fluid;

  /// Interface base class for FluidEleCalc
  /*!
    This class exists to provide a common interface for all template
    versions of FluidEleCalc. The only function this class actually defines
    is Ele, which returns a pointer to the appropriate version of FluidEleCalc.
   */
  class FluidEleInterface : public DRT::SingletonDestruction
  {
  public:

    /// Empty constructor
    FluidEleInterface() {}

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function; it
      must be defined in FluidEleCalc.
     */
    virtual int Evaluate(
        DRT::ELEMENTS::Fluid*        ele,
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
    ) = 0;

    /// evaluate element at specified Gauss points
    virtual int Evaluate(
        DRT::ELEMENTS::Fluid*               ele,
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
        bool                                 offdiag = false
    ) = 0;

    /// Evaluate the XFEM cut element
    virtual int EvaluateXFEM(
        DRT::ELEMENTS::Fluid*                             ele,
        DRT::Discretization &                             discretization,
        const std::vector<int> &                          lm,
        Teuchos::ParameterList&                           params,
        Teuchos::RCP<MAT::Material> &                     mat,
        Epetra_SerialDenseMatrix&                         elemat1_epetra,
        Epetra_SerialDenseMatrix&                         elemat2_epetra,
        Epetra_SerialDenseVector&                         elevec1_epetra,
        Epetra_SerialDenseVector&                         elevec2_epetra,
        Epetra_SerialDenseVector&                         elevec3_epetra,
        const std::vector<DRT::UTILS::GaussIntegration> & intpoints,
        const GEO::CUT::plain_volumecell_set &            cells,
        bool                                              offdiag = false
    ) = 0;

    virtual int IntegrateShapeFunction(
        DRT::ELEMENTS::Fluid*                ele,
        DRT::Discretization&                 discretization,
        const std::vector<int>&              lm,
        Epetra_SerialDenseVector&            elevec1,
        const DRT::UTILS::GaussIntegration & intpoints
    ) = 0;

    virtual int IntegrateShapeFunctionXFEM(
        DRT::ELEMENTS::Fluid*                             ele,
        DRT::Discretization&                              discretization,
        const std::vector<int>&                           lm,
        Epetra_SerialDenseVector&                         elevec1,
        const std::vector<DRT::UTILS::GaussIntegration> & intpoints,
        const GEO::CUT::plain_volumecell_set &            cells
    ) = 0;

    /// Evaluate supporting methods of the element
    virtual int EvaluateService(
        DRT::ELEMENTS::Fluid*     ele,
        Teuchos::ParameterList&   params,
        Teuchos::RCP<MAT::Material> & mat,
        DRT::Discretization&      discretization,
        std::vector<int>&         lm,
        Epetra_SerialDenseMatrix& elemat1,
        Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1,
        Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3
    ) = 0;

    virtual int ComputeError(
        DRT::ELEMENTS::Fluid*                ele,
        Teuchos::ParameterList&              params,
        Teuchos::RCP<MAT::Material>&         mat,
        DRT::Discretization&                 discretization,
        std::vector<int>&                    lm,
        Epetra_SerialDenseVector&            elevec1,
        const DRT::UTILS::GaussIntegration & intpoints2
    ) = 0;

    virtual int ComputeErrorInterface(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
        DRT::Discretization &                                               dis,               ///< background discretization
        const std::vector<int> &                                            lm,                ///< element local map
        Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
        Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
        DRT::Discretization &                                               cutdis,            ///< cut discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
        Teuchos::ParameterList&                                             params,            ///< parameter list
        const GEO::CUT::plain_volumecell_set&                               vcSet              ///< volumecell sets in this element
    ) = 0;

    virtual int ComputeErrorInterfaceXFluidFluid(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
        DRT::Discretization &                                               dis,               ///< background discretization
        const std::vector<int> &                                            lm,                ///< element local map
        Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
        Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
        DRT::Discretization &                                               cutdis,            ///< cut discretization
        DRT::Discretization &                                               embdis,            ///< embedded discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
        Teuchos::ParameterList&                                             params,            ///< parameter list
        const GEO::CUT::plain_volumecell_set&                               vcSet,              ///< volumecell sets in this element
        std::map<int,int> &                                                 boundary_emb_gid_map
    ) = 0;

    virtual void ElementXfemInterfaceHybridLM(
        DRT::ELEMENTS::Fluid *                                            ele,                      ///< fluid element
        DRT::Discretization &                                             dis,                      ///< background discretization
        const std::vector<int> &                                          lm,                       ///< element local map
        const std::vector<DRT::UTILS::GaussIntegration> &                 intpoints,                ///< element gauss points
        DRT::Discretization &                                             cutdis,                   ///< boundary discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> >&       bcells,                   ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,               ///< boundary integration points
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,            ///< side coupling matrices
        Teuchos::ParameterList&                                           params,                   ///< parameter list
        Epetra_SerialDenseMatrix&                                         elemat1_epetra,           ///< local system matrix of intersected element
        Epetra_SerialDenseVector&                                         elevec1_epetra,           ///< local element vector of intersected element
        Epetra_SerialDenseMatrix&                                         Cuiui,                    ///< coupling matrix of a side with itself
        const GEO::CUT::plain_volumecell_set &                            vcSet,                    ///< set of plain volume cells
        bool                                                              fluidfluidcoupling = false///< indicates fluid-fluid coupling context
    ) = 0;

    /// add interface condition at cut to element matrix and rhs (Nitsche coupling)
    virtual void ElementXfemInterfaceNIT(
        DRT::ELEMENTS::Fluid *                                            ele,
        DRT::Discretization &                                             dis,
        const std::vector<int> &                                          lm,
        DRT::Discretization &                                             cutdis,
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &      bcells,
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,
        Teuchos::ParameterList&                                           params,
        Epetra_SerialDenseMatrix&                                         elemat1_epetra,
        Epetra_SerialDenseVector&                                         elevec1_epetra,
        Epetra_SerialDenseMatrix&                                         Cuiui,
        const GEO::CUT::plain_volumecell_set &                            cells,
        bool                                                              fluidfluidcoupling = false) = 0;

    /// add interface condition at cut to element matrix and rhs (two-sided Nitsche coupling)
    virtual void ElementXfemInterfaceNIT2(
        DRT::ELEMENTS::Fluid *                                            ele,
        DRT::Discretization &                                             dis,
        const std::vector<int> &                                          lm,
        DRT::Discretization &                                             cutdis,
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &      bcells,
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,
        Teuchos::ParameterList&                                           params,
        DRT::Discretization &                                             alediscret,
        std::map<int,int> &                                               boundary_emb_gid_map,
        Epetra_SerialDenseMatrix&                                         elemat1_epetra,
        Epetra_SerialDenseVector&                                         elevec1_epetra,
        Epetra_SerialDenseMatrix&                                         Cuiui,
        const GEO::CUT::plain_volumecell_set &                            cells) = 0;

    virtual void CalculateContinuityXFEM(
        DRT::ELEMENTS::Fluid *              ele,
        DRT::Discretization &                dis,
        const std::vector<int> &             lm,
        Epetra_SerialDenseVector&            elevec1_epetra,
        const DRT::UTILS::GaussIntegration & intpoints
    ) = 0;

    virtual void CalculateContinuityXFEM(
        DRT::ELEMENTS::Fluid *    ele,
        DRT::Discretization &      dis,
        const std::vector<int> &   lm,
        Epetra_SerialDenseVector&  elevec1_epetra
    ) = 0;

  };

} //namespace Elements
} //namespace DRT

#endif

