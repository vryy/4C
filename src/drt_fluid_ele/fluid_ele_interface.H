/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_interface.H

\brief Interface of fluid elements

<pre>
Maintainer: Volker Gravemeier & Andreas Ehrl
            {vgravem,ehrl}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15245/15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_INTERFACE_H
#define FLUID_ELE_INTERFACE_H

#include <vector>
#include <map>
#include "Teuchos_RefCountPtr.hpp"
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
    virtual int Evaluate(DRT::ELEMENTS::Fluid*        ele,
                         DRT::Discretization &         discretization,
                         const std::vector<int> &      lm,
                         Teuchos::ParameterList&       params,
                         Teuchos::RCP<MAT::Material> & mat,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra,
                         bool                          offdiag = false) = 0;

    /// evaluate element at specified Gauss points
    virtual int Evaluate(DRT::ELEMENTS::Fluid*               ele,
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
                         bool                                 offdiag = false) = 0;

    /// Evaluate the XFEM cut element
    virtual int EvaluateXFEM(DRT::ELEMENTS::Fluid*                             ele,
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
                             std::string&                                      VCellGaussPts,
                             const GEO::CUT::plain_volumecell_set &            cells,
                             bool                                              offdiag = false) = 0;

    /// Integrate shape function
    /*!
      This class does not provide a definition for this function; it
      must be defined in FluidEleCalc.
     */
    virtual int IntegrateShapeFunction(DRT::ELEMENTS::Fluid*    ele,
                                       DRT::Discretization&      discretization,
                                       std::vector<int>&         lm,
                                       Epetra_SerialDenseVector& elevec1) = 0;


    /*! \brief Calculate a integrated divergence operator in vector form
     *
     *   The vector valued operator \f$B\f$ is constructed such that
     *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
     */
    virtual int CalcDivOp(DRT::ELEMENTS::Fluid*     ele,
                          DRT::Discretization&      discretization,
                          std::vector<int>&         lm,
                          Epetra_SerialDenseVector& elevec1) = 0;

    /// error computation
    virtual int ComputeError(DRT::ELEMENTS::Fluid*        ele,
                             Teuchos::ParameterList&       params,
                             Teuchos::RCP<MAT::Material>&  mat,
                             DRT::Discretization&          discretization,
                             std::vector<int>&             lm,
                             Epetra_SerialDenseVector&     elevec1) = 0;

    virtual int ComputeError(DRT::ELEMENTS::Fluid*               ele,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<MAT::Material>&         mat,
                             DRT::Discretization&                 discretization,
                             std::vector<int>&                    lm,
                             Epetra_SerialDenseVector&            elevec1,
                             const DRT::UTILS::GaussIntegration & intpoints2) = 0;

    /// calculation of dissipation
    virtual int CalcDissipation(Fluid*                     ele,
                                Teuchos::ParameterList&     params,
                                DRT::Discretization&        discretization,
                                std::vector<int>&           lm,
                                RefCountPtr<MAT::Material>  mat) = 0;

    virtual int ComputeErrorInterface(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
        DRT::Discretization &                                               dis,               ///< background discretization
        const std::vector<int> &                                            lm,                ///< element local map
        Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
        Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
        DRT::Discretization &                                               cutdis,            ///< cut discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
        Teuchos::ParameterList&                                             params,            ///< parameter list
        const GEO::CUT::plain_volumecell_set&                               vcSet              ///< volumecell sets in this element
    ) = 0;

    /// add interface condition at cut to element matrix and rhs (mixed/stress/hybrid coupling)
    virtual void ElementXfemInterfaceMSH(
        DRT::ELEMENTS::Fluid *                                            ele,
        DRT::Discretization &                                             dis,
        const std::vector<int> &                                          lm,
        const std::vector<DRT::UTILS::GaussIntegration> &                 intpoints,
        DRT::Discretization &                                             cutdis,
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &      bcells,
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,
        Teuchos::ParameterList&                                           params,
        Epetra_SerialDenseMatrix&                                         elemat1_epetra,
        Epetra_SerialDenseVector&                                         elevec1_epetra,
        Epetra_SerialDenseMatrix&                                         Cuiui,
        std::string&                                                      VCellGaussPts,
        const GEO::CUT::plain_volumecell_set &                            cells) = 0;

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
        const GEO::CUT::plain_volumecell_set &                            cells) = 0;

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
        map<int,int> &                                                    boundary_emb_gid_map,
        Epetra_SerialDenseMatrix&                                         elemat1_epetra,
        Epetra_SerialDenseVector&                                         elevec1_epetra,
        Epetra_SerialDenseMatrix&                                         Cuiui,
        const GEO::CUT::plain_volumecell_set &                            cells) = 0;

    virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid *              ele,
                                         DRT::Discretization &                dis,
                                         const std::vector<int> &             lm,
                                         Epetra_SerialDenseVector&            elevec1_epetra,
                                         const DRT::UTILS::GaussIntegration & intpoints) = 0;

    virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid *    ele,
                                         DRT::Discretization &      dis,
                                         const std::vector<int> &   lm,
                                         Epetra_SerialDenseVector&  elevec1_epetra) = 0;

  };

} //namespace Elements
} //namespace DRT

#endif

