/*----------------------------------------------------------------------*/
/*!
\file fluid3_interface.H

\brief Internal implementation of fluid3 elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID3_INTERFACE_H
#define FLUID3_INTERFACE_H

#ifdef D_FLUID3

#include "../drt_f3/fluid3.H"
#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

namespace GEO
{
  namespace CUT
  {
    class BoundaryCell;
  }
}

namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid3Impl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid3Impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid3Impl.
   */
  class Fluid3ImplInterface : public DRT::SingletonDestruction
  {
  public:

    /// Empty constructor
    Fluid3ImplInterface() {}

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int Evaluate(DRT::ELEMENTS::Fluid3*    ele,
                         DRT::Discretization & discretization,
                         const std::vector<int> & lm,
                         Teuchos::ParameterList&    params,
                         Teuchos::RCP<MAT::Material> & mat,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra
      ) = 0;

    /// Evaluate the element at specified gauss points
    virtual int Evaluate(DRT::ELEMENTS::Fluid3*    ele,
                         DRT::Discretization & discretization,
                         const std::vector<int> & lm,
                         Teuchos::ParameterList&    params,
                         Teuchos::RCP<MAT::Material> & mat,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra,
                         const DRT::UTILS::GaussIntegration & intpoints ) = 0;

    /// Integrate shape function
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int IntegrateShapeFunction(
        DRT::ELEMENTS::Fluid3*    ele,
        DRT::Discretization&      discretization,
        vector<int>&              lm            ,
        Epetra_SerialDenseVector& elevec1       ) = 0;

    virtual int ComputeError(
        DRT::ELEMENTS::Fluid3*          ele,
        ParameterList&                  params,
        Teuchos::RCP<MAT::Material>&    mat,
        DRT::Discretization&            discretization,
        vector<int>&                    lm,
        Epetra_SerialDenseVector&       elevec1) = 0;

    virtual int CalcDissipation(Fluid3*                    ele,
                                ParameterList&             params,
                                DRT::Discretization&       discretization,
                                vector<int>&               lm,
                                RefCountPtr<MAT::Material> mat) = 0;

    /// add interface condition at cut to element matrix and rhs
    virtual void ElementXfemInterface( DRT::ELEMENTS::Fluid3 * ele,
                                       DRT::Discretization & dis,
                                       const std::vector<int> & lm,
                                       const DRT::UTILS::GaussIntegration & intpoints,
                                       DRT::Discretization & cutdis,
                                       const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                       const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
                                       Teuchos::ParameterList&    params,
                                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                                       Epetra_SerialDenseVector&  elevec1_epetra
      ) = 0;

    /// Internal implementation class for fluid element
    static Fluid3ImplInterface* Impl(DRT::Element::DiscretizationType distype);

  };

}
}

#endif
#endif
