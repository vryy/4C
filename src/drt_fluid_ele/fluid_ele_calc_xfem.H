/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem.H

\brief Internal implementation of XFluid element interface coupling

<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_H
#define FLUID_ELE_CALC_XFEM_H

#include "fluid_ele_calc.H"


namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcXFEM : public FluidEleCalc<distype>
{
    /// private Constructor since we are a Singleton.
    FluidEleCalcXFEM();

    typedef FluidEleCalc<distype> my;

public:
    virtual ~FluidEleCalcXFEM() {}

    /// Singleton access method
    static FluidEleCalcXFEM<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();

    /// add interface condition at cut to element matrix and rhs
    virtual void ElementXfemInterface( DRT::ELEMENTS::Fluid3 * ele,
                                       DRT::Discretization & dis,
                                       const std::vector<int> & lm,
                                       const DRT::UTILS::GaussIntegration & intpoints,
                                       DRT::Discretization & cutdis,
                                       const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                       const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
                                       std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
                                       //std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
                                       Teuchos::ParameterList&    params,
                                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                                       Epetra_SerialDenseVector&  elevec1_epetra,
                                       Epetra_SerialDenseMatrix&  Cuiui
      );

    /// add interface condition at cut to element matrix and rhs (Nitsche coupling)
    virtual void ElementXfemInterfaceNitsche( DRT::ELEMENTS::Fluid3 * ele,
                                       DRT::Discretization & dis,
                                       const std::vector<int> & lm,
                                       const DRT::UTILS::GaussIntegration & intpoints,
                                       DRT::Discretization & cutdis,
                                       const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                       const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
                                       std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
                                       //std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
                                       Teuchos::ParameterList&    params,
                                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                                       Epetra_SerialDenseVector&  elevec1_epetra,
                                       Epetra_SerialDenseMatrix&  Cuiui
      );

    virtual void ElementXfemInterfaceNitscheTwoSided( DRT::ELEMENTS::Fluid3 * ele,
                                       DRT::Discretization & dis,
                                       const std::vector<int> & lm,
                                       const DRT::UTILS::GaussIntegration & intpoints,
                                       DRT::Discretization & cutdis,
                                       const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                       const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & cutintpoints,
                                       std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
                                       //std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
                                       Teuchos::ParameterList&    params,
                                       DRT::Discretization &  alediscret,
                                       map<int,int> & boundary_emb_gid_map,
                                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                                       Epetra_SerialDenseVector&  elevec1_epetra,
                                       Epetra_SerialDenseMatrix&  Cuiui
      );

    virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid3 * ele,
                                         DRT::Discretization & dis,
                                         const std::vector<int> & lm,
                                         Epetra_SerialDenseVector&  elevec1_epetra,
                                         const DRT::UTILS::GaussIntegration & intpoints
      );

    virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid3 * ele,
                                         DRT::Discretization & dis,
                                         const std::vector<int> & lm,
                                         Epetra_SerialDenseVector&  elevec1_epetra
      );

private:

};
}
}

#endif

