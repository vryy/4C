
/*----------------------------------------------------------------------*/
/*!
\file inter_acinar_dep_impl.H

\brief Internal implementation of inter_acinar_dep_impl element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef INTER_ACINAR_DEP_IMPL_H
#define INTER_ACINAR_DEP_IMPL_H

#include "red_airway.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
namespace ELEMENTS
{
  /// Interface base class for inter_acinar_dep_impl
  /*!
    This class exists to provide a common interface for all template
    versions of inter_acinar_dep_impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of inter_acinar_dep_impl.
   */
  class RedInterAcinarDepImplInterface
  {
  public:
    /// Empty constructor
    RedInterAcinarDepImplInterface() {}
    /// Empty destructor
    virtual ~RedInterAcinarDepImplInterface() {}
    /// Evaluate the element
    virtual int Evaluate(RedInterAcinarDep*        ele,
                         Teuchos::ParameterList&   params,
                         DRT::Discretization&      discretization,
                         std::vector<int>&         lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         Teuchos::RCP<MAT::Material> mat) = 0;

    virtual void Initial(RedInterAcinarDep*                       ele,
                         Teuchos::ParameterList&                           params,
                         DRT::Discretization&                     discretization,
                         std::vector<int>&                        lm,
                         Teuchos::RCP<const MAT::Material>        material)=0;

    virtual void EvaluateTerminalBC(RedInterAcinarDep*            ele,
                                    Teuchos::ParameterList&       params,
                                    DRT::Discretization&          discretization,
                                    std::vector<int>&             lm,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Teuchos::RCP<MAT::Material>    mat) = 0;

    virtual void CalcFlowRates(RedInterAcinarDep*            ele,
                               Teuchos::ParameterList&       params,
                               DRT::Discretization&          discretization,
                               Epetra_SerialDenseVector&     a_volumen,
                               Epetra_SerialDenseVector&     a_volumenp,
                               std::vector<int>&             lm,
                               Teuchos::RCP<MAT::Material>    mat) = 0;

    virtual void GetCoupledValues(RedInterAcinarDep*           ele,
                                  Teuchos::ParameterList&      params,
                                  DRT::Discretization&         discretization,
                                  std::vector<int>&            lm,
                                  Teuchos::RCP<MAT::Material>   material) = 0;

    /// Internal implementation class for acinus element
    static RedInterAcinarDepImplInterface* Impl(DRT::ELEMENTS::RedInterAcinarDep* acinus);

  };


  /// Internal acinus implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the acinus element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    \author ismail
    \date 01/09
  */

  template <DRT::Element::DiscretizationType distype>
  class InterAcinarDepImpl: public RedInterAcinarDepImplInterface
  {
  public:
    /// Constructor
    explicit InterAcinarDepImpl();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    

    /// Evaluate
    /*!
      The evaluate function for the general acinus case.
     */
    virtual int Evaluate(RedInterAcinarDep*        ele,
                         Teuchos::ParameterList&   params,
                         DRT::Discretization&      discretization,
                         std::vector<int>&         lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         Teuchos::RCP<MAT::Material> mat);

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param eqnp             (i) nodal volumetric flow rate at n+1
      \param evelnp           (i) nodal velocity at n+1
      \param eareanp          (i) nodal cross-sectional area at n+1
      \param eprenp           (i) nodal pressure at n+1
      \param estif            (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) acinus material/dimesion
      \param time             (i) current simulation time
      \param dt               (i) timestep
      */
    void Sysmat(RedInterAcinarDep*                       ele,
                Epetra_SerialDenseVector&                epnp,
                Epetra_SerialDenseVector&                epn,
                Epetra_SerialDenseVector&                epnm,
                Epetra_SerialDenseMatrix&                estif,
                Epetra_SerialDenseVector&                eforce,
                Teuchos::RCP<const MAT::Material>        material,
                Teuchos::ParameterList &                          params,
                double                                   time,
                double                                   dt);


    virtual void EvaluateTerminalBC(RedInterAcinarDep *          ele,
                                    Teuchos::ParameterList&      params,
                                    DRT::Discretization&         disctretization,
                                    std::vector<int>&            lm,
                                    Epetra_SerialDenseVector&    elevec1_epetra,
                                    Teuchos::RCP<MAT::Material>   mat);

    /*!
      \brief get the initial values of the degrees of freedome at the node

      \param ele              (i) the element those matrix is calculated
      \param eqnp             (i) nodal volumetric flow rate at n+1
      \param evelnp           (i) nodal velocity at n+1
      \param eareanp          (i) nodal cross-sectional area at n+1
      \param eprenp           (i) nodal pressure at n+1
      \param estif            (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) acinus material/dimesion
      \param time             (i) current simulation time
      \param dt               (i) timestep
      */
    virtual void Initial(RedInterAcinarDep*                       ele,
                         Teuchos::ParameterList&                           params,
                         DRT::Discretization&                     discretization,
                         std::vector<int>&                        lm,
                         Teuchos::RCP<const MAT::Material>        material);

    /*!
     \Essential functions to compute the results of essentail matrices 
    */
    virtual void CalcFlowRates(RedInterAcinarDep*            ele,
                               Teuchos::ParameterList&       params,
                               DRT::Discretization&          discretization,
                               Epetra_SerialDenseVector&     a_volumen_strain_np,
                               Epetra_SerialDenseVector&     a_volumenp,
                               std::vector<int>&             lm,
                               Teuchos::RCP<MAT::Material>    mat);

    /*!
     \Essential functions to evaluate the coupled results
    */
    virtual void GetCoupledValues(RedInterAcinarDep*           ele,
                                  Teuchos::ParameterList&      params,
                                  DRT::Discretization&         discretization,
                                  std::vector<int>&            lm,
                                  Teuchos::RCP<MAT::Material>   material);

  private:
  };

}
}

#endif

