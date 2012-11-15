
/*----------------------------------------------------------------------*/
/*!
\file airway_impl.H

\brief Internal implementation of airway_impl element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef AIRWAY_IMPL_H
#define AIRWAY_IMPL_H

#include "red_airway.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
namespace ELEMENTS
{
  /// Interface base class for airway_impl
  /*!
    This class exists to provide a common interface for all template
    versions of airway_impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of airway_impl.
   */
  class RedAirwayImplInterface
  {
  public:
    /// Empty constructor
    RedAirwayImplInterface() {}
    /// Empty destructor
    virtual ~RedAirwayImplInterface() {}
    /// Evaluate the element
    virtual int Evaluate(RedAirway*                ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RCP<MAT::Material> mat) = 0;

    virtual void Initial(RedAirway*                               ele,
                         ParameterList&                           params,
                         DRT::Discretization&                     discretization,
                         vector<int>&                             lm,
                         Teuchos::RCP<const MAT::Material>        material)=0;

    virtual void EvaluateTerminalBC(RedAirway*                    ele,
                                    ParameterList&                params,
                                    DRT::Discretization&          discretization,
                                    vector<int>&                  lm,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    RCP<MAT::Material>    mat) = 0;

    virtual void CalcFlowRates(RedAirway*                    ele,
                               ParameterList&                params,
                               DRT::Discretization&          discretization,
                               Epetra_SerialDenseVector&     a_volumen,
                               Epetra_SerialDenseVector&     a_volumenp,
                               vector<int>&                  lm,
                               RCP<MAT::Material>    mat) = 0;

    virtual void GetCoupledValues(RedAirway*                   ele,
                                  ParameterList&               params,
                                  DRT::Discretization&         discretization,
                                  vector<int>&                 lm,
                                  RCP<MAT::Material>   material) = 0;

    /// Internal implementation class for airway element
    static RedAirwayImplInterface* Impl(DRT::ELEMENTS::RedAirway* airway);

  };


  /// Internal airway implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the airway element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    \author ismail
    \date 01/09
  */

  template <DRT::Element::DiscretizationType distype>
  class AirwayImpl: public RedAirwayImplInterface
  {
  public:
    /// Constructor
    explicit AirwayImpl();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    

    /// Evaluate
    /*!
      The evaluate function for the general airway case.
     */
    virtual int Evaluate(RedAirway*                ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RCP<MAT::Material> mat);

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param eqnp             (i) nodal volumetric flow rate at n+1
      \param evelnp           (i) nodal velocity at n+1
      \param eareanp          (i) nodal cross-sectional area at n+1
      \param eprenp           (i) nodal pressure at n+1
      \param estif            (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) airway material/dimesion
      \param time             (i) current simulation time
      \param dt               (i) timestep
      */
    void Sysmat(RedAirway*                               ele,
                Epetra_SerialDenseVector&                epnp,
                Epetra_SerialDenseVector&                epn,
                Epetra_SerialDenseVector&                epnm,
                Epetra_SerialDenseMatrix&                estif,
                Epetra_SerialDenseVector&                eforce,
                Teuchos::RCP<const MAT::Material>        material,
                ParameterList &                          params,
                double                                   time,
                double                                   dt);


    virtual void EvaluateTerminalBC(RedAirway *                  ele,
                                    ParameterList&               params,
                                    DRT::Discretization&         disctretization,
                                    vector<int>&                 lm,
                                    Epetra_SerialDenseVector&    elevec1_epetra,
                                    RCP<MAT::Material>   mat);

    /*!
      \brief get the initial values of the degrees of freedome at the node

      \param ele              (i) the element those matrix is calculated
      \param eqnp             (i) nodal volumetric flow rate at n+1
      \param evelnp           (i) nodal velocity at n+1
      \param eareanp          (i) nodal cross-sectional area at n+1
      \param eprenp           (i) nodal pressure at n+1
      \param estif            (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) airway material/dimesion
      \param time             (i) current simulation time
      \param dt               (i) timestep
      */
    virtual void Initial(RedAirway*                               ele,
                         ParameterList&                           params,
                         DRT::Discretization&                     discretization,
                         vector<int>&                             lm,
                         Teuchos::RCP<const MAT::Material>        material);

    /*!
     \Essential functions to compute the results of essentail matrices 
    */
    virtual void CalcFlowRates(RedAirway*                    ele,
                               ParameterList&                params,
                               DRT::Discretization&          discretization,
                               Epetra_SerialDenseVector&     a_volumen_strain_np,
                               Epetra_SerialDenseVector&     a_volumenp,
                               vector<int>&                  lm,
                               RCP<MAT::Material>    mat);

    /*!
     \Essential functions to evaluate the coupled results
    */
    virtual void GetCoupledValues(RedAirway*                   ele,
                                  ParameterList&               params,
                                  DRT::Discretization&         discretization,
                                  vector<int>&                 lm,
                                  RCP<MAT::Material>   material);

  private:
  };

}
}

#endif

