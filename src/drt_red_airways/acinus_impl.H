
/*----------------------------------------------------------------------*/
/*!
\file acinus_impl.H

\brief Internal implementation of acinus_impl element

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef ACINUS_IMPL_H
#define ACINUS_IMPL_H

#include "red_airway.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for acinus_impl
    /*!
      This class exists to provide a common interface for all template
      versions of acinus_impl. The only function this class actually
      defines is Impl, which returns a pointer to the appropriate version
      of acinus_impl.
    */
    class RedAcinusImplInterface
    {
     public:
      /// Empty constructor
      RedAcinusImplInterface() {}
      /// Empty destructor
      virtual ~RedAcinusImplInterface() {}
      /// Evaluate the element
      virtual int Evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void Initial(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcElemVolume(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void GetCoupledValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void GetJunctionVolumeMix(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatra_np,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatraBifurcations(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatranp,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateElem12Scatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalPO2FromScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalNodalEssentialValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& nodal_surface,
          Epetra_SerialDenseVector& nodal_volume, Epetra_SerialDenseVector& nodal_flow,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) = 0;

      /// Internal implementation class for acinus element
      static RedAcinusImplInterface* Impl(DRT::ELEMENTS::RedAcinus* acinus);
    };


    /// Internal acinus implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the acinus element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/13
    */

    template <DRT::Element::DiscretizationType distype>
    class AcinusImpl : public RedAcinusImplInterface
    {
     public:
      /// Constructor
      explicit AcinusImpl();

      //! number of nodes
      static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;


      /// Evaluate
      /*!
        The evaluate function for the general acinus case.
      */
      virtual int Evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat);

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
      void Sysmat(RedAcinus* ele, Epetra_SerialDenseVector& epnp, Epetra_SerialDenseVector& epn,
          Epetra_SerialDenseVector& epnm, Epetra_SerialDenseMatrix& estif,
          Epetra_SerialDenseVector& eforce, Teuchos::RCP<const MAT::Material> material,
          Teuchos::ParameterList& params, double time, double dt);


      virtual void EvaluateTerminalBC(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& disctretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat);

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
      virtual void Initial(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material);

      /*!
        \Essential functions to compute the results of essentail matrices
      */
      virtual void CalcFlowRates(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      /*!
        \Essential functions to compute the volume of an element
      */
      virtual void CalcElemVolume(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat);

      /*!
        \Essential functions to evaluate the coupled results
      */
      virtual void GetCoupledValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      /*!
        \Essential functions to evaluate mixed volume flowing into a junction
      */
      virtual void GetJunctionVolumeMix(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material);

      /*!
        \Essential functions to solve the forward scatra flow
      */
      virtual void SolveScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatra_np,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void SolveScatraBifurcations(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& scatranp,
          Epetra_SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void UpdateScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void UpdateElem12Scatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void EvalPO2FromScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material);

      virtual void EvalNodalEssentialValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& nodal_surface,
          Epetra_SerialDenseVector& nodal_volume, Epetra_SerialDenseVector& nodal_flow,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material);


     private:
    };
  }  // namespace ELEMENTS
}  // namespace DRT


#endif
