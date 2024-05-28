/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of LinExp artery element

\level 3


*----------------------------------------------------------------------*/


#ifndef FOUR_C_ART_NET_ARTERY_ELE_CALC_LIN_EXP_HPP
#define FOUR_C_ART_NET_ARTERY_ELE_CALC_LIN_EXP_HPP

#include "4C_config.hpp"

#include "4C_art_net_artery.hpp"
#include "4C_art_net_artery_ele_action.hpp"
#include "4C_art_net_artery_ele_calc.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    /// Internal artery implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the artery element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/09
    */

    template <CORE::FE::CellType distype>
    class ArteryEleCalcLinExp : public ArteryEleCalc<distype>
    {
     private:
      typedef ArteryEleCalc<distype> my;

      /// private constructor, since we are a Singleton.
      ArteryEleCalcLinExp(const int numdofpernode, const std::string& disname);

     public:
      //! Singleton access method
      static ArteryEleCalcLinExp<distype>* Instance(
          const int numdofpernode, const std::string& disname);

      int Evaluate(Artery* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) override;

      int ScatraEvaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) override;

      int EvaluateService(Artery* ele, const ARTERY::Action action, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) override;

      /*!
        \brief calculate element matrix and rhs

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) artery material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      void sysmat(Artery* ele, const CORE::LINALG::Matrix<my::iel_, 1>& eqnp,
          const CORE::LINALG::Matrix<my::iel_, 1>& eareanp,
          CORE::LINALG::Matrix<2 * my::iel_, 2 * my::iel_>& sysmat,
          CORE::LINALG::Matrix<2 * my::iel_, 1>& rhs,
          Teuchos::RCP<const CORE::MAT::Material> material, double dt);

      void ScatraSysmat(Artery* ele, const CORE::LINALG::Matrix<2 * my::iel_, 1>& escatran,
          const CORE::LINALG::Matrix<my::iel_, 1>& ewfnp,
          const CORE::LINALG::Matrix<my::iel_, 1>& ewbnp,
          const CORE::LINALG::Matrix<my::iel_, 1>& eareanp,
          const CORE::LINALG::Matrix<my::iel_, 1>& earean,
          CORE::LINALG::Matrix<2 * my::iel_, 2 * my::iel_>& sysmat,
          CORE::LINALG::Matrix<2 * my::iel_, 1>& rhs,
          Teuchos::RCP<const CORE::MAT::Material> material, double dt);

      virtual bool SolveRiemann(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const CORE::MAT::Material> mat);

      virtual void EvaluateTerminalBC(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> mat);

      virtual void EvaluateScatraBC(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& disctretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material);

      virtual void calc_postprocessing_values(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> mat);

      virtual void calc_scatra_from_scatra_fw(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material);

      virtual void EvaluateWfAndWb(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material);

      virtual void solve_scatra_analytically(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material);

      /*!
        \brief get the initial values of the degrees of freedome at the node

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) artery material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      virtual void Initial(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const CORE::MAT::Material> material);

      /*!
       \Essential functions to compute the results of essential matrices
      */
     private:
      //! nodal volumetric flow rate at time step "n"
      CORE::LINALG::Matrix<my::iel_, 1> qn_;
      //! nodal cross-sectional area at time step "n"
      CORE::LINALG::Matrix<my::iel_, 1> an_;
      //! vector containing the initial cross-sectional area at the element nodes
      CORE::LINALG::Matrix<my::iel_, 1> area0_;
      //! vector containing the initial thickness at the element nodes
      CORE::LINALG::Matrix<my::iel_, 1> th_;
      //! vector containing the initial Youngs modulus at the element nodes
      CORE::LINALG::Matrix<my::iel_, 1> young_;
      //! vector containing the fixed external pressure
      CORE::LINALG::Matrix<my::iel_, 1> pext_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
