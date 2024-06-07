/*----------------------------------------------------------------------*/
/*! \file

\brief integration routines for the volmortar framework

\level 1


*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 01/14 |
 *---------------------------------------------------------------------*/
#ifndef FOUR_C_COUPLING_VOLMORTAR_INTEGRATOR_HPP
#define FOUR_C_COUPLING_VOLMORTAR_INTEGRATOR_HPP

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_coupling_volmortar.hpp"
#include "4C_cut_utils.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace Core::LinAlg
{
  class SerialDenseVector;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Mortar
{
  class IntCell;
}

namespace Core::VolMortar
{
  class Cell;

  template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
  class VolMortarIntegrator
  {
   public:
    // constructor
    VolMortarIntegrator(Teuchos::ParameterList& params);

    // destructor
    virtual ~VolMortarIntegrator() = default;

    //! ns_: number of slave element nodes
    static constexpr int ns_ = Core::FE::num_nodes<distypeS>;

    //! nm_: number of master element nodes
    static constexpr int nm_ = Core::FE::num_nodes<distypeM>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = Core::FE::dim<distypeS>;

    /*!
    \brief Initialize Gauss rule (points, weights)

    */
    void initialize_gp(bool integrateele = false, int domain = 0,
        Core::FE::CellType shape = Core::FE::CellType::tet4);

    /*!
    \brief Integrate cell for 2D problems

    */
    void IntegrateCells2D(Core::Elements::Element& sele, Core::Elements::Element& mele,
        Teuchos::RCP<Mortar::IntCell> cell, Core::LinAlg::SparseMatrix& dmatrix,
        Core::LinAlg::SparseMatrix& mmatrix, Teuchos::RCP<const Discret::Discretization> slavedis,
        Teuchos::RCP<const Discret::Discretization> masterdis, int sdofset, int mdofset);

    /*!
    \brief Integrate cell for 3D problems

    */
    void IntegrateCells3D(Core::Elements::Element& Aele, Core::Elements::Element& Bele,
        Teuchos::RCP<Core::VolMortar::Cell> cell, Core::LinAlg::SparseMatrix& dmatrix_A,
        Core::LinAlg::SparseMatrix& mmatrix_A, Core::LinAlg::SparseMatrix& dmatrix_B,
        Core::LinAlg::SparseMatrix& mmatrix_B, Teuchos::RCP<const Discret::Discretization> Adis,
        Teuchos::RCP<const Discret::Discretization> Bdis, int sdofset_A, int mdofset_A,
        int sdofset_B, int mdofset_B);

    /*!
    \brief Integrate cell for 3D problems

    */
    void integrate_cells3_d_direct_diveregence(Core::Elements::Element& Aele,
        Core::Elements::Element& Bele, Core::Geo::Cut::VolumeCell& vc,
        Teuchos::RCP<Core::FE::GaussPoints> intpoints, bool switched_conf,
        Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
        Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
        Teuchos::RCP<const Discret::Discretization> Adis,
        Teuchos::RCP<const Discret::Discretization> Bdis, int sdofset_A, int mdofset_A,
        int sdofset_B, int mdofset_B);

    /*!
    \brief Integrate ele for 3D problems

    */
    void IntegrateEle3D(int domain, Core::Elements::Element& Aele, Core::Elements::Element& Bele,
        Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
        Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
        Teuchos::RCP<const Discret::Discretization> Adis,
        Teuchos::RCP<const Discret::Discretization> Bdis, int sdofset_A, int mdofset_A,
        int sdofset_B, int mdofset_B);

    /*!
    \brief Integrate ele for 3D problems

    */
    void integrate_ele_based3_d_a_dis(Core::Elements::Element& Aele, std::vector<int>& foundeles,
        Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
        Teuchos::RCP<const Discret::Discretization> Adiscret,
        Teuchos::RCP<const Discret::Discretization> Bdiscret, int dofsetA, int dofsetB);

    /*!
    \brief Integrate ele for 3D problems

    */
    void integrate_ele_based3_d_b_dis(Core::Elements::Element& Bele, std::vector<int>& foundeles,
        Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
        Teuchos::RCP<const Discret::Discretization> Adiscret,
        Teuchos::RCP<const Discret::Discretization> Bdiscret, int dofsetA, int dofsetB);

   protected:
    /*!
    \brief Check integration point mapping (2D)

    */
    bool check_mapping2_d(
        Core::Elements::Element& sele, Core::Elements::Element& mele, double* sxi, double* mxi);

    /*!
    \brief Check integration point mapping (3D)

    */
    bool check_mapping3_d(
        Core::Elements::Element& sele, Core::Elements::Element& mele, double* sxi, double* mxi);

    //@}
    int ngp_;                                 // number of Gauss points
    Core::LinAlg::SerialDenseMatrix coords_;  // Gauss point coordinates
    std::vector<double> weights_;             // Gauss point weights

    // input parameter
    DualQuad dualquad_;  // type of quadratic weighting interpolation
    Shapefcn shape_;     // type of test function
  };

  template <Core::FE::CellType distypeS>
  class VolMortarIntegratorEleBased
  {
   public:
    // constructor
    VolMortarIntegratorEleBased(Teuchos::ParameterList& params);

    // destructor
    virtual ~VolMortarIntegratorEleBased() = default;

    //! ns_: number of slave element nodes
    static constexpr int ns_ = Core::FE::num_nodes<distypeS>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = Core::FE::dim<distypeS>;

    /*!
    \brief Initialize Gauss rule (points, weights)

    */
    void initialize_gp();

    /*!
    \brief Integrate ele for 3D problems

    */
    void IntegrateEleBased3D(Core::Elements::Element& sele, std::vector<int>& foundeles,
        Core::LinAlg::SparseMatrix& dmatrixA, Core::LinAlg::SparseMatrix& mmatrixA,
        Teuchos::RCP<const Discret::Discretization> Adiscret,
        Teuchos::RCP<const Discret::Discretization> Bdiscret, int dofseta, int dofsetb,
        const Teuchos::RCP<const Epetra_Map>& PAB_dofrowmap,
        const Teuchos::RCP<const Epetra_Map>& PAB_dofcolmap);

   protected:
    /*!
    \brief Check integration point mapping (3D)

    */
    bool check_mapping3_d(
        Core::Elements::Element& sele, Core::Elements::Element& mele, double* sxi, double* mxi);

    //@}
    int ngp_;                                 // number of Gauss points
    Core::LinAlg::SerialDenseMatrix coords_;  // Gauss point coordinates
    std::vector<double> weights_;             // Gauss point weights

    // input parameter
    DualQuad dualquad_;  // type of quadratic weighting interpolation
    Shapefcn shape_;     // type of test function
  };

  /*----------------------------------------------------------------------*
   |  Check paraspace mapping                                  farah 02/15|
   *----------------------------------------------------------------------*/
  template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
  bool CheckMapping(
      Core::Elements::Element& sele, Core::Elements::Element& mele, double* sxi, double* mxi)
  {
    // check GP projection (SLAVE)
    double tol = 1e-5;
    if (distypeS == Core::FE::CellType::hex8 || distypeS == Core::FE::CellType::hex20 ||
        distypeS == Core::FE::CellType::hex27)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[2] < -1.0 - tol || sxi[0] > 1.0 + tol ||
          sxi[1] > 1.0 + tol || sxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == Core::FE::CellType::pyramid5)
    {
      if (sxi[2] < 0.0 - tol || -sxi[0] + sxi[2] > 1.0 + tol || sxi[0] + sxi[2] > 1.0 + tol ||
          -sxi[1] + sxi[2] > 1.0 + tol || sxi[1] + sxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == Core::FE::CellType::tet4 || distypeS == Core::FE::CellType::tet10)
    {
      if (sxi[0] < 0.0 - tol || sxi[1] < 0.0 - tol || sxi[2] < 0.0 - tol ||
          (sxi[0] + sxi[1] + sxi[2]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == Core::FE::CellType::tri3 || distypeS == Core::FE::CellType::tri6)
    {
      if (sxi[0] < 0.0 - tol || sxi[1] < 0.0 - tol || (sxi[0] + sxi[1]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == Core::FE::CellType::quad4 || distypeS == Core::FE::CellType::quad8 ||
             distypeS == Core::FE::CellType::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        return false;
      }
    }
    else
      FOUR_C_THROW("Wrong element type!");

    // check GP projection (MASTER)
    if (distypeM == Core::FE::CellType::hex8 || distypeM == Core::FE::CellType::hex20 ||
        distypeM == Core::FE::CellType::hex27)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[2] < -1.0 - tol || mxi[0] > 1.0 + tol ||
          mxi[1] > 1.0 + tol || mxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == Core::FE::CellType::pyramid5)
    {
      if (mxi[2] < 0.0 - tol || -mxi[0] + mxi[2] > 1.0 + tol || mxi[0] + mxi[2] > 1.0 + tol ||
          -mxi[1] + mxi[2] > 1.0 + tol || mxi[1] + mxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == Core::FE::CellType::tet4 || distypeM == Core::FE::CellType::tet10)
    {
      if (mxi[0] < 0.0 - tol || mxi[1] < 0.0 - tol || mxi[2] < 0.0 - tol ||
          (mxi[0] + mxi[1] + mxi[2]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == Core::FE::CellType::tri3 || distypeM == Core::FE::CellType::tri6)
    {
      if (mxi[0] < 0.0 - tol || mxi[1] < 0.0 - tol || (mxi[0] + mxi[1]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == Core::FE::CellType::quad4 || distypeM == Core::FE::CellType::quad8 ||
             distypeM == Core::FE::CellType::quad9)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
      {
        return false;
      }
    }
    else
      FOUR_C_THROW("Wrong element type!");

    return true;
  };

  /*----------------------------------------------------------------------*
   |  Check mapping                                            farah 06/14|
   *----------------------------------------------------------------------*/
  template <Core::FE::CellType distype>
  bool CheckMapping(Core::Elements::Element& ele, double* xi)
  {
    // check node projection
    double tol = 1e-5;
    if (distype == Core::FE::CellType::hex8 || distype == Core::FE::CellType::hex20 ||
        distype == Core::FE::CellType::hex27)
    {
      if (xi[0] < -1.0 - tol || xi[1] < -1.0 - tol || xi[2] < -1.0 - tol || xi[0] > 1.0 + tol ||
          xi[1] > 1.0 + tol || xi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == Core::FE::CellType::pyramid5)
    {
      if (xi[2] < 0.0 - tol || -xi[0] + xi[2] > 1.0 + tol || xi[0] + xi[2] > 1.0 + tol ||
          -xi[1] + xi[2] > 1.0 + tol || xi[1] + xi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == Core::FE::CellType::tet4 || distype == Core::FE::CellType::tet10)
    {
      if (xi[0] < 0.0 - tol || xi[1] < 0.0 - tol || xi[2] < 0.0 - tol ||
          (xi[0] + xi[1] + xi[2]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == Core::FE::CellType::quad4 || distype == Core::FE::CellType::quad8 ||
             distype == Core::FE::CellType::quad9)
    {
      if (xi[0] < -1.0 - tol || xi[1] < -1.0 - tol || xi[0] > 1.0 + tol || xi[1] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == Core::FE::CellType::tri3 || distype == Core::FE::CellType::tri6)
    {
      if (xi[0] < -tol || xi[1] < -tol || xi[0] > 1.0 + tol || xi[1] > 1.0 + tol ||
          xi[0] + xi[1] > 1.0 + 2 * tol)
      {
        return false;
      }
    }
    else
    {
      FOUR_C_THROW("ERROR: Wrong element type!");
    }

    return true;
  }

  //===================================
  template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
  bool VolMortarEleBasedGP(Core::Elements::Element& sele, Core::Elements::Element* mele,
      std::vector<int>& foundeles, int& found, int& gpid, double& jac, double& wgt, double& gpdist,
      double* Axi, double* AuxXi, double* globgp, DualQuad& dq, Shapefcn& shape,
      Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
      Teuchos::RCP<const Discret::Discretization> Adis,
      Teuchos::RCP<const Discret::Discretization> Bdis, int dofseta, int dofsetb,
      const Teuchos::RCP<const Epetra_Map>& PAB_dofrowmap,
      const Teuchos::RCP<const Epetra_Map>& PAB_dofcolmap);

  // evaluation of nts approach
  template <Core::FE::CellType distype>
  bool ConsInterpolatorEval(Core::Nodes::Node* node, Core::Elements::Element* ele,
      Core::LinAlg::SparseMatrix& pmatrix, Teuchos::RCP<const Discret::Discretization> nodediscret,
      Teuchos::RCP<const Discret::Discretization> elediscret, std::vector<int>& foundeles,
      int& found, int& eleid, double& dist, double* AuxXi, double* nodepos,
      std::pair<int, int>& dofset, const Teuchos::RCP<const Epetra_Map>& P_dofrowmap,
      const Teuchos::RCP<const Epetra_Map>& P_dofcolmap);

  //===================================
  // Alternative to VolMortar coupling:
  class ConsInterpolator
  {
   public:
    // constructor
    ConsInterpolator();

    // destructor
    virtual ~ConsInterpolator() = default;

    /*!
    \brief interpolation functionality

    */
    void Interpolate(Core::Nodes::Node* node, Core::LinAlg::SparseMatrix& pmatrix_,
        Teuchos::RCP<const Discret::Discretization> nodediscret,
        Teuchos::RCP<const Discret::Discretization> elediscret, std::vector<int>& foundeles,
        std::pair<int, int>& dofset, const Teuchos::RCP<const Epetra_Map>& P_dofrowmap,
        const Teuchos::RCP<const Epetra_Map>& P_dofcolmap);
  };

}  // namespace Core::VolMortar

FOUR_C_NAMESPACE_CLOSE

#endif
