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
#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace MORTAR
{
  class IntCell;
}

namespace CORE::VOLMORTAR
{
  class Cell;

  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class VolMortarIntegrator
  {
   public:
    // constructor
    VolMortarIntegrator(Teuchos::ParameterList& params);

    // destructor
    virtual ~VolMortarIntegrator() = default;

    //! ns_: number of slave element nodes
    static constexpr int ns_ = CORE::FE::num_nodes<distypeS>;

    //! nm_: number of master element nodes
    static constexpr int nm_ = CORE::FE::num_nodes<distypeM>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = CORE::FE::dim<distypeS>;

    /*!
    \brief Initialize Gauss rule (points, weights)

    */
    void InitializeGP(bool integrateele = false, int domain = 0,
        CORE::FE::CellType shape = CORE::FE::CellType::tet4);

    /*!
    \brief Integrate cell for 2D problems

    */
    void IntegrateCells2D(DRT::Element& sele, DRT::Element& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, CORE::LINALG::SparseMatrix& dmatrix,
        CORE::LINALG::SparseMatrix& mmatrix, Teuchos::RCP<const DRT::Discretization> slavedis,
        Teuchos::RCP<const DRT::Discretization> masterdis, int sdofset, int mdofset);

    /*!
    \brief Integrate cell for 3D problems

    */
    void IntegrateCells3D(DRT::Element& Aele, DRT::Element& Bele,
        Teuchos::RCP<CORE::VOLMORTAR::Cell> cell, CORE::LINALG::SparseMatrix& dmatrix_A,
        CORE::LINALG::SparseMatrix& mmatrix_A, CORE::LINALG::SparseMatrix& dmatrix_B,
        CORE::LINALG::SparseMatrix& mmatrix_B, Teuchos::RCP<const DRT::Discretization> Adis,
        Teuchos::RCP<const DRT::Discretization> Bdis, int sdofset_A, int mdofset_A, int sdofset_B,
        int mdofset_B);

    /*!
    \brief Integrate cell for 3D problems

    */
    void IntegrateCells3D_DirectDiveregence(DRT::Element& Aele, DRT::Element& Bele,
        CORE::GEO::CUT::VolumeCell& vc, Teuchos::RCP<CORE::FE::GaussPoints> intpoints,
        bool switched_conf, CORE::LINALG::SparseMatrix& dmatrix_A,
        CORE::LINALG::SparseMatrix& mmatrix_A, CORE::LINALG::SparseMatrix& dmatrix_B,
        CORE::LINALG::SparseMatrix& mmatrix_B, Teuchos::RCP<const DRT::Discretization> Adis,
        Teuchos::RCP<const DRT::Discretization> Bdis, int sdofset_A, int mdofset_A, int sdofset_B,
        int mdofset_B);

    /*!
    \brief Integrate ele for 3D problems

    */
    void IntegrateEle3D(int domain, DRT::Element& Aele, DRT::Element& Bele,
        CORE::LINALG::SparseMatrix& dmatrix_A, CORE::LINALG::SparseMatrix& mmatrix_A,
        CORE::LINALG::SparseMatrix& dmatrix_B, CORE::LINALG::SparseMatrix& mmatrix_B,
        Teuchos::RCP<const DRT::Discretization> Adis, Teuchos::RCP<const DRT::Discretization> Bdis,
        int sdofset_A, int mdofset_A, int sdofset_B, int mdofset_B);

    /*!
    \brief Integrate ele for 3D problems

    */
    void IntegrateEleBased3D_ADis(DRT::Element& Aele, std::vector<int>& foundeles,
        CORE::LINALG::SparseMatrix& dmatrix_A, CORE::LINALG::SparseMatrix& mmatrix_A,
        Teuchos::RCP<const DRT::Discretization> Adiscret,
        Teuchos::RCP<const DRT::Discretization> Bdiscret, int dofsetA, int dofsetB);

    /*!
    \brief Integrate ele for 3D problems

    */
    void IntegrateEleBased3D_BDis(DRT::Element& Bele, std::vector<int>& foundeles,
        CORE::LINALG::SparseMatrix& dmatrix_B, CORE::LINALG::SparseMatrix& mmatrix_B,
        Teuchos::RCP<const DRT::Discretization> Adiscret,
        Teuchos::RCP<const DRT::Discretization> Bdiscret, int dofsetA, int dofsetB);

   protected:
    /*!
    \brief Check integration point mapping (2D)

    */
    bool CheckMapping2D(DRT::Element& sele, DRT::Element& mele, double* sxi, double* mxi);

    /*!
    \brief Check integration point mapping (3D)

    */
    bool CheckMapping3D(DRT::Element& sele, DRT::Element& mele, double* sxi, double* mxi);

    //@}
    int ngp_;                                 // number of Gauss points
    CORE::LINALG::SerialDenseMatrix coords_;  // Gauss point coordinates
    std::vector<double> weights_;             // Gauss point weights

    // input parameter
    DualQuad dualquad_;  // type of quadratic weighting interpolation
    Shapefcn shape_;     // type of test function
  };

  template <CORE::FE::CellType distypeS>
  class VolMortarIntegratorEleBased
  {
   public:
    // constructor
    VolMortarIntegratorEleBased(Teuchos::ParameterList& params);

    // destructor
    virtual ~VolMortarIntegratorEleBased() = default;

    //! ns_: number of slave element nodes
    static constexpr int ns_ = CORE::FE::num_nodes<distypeS>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = CORE::FE::dim<distypeS>;

    /*!
    \brief Initialize Gauss rule (points, weights)

    */
    void InitializeGP();

    /*!
    \brief Integrate ele for 3D problems

    */
    void IntegrateEleBased3D(DRT::Element& sele, std::vector<int>& foundeles,
        CORE::LINALG::SparseMatrix& dmatrixA, CORE::LINALG::SparseMatrix& mmatrixA,
        Teuchos::RCP<const DRT::Discretization> Adiscret,
        Teuchos::RCP<const DRT::Discretization> Bdiscret, int dofseta, int dofsetb,
        const Teuchos::RCP<const Epetra_Map>& PAB_dofrowmap,
        const Teuchos::RCP<const Epetra_Map>& PAB_dofcolmap);

   protected:
    /*!
    \brief Check integration point mapping (3D)

    */
    bool CheckMapping3D(DRT::Element& sele, DRT::Element& mele, double* sxi, double* mxi);

    //@}
    int ngp_;                                 // number of Gauss points
    CORE::LINALG::SerialDenseMatrix coords_;  // Gauss point coordinates
    std::vector<double> weights_;             // Gauss point weights

    // input parameter
    DualQuad dualquad_;  // type of quadratic weighting interpolation
    Shapefcn shape_;     // type of test function
  };

  /*----------------------------------------------------------------------*
   |  Check paraspace mapping                                  farah 02/15|
   *----------------------------------------------------------------------*/
  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  bool CheckMapping(DRT::Element& sele, DRT::Element& mele, double* sxi, double* mxi)
  {
    // check GP projection (SLAVE)
    double tol = 1e-5;
    if (distypeS == CORE::FE::CellType::hex8 || distypeS == CORE::FE::CellType::hex20 ||
        distypeS == CORE::FE::CellType::hex27)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[2] < -1.0 - tol || sxi[0] > 1.0 + tol ||
          sxi[1] > 1.0 + tol || sxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == CORE::FE::CellType::pyramid5)
    {
      if (sxi[2] < 0.0 - tol || -sxi[0] + sxi[2] > 1.0 + tol || sxi[0] + sxi[2] > 1.0 + tol ||
          -sxi[1] + sxi[2] > 1.0 + tol || sxi[1] + sxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == CORE::FE::CellType::tet4 || distypeS == CORE::FE::CellType::tet10)
    {
      if (sxi[0] < 0.0 - tol || sxi[1] < 0.0 - tol || sxi[2] < 0.0 - tol ||
          (sxi[0] + sxi[1] + sxi[2]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == CORE::FE::CellType::tri3 || distypeS == CORE::FE::CellType::tri6)
    {
      if (sxi[0] < 0.0 - tol || sxi[1] < 0.0 - tol || (sxi[0] + sxi[1]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeS == CORE::FE::CellType::quad4 || distypeS == CORE::FE::CellType::quad8 ||
             distypeS == CORE::FE::CellType::quad9)
    {
      if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
      {
        return false;
      }
    }
    else
      FOUR_C_THROW("Wrong element type!");

    // check GP projection (MASTER)
    if (distypeM == CORE::FE::CellType::hex8 || distypeM == CORE::FE::CellType::hex20 ||
        distypeM == CORE::FE::CellType::hex27)
    {
      if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[2] < -1.0 - tol || mxi[0] > 1.0 + tol ||
          mxi[1] > 1.0 + tol || mxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == CORE::FE::CellType::pyramid5)
    {
      if (mxi[2] < 0.0 - tol || -mxi[0] + mxi[2] > 1.0 + tol || mxi[0] + mxi[2] > 1.0 + tol ||
          -mxi[1] + mxi[2] > 1.0 + tol || mxi[1] + mxi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == CORE::FE::CellType::tet4 || distypeM == CORE::FE::CellType::tet10)
    {
      if (mxi[0] < 0.0 - tol || mxi[1] < 0.0 - tol || mxi[2] < 0.0 - tol ||
          (mxi[0] + mxi[1] + mxi[2]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == CORE::FE::CellType::tri3 || distypeM == CORE::FE::CellType::tri6)
    {
      if (mxi[0] < 0.0 - tol || mxi[1] < 0.0 - tol || (mxi[0] + mxi[1]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distypeM == CORE::FE::CellType::quad4 || distypeM == CORE::FE::CellType::quad8 ||
             distypeM == CORE::FE::CellType::quad9)
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
  template <CORE::FE::CellType distype>
  bool CheckMapping(DRT::Element& ele, double* xi)
  {
    // check node projection
    double tol = 1e-5;
    if (distype == CORE::FE::CellType::hex8 || distype == CORE::FE::CellType::hex20 ||
        distype == CORE::FE::CellType::hex27)
    {
      if (xi[0] < -1.0 - tol || xi[1] < -1.0 - tol || xi[2] < -1.0 - tol || xi[0] > 1.0 + tol ||
          xi[1] > 1.0 + tol || xi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == CORE::FE::CellType::pyramid5)
    {
      if (xi[2] < 0.0 - tol || -xi[0] + xi[2] > 1.0 + tol || xi[0] + xi[2] > 1.0 + tol ||
          -xi[1] + xi[2] > 1.0 + tol || xi[1] + xi[2] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == CORE::FE::CellType::tet4 || distype == CORE::FE::CellType::tet10)
    {
      if (xi[0] < 0.0 - tol || xi[1] < 0.0 - tol || xi[2] < 0.0 - tol ||
          (xi[0] + xi[1] + xi[2]) > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == CORE::FE::CellType::quad4 || distype == CORE::FE::CellType::quad8 ||
             distype == CORE::FE::CellType::quad9)
    {
      if (xi[0] < -1.0 - tol || xi[1] < -1.0 - tol || xi[0] > 1.0 + tol || xi[1] > 1.0 + tol)
      {
        return false;
      }
    }
    else if (distype == CORE::FE::CellType::tri3 || distype == CORE::FE::CellType::tri6)
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
  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  bool VolMortarEleBasedGP(DRT::Element& sele, DRT::Element* mele, std::vector<int>& foundeles,
      int& found, int& gpid, double& jac, double& wgt, double& gpdist, double* Axi, double* AuxXi,
      double* globgp, DualQuad& dq, Shapefcn& shape, CORE::LINALG::SparseMatrix& dmatrix_A,
      CORE::LINALG::SparseMatrix& mmatrix_A, Teuchos::RCP<const DRT::Discretization> Adis,
      Teuchos::RCP<const DRT::Discretization> Bdis, int dofseta, int dofsetb,
      const Teuchos::RCP<const Epetra_Map>& PAB_dofrowmap,
      const Teuchos::RCP<const Epetra_Map>& PAB_dofcolmap);

  // evaluation of nts approach
  template <CORE::FE::CellType distype>
  bool ConsInterpolatorEval(DRT::Node* node, DRT::Element* ele, CORE::LINALG::SparseMatrix& pmatrix,
      Teuchos::RCP<const DRT::Discretization> nodediscret,
      Teuchos::RCP<const DRT::Discretization> elediscret, std::vector<int>& foundeles, int& found,
      int& eleid, double& dist, double* AuxXi, double* nodepos, std::pair<int, int>& dofset,
      const Teuchos::RCP<const Epetra_Map>& P_dofrowmap,
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
    void Interpolate(DRT::Node* node, CORE::LINALG::SparseMatrix& pmatrix_,
        Teuchos::RCP<const DRT::Discretization> nodediscret,
        Teuchos::RCP<const DRT::Discretization> elediscret, std::vector<int>& foundeles,
        std::pair<int, int>& dofset, const Teuchos::RCP<const Epetra_Map>& P_dofrowmap,
        const Teuchos::RCP<const Epetra_Map>& P_dofcolmap);
  };

}  // namespace CORE::VOLMORTAR

FOUR_C_NAMESPACE_CLOSE

#endif
