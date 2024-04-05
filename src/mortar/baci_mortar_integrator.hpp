/*-----------------------------------------------------------------------*/
/*! \file
\brief integrate mortar terms
\level 1
*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_INTEGRATOR_HPP
#define FOUR_C_MORTAR_INTEGRATOR_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_inpar_mortar.hpp"
#include "baci_utils_singleton_owner.hpp"

#include <Epetra_Comm.h>

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace DRT
{
  class Element;
}

namespace MORTAR
{
  // forward declarations
  class Element;
  class IntElement;
  class IntCell;


  /*!
  \brief A class to implement MORTAR::IntegratorCalc

  */
  class Integrator
  {
   public:
    Integrator(){};

    virtual ~Integrator() = default;
    //! @name Access methods
    /// Internal implementation class
    static Integrator* Impl(
        MORTAR::Element& sele, MORTAR::Element& mele, Teuchos::ParameterList& params);

    //! @ pure virtual functions --> access per MORTAR::IntegratorCalc
    virtual void IntegrateEleBased2D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, const Epetra_Comm& comm) = 0;

    virtual void IntegrateSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
        MORTAR::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm) = 0;

    virtual Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> IntegrateMmod2D(MORTAR::Element& sele,
        double& sxia, double& sxib, MORTAR::Element& mele, double& mxia, double& mxib) = 0;

    virtual void IntegrateEleBased3D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, const Epetra_Comm& comm) = 0;

    virtual void IntegrateCell3DAuxPlane(MORTAR::Element& sele, MORTAR::Element& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm) = 0;

    virtual void IntegrateCell3DAuxPlaneQuad(MORTAR::Element& sele, MORTAR::Element& mele,
        MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn) = 0;

    virtual int nGP() = 0;

    virtual double Coordinate(int& gp, int dir) = 0;

    virtual double Weight(int& gp) = 0;
  };


  /*!
  \brief A class to perform Gaussian integration and assembly of Mortar
         matrices on the overlap of two MORTAR::Elements (1 Slave, 1 Master)
         in 1D (which is equivalent to a 2D coupling problem) and in 2D
         (which is equivalent to a 3D coupling problem)

  */
  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class IntegratorCalc : public Integrator
  {
   public:
    /*!
    \brief Constructor  with shape function specification

    Constructs an instance of this class using a specific type of shape functions.<br>
    Note that this is \b not a collective call as overlaps are
    integrated in parallel by individual processes.<br>
    Note also that this constructor relies heavily on the
    CORE::FE::IntegrationPoints structs to get Gauss points
    and corresponding weights.

    */
    IntegratorCalc(const Teuchos::ParameterList& params);


    /// Singleton access method
    static IntegratorCalc<distypeS, distypeM>* Instance(
        CORE::UTILS::SingletonAction action, const Teuchos::ParameterList& params);

    //! ns_: number of slave element nodes
    static constexpr int ns_ = CORE::FE::num_nodes<distypeS>;

    //! nm_: number of master element nodes
    static constexpr int nm_ = CORE::FE::num_nodes<distypeM>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = CORE::FE::dim<distypeS> + 1;

    //@}
    //! @name 2D and 3D integration methods
    /*!
    \brief Perform mortar-integration without previous segmentation -- 2D

    */
    void IntegrateEleBased2D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, const Epetra_Comm& comm) override;

    /*!
    \brief Build all integrals and linearizations on a 1D slave /
           master overlap (i.e. D, M, g, LindD, LinM, Ling)

    */
    void IntegrateSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
        MORTAR::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm) override;

    /*!
    \brief Integrate modification Mmod on a 1D slave / master overlap

    This modification is based on a paper by Puso/Wohlmuth, 2005.<br>
    It is necessary in the case of linear slave side elements
    and dual shape functions for the Lagrange multipliers, when
    the interface is curved (but only for mesh tying)!

    */
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> IntegrateMmod2D(MORTAR::Element& sele,
        double& sxia, double& sxib, MORTAR::Element& mele, double& mxia, double& mxib) override;

    /*!
    \brief Build all integrals and linearizations without segmentation
           (i.e. M, g, LinM, Ling and possibly D, LinD)

    */
    void IntegrateEleBased3D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
        bool* boundary_ele, const Epetra_Comm& comm) override;

    /*!
    \brief Build all integrals and linearizations on a 2D slave /
           master integration cell (i.e. D, M, g, LinD, LinM, Ling)
           using a so-called auxiliary plane

    */
    void IntegrateCell3DAuxPlane(MORTAR::Element& sele, MORTAR::Element& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm) override;

    /*!
    \brief Build all integrals and linearizations on a 2D slave /
           master integration cell  (i.e. D, M, g, LinD, LinM, Ling)
           using a so-called auxiliary plane with quadratic interpolation

    */
    void IntegrateCell3DAuxPlaneQuad(MORTAR::Element& sele, MORTAR::Element& mele,
        MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn) override;

    // protected:
    /*!
    \brief Initialize Gauss rule (points, weights) for this MORTAR::Integrator

    */
    void InitializeGP();

    //! @name Access methods
    /*!
    \brief Return number of Gauss points for this instance

    */
    int nGP() override { return ngp_; }

    /*!
    \brief Return coordinates of a specific GP in 1D/2D CElement

    */
    double Coordinate(int& gp, int dir) override { return coords_(gp, dir); }

    /*!
    \brief Return weight of a specific GP in 1D/2D CElement

    */
    double Weight(int& gp) override { return weights_[gp]; }

   private:
    //----------------- GP EVALUATIONS ---------------
    /*!
    \brief evaluate D/M-matrix entries at GP

    */
    void inline GP_DM(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::Matrix<ns_, 1>& lmval, CORE::LINALG::Matrix<ns_, 1>& sval,
        CORE::LINALG::Matrix<nm_, 1>& mval, double& jac, double& wgt, int& nrow, int& ncol,
        int& ndof, bool& bound, const Epetra_Comm& comm);

    /*!
    \brief evaluate D/M-matrix entries at GP (3D and quadratic)

    */
    void inline GP_3D_DM_Quad(MORTAR::Element& sele, MORTAR::Element& mele,
        MORTAR::IntElement& sintele, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& lmintval, CORE::LINALG::Matrix<ns_, 1>& sval,
        CORE::LINALG::Matrix<nm_, 1>& mval, double& jac, double& wgt, int& nrow, int& nintrow,
        int& ncol, int& ndof, bool& bound);
    //@}

    Teuchos::ParameterList imortar_;          // merged parameter list
    INPAR::MORTAR::ShapeFcn shapefcn_;        // lm shape function type
    INPAR::MORTAR::LagMultQuad lmquadtype_;   // type of quadratic lm interpolation
    int ngp_;                                 // number of Gauss points
    CORE::LINALG::SerialDenseMatrix coords_;  // Gauss point coordinates
    std::vector<double> weights_;             // Gauss point weights
  };                                          // class MORTAR::Integrator
}  // namespace MORTAR


BACI_NAMESPACE_CLOSE

#endif  // MORTAR_INTEGRATOR_H
