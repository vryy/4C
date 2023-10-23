/*-----------------------------------------------------------------------*/
/*! \file
\brief integrate mortar terms
\level 1
*/
/*---------------------------------------------------------------------*/
#ifndef BACI_MORTAR_INTEGRATOR_H
#define BACI_MORTAR_INTEGRATOR_H

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_inpar_contact.H"
#include "baci_inpar_mortar.H"
#include "baci_utils_singleton_owner.H"

#include <Epetra_Comm.h>

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
  class MortarElement;
  class IntElement;
  class IntCell;


  /*!
  \brief A class to implement MortarIntegratorCalc

  */
  class MortarIntegrator
  {
   public:
    MortarIntegrator(){};

    virtual ~MortarIntegrator() = default;
    //! @name Access methods
    /// Internal implementation class
    static MortarIntegrator* Impl(
        MortarElement& sele, MortarElement& mele, Teuchos::ParameterList& params);

    //! @ pure virtual functions --> access per MortarIntegratorCalc
    virtual void IntegrateEleBased2D(MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, const Epetra_Comm& comm) = 0;

    virtual void IntegrateSegment2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
        MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm) = 0;

    virtual Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> IntegrateMmod2D(
        MORTAR::MortarElement& sele, double& sxia, double& sxib, MORTAR::MortarElement& mele,
        double& mxia, double& mxib) = 0;

    virtual void IntegrateEleBased3D(MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, const Epetra_Comm& comm) = 0;

    virtual void IntegrateCell3DAuxPlane(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm) = 0;

    virtual void IntegrateCell3DAuxPlaneQuad(MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele, MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn) = 0;

    virtual int nGP() = 0;

    virtual double Coordinate(int& gp, int dir) = 0;

    virtual double Weight(int& gp) = 0;
  };


  /*!
  \brief A class to perform Gaussian integration and assembly of Mortar
         matrices on the overlap of two MortarElements (1 Slave, 1 Master)
         in 1D (which is equivalent to a 2D coupling problem) and in 2D
         (which is equivalent to a 3D coupling problem)

  */
  template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
  class MortarIntegratorCalc : public MortarIntegrator
  {
   public:
    /*!
    \brief Constructor  with shape function specification

    Constructs an instance of this class using a specific type of shape functions.<br>
    Note that this is \b not a collective call as overlaps are
    integrated in parallel by individual processes.<br>
    Note also that this constructor relies heavily on the
    CORE::DRT::UTILS::IntegrationPoints structs to get Gauss points
    and corresponding weights.

    */
    MortarIntegratorCalc(const Teuchos::ParameterList& params);


    /// Singleton access method
    static MortarIntegratorCalc<distypeS, distypeM>* Instance(
        CORE::UTILS::SingletonAction action, const Teuchos::ParameterList& params);

    //! ns_: number of slave element nodes
    static const int ns_ = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;

    //! nm_: number of master element nodes
    static const int nm_ = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static const int ndim_ = CORE::DRT::UTILS::DisTypeToDim<distypeS>::dim + 1;

    //@}
    //! @name 2D and 3D integration methods
    /*!
    \brief Perform mortar-integration without previous segmentation -- 2D

    */
    void IntegrateEleBased2D(MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
        bool* boundary_ele, const Epetra_Comm& comm) override;

    /*!
    \brief Build all integrals and linearizations on a 1D slave /
           master overlap (i.e. D, M, g, LindD, LinM, Ling)

    */
    void IntegrateSegment2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
        MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm) override;

    /*!
    \brief Integrate modification Mmod on a 1D slave / master overlap

    This modification is based on a paper by Puso/Wohlmuth, 2005.<br>
    It is necessary in the case of linear slave side elements
    and dual shape functions for the Lagrange multipliers, when
    the interface is curved (but only for mesh tying)!

    */
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> IntegrateMmod2D(MORTAR::MortarElement& sele,
        double& sxia, double& sxib, MORTAR::MortarElement& mele, double& mxia,
        double& mxib) override;

    /*!
    \brief Build all integrals and linearizations without segmentation
           (i.e. M, g, LinM, Ling and possibly D, LinD)

    */
    void IntegrateEleBased3D(MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
        bool* boundary_ele, const Epetra_Comm& comm) override;

    /*!
    \brief Build all integrals and linearizations on a 2D slave /
           master integration cell (i.e. D, M, g, LinD, LinM, Ling)
           using a so-called auxiliary plane

    */
    void IntegrateCell3DAuxPlane(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm) override;

    /*!
    \brief Build all integrals and linearizations on a 2D slave /
           master integration cell  (i.e. D, M, g, LinD, LinM, Ling)
           using a so-called auxiliary plane with quadratic interpolation

    */
    void IntegrateCell3DAuxPlaneQuad(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn) override;

    // protected:
    /*!
    \brief Initialize Gauss rule (points, weights) for this MortarIntegrator

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
    void inline GP_DM(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        CORE::LINALG::Matrix<ns_, 1>& lmval, CORE::LINALG::Matrix<ns_, 1>& sval,
        CORE::LINALG::Matrix<nm_, 1>& mval, double& jac, double& wgt, int& nrow, int& ncol,
        int& ndof, bool& bound, const Epetra_Comm& comm);

    /*!
    \brief evaluate D/M-matrix entries at GP (3D and quadratic)

    */
    void inline GP_3D_DM_Quad(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
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
  };                                          // class MortarIntegrator
}  // namespace MORTAR


#endif  // MORTAR_INTEGRATOR_H
