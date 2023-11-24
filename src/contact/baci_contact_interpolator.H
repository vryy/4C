/*----------------------------------------------------------------------*/
/*! \file

\level 2


*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 09/14 |
 *---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_INTERPOLATOR_H
#define BACI_CONTACT_INTERPOLATOR_H

/*---------------------------------------------------------------------*
 | headers                                                 farah 09/14 |
 *---------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_inpar_wear.H"
#include "baci_utils_pairedvector.H"
#include "baci_utils_singleton_owner.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 09/14 |
 *---------------------------------------------------------------------*/
namespace MORTAR
{
  class MortarElement;
  class MortarNode;
}  // namespace MORTAR

namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace CONTACT
{
  class CoNode;
}

namespace NTS
{
  class CoInterpolator
  {
   public:
    /*!
     \brief Constructor

     */
    CoInterpolator(Teuchos::ParameterList& params, const int& dim);

    /*!
     \brief Destructor

     */
    virtual ~CoInterpolator() = default;

    /*!
     \brief Interpolate for nts algorithm

     */
    bool Interpolate(MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles);

    /*!
     \brief Interpolate temperature of master side at a slave node
     for 3D problems

     */
    void InterpolateMasterTemp3D(
        MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles);

    /*!
     \brief lin 3D projection

     */
    void DerivXiGP3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, double* sxigp,
        double* mxigp, const std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi, double& alpha);

    /*!
     \brief node-wise gap calculation for 3D problems

     */
    void nwGap3D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, double* gpn);

   private:
    /*!
     \brief Interpolate for 2D problems

     */
    void Interpolate2D(MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles);

    /*!
     \brief Interpolate for 3D problems

     */
    bool Interpolate3D(MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles);

    /*!
     \brief lin 2D projection

     */
    void DerivXiGP2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, double& sxigp,
        double& mxigp, const CORE::GEN::pairedvector<int, double>& derivsxi,
        CORE::GEN::pairedvector<int, double>& derivmxi, int& linsize);

    /*!
     \brief node-wise D/M calculation

     */
    void nwDM2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::GEN::pairedvector<int, double>& dmxi);

    /*!
     \brief node-wise D/M calculation for 3D problems

     */
    void nwDM3D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        std::vector<CORE::GEN::pairedvector<int, double>>& dmxi);

    /*!
     \brief node-wise gap calculation

     */
    void nwGap2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::GEN::pairedvector<int, double>& dmxi, double* gpn);

    /*!
     \brief node-wise master temperature calculation for 3D problems

     */
    void nwMasterTemp(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi);

    /*!
     \brief node-wise slip calculation

     */
    void nwSlip2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseMatrix& scoord, CORE::LINALG::SerialDenseMatrix& mcoord,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> scoordold,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> mcoordold, int& snodes, int& linsize,
        CORE::GEN::pairedvector<int, double>& dmxi);

    /*!
     \brief node-wise wear calculation (internal state var.)

     */
    void nwWear2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& mderiv,
        CORE::LINALG::SerialDenseMatrix& scoord, CORE::LINALG::SerialDenseMatrix& mcoord,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> scoordold,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> mcoordold,
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> lagmult, int& snodes, int& linsize,
        double& jumpval, double& area, double* gpn, CORE::GEN::pairedvector<int, double>& dmxi,
        CORE::GEN::pairedvector<int, double>& dslipmatrix,
        CORE::GEN::pairedvector<int, double>& dwear);

    /*!
     \brief node-wise wear calculation (primary variable)

     */
    void nwTE2D(CONTACT::CoNode& mynode, double& area, double& jumpval,
        CORE::GEN::pairedvector<int, double>& dslipmatrix);

    Teuchos::ParameterList& iparams_;  //< containing contact input parameters
    int dim_;                          //< problem dimension
    bool pwslip_;                      //< point-wise evaluated slip increment

    // wear inputs from parameter list
    INPAR::WEAR::WearLaw wearlaw_;         //< type of wear law
    bool wearimpl_;                        //< flag for implicit wear algorithm
    INPAR::WEAR::WearSide wearside_;       //< definition of wear surface
    INPAR::WEAR::WearType weartype_;       //< definition of contact wear algorithm
    INPAR::WEAR::WearShape wearshapefcn_;  //< type of wear shape function
    double wearcoeff_;                     //< wear coefficient
    double wearcoeffm_;                    //< wear coefficient master
    bool sswear_;                          //< flag for steady state wear
    double ssslip_;                        //< fixed slip for steady state wear
  };


  /*!
  \brief A class to implement MTInterpolator

  \author farah
  */
  class MTInterpolator
  {
   public:
    MTInterpolator(){};

    // destructor
    virtual ~MTInterpolator() = default;
    //! @name Access methods
    /// Internal implementation class
    static MTInterpolator* Impl(std::vector<MORTAR::MortarElement*> meles);

    /*!
     \brief Interpolate for nts algorithm

     */
    virtual void Interpolate(
        MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles) = 0;
  };


  /*!
  \author farah
  */
  template <CORE::FE::CellType distypeM>
  class MTInterpolatorCalc : public MTInterpolator
  {
   public:
    MTInterpolatorCalc();

    /// Singleton access method
    static MTInterpolatorCalc<distypeM>* Instance(CORE::UTILS::SingletonAction action);

    //! nm_: number of master element nodes
    static const int nm_ = CORE::FE::num_nodes<distypeM>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static const int ndim_ = CORE::DRT::UTILS::DisTypeToDim<distypeM>::dim + 1;

    /*!
     \brief Interpolate for nts problems

     */
    void Interpolate(MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles) override;

   private:
    /*!
     \brief Interpolate for 2D problems

     */
    virtual void Interpolate2D(
        MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles);

    /*!
     \brief Interpolate for 3D problems

     */
    virtual void Interpolate3D(
        MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles);
  };

}  // namespace NTS

#endif  // CONTACT_INTERPOLATOR_H
