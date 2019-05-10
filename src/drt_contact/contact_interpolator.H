/*!----------------------------------------------------------------------
\file contact_interpolator.H

\maintainer Alexander Popp

\level 2

</pre>

*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 09/14 |
 *---------------------------------------------------------------------*/
#ifndef CONTACT_INTERPOLATOR_H_
#define CONTACT_INTERPOLATOR_H_

/*---------------------------------------------------------------------*
 | headers                                                 farah 09/14 |
 *---------------------------------------------------------------------*/
#include "../headers/pairedvector.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_inpar/inpar_wear.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 09/14 |
 *---------------------------------------------------------------------*/
namespace MORTAR
{
  class MortarElement;
  class MortarNode;
}  // namespace MORTAR

namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG

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
    virtual ~CoInterpolator(){};

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
        double* mxigp, const std::vector<GEN::pairedvector<int, double>>& derivsxi,
        std::vector<GEN::pairedvector<int, double>>& derivmxi, double& alpha);

    /*!
     \brief node-wise gap calculation for 3D problems

     */
    void nwGap3D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
        std::vector<GEN::pairedvector<int, double>>& dmxi, double* gpn);

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
        double& mxigp, const GEN::pairedvector<int, double>& derivsxi,
        GEN::pairedvector<int, double>& derivmxi, int& linsize);

    /*!
     \brief node-wise D/M calculation

     */
    void nwDM2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
        GEN::pairedvector<int, double>& dmxi);

    /*!
     \brief node-wise D/M calculation for 3D problems

     */
    void nwDM3D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
        std::vector<GEN::pairedvector<int, double>>& dmxi);

    /*!
     \brief node-wise gap calculation

     */
    void nwGap2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
        GEN::pairedvector<int, double>& dmxi, double* gpn);

    /*!
     \brief node-wise master temperature calculation for 3D problems

     */
    void nwMasterTemp(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        const LINALG::SerialDenseVector& mval, const LINALG::SerialDenseMatrix& mderiv,
        const std::vector<GEN::pairedvector<int, double>>& dmxi);

    /*!
     \brief node-wise slip calculation

     */
    void nwSlip2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
        LINALG::SerialDenseMatrix& scoord, LINALG::SerialDenseMatrix& mcoord,
        Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
        Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold, int& snodes, int& linsize,
        GEN::pairedvector<int, double>& dmxi);

    /*!
     \brief node-wise wear calculation (internal state var.)

     */
    void nwWear2D(CONTACT::CoNode& mynode, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& mderiv,
        LINALG::SerialDenseMatrix& scoord, LINALG::SerialDenseMatrix& mcoord,
        Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
        Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
        Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult, int& snodes, int& linsize, double& jumpval,
        double& area, double* gpn, GEN::pairedvector<int, double>& dmxi,
        GEN::pairedvector<int, double>& dslipmatrix, GEN::pairedvector<int, double>& dwear);

    /*!
     \brief node-wise wear calculation (primary variable)

     */
    void nwTE2D(CONTACT::CoNode& mynode, double& area, double& jumpval,
        GEN::pairedvector<int, double>& dslipmatrix);

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

  \author farah (farah@lnm.mw.tum.de)
  */
  class MTInterpolator : DRT::SingletonDestruction
  {
   public:
    MTInterpolator(){};

    // destructor
    virtual ~MTInterpolator() {}

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
  \author farah (farah@lnm.mw.tum.de)
  */
  template <DRT::Element::DiscretizationType distypeM>
  class MTInterpolatorCalc : public MTInterpolator
  {
   public:
    MTInterpolatorCalc();

    // destructor
    virtual ~MTInterpolatorCalc(){};

    /// called upon destruction
    virtual void Done();

    /// Singleton access method
    static MTInterpolatorCalc<distypeM>* Instance(bool create);

    //! nm_: number of master element nodes
    static const int nm_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static const int ndim_ = DRT::UTILS::DisTypeToDim<distypeM>::dim + 1;

    /*!
     \brief Interpolate for nts problems

     */
    virtual void Interpolate(MORTAR::MortarNode& snode, std::vector<MORTAR::MortarElement*> meles);

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

#endif /* CONTACT_INTERPOLATOR_H_ */
