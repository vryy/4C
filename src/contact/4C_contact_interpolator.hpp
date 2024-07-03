/*----------------------------------------------------------------------*/
/*! \file

\level 2


*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 09/14 |
 *---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_INTERPOLATOR_HPP
#define FOUR_C_CONTACT_INTERPOLATOR_HPP

/*---------------------------------------------------------------------*
 | headers                                                 farah 09/14 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_utils_pairedvector.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 09/14 |
 *---------------------------------------------------------------------*/
namespace Mortar
{
  class Element;
  class Node;
}  // namespace Mortar

namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace CONTACT
{
  class Node;
}

namespace NTS
{
  class Interpolator
  {
   public:
    /*!
     \brief Constructor

     */
    Interpolator(Teuchos::ParameterList& params, const int& dim);

    /*!
     \brief Destructor

     */
    virtual ~Interpolator() = default;

    /*!
     \brief Interpolate for nts algorithm

     */
    bool Interpolate(Mortar::Node& snode, std::vector<Mortar::Element*> meles);

    /*!
     \brief Interpolate temperature of master side at a slave node
     for 3D problems

     */
    void interpolate_master_temp3_d(Mortar::Element& sele, std::vector<Mortar::Element*> meles);

    /*!
     \brief lin 3D projection

     */
    void DerivXiGP3D(Mortar::Element& sele, Mortar::Element& mele, double* sxigp, double* mxigp,
        const std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
        std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi, double& alpha);

    /*!
     \brief node-wise gap calculation for 3D problems

     */
    void nwGap3D(CONTACT::Node& mynode, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
        std::vector<Core::Gen::Pairedvector<int, double>>& dmxi, double* gpn);

   private:
    /*!
     \brief Interpolate for 2D problems

     */
    void interpolate2_d(Mortar::Node& snode, std::vector<Mortar::Element*> meles);

    /*!
     \brief Interpolate for 3D problems

     */
    bool interpolate3_d(Mortar::Node& snode, std::vector<Mortar::Element*> meles);

    /*!
     \brief lin 2D projection

     */
    void deriv_xi_g_p2_d(Mortar::Element& sele, Mortar::Element& mele, double& sxigp, double& mxigp,
        const Core::Gen::Pairedvector<int, double>& derivsxi,
        Core::Gen::Pairedvector<int, double>& derivmxi, int& linsize);

    /*!
     \brief node-wise D/M calculation

     */
    void nw_d_m2_d(CONTACT::Node& mynode, Mortar::Element& sele, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
        Core::Gen::Pairedvector<int, double>& dmxi);

    /*!
     \brief node-wise D/M calculation for 3D problems

     */
    void nw_d_m3_d(CONTACT::Node& mynode, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
        std::vector<Core::Gen::Pairedvector<int, double>>& dmxi);

    /*!
     \brief node-wise gap calculation

     */
    void nw_gap2_d(CONTACT::Node& mynode, Mortar::Element& sele, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
        Core::Gen::Pairedvector<int, double>& dmxi, double* gpn);

    /*!
     \brief node-wise master temperature calculation for 3D problems

     */
    void nw_master_temp(CONTACT::Node& mynode, Mortar::Element& mele,
        const Core::LinAlg::SerialDenseVector& mval, const Core::LinAlg::SerialDenseMatrix& mderiv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& dmxi);

    /*!
     \brief node-wise slip calculation

     */
    void nw_slip2_d(CONTACT::Node& mynode, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
        Core::LinAlg::SerialDenseMatrix& scoord, Core::LinAlg::SerialDenseMatrix& mcoord,
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> scoordold,
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> mcoordold, int& snodes, int& linsize,
        Core::Gen::Pairedvector<int, double>& dmxi);

    /*!
     \brief node-wise wear calculation (internal state var.)

     */
    void nw_wear2_d(CONTACT::Node& mynode, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& mderiv,
        Core::LinAlg::SerialDenseMatrix& scoord, Core::LinAlg::SerialDenseMatrix& mcoord,
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> scoordold,
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> mcoordold,
        Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> lagmult, int& snodes, int& linsize,
        double& jumpval, double& area, double* gpn, Core::Gen::Pairedvector<int, double>& dmxi,
        Core::Gen::Pairedvector<int, double>& dslipmatrix,
        Core::Gen::Pairedvector<int, double>& dwear);

    /*!
     \brief node-wise wear calculation (primary variable)

     */
    void nw_t_e2_d(CONTACT::Node& mynode, double& area, double& jumpval,
        Core::Gen::Pairedvector<int, double>& dslipmatrix);

    Teuchos::ParameterList& iparams_;  //< containing contact input parameters
    int dim_;                          //< problem dimension
    bool pwslip_;                      //< point-wise evaluated slip increment

    // wear inputs from parameter list
    Inpar::Wear::WearLaw wearlaw_;         //< type of wear law
    bool wearimpl_;                        //< flag for implicit wear algorithm
    Inpar::Wear::WearSide wearside_;       //< definition of wear surface
    Inpar::Wear::WearType weartype_;       //< definition of contact wear algorithm
    Inpar::Wear::WearShape wearshapefcn_;  //< type of wear shape function
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
    static MTInterpolator* Impl(std::vector<Mortar::Element*> meles);

    /*!
     \brief Interpolate for nts algorithm

     */
    virtual void Interpolate(Mortar::Node& snode, std::vector<Mortar::Element*> meles) = 0;
  };


  /*!
  \author farah
  */
  template <Core::FE::CellType distype_m>
  class MTInterpolatorCalc : public MTInterpolator
  {
   public:
    MTInterpolatorCalc();

    /// Singleton access method
    static MTInterpolatorCalc<distype_m>* Instance(Core::UTILS::SingletonAction action);

    //! nm_: number of master element nodes
    static constexpr int nm_ = Core::FE::num_nodes<distype_m>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = Core::FE::dim<distype_m> + 1;

    /*!
     \brief Interpolate for nts problems

     */
    void Interpolate(Mortar::Node& snode, std::vector<Mortar::Element*> meles) override;

   private:
    /*!
     \brief Interpolate for 2D problems

     */
    virtual void interpolate2_d(Mortar::Node& snode, std::vector<Mortar::Element*> meles);

    /*!
     \brief Interpolate for 3D problems

     */
    virtual void interpolate3_d(Mortar::Node& snode, std::vector<Mortar::Element*> meles);
  };

}  // namespace NTS

FOUR_C_NAMESPACE_CLOSE

#endif
