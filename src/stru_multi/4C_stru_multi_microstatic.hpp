/*---------------------------------------------------------------------*/
/*! \file

\brief Quasi-static control for microstructural analysis


\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_STRU_MULTI_MICROSTATIC_HPP
#define FOUR_C_STRU_MULTI_MICROSTATIC_HPP



#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace MultiScale
{
  /*!
  \brief Stop nested parallelism support by sending a message to the
  supporting procs

  */
  void stop_np_multiscale();

  /*!
  \brief Quasi-static control for microstructural analysis
  in case of multi-scale problems

  Note that implementation currently only holds for imr-like generalized
  alpha time integration. Corresponding functions (e.g. UpdateNewTimeStep,
  but also calls to SurfaceStressManager!) need to be adapted accordingly
  if usage of other time integration schemes should be enabled.

  */

  class MicroStatic
  {
   public:
    /*!
    \brief Standard Constructor

    */
    MicroStatic(const int microdisnum, const double V0);

    /*!
    \brief Destructor

    */
    virtual ~MicroStatic() = default;

    /*!
    \brief Read restart

    */
    void read_restart(int step, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> lastalpha,
        std::string name);

    /*!
    \brief Return time from parameter list

    */
    double TimeOld() { return time_; }

    /*!
    \brief Predictor step

    */
    void Predictor(Core::LinAlg::Matrix<3, 3>* defgrd);

    /*!
    \brief Predictor step

    */
    void PredictConstDis(Core::LinAlg::Matrix<3, 3>* defgrd);

    /*!
   \brief Predictor step

   */
    void PredictTangDis(Core::LinAlg::Matrix<3, 3>* defgrd);

    /*!
    \brief Full Newton iteration

    */
    void FullNewton();

    /*!
    \brief Calculate stresses and strains

    */
    void prepare_output();

    /*!
    \brief Write output and (possibly) restart

    */
    void Output(Teuchos::RCP<Core::IO::DiscretizationWriter> output, const double time,
        const int istep, const double dt);

    /*!
    \brief Determine toggle vector identifying prescribed boundary dofs

    */
    void DetermineToggle();

    /*!
    \brief Evaluate microscale boundary displacement according to
    associated macroscale deformation gradient

    */
    void EvaluateMicroBC(Core::LinAlg::Matrix<3, 3>* defgrd, Teuchos::RCP<Epetra_Vector> disp);

    /*!
    \brief Set old state given from micromaterialgp

    */
    void set_state(Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<std::vector<char>> stress, Teuchos::RCP<std::vector<char>> strain,
        Teuchos::RCP<std::vector<char>> plstrain,
        Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> lastalpha,
        Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldalpha,
        Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldfeas,
        Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldKaainv,
        Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldKda);

    /*!
    \brief Set time and step

    */
    void set_time(
        const double time, const double timen, const double dt, const int step, const int stepn);

    /*!
    \brief Clear all displacement states

    */
    void ClearState();

    /*!
    \brief Set up everything for homogenization
    (e.g. calculation of matrix D containing reference boundary coordinates)

    */
    void SetUpHomogenization();

    /*!
    \brief Perform homogenization, i.e. calculate second Piola-Kirchhoff
    stresses and constitutive tensor by averaging over RVE

    */
    void static_homogenization(Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat,
        Core::LinAlg::Matrix<3, 3>* defgrd, const bool mod_newton, bool& build_stiff);

    /*!
    \brief Convert constitutive tensor relating first Piola-Kirchhoff
    stresses and deformation gradient to tensor relating second
    Piola-Kirchhoff stresses and Green-Lagrange strains

    For details cf.

    Marsden and Hughes, Mathematical Foundations of Elasticity,
    Dover, pg. 215
    */
    void ConvertMat(const Epetra_MultiVector& cmatpf, const Core::LinAlg::Matrix<3, 3>& F_inv,
        const Core::LinAlg::Matrix<6, 1>& S, Core::LinAlg::Matrix<6, 6>& cmat);


    /*!
    \brief Check for Newton convergence
    */

    bool Converged();

    /*!
    \brief Calculate reference norms for relative convergence checks

    */
    void CalcRefNorms();

    /*!
    \brief Output of Newton details

    Note that this is currently disabled for the sake of clearness
    */
    void PrintNewton(bool print_unconv, Teuchos::Time timer);

    /*!
    \brief Output of predictor details

    Note that this is currently disabled for the sake of clearness
    */
    void print_predictor();

    /*!
    \brief Set EAS internal data if necessary

    */
    void SetEASData();

    double Density() const { return density_; };

   protected:
    // don't want = operator and cctor
    MicroStatic operator=(const MicroStatic& old);
    MicroStatic(const MicroStatic& old);

    Teuchos::RCP<Core::FE::Discretization> discret_;
    Teuchos::RCP<Core::LinAlg::Solver> solver_;
    int myrank_;
    int maxentriesperrow_;

    double dt_;
    double time_;
    double timen_;

    Inpar::STR::PredEnum pred_;  //!< predictor

    bool isadapttol_;
    double adaptolbetter_;

    int maxiter_;
    int numiter_;
    int numstep_;
    int step_;
    int stepn_;

    bool iodisp_;
    int resevrydisp_;
    Inpar::STR::StressType iostress_;
    int resevrystrs_;
    Inpar::STR::StrainType iostrain_;
    Inpar::STR::StrainType ioplstrain_;
    bool iosurfactant_;
    int restart_;
    int restartevry_;
    int printscreen_;

    Inpar::STR::VectorNorm iternorm_;
    double tolfres_;
    double toldisi_;


    enum Inpar::STR::BinaryOp combdisifres_;  //!< binary operator to
                                              // combine displacement and forces
    enum Inpar::STR::ConvNorm normtypedisi_;  //!< convergence check for residual displacements
    enum Inpar::STR::ConvNorm normtypefres_;  //!< convergence check for residual forces
    double normcharforce_;
    double normfres_;
    double normchardis_;
    double normdisi_;

    Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_dirich_;

    Teuchos::RCP<Epetra_Vector> dirichtoggle_;
    Teuchos::RCP<Epetra_Vector> invtoggle_;
    Teuchos::RCP<Epetra_Vector> zeros_;
    Teuchos::RCP<Epetra_Vector>
        dis_;  //!< displacements at t_{n} (needed for convergence check only)
    Teuchos::RCP<Epetra_Vector> disn_;  //!< displacements at t_{n+1}
    Teuchos::RCP<Epetra_Vector> disi_;
    Teuchos::RCP<Epetra_Vector> fintn_;
    Teuchos::RCP<Epetra_Vector> fresn_;
    Teuchos::RCP<Epetra_Vector> freactn_;

    Teuchos::RCP<std::vector<char>> stress_;
    Teuchos::RCP<std::vector<char>> strain_;
    Teuchos::RCP<std::vector<char>> plstrain_;

    // EAS history data
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> lastalpha_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldalpha_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldfeas_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldKaainv_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldKda_;

    Teuchos::RCP<Epetra_MultiVector> D_;    //!< D Matrix following Miehe et al., 2002
    Teuchos::RCP<Epetra_MultiVector> rhs_;  //!< exported transpose of D (pdof -> dofrowmap)

    int microdisnum_;  //!< number of RVE

    double V0_;       //!< initial volume of RVE
    double density_;  //!< initial density of RVE

    int ndof_;                        //!< number of dofs overall
    int np_;                          //!< number of boundary dofs
    Teuchos::RCP<Epetra_Vector> Xp_;  //!< vector containing material
                                      //!< coordinates of boundary nodes
    Teuchos::RCP<Epetra_Map> pdof_;   //!< prescribed dofs
    Teuchos::RCP<Epetra_Map> fdof_;   //!< free dofs
    Teuchos::RCP<Epetra_Import> importp_;
    Teuchos::RCP<Epetra_Import> importf_;
  };


  class MicroStaticParObjectType : public Core::Communication::ParObjectType
  {
   public:
    [[nodiscard]] std::string Name() const override { return "MicroStaticParObjectType"; }

    static MicroStaticParObjectType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MicroStaticParObjectType instance_;
  };

  class MicroStaticParObject : public Core::Communication::ParObject
  {
   public:
    [[nodiscard]] inline int UniqueParObjectId() const override
    {
      return MultiScale::MicroStaticParObjectType::Instance().UniqueParObjectId();
    };

    void pack(Core::Communication::PackBuffer& data) const override;

    void unpack(const std::vector<char>& data) override;

    struct MicroStaticData
    {
      int gp_{};
      int microdisnum_{};
      int eleowner_{};
      double V0_{};
      Core::LinAlg::SerialDenseMatrix defgrd_;
      Core::LinAlg::SerialDenseMatrix stress_;
      Core::LinAlg::SerialDenseMatrix cmat_;
    };

    [[nodiscard]] inline const MicroStaticData* get_micro_static_data_ptr() const
    {
      return std::addressof(microstatic_data_);
    };

    inline void SetMicroStaticData(MicroStaticData& micro_data) { microstatic_data_ = micro_data; };

   private:
    MicroStaticData microstatic_data_{};
  };
}  // namespace MultiScale
FOUR_C_NAMESPACE_CLOSE

#endif
