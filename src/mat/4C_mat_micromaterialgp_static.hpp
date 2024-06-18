/*----------------------------------------------------------------------*/
/*! \file
\brief class for handling of micro-macro transitions

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MICROMATERIALGP_STATIC_HPP
#define FOUR_C_MAT_MICROMATERIALGP_STATIC_HPP



#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations

namespace MultiScale
{
  class MicroStatic;
}

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Mat
{
  /// one Gauss point of the micro material
  class MicroMaterialGP
  {
   public:
    /// construct an instance of MicroMaterial for a given Gauss point
    /// and microscale discretization
    MicroMaterialGP(const int gp, const int ele_ID, const bool eleowner, const int microdisnum,
        const double V0);

    /// destructor
    ~MicroMaterialGP();

    /// Read restart
    void read_restart();

    /// New result file
    void new_result_file(bool eleowner, std::string& newfilename);

    //! create path of new result file
    std::string new_result_file_path(const std::string& newprefix);

    /// Perform microscale simulation
    void perform_micro_simulation(Core::LinAlg::Matrix<3, 3>* defgrd,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat);

    void update();

    /// Calculate stresses and strains on the micro-scale
    void prepare_output();

    void Output();

    /// Create and initialize "empty" EAS history map
    void eas_init();

    /// Reset global time and step number (needed for multi-scale inverse
    /// analyses with multiple runs)
    void ResetTimeAndStep();

    /// get ele id
    int eleID() { return ele_id_; }

    /// get density
    double Density() const { return density_; }


   private:
    /// corresponding macroscale Gauss point
    const int gp_;

    /// corresponding macroscale element
    const int ele_id_;

    /// corresponding microstructure discretization number
    const int microdisnum_;

    /// microstructure "time integration" classes (one for each micro-discretization)
    static std::map<int, Teuchos::RCP<MultiScale::MicroStatic>> microstaticmap_;

    static std::map<int, int> microstaticcounter_;

    /// microstructure discretization writer
    Teuchos::RCP<Core::IO::DiscretizationWriter> micro_output_;

    /// homogenized density
    double density_;

    /// my vector of old displacements
    Teuchos::RCP<Epetra_Vector> dis_;

    /// my vector of new displacements
    Teuchos::RCP<Epetra_Vector> disn_;

    // my EAS history data -> note that microstructure is not parallel
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> lastalpha_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldalpha_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> oldfeas_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> old_kaainv_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>>> old_kda_;

    /// my stresses and strains
    Teuchos::RCP<std::vector<char>> stress_;
    Teuchos::RCP<std::vector<char>> strain_;
    Teuchos::RCP<std::vector<char>> plstrain_;

    /// old absolute time
    double time_;

    /// current absolute time
    double timen_;

    /// old step
    int step_;

    /// current step
    int stepn_;

    /// timestep size
    double dt_;

    /// restart name
    std::string restartname_;

    /// flag for modified Newton on macroscale
    bool mod_newton_;

    /// flag for build of stiffness matrix
    bool build_stiff_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
