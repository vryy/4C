/*----------------------------------------------------------------------------*/
/*! \file
\brief unite all necessary methods to generate the data for external plots in
MATLAB, PGFPlot or other tools.

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_PLOT_HPP
#define FOUR_C_CONTACT_AUG_PLOT_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_serialdensematrix.hpp"

#include <NOX_Observer.hpp>
#include <Teuchos_RCP.hpp>

class Epetra_Vector;
class Epetra_Map;
namespace Teuchos
{
  class ParameterList;
}
namespace NOX
{
  namespace Solver
  {
    class Generic;
  }  // namespace Solver
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT
namespace STR
{
  namespace MODELEVALUATOR
  {
    class Contact;
  }  // namespace MODELEVALUATOR
}  // namespace STR
namespace NOX
{
  namespace NLN
  {
    namespace CONSTRAINT
    {
      class Group;
    }  // namespace CONSTRAINT
    namespace MeritFunction
    {
      enum MeritFctName : int;
    }  // namespace MeritFunction
  }    // namespace NLN
}  // namespace NOX
namespace INPAR
{
  namespace CONTACT
  {
    enum class PlotMode : int;
    enum class PlotType : int;
    enum class PlotFileFormat : char;
    enum class PlotFuncName : char;
    enum class PlotReferenceType : char;
    enum class PlotSupportType : char;
    enum class PlotDirection : char;
    enum class PlotDirectionSplit : char;
  }  // namespace CONTACT
}  // namespace INPAR
namespace CONTACT
{
  class AbstractStrategy;
  namespace AUG
  {
    enum class WGapGradientType : char;
    enum class SideType : char;

    class Strategy;

    class Plot
    {
     public:
      static void Create(Teuchos::ParameterList& nox_params,
          const Teuchos::ParameterList& plot_params, const CONTACT::AbstractStrategy* strat);

      void Do(const ::NOX::Solver::Generic& solver);

      /// extra wrapper function to enable the plot of the predictor state as well
      void DoPredictor(const ::NOX::Solver::Generic& solver);

     private:
      static bool Activated(const Teuchos::ParameterList& plot_params);

      Plot();

      void Init(const Teuchos::ParameterList& plot_params, const CONTACT::AbstractStrategy* strat);

      void Setup();

      void Execute(const ::NOX::Solver::Generic& solver);

      void LinSpace(
          const double a, const double b, const unsigned n, std::vector<double>& res) const;

      const CONTACT::AUG::Strategy& Strategy() const;

      void GetSupportPoints(enum INPAR::CONTACT::PlotSupportType stype,
          CORE::LINALG::SerialDenseMatrix& support_mat_x);

      void ComputeAnglePosition();

      void ComputeDistancePosition();

      void PlotScalar(const NOX::NLN::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::NLN::CONSTRAINT::Group& plot_grp);

      void PlotLine(const NOX::NLN::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::NLN::CONSTRAINT::Group& plot_grp);

      void PlotSurface(const NOX::NLN::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::NLN::CONSTRAINT::Group& plot_grp);

      void PlotVectorField2D(const NOX::NLN::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::NLN::CONSTRAINT::Group& plot_grp);

      void WriteSurfaceDataToFile() const;

      void WriteLineDataToFile() const;

      void WriteVectorFieldToFile() const;

      void AddFileNameToPath();

      enum NOX::NLN::MeritFunction::MeritFctName ConvertPlotFuncName2MeritFuncName(
          const enum INPAR::CONTACT::PlotFuncName pfunc_name) const;

      enum CONTACT::AUG::WGapGradientType ConvertPlotFuncName2WGapGradientType(
          const enum INPAR::CONTACT::PlotFuncName pfunc_name) const;

      const NOX::NLN::CONSTRAINT::Group* GetReferenceGroup(
          const ::NOX::Solver::Generic& solver) const;

      double GetValue(const enum INPAR::CONTACT::PlotFuncName functype,
          NOX::NLN::CONSTRAINT::Group& plot_grp, const double* curr_xy = nullptr,
          const Epetra_Vector* curr_dir = nullptr) const;

      double GetNodalErrorAtPosition(
          const double* pos, const std::vector<std::pair<int, double>>& nodal_error) const;

      void GetVectorValues(const enum INPAR::CONTACT::PlotFuncName functype,
          NOX::NLN::CONSTRAINT::Group& plot_grp, const std::vector<const Epetra_Vector*>& dirs,
          std::vector<double>& vec_vals) const;

      void GetWGapDirectionGradients(const enum CONTACT::AUG::WGapGradientType wgap_type,
          const std::vector<const Epetra_Vector*>& dirs, std::vector<double>& grad_vals) const;

      void GetEnergyDirectionGradients(
          const std::vector<const Epetra_Vector*>& dirs, std::vector<double>& grad_vals) const;

      int MapSlNodeGID2NDofGID(int node_gid) const;

      double CharacteristicInterfaceElementLength(const enum CONTACT::AUG::SideType stype) const;

      void ModifyStepLength(const INPAR::CONTACT::PlotSupportType stype, const double alpha,
          const Epetra_Vector& full_x_dir, Epetra_Vector& mod_step) const;

      void ReadRefPoints(const Teuchos::ParameterList& plot_params);

      void ReadRefPoint(const Teuchos::ParameterList& plot_params, const std::string& param_name,
          double* coords) const;

     private:
      struct DoPlot
      {
        DoPlot() : step_(0), iter_(0){};

        int step_;
        int iter_;
      };

      class Direction
      {
       public:
        Direction(const Plot& plot);

        void ReadInput(const Teuchos::ParameterList& pp);

        void SplitIntoSurfaceDirections(const Epetra_Vector& dir,
            Teuchos::RCP<Epetra_Vector>& x_dir_ptr, Teuchos::RCP<Epetra_Vector>& y_dir_ptr) const;

        Teuchos::RCP<const Epetra_Vector> Get(const ::NOX::Solver::Generic& solver) const;

        enum INPAR::CONTACT::PlotDirection type_;
        enum INPAR::CONTACT::PlotDirectionSplit split_;
        Teuchos::RCP<Epetra_Vector> from_file_;

       private:
        /// call-back
        const Plot& plot_;

        std::string GetFullFilePath(
            const std::string& input_file, const std::string& dir_file) const;

        void SplitIntoSlaveMasterBody(const Epetra_Vector& dir,
            Teuchos::RCP<Epetra_Vector>& x_dir_ptr, Teuchos::RCP<Epetra_Vector>& y_dir_ptr) const;

        Teuchos::RCP<Epetra_Vector> ReadSparseVectorFromMatlab(const std::string& dir_file) const;

        bool ExtendFileName(std::string& file_name, const std::string& file_path) const;

        Teuchos::RCP<Epetra_Map> FindConnectedDofs(
            const DRT::Node* node, const DRT::Discretization& discret) const;
      };

      struct Options
      {
        Options()
            : resolution_x_(0), resolution_y_(0), min_x_(0.0), max_x_(0.0), min_y_(0.0), max_y_(0.0)
        {
        }

        unsigned output_precision_ = 0;

        unsigned resolution_x_;
        unsigned resolution_y_;

        double min_x_;
        double max_x_;

        double min_y_;
        double max_y_;
      };

      DoPlot do_plot_;

      Direction dir_;

      Options opt_;

      std::string filepath_;

      const CONTACT::AUG::Strategy* strat_ = nullptr;

      /// full discretization
      const DRT::Discretization* discret_ = nullptr;

      STR::MODELEVALUATOR::Contact* model_ = nullptr;

      const int* curr_step_np_ = nullptr;

      std::ios_base::openmode file_open_mode_;

      INPAR::CONTACT::PlotMode mode_;

      INPAR::CONTACT::PlotFuncName func_type_;

      INPAR::CONTACT::PlotType type_;

      INPAR::CONTACT::PlotReferenceType reference_type_;

      INPAR::CONTACT::PlotFileFormat format_;

      INPAR::CONTACT::PlotSupportType x_type_;
      INPAR::CONTACT::PlotSupportType y_type_;

      int wgap_node_gid_ = -1;

      std::map<double, int> position_node_id_map_;

      std::vector<CORE::LINALG::Matrix<3, 1>> ref_points_;

      CORE::LINALG::SerialDenseMatrix x_;
      CORE::LINALG::SerialDenseMatrix y_;
      std::vector<CORE::LINALG::SerialDenseMatrix> z_;
    };  // class Plot

    template <typename T>
    void WriteMatrixToFile(std::ofstream& outputfile, const T& mat, const unsigned precision);

    template <typename T>
    void WriteColumnDataToFile(std::ofstream& outputfile, const std::vector<const T*>& columndata,
        const unsigned precision);
  }  // namespace AUG
}  // namespace CONTACT

namespace NOX
{
  namespace NLN
  {
    namespace Solver
    {
      namespace PrePostOp
      {
        namespace CONTACT
        {
          class Plot : public ::NOX::Observer
          {
           public:
            Plot(const Teuchos::RCP<FourC::CONTACT::AUG::Plot>& plot_ptr);

            void runPreIterate(const ::NOX::Solver::Generic& solver) override;
            void runPostIterate(const ::NOX::Solver::Generic& solver) override;

           private:
            Teuchos::RCP<FourC::CONTACT::AUG::Plot> plot_ptr_;
            FourC::CONTACT::AUG::Plot& plot_;
          };  // class Plot
        }     // namespace CONTACT
      }       // namespace PrePostOp
    }         // namespace Solver
  }           // namespace NLN
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
