/*----------------------------------------------------------------------------*/
/*! \file
\brief unite all necessary methods to generate the data for external plots in
MATLAB, PGFPlot or other tools.

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_PLOT_HPP
#define FOUR_C_CONTACT_AUG_PLOT_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

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

namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::Nodes
{
  class Node;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class Contact;
  }  // namespace MODELEVALUATOR
}  // namespace STR
namespace NOX
{
  namespace Nln
  {
    namespace CONSTRAINT
    {
      class Group;
    }  // namespace CONSTRAINT
    namespace MeritFunction
    {
      enum MeritFctName : int;
    }  // namespace MeritFunction
  }    // namespace Nln
}  // namespace NOX
namespace Inpar
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
}  // namespace Inpar
namespace CONTACT
{
  class AbstractStrategy;
  namespace Aug
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
      static bool activated(const Teuchos::ParameterList& plot_params);

      Plot();

      void init(const Teuchos::ParameterList& plot_params, const CONTACT::AbstractStrategy* strat);

      void setup();

      void execute(const ::NOX::Solver::Generic& solver);

      void lin_space(
          const double a, const double b, const unsigned n, std::vector<double>& res) const;

      const CONTACT::Aug::Strategy& strategy() const;

      void get_support_points(enum Inpar::CONTACT::PlotSupportType stype,
          Core::LinAlg::SerialDenseMatrix& support_mat_x);

      void compute_angle_position();

      void compute_distance_position();

      void plot_scalar(const NOX::Nln::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::Nln::CONSTRAINT::Group& plot_grp);

      void plot_line(const NOX::Nln::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::Nln::CONSTRAINT::Group& plot_grp);

      void plot_surface(const NOX::Nln::CONSTRAINT::Group& ref_grp, const Epetra_Vector& dir,
          NOX::Nln::CONSTRAINT::Group& plot_grp);

      void plot_vector_field2_d(const NOX::Nln::CONSTRAINT::Group& ref_grp,
          const Epetra_Vector& dir, NOX::Nln::CONSTRAINT::Group& plot_grp);

      void write_surface_data_to_file() const;

      void write_line_data_to_file() const;

      void write_vector_field_to_file() const;

      void add_file_name_to_path();

      enum NOX::Nln::MeritFunction::MeritFctName convert_plot_func_name2_merit_func_name(
          const enum Inpar::CONTACT::PlotFuncName pfunc_name) const;

      enum CONTACT::Aug::WGapGradientType convert_plot_func_name2_w_gap_gradient_type(
          const enum Inpar::CONTACT::PlotFuncName pfunc_name) const;

      const NOX::Nln::CONSTRAINT::Group* get_reference_group(
          const ::NOX::Solver::Generic& solver) const;

      double get_value(const enum Inpar::CONTACT::PlotFuncName functype,
          NOX::Nln::CONSTRAINT::Group& plot_grp, const double* curr_xy = nullptr,
          const Epetra_Vector* curr_dir = nullptr) const;

      double get_nodal_error_at_position(
          const double* pos, const std::vector<std::pair<int, double>>& nodal_error) const;

      void get_vector_values(const enum Inpar::CONTACT::PlotFuncName functype,
          NOX::Nln::CONSTRAINT::Group& plot_grp, const std::vector<const Epetra_Vector*>& dirs,
          std::vector<double>& vec_vals) const;

      void get_w_gap_direction_gradients(const enum CONTACT::Aug::WGapGradientType wgap_type,
          const std::vector<const Epetra_Vector*>& dirs, std::vector<double>& grad_vals) const;

      void get_energy_direction_gradients(
          const std::vector<const Epetra_Vector*>& dirs, std::vector<double>& grad_vals) const;

      int map_sl_node_gi_d2_n_dof_gid(int node_gid) const;

      double characteristic_interface_element_length(const enum CONTACT::Aug::SideType stype) const;

      void modify_step_length(const Inpar::CONTACT::PlotSupportType stype, const double alpha,
          const Epetra_Vector& full_x_dir, Epetra_Vector& mod_step) const;

      void read_ref_points(const Teuchos::ParameterList& plot_params);

      void read_ref_point(const Teuchos::ParameterList& plot_params, const std::string& param_name,
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

        void split_into_surface_directions(const Epetra_Vector& dir,
            Teuchos::RCP<Epetra_Vector>& x_dir_ptr, Teuchos::RCP<Epetra_Vector>& y_dir_ptr) const;

        Teuchos::RCP<const Epetra_Vector> Get(const ::NOX::Solver::Generic& solver) const;

        enum Inpar::CONTACT::PlotDirection type_;
        enum Inpar::CONTACT::PlotDirectionSplit split_;
        Teuchos::RCP<Epetra_Vector> from_file_;

       private:
        /// call-back
        const Plot& plot_;

        std::string get_full_file_path(
            const std::string& input_file, const std::string& dir_file) const;

        void split_into_slave_master_body(const Epetra_Vector& dir,
            Teuchos::RCP<Epetra_Vector>& x_dir_ptr, Teuchos::RCP<Epetra_Vector>& y_dir_ptr) const;

        Teuchos::RCP<Epetra_Vector> read_sparse_vector_from_matlab(
            const std::string& dir_file) const;

        bool extend_file_name(std::string& file_name, const std::string& file_path) const;

        Teuchos::RCP<Epetra_Map> find_connected_dofs(
            const Core::Nodes::Node* node, const Discret::Discretization& discret) const;
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

      const CONTACT::Aug::Strategy* strat_ = nullptr;

      /// full discretization
      const Discret::Discretization* discret_ = nullptr;

      STR::MODELEVALUATOR::Contact* model_ = nullptr;

      const int* curr_step_np_ = nullptr;

      std::ios_base::openmode file_open_mode_;

      Inpar::CONTACT::PlotMode mode_;

      Inpar::CONTACT::PlotFuncName func_type_;

      Inpar::CONTACT::PlotType type_;

      Inpar::CONTACT::PlotReferenceType reference_type_;

      Inpar::CONTACT::PlotFileFormat format_;

      Inpar::CONTACT::PlotSupportType x_type_;
      Inpar::CONTACT::PlotSupportType y_type_;

      int wgap_node_gid_ = -1;

      std::map<double, int> position_node_id_map_;

      std::vector<Core::LinAlg::Matrix<3, 1>> ref_points_;

      Core::LinAlg::SerialDenseMatrix x_;
      Core::LinAlg::SerialDenseMatrix y_;
      std::vector<Core::LinAlg::SerialDenseMatrix> z_;
    };  // class Plot

    template <typename T>
    void WriteMatrixToFile(std::ofstream& outputfile, const T& mat, const unsigned precision);

    template <typename T>
    void WriteColumnDataToFile(std::ofstream& outputfile, const std::vector<const T*>& columndata,
        const unsigned precision);
  }  // namespace Aug
}  // namespace CONTACT

namespace NOX
{
  namespace Nln
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
            Plot(const Teuchos::RCP<FourC::CONTACT::Aug::Plot>& plot_ptr);

            void runPreIterate(const ::NOX::Solver::Generic& solver) override;
            void runPostIterate(const ::NOX::Solver::Generic& solver) override;

           private:
            Teuchos::RCP<FourC::CONTACT::Aug::Plot> plot_ptr_;
            FourC::CONTACT::Aug::Plot& plot_;
          };  // class Plot
        }     // namespace CONTACT
      }       // namespace PrePostOp
    }         // namespace Solver
  }           // namespace Nln
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
