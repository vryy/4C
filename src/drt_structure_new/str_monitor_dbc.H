/*----------------------------------------------------------------------------*/
/*!
\file str_monitor_dbc.H

\brief Monitor tagged Dirichlet boundary conditions

\level 3

\maintainer Michael Hiermeier
\date Jan 2, 2018

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_MONITOR_DBC_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MONITOR_DBC_H_

#include "../drt_lib/drt_condition.H"

namespace IO
{
  class DiscretizationWriter;
}
namespace DRT
{
  class DiscretizationInterface;
}
namespace STR
{
  class Dbc;
  namespace TIMINT
  {
    class BaseDataGlobalState;
  }

  /** \brief Monitor Dirichlet boundary conditions
   *
   *  This class can be used to monitor e.g. the reaction forces and the area
   *  change of a tagged Dirichlet condition during a simulation. To tag a Dirichlet
   *  condition just add the corresponding TAG, e.g. \"monitor_reaction\"
   *
   *  E 1 - NUMDOF 3 ONOFF 1 0 0 VAL 0.0 0.0 0.0 FUNCT 0 0 0 TAG monitor_reaction
   *
   *  If the TAG can be found for any Dirichlet condition the reaction force as
   *  well as the reference and current area will be stored in a text file
   *  located at
   *
   *  <OUTPUT_PATH>/<OUTPUT_FILE_NAME>_monitor_dbc/<ID>_monitor_dbc.data
   *
   *  If no tag is found nothing is happening.
   *
   *  \author hiermeier \date 01/18 */
  class MonitorDbc
  {
    const static unsigned DIM = 3;

    /// constants for the FILE output
    const static unsigned OF_WIDTH = 24;

    /// constants for the SCREEN output
    const static unsigned OS_WIDTH = 14;

   public:
    MonitorDbc() = default;

    ~MonitorDbc() = default;

    /// initialize class members
    void Init(DRT::DiscretizationInterface& discret, STR::TIMINT::BaseDataGlobalState& gstate,
        STR::Dbc& dbc);

    /// setup new class members
    void Setup();

    /// monitor the tensile test results and write them to a text file
    void Execute(IO::DiscretizationWriter& writer);

   private:
    int GetUniqueId(int tagged_id, DRT::Condition::GeometryType gtype) const;

    void CreateReactionForceCondition(
        const DRT::Condition& tagged_cond, DRT::DiscretizationInterface& discret) const;

    void GetTaggedCondition(std::vector<const DRT::Condition*>& tagged_conds,
        const std::string& cond_name, const std::string& tag_name,
        const DRT::DiscretizationInterface& discret) const;

    void CreateReactionMaps(const DRT::DiscretizationInterface& discret,
        const DRT::Condition& rcond, Teuchos::RCP<Epetra_Map>* react_maps) const;

    void ReadResultsPriorRestartStepAndWriteToFile(
        const std::vector<std::string>& full_restart_filepaths, int restart_step) const;

    void GetArea(double area_ref[], const DRT::Condition* rcond) const;

    double GetReactionForce(
        LINALG::Matrix<3, 1>& rforce_xyz, const Teuchos::RCP<Epetra_Map>* react_maps) const;

    std::vector<std::string> CreateFilePaths(
        const std::vector<Teuchos::RCP<DRT::Condition>>& rconds, const std::string& full_dirpath,
        const std::string& filename_only_prefix, const std::string& file_type) const;

    void ClearFilesAndWriteHeader(const std::vector<Teuchos::RCP<DRT::Condition>>& rconds,
        std::vector<std::string>& full_filepaths, bool write_condition_header);

    void WriteConditionHeader(
        std::ostream& os, const int col_width, const DRT::Condition* cond = NULL) const;

    void WriteColumnHeader(std::ostream& os, const int col_width) const;

    void WriteResultsToFile(const std::string& full_filepath, const LINALG::Matrix<DIM, 1>& rforce,
        const double& area_ref, const double& area_curr) const;

    void WriteResultsToScreen(const Teuchos::RCP<DRT::Condition>& rcond_ptr,
        const LINALG::Matrix<DIM, 1>& rforce, const double& area_ref,
        const double& area_curr) const;

    void WriteResults(std::ostream& os, const int col_width, const int precision,
        const unsigned step, const double time, const LINALG::Matrix<DIM, 1>& rforce,
        const double& area_ref, const double& area_cur) const;

    inline const Epetra_Comm& Comm() const;

    inline void ThrowIfNotInit() const
    {
      if (not isinit_) dserror("Call Init() first!");
    }

    inline void ThrowIfNotSetup() const
    {
      if (not issetup_) dserror("Call Setup() first!");
    }

   private:
    DRT::DiscretizationInterface* discret_ptr_ = NULL;
    STR::TIMINT::BaseDataGlobalState* gstate_ptr_ = NULL;
    STR::Dbc* dbc_ptr_ = NULL;

    std::vector<std::string> full_filepaths_ = std::vector<std::string>();

    /// extract the dofs of the reaction forces which shall be monitored
    std::map<int, std::vector<Teuchos::RCP<Epetra_Map>>> react_maps_;
    DRT::Condition::GeometryType geo_type_ = DRT::Condition::NoGeom;

    unsigned of_precision_ = -1;
    unsigned os_precision_ = -1;

    bool isempty_ = true;
    bool isinit_ = false;
    bool issetup_ = false;
  };

}  // namespace STR

#endif /* SRC_DRT_STRUCTURE_NEW_STR_MONITOR_DBC_H_ */
