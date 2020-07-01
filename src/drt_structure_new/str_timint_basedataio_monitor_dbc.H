/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related monitoring reaction forces for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_MONITOR_DBC_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_MONITOR_DBC_H_

#include <string>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

namespace STR
{
  namespace TIMINT
  {
    /** \brief Input data container for monitoring reaction forces for structural (time) integration
     *
     * \author Jonas Eichinger */
    class ParamsMonitorDBC
    {
     public:
      /// constructor
      ParamsMonitorDBC();

      /// destructor
      virtual ~ParamsMonitorDBC(){};

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtk_structure_paramslist);

      /// setup new class variables
      void Setup();


      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int OutputIntervalInSteps() const
      {
        CheckInitSetup();
        return output_interval_steps_;
      };

      /// precision for file output
      int FilePrecision() const
      {
        CheckInitSetup();
        return of_precision_;
      };

      /// precision for screen output
      int ScreenPrecision() const
      {
        CheckInitSetup();
        return os_precision_;
      };

      /// file tpye ending
      std::string const& FileType() const
      {
        CheckInitSetup();
        return file_type_;
      };

      /// whether to write header in csv files
      bool WriteHeader() const
      {
        CheckInitSetup();
        return write_header_;
      }


     private:
      /// get the init indicator status
      const bool& IsInit() const { return isinit_; };

      /// get the setup indicator status
      const bool& IsSetup() const { return issetup_; };

      /// Check if Init() and Setup() have been called, yet.
      void CheckInitSetup() const;


     private:
      /// @name variables for internal use only
      /// @{
      ///
      bool isinit_;

      bool issetup_;
      /// @}

      /// @name variables controlling output
      /// @{

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int output_interval_steps_;

      /// precision for file output
      unsigned of_precision_;

      /// precision for screen output
      unsigned os_precision_;

      /// file type
      std::string file_type_;

      /// write header in csv files
      bool write_header_;

      /// @}
    };

  }  // namespace TIMINT
}  // namespace STR


#endif
