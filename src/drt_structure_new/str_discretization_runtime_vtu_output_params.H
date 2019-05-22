/*-----------------------------------------------------------------------------------------------*/
/*!

\brief input parameters related to VTU output at runtime for beams

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef STR_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_
#define STR_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_


#include "../drt_inpar/inpar_IO_runtime_vtk_output_structure_beams.H"

namespace Teuchos
{
  class ParameterList;
}



namespace DRT
{
  namespace ELEMENTS
  {
    /** \brief Input data container for VTU output at runtime for structure*/
    class StructureRuntimeVtuOutputParams
    {
     public:
      /// constructor
      StructureRuntimeVtuOutputParams();

      /// destructor
      virtual ~StructureRuntimeVtuOutputParams(){};

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtk_structure_paramslist);

      /// setup new class variables
      void Setup();


      /// whether to write displacements
      bool OutputDisplacementState() const
      {
        CheckInitSetup();
        return output_displacement_state_;
      };

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

      /// whether to write displacement output
      bool output_displacement_state_;

      //@}
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
