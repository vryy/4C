/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_solid_volume_meshtying_evaluation_data.H

\brief data container holding all beam to solid volume meshtying
       evaluation data, that can not be stored pairwise.

\level 3

\maintainer Alexander Popp
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAMTOSOLIDVOLUMEMESHTYINGEVALUATIONDATA_H
#define BEAMTOSOLIDVOLUMEMESHTYINGEVALUATIONDATA_H

#include <Teuchos_RCP.hpp>

#include "../drt_lib/drt_dserror.H"


namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingEvaluationData
  {
   public:
    //! constructor
    BeamToSolidVolumeMeshtyingEvaluationData();

    //! destructor
    virtual ~BeamToSolidVolumeMeshtyingEvaluationData(){};

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit()) dserror("Init() has not been called, yet!");
    }

    inline std::map<int, std::vector<bool>>* BTSVOLMTGaussPointProjectionTracker() const
    {
      return BTSVOLMT_GaussPointProjectionTracker_ptr_;
    }

    inline std::map<int, std::vector<bool>>* BTSVOLMTGaussPointEvaluateTracker() const
    {
      return BTSVOLMT_GaussPointEvaluateTracker_ptr_;
    }


   private:
    bool isinit_;

    bool issetup_;

    /*!
    \brief map that keeps track which Gauss points of a beam element have a valid projection
           int --> global element id of the beam element
           vector of bools --> one bool for each Gauss point of the beam element
    */
    std::map<int, std::vector<bool>>* BTSVOLMT_GaussPointProjectionTracker_ptr_;

    /*!
    \brief map that keeps track which Gauss points of a beam element have a valid projection
           int --> global element id of the beam element
           vector of bools --> one bool for each Gauss point of the beam element
    */
    std::map<int, std::vector<bool>>* BTSVOLMT_GaussPointEvaluateTracker_ptr_;
  };

}  // namespace BEAMINTERACTION

#endif
