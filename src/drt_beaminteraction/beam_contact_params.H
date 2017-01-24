/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_contact_params.H

\brief data container holding all beam contact input parameters

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMCONTACTPARAMS_H
#define BEAMCONTACTPARAMS_H

#include "../drt_inpar/inpar_beamcontact.H"

#include "../drt_lib/drt_dserror.H"

namespace BEAMINTERACTION {

  /*!
   *  */
  class BeamContactParams
  {
  public:
    //! constructor
    BeamContactParams();

    //! destructor
    virtual ~BeamContactParams() {};

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const
    { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

    inline enum INPAR::BEAMCONTACT::Strategy Strategy() const
    {
      return strategy_;
    }

    inline enum INPAR::BEAMCONTACT::PenaltyLaw PenaltyLaw() const
    {
      return penalty_law_;
    }

    inline double BeamToBeamPenaltyLawRegularizationG0() const
    {
      return BTB_penalty_law_regularization_G0_;
    }

    inline double BeamToBeamPenaltyLawRegularizationF0() const
    {
      return BTB_penalty_law_regularization_F0_;
    }

    inline double BeamToBeamPenaltyLawRegularizationC0() const
    {
      return BTB_penalty_law_regularization_C0_;
    }

    inline double GapShift() const
    {
      return gap_shift_;
    }

    inline double BeamToBeamPointPenaltyParam() const
    {
      return BTB_point_penalty_param_;
    }

    inline double BeamToBeamLinePenaltyParam() const
    {
      return BTB_line_penalty_param_;
    }

    inline double BeamToBeamPerpShiftingAngle1() const
    {
      return BTB_perp_shifting_angle1_;
    }

    inline double BeamToBeamPerpShiftingAngle2() const
    {
      return BTB_perp_shifting_angle2_;
    }

    inline double BeamToBeamParallelShiftingAngle1() const
    {
      return BTB_parallel_shifting_angle1_;
    }

    inline double BeamToBeamParallelShiftingAngle2() const
    {
      return BTB_parallel_shifting_angle2_;
    }

    inline double SegmentationAngle() const
    {
      return segangle_;
    }

    inline int NumIntegrationIntervals() const
    {
      return num_integration_intervals_;
    }

    inline double BasicStiffGap() const
    {
      return BTB_basicstiff_gap_;
    }

    inline bool EndPointPenalty() const
    {
      return BTB_endpoint_penalty_;
    }

    // ********** Fixme: this needs to go to BeamToSphereContactParams ******************
    inline double BeamToSpherePenaltyParam() const
    {
      return BTSPH_penalty_param_;
    }
    // ***********************************************************************************

  private:

    bool isinit_;

    bool issetup_;

    //! strategy
    enum INPAR::BEAMCONTACT::Strategy strategy_;

    //! penalty law
    enum INPAR::BEAMCONTACT::PenaltyLaw penalty_law_;

    //! regularization parameters for penalty law
    double BTB_penalty_law_regularization_G0_;
    double BTB_penalty_law_regularization_F0_;
    double BTB_penalty_law_regularization_C0_;

    //! Todo understand and check usage of this parameter
    double gap_shift_;

    //! beam-to-beam point penalty parameter
    double BTB_point_penalty_param_;

    //! beam-to-beam line penalty parameter
    double BTB_line_penalty_param_;

    //! shifting angles [radians] for point contact (near perpendicular configurations) fade
    double BTB_perp_shifting_angle1_;
    double BTB_perp_shifting_angle2_;

    //! shifting angles [radians] for line contact (near parallel configurations) fade
    double BTB_parallel_shifting_angle1_;
    double BTB_parallel_shifting_angle2_;

    //! maximum difference in tangent orientation between the endpoints of one created segment
    //  if the angle is larger => subdivide and create more segments
    double segangle_;

    //! number of integration intervals
    int num_integration_intervals_;

    //! gap value, from which on only the basic part of the stiffness contribution is applied
    //  this should accelerate convergence
    double BTB_basicstiff_gap_;

    //! flag indicating if the integration should take special care of physical
    //  end points of beams in order to avoid strong discontinuities
    bool BTB_endpoint_penalty_;

    // ********** Fixme: this needs to go to BeamToSphereContactParams ******************

    //! beam-to-sphere penalty parameter
    double BTSPH_penalty_param_;

    // ***********************************************************************************

  };

}

#endif
