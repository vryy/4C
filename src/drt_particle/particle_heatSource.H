/*----------------------------------------------------------------------*/
/*!
\file particle_heatSource.H

\brief HeatSource object

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_HEATSOURCE_H_
#define SRC_DRT_PARTICLE_PARTICLE_HEATSOURCE_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"


namespace PARTICLE
{
  class HeatSource
  {
  public:

  /// constructor

  HeatSource(
    bool active,                            ///< flag that simplifies the connectivity algorithm
    const int id,                           ///< unique id over all heat sources conditions
    const std::vector<double> minVerZone,   ///< heat zone, vertex 0
    const std::vector<double> maxVerZone,   ///< heat zone, vertex 1
    const double QDot,                      ///< intake q/s
    const double Tstart,                    ///< t start
    const double Tend) :                    ///< t end
       active_(active),
       id_(id),
       minVerZone_(minVerZone),
       maxVerZone_(maxVerZone),
       QDot_(QDot),
       Tstart_(Tstart),
       Tend_(Tend) {}

  /// access functions (for the Less comparator)

  int Id() { return id_; }

  /// internal variables

  bool active_;                             ///< flag that simplifies the connectivity algorithm
  const int id_;                            ///< unique id over all heat sources conditions
  const std::vector<double> minVerZone_;    ///< heat zone, vertex 0
  const std::vector<double> maxVerZone_;    ///< heat zone, vertex 1
  const double QDot_;                       ///< intake heat per second per unit volume ( q/(s * V) )
  const double Tstart_;                     ///< t start
  const double Tend_;                       ///< t end

  }; // HeatSource
}

#endif /* SRC_DRT_PARTICLE_PARTICLE_HEATSOURCE_H_ */
