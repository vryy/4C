/*----------------------------------------------------------------------*/
/*!
\file particle_interData.H

\brief Interaction data objects: particle vs particle and particle vs wall

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_INTERDATA_H_
#define SRC_DRT_PARTICLE_PARTICLE_INTERDATA_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"

namespace PARTICLE
{
// -------------------------------------------------------------------- //
  /// small interaction data object for pvp
  class InterDataPvP
  {
  public:

    /// constructors

    InterDataPvP(
    const LINALG::Matrix<3,1> rRelVersor_ij,
    const double rRelNorm2,
    const double step,
    const double w_ij,
    const double w_ji,
    const double dw_ij,
    const double dw_ji,
    const double ddw_ij,
    const double ddw_ji
    ) :
      rRelVersor_ij_(rRelVersor_ij),
      rRelNorm2_(rRelNorm2),
      step_(step),
      w_ij_(w_ij),
      w_ji_(w_ji),
      dw_ij_(dw_ij),
      dw_ji_(dw_ji),
      ddw_ij_(ddw_ij),
      ddw_ji_(ddw_ji) {}

    InterDataPvP() :
      rRelNorm2_(0.0),
      step_(-1),
      w_ij_(0.0),
      w_ji_(0.0),
      dw_ij_(0.0),
      dw_ji_(0.0),
      ddw_ij_(0.0),
      ddw_ji_(0.0) {}

    ~InterDataPvP() {}

    /// internal variables
    LINALG::Matrix<3,1> rRelVersor_ij_;       ///< versor i->j
    double rRelNorm2_;                        ///< Norm2
    int step_;                                ///< time step of the interaction
    double w_ij_;                             ///< weight i->j
    double w_ji_;                             ///< weight j->i
    double dw_ij_;                            ///< weight 1 derivative i->j
    double dw_ji_;                            ///< weight 1 derivative j->i
    double ddw_ij_;                           ///< weight 2 derivative i->j
    double ddw_ji_;                           ///< weight 2 derivative j->i

  }; // InterDataPvP

}

#endif /* SRC_DRT_PARTICLE_PARTICLE_INTERDATA_H_ */
