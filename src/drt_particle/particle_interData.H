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

// -------------------------------------------------------------------- //
  /// small interaction data object for pvw
  class InterDataPvW
  {
  public:

  /// constructors

    InterDataPvW(
    const LINALG::Matrix<3,1> rRelVersor,
    const double rRelNorm2,
    const int step,
    const double w,
    const double dw,
    const double ddw
    ) :
      rRelVersor_(rRelVersor),
      rRelNorm2_(rRelNorm2),
      step_(step),
      w_(w),
      dw_(dw),
      ddw_(ddw) {}

    InterDataPvW() :
      rRelNorm2_(0.0),
      step_(-1),
      w_(0.0),
      dw_(0.0),
      ddw_(0.0) {}

    ~InterDataPvW() {}

    /// internal variables
    LINALG::Matrix<3,1> rRelVersor_;       ///< versor i->j
    double rRelNorm2_;                     ///< Norm2
    int step_;                             ///< time step of the interaction
    double w_;                             ///< weight
    double dw_;                            ///< weight 1 derivative
    double ddw_;                           ///< weight 2 derivative

  }; // InterDataPvW

  // -------------------------------------------------------------------- //
    /// struct to store wall contact point and additional information
    /// about the corresponding wall element
  class WallInteractionPoint
  {
  public :

    // constructors

  WallInteractionPoint(
      const int elemId,
      const std::vector<int> lm,
      const std::vector<int> lmowner,
      const double penetration,
      const LINALG::Matrix<3,1> point,
      const std::map<int,LINALG::Matrix<3,1> > nodalCoordinates
      ) :
        elemId_(elemId),
        lm_(lm),
        lmowner_(lmowner),
        penetration_(penetration),
        point_(point),
        nodalCoordinates_(nodalCoordinates) {}

  WallInteractionPoint(
      const int elemId,
      const std::vector<int> lm,
      const std::vector<int> lmowner
      ) :
        elemId_(elemId),
        lm_(lm),
        lmowner_(lmowner),
        penetration_(-1) {}

  WallInteractionPoint() :
    elemId_(-1),
    penetration_(-1) {}

  // internal variables

  int elemId_;
  std::vector<int> lm_;
  std::vector<int> lmowner_;
  double penetration_;
  LINALG::Matrix<3,1> point_;
  std::map<int,LINALG::Matrix<3,1> > nodalCoordinates_;

  };

}

#endif /* SRC_DRT_PARTICLE_PARTICLE_INTERDATA_H_ */
