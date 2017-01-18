/*---------------------------------------------------------------------*/
/*!
\file particle_interData.H

\brief Interaction data objects: particle vs particle and particle vs wall

\level 3

\maintainer Alessandro Cattabiani
*----------------------------------------------------------------------*/

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
    const double weightDerivative_ij,
    const double weightDerivative_ji,
    const double weight_ij,
    const double weight_ji,
    const double rRelNorm2,
    const double step
    ) :
      rRelVersor_ij_(rRelVersor_ij),
      weightDerivative_ij_(weightDerivative_ij),
      weightDerivative_ji_(weightDerivative_ji),
      weight_ij_(weight_ij),
      weight_ji_(weight_ji),
      rRelNorm2_(rRelNorm2),
      step_(step) {}

    InterDataPvP() :
      weightDerivative_ij_(0.0),
      weightDerivative_ji_(0.0),
      weight_ij_(0.0),
      weight_ji_(0.0),
      rRelNorm2_(0.0),
      step_(-1) {}

    ~InterDataPvP() {}

    /// internal variables
    LINALG::Matrix<3,1> rRelVersor_ij_;        ///< versor i->j
    double weightDerivative_ij_;               ///< weightDerivative i->j
    double weightDerivative_ji_;               ///< weightDerivative j->i
    double weight_ij_;                         ///< weight i->j
    double weight_ji_;                         ///< weight j->i
    double rRelNorm2_;                         ///< Norm2
    int step_;                                 ///< time step of the interaction

  }; // InterDataPvP

// -------------------------------------------------------------------- //
  /// small interaction data object for pvw
  class InterDataPvW
  {
  public:

  /// constructors

    InterDataPvW(
    const LINALG::Matrix<3,1> rRelVersor,
    const double weightDerivative,
    const double rRelNorm2,
    const int step
    ) :
      rRelVersor_(rRelVersor),
      weightDerivative_(weightDerivative),
      rRelNorm2_(rRelNorm2),
      step_(step) {}

    InterDataPvW() :
       weightDerivative_(0.0),
       rRelNorm2_(0.0),
       step_(-1) {}

    ~InterDataPvW() {}

    /// internal variables
    LINALG::Matrix<3,1> rRelVersor_;         ///< versor toward the particle
    double weightDerivative_;                ///< weightDerivative toward the particle
    double rRelNorm2_;                       ///< Norm2
    int step_;                               ///< time step of the interaction

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
