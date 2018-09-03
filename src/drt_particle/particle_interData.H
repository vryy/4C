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

    InterDataPvP(const LINALG::Matrix<3, 1> rRelVersor_ij, const double rRelNorm2,
        const double w_ij, const double dw_ij, const double ddw_ij)
        : rRelVersor_ij_(rRelVersor_ij),
          rRelNorm2_(rRelNorm2),
          w_ij_(w_ij),
          dw_ij_(dw_ij),
          ddw_ij_(ddw_ij)
    {
    }

    InterDataPvP() : rRelNorm2_(0.0), w_ij_(0.0), dw_ij_(0.0), ddw_ij_(0.0) {}

    ~InterDataPvP() {}

    /// internal variables
    LINALG::Matrix<3, 1> rRelVersor_ij_;  ///< versor i->j
    double rRelNorm2_;                    ///< Norm2
    double w_ij_;                         ///< weight i->j
    double dw_ij_;                        ///< weight 1 derivative i->j
    double ddw_ij_;                       ///< weight 2 derivative i->j

  };  // InterDataPvP

}  // namespace PARTICLE

#endif /* SRC_DRT_PARTICLE_PARTICLE_INTERDATA_H_ */
