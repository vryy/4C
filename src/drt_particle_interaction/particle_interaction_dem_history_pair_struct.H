/*---------------------------------------------------------------------------*/
/*!
\brief history pair struct for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_HISTORY_PAIR_STRUCT_H
#define PARTICLE_INTERACTION_DEM_HISTORY_PAIR_STRUCT_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
#include "../drt_lib/drt_parobject.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! struct to store tangential contact history of interacting particles
  struct DEMHistoryPairTangential
  {
    //! tangential stick flag
    bool stick_ = true;

    //! tangential gap
    double gap_t_[3] = {0.0};

    //! pack history pair data
    void Pack(DRT::PackBuffer& data) const
    {
      data.AddtoPack(stick_);

      for (int i = 0; i < 3; ++i) data.AddtoPack(gap_t_[i]);
    }

    //! unpack history pair data
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      DRT::ParObject::ExtractfromPack(position, data, stick_);

      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, data, gap_t_[i]);
    }
  };

  //! struct to store rolling contact history of interacting particles
  struct DEMHistoryPairRolling
  {
    //! rolling stick flag
    bool stick_ = true;

    //! rolling gap
    double gap_r_[3] = {0.0};

    //! pack history pair data
    void Pack(DRT::PackBuffer& data) const
    {
      data.AddtoPack(stick_);

      for (int i = 0; i < 3; ++i) data.AddtoPack(gap_r_[i]);
    }

    //! unpack history pair data
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      DRT::ParObject::ExtractfromPack(position, data, stick_);

      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, data, gap_r_[i]);
    }
  };

  //! struct to store adhesion history of interacting particles
  struct DEMHistoryPairAdhesion
  {
    //! surface energy
    double surface_energy_ = 0.0;

    //! adhesion force
    double adhesion_force_ = 0.0;

    //! pack history pair data
    void Pack(DRT::PackBuffer& data) const
    {
      data.AddtoPack(surface_energy_);

      data.AddtoPack(adhesion_force_);
    }

    //! unpack history pair data
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      DRT::ParObject::ExtractfromPack(position, data, surface_energy_);

      DRT::ParObject::ExtractfromPack(position, data, adhesion_force_);
    }
  };
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
