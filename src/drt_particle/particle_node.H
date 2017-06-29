/*----------------------------------------------------------------------*/
/*!
\file particle_node.H

\brief A particle is a DRT::Node with additional knowledge of its collision status

\level 2

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_NODE_H
#define PARTICLE_NODE_H

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"


/*!
\brief namespace PARTICLE: Particle simulations

*/
namespace PARTICLE
{

class ParticleNodeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ParticleNodeType"; }

  static ParticleNodeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ParticleNodeType instance_;
};


//! @name Struct to store collision info

struct Collision {
  bool stick;
  double g_t[3];
};

//@}


/*!
\brief A class for a particle derived from DRT::Node

This class represents a single particle involved in a particle collision simulation.

\author hammerl (hammerl@lnm.mw.tum.de)
*/
class ParticleNode : public DRT::Node
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of ParticleNode
  */
//  friend class DRT::Discretization;
  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id     (in): A globally unique node id
  \param coords (in): vector of nodal coordinates, length 3
  \param owner  (in): Owner of this node.

  */
  ParticleNode(int id, const double* coords, const int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a ParticleNode

  */
  ParticleNode(const PARTICLE::ParticleNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual PARTICLE::ParticleNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~ParticleNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.

  */
  virtual int UniqueParObjectId() const { return ParticleNodeType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Unpack(const std::vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Print this ParticleNode
  */
  virtual void Print(std::ostream& os) const ;

  /*!
  \brief Return collision history of this ParticleNode with particles
  */
  virtual std::map<int, Collision>& Get_history_particle(){ return  history_particle_; }

  /*!
  \brief Return collision history of this ParticleNode with walls
  */
  virtual std::map<int, Collision>& Get_history_wall(){ return  history_wall_; }

  /*!
  \brief Return flag for boundary particles
  */
  virtual bool Is_bdry_particle(){ return  is_bdry_particle_; }

  /*!
  \brief Set flag for boundary particles
  */
  virtual void Set_bdry_particle(bool is_bdry_particle)
  {
    is_bdry_particle_ = is_bdry_particle;
    return;
  }

  //@}

protected:

  std::map<int, Collision> history_particle_;
  std::map<int, Collision> history_wall_;

  bool is_bdry_particle_;

}; // class ParticleNode
} // namespace PARTICLE


// << operator
std::ostream& operator << (std::ostream& os, const PARTICLE::ParticleNode& particle);


#endif  // #ifndef PARTICLE_NODE_H
