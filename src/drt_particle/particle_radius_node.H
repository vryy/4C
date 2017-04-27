/*--------------------------------------------------------------------------*/
/*!
\file particle_radius_node.H

\brief particle node with specified radius

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef PARTICLE_RADIUS_NODE_H
#define PARTICLE_RADIUS_NODE_H

#include "../drt_particle/particle_node.H"

namespace PARTICLE
{
  class ParticleRadiusNodeType : public DRT::ParObjectType
  {
    public:

      //! return name
      std::string Name() const {return "ParticleRadiusNodeType";}

      //! return instance
      static ParticleRadiusNodeType& Instance(){return instance_;};

      //! create particle node with specified radius
      DRT::ParObject* Create(const std::vector<char>& data);

    private:

      //! instance
      static ParticleRadiusNodeType instance_;
  };


  /*!
  \brief particle node with specified radius

  This class represents a single particle with a specified radius in a particle simulation.

  \author fang (fang@lnm.mw.tum.de)

  \date 04/17
  */
  class ParticleRadiusNode : public ParticleNode
  {
    public:

      //! @name Constructors and destructors and related methods
      /*!
      \brief Standard Constructor

      \param id     (in): globally unique node ID
      \param coords (in): vector of nodal coordinates (length 3)
      \param radius (in): particle radius
      \param owner  (in): particle owner
      */
      ParticleRadiusNode(int id, const double* coords, const double radius, const int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a ParticleRadiusNode.
      */
      ParticleRadiusNode(const ParticleRadiusNode& old);

      /*!
      \brief Deep copy the derived class and return pointer to it.
      */
      PARTICLE::ParticleRadiusNode* Clone() const;

      /*!
      \brief Destructor
      */
      virtual ~ParticleRadiusNode(){return;};

      /*!
      \brief Return unique ParObject ID.

      Every class implementing ParObject needs a unique ID defined at the top of drt_lib/drt_parobject.H.
      */
      int UniqueParObjectId() const {return ParticleRadiusNodeType::Instance().UniqueParObjectId();};

      /*!
      \brief Pack this class so it can be communicated.

      \ref Pack and \ref Unpack are used to communicate this node.
      */
      void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class.

      \ref Pack and \ref Unpack are used to communicate this node.
      */
      void Unpack(const std::vector<char>& data);
      //@}

      //! @name Access methods
      /*!
      \brief Print this ParticleRadiusNode.
      */
      void Print(std::ostream& os) const;

      /*!
      \brief Access particle radius.
      */
      double Radius() const {return radius_;};
      //@}

    private:

      //! particle radius
      double radius_;
  }; // class ParticleRadiusNode
} // namespace PARTICLE
#endif
