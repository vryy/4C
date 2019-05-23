/*--------------------------------------------------------------------------*/
/*!
\file particle_ellipsoid_node.H

\brief node representing ellipsoidal particle with possibly specified semi-axes and orientation

\level 2

\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15262
*/
/*--------------------------------------------------------------------------*/
#ifndef PARTICLE_ELLIPSOID_NODE_H
#define PARTICLE_ELLIPSOID_NODE_H

#include "../drt_particle/particle_node.H"

namespace PARTICLE
{
  class ParticleEllipsoidNodeType : public DRT::ParObjectType
  {
   public:
    //! return name
    std::string Name() const { return "ParticleEllipsoidNodeType"; }

    //! return instance
    static ParticleEllipsoidNodeType& Instance() { return instance_; };

    //! create node representing ellipsoidal particle
    DRT::ParObject* Create(const std::vector<char>& data);

   private:
    //! instance
    static ParticleEllipsoidNodeType instance_;
  };


  /*!
  \brief node representing ellipsoidal particle with possibly specified semi-axes and orientation

  This class represents a single ellipsoidal particle with possibly specified semi-axes and
  orientation.

  \author fang (fang@lnm.mw.tum.de)

  \date 10/17
  */
  class ParticleEllipsoidNode : public ParticleNode
  {
   public:
    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor

    \param id       (in): globally unique node ID
    \param coords   (in): vector of nodal coordinates (length 3)
    \param semiaxes (in): vector of semi-axes (length 3, zero if not specified)
    \param orient   (in): vector of orientation (length 3, zero if not specified)
    \param owner    (in): particle owner
    */
    ParticleEllipsoidNode(int id, const double* const coords, const double* const semiaxes,
        const double* const orient, const int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a ParticleEllipsoidNode.
    */
    ParticleEllipsoidNode(const ParticleEllipsoidNode& old);

    /*!
    \brief Deep copy the derived class and return pointer to it.
    */
    PARTICLE::ParticleEllipsoidNode* Clone() const;

    /*!
    \brief Destructor
    */
    virtual ~ParticleEllipsoidNode() { return; };

    /*!
    \brief Return unique ParObject ID.

    Every class implementing ParObject needs a unique ID defined at the top of
    drt_lib/drt_parobject.H.
    */
    int UniqueParObjectId() const
    {
      return ParticleEllipsoidNodeType::Instance().UniqueParObjectId();
    };

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
    \brief Access particle orientation.
    */
    const LINALG::Matrix<3, 1>& Orientation() const { return orientation_; };

    /*!
    \brief Print this ParticleEllipsoidNode.
    */
    void Print(std::ostream& os) const;

    /*!
    \brief Access particle semi-axes.
    */
    const LINALG::Matrix<3, 1>& SemiAxes() const { return semiaxes_; };
    //@}

   private:
    //! particle orientation
    LINALG::Matrix<3, 1> orientation_;

    //! particle semi-axes
    LINALG::Matrix<3, 1> semiaxes_;
  };  // class ParticleEllipsoidNode
}  // namespace PARTICLE
#endif
