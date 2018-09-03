/*--------------------------------------------------------------------------*/
/*!
\file particle_ellipsoid_node.cpp

\brief node representing ellipsoidal particle with possibly specified semi-axes and orientation

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "particle_ellipsoid_node.H"

// instantiation
PARTICLE::ParticleEllipsoidNodeType PARTICLE::ParticleEllipsoidNodeType::instance_;

/*----------------------------------------------------------------------*
 | create particle node with specified radius                fang 10/17 |
 *----------------------------------------------------------------------*/
DRT::ParObject* PARTICLE::ParticleEllipsoidNodeType::Create(const std::vector<char>& data)
{
  const double dummy[3] = {0., 0., 0.};
  DRT::Node* const particle = new PARTICLE::ParticleEllipsoidNode(-1, dummy, dummy, dummy, -1);
  particle->Unpack(data);
  return particle;
}


/*----------------------------------------------------------------------*
 | standard constructor                                      fang 10/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleEllipsoidNode::ParticleEllipsoidNode(int id, const double* const coords,
    const double* const semiaxes, const double* const orient, const int owner)
    : ParticleNode(id, coords, owner), orientation_(orient), semiaxes_(semiaxes)
{
  return;
}


/*----------------------------------------------------------------------*
 | copy constructor                                          fang 10/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleEllipsoidNode::ParticleEllipsoidNode(const ParticleEllipsoidNode& old)
    : ParticleNode(old), orientation_(old.orientation_), semiaxes_(old.semiaxes_)
{
  dserror("Not yet needed!");
  return;
}


/*----------------------------------------------------------------------*
 | deep copy the derived class and return pointer to it      fang 10/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleEllipsoidNode* PARTICLE::ParticleEllipsoidNode::Clone() const
{
  return new PARTICLE::ParticleEllipsoidNode(*this);
}


/*----------------------------------------------------------------------*
 | pack this class so it can be communicated                 fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleEllipsoidNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  AddtoPack(data, UniqueParObjectId());

  // add base class
  ParticleNode::Pack(data);

  // add particle orientation
  AddtoPack(data, orientation_);

  // add particle semi-axes
  AddtoPack(data, semiaxes_);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector into this class            fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleEllipsoidNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position(0);

  // extract type
  int type(0);
  ExtractfromPack(position, data, type);

  // safety check
  if (type != UniqueParObjectId()) dserror("Wrong instance type!");

  // extract base class
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ParticleNode::Unpack(basedata);

  // extract particle orientation
  ExtractfromPack(position, data, orientation_);

  // extract particle semi-axes
  ExtractfromPack(position, data, semiaxes_);

  // safety check
  if (position != data.size())
    dserror("Mismatch in size of data: %d <-> %d!", (int)data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 | print this ParticleEllipsoidNode                          fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleEllipsoidNode::Print(std::ostream& os) const
{
  // print base class information
  ParticleNode::Print(os);

  // print particle orientation and semi-axes
  os << "  with orientation  " << orientation_ << "  and semi-axes  " << semiaxes_ << std::endl;

  return;
}
