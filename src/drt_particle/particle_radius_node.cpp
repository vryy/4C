/*--------------------------------------------------------------------------*/
/*!

\brief particle node with specified radius

\level 2

\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15262
*/
/*--------------------------------------------------------------------------*/
#include "particle_radius_node.H"

// instantiation
PARTICLE::ParticleRadiusNodeType PARTICLE::ParticleRadiusNodeType::instance_;

/*----------------------------------------------------------------------*
 | create particle node with specified radius                fang 04/17 |
 *----------------------------------------------------------------------*/
DRT::ParObject* PARTICLE::ParticleRadiusNodeType::Create(const std::vector<char>& data)
{
  double dummycoords[3] = {0., 0., 0.};
  DRT::Node* const particle = new PARTICLE::ParticleRadiusNode(-1, dummycoords, 0., -1);
  particle->Unpack(data);
  return particle;
}


/*----------------------------------------------------------------------*
 | standard constructor                                      fang 04/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleRadiusNode::ParticleRadiusNode(
    int id, const double* coords, const double radius, const int owner)
    : ParticleNode(id, coords, owner), radius_(radius)
{
  return;
}


/*----------------------------------------------------------------------*
 | copy constructor                                          fang 04/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleRadiusNode::ParticleRadiusNode(const ParticleRadiusNode& old)
    : ParticleNode(old), radius_(old.radius_)
{
  dserror("Not yet needed!");
  return;
}


/*----------------------------------------------------------------------*
 | deep copy the derived class and return pointer to it      fang 04/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleRadiusNode* PARTICLE::ParticleRadiusNode::Clone() const
{
  return new PARTICLE::ParticleRadiusNode(*this);
}


/*----------------------------------------------------------------------*
 | pack this class so it can be communicated                 fang 04/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleRadiusNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  AddtoPack(data, UniqueParObjectId());

  // add base class
  ParticleNode::Pack(data);

  // add particle radius
  AddtoPack(data, radius_);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector into this class            fang 04/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleRadiusNode::Unpack(const std::vector<char>& data)
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

  // extract particle radius
  ExtractfromPack(position, data, radius_);

  // safety check
  if (position != data.size())
    dserror("Mismatch in size of data: %d <-> %d!", (int)data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 | print this ParticleRadiusNode                             fang 04/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleRadiusNode::Print(std::ostream& os) const
{
  // print base class information
  ParticleNode::Print(os);

  // print particle radius
  os << "  with radius  " << radius_ << std::endl;

  return;
}
