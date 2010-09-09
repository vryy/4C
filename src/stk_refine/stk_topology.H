#ifdef STKADAPTIVE

#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "stk_types.H"

namespace STK
{

class UnrefineSet;
class RefineSet;
class Mesh;

template <typename topology>
class TopologyType;

template <>
class TopologyType<void>
{
public:
  virtual ~TopologyType() {}

  virtual int NumChildNodes() const = 0;

  virtual int NumChildElements() const = 0;

  virtual void FindChildNodes(RefineSet& hns,
                              const stk::mesh::Entity* e) const = 0;

  virtual void HangingNodes( const RefineSet* rs,
                             Mesh* mesh,
                             const std::vector<stk::mesh::Entity*> & nodes,
                             std::vector<stk::mesh::Entity*> & hanging ) const = 0;

  virtual void CreateElements( Mesh* mesh,
                               stk::mesh::EntityRank rank,
                               const std::pair<const unsigned *, const unsigned *> & part_ord,
                               std::vector<stk::mesh::Entity*> & nodes,
                               std::vector<stk::mesh::Entity*> & hanging,
                               std::vector<stk::mesh::Entity *>::iterator & elements_iter ) const = 0;

  static const TopologyType* GetTopologyType(const stk::mesh::Entity* e);
};


template <>
class TopologyType<shards::Line<2> > : public TopologyType<void>
{
  virtual int NumChildNodes() const { return 1; }

  virtual int NumChildElements() const { return 2; }

  virtual void FindChildNodes(RefineSet& hns,
                              const stk::mesh::Entity* e) const;

  virtual void HangingNodes( const RefineSet* rs,
                             Mesh* mesh,
                             const std::vector<stk::mesh::Entity*> & nodes,
                             std::vector<stk::mesh::Entity*> & hanging ) const;

  virtual void CreateElements( Mesh* mesh,
                               stk::mesh::EntityRank rank,
                               const std::pair<const unsigned *, const unsigned *> & part_ord,
                               std::vector<stk::mesh::Entity*> & nodes,
                               std::vector<stk::mesh::Entity*> & hanging,
                               std::vector<stk::mesh::Entity *>::iterator & elements_iter ) const;

};


template <>
class TopologyType<shards::Quadrilateral<4> > : public TopologyType<void>
{
  virtual int NumChildNodes() const { return 5; }

  virtual int NumChildElements() const { return 4; }

  virtual void FindChildNodes(RefineSet& hns,
                              const stk::mesh::Entity* e) const;

  virtual void HangingNodes( const RefineSet* rs,
                             Mesh* mesh,
                             const std::vector<stk::mesh::Entity*> & nodes,
                             std::vector<stk::mesh::Entity*> & hanging ) const;

  virtual void CreateElements( Mesh* mesh,
                               stk::mesh::EntityRank rank,
                               const std::pair<const unsigned *, const unsigned *> & part_ord,
                               std::vector<stk::mesh::Entity*> & nodes,
                               std::vector<stk::mesh::Entity*> & hanging,
                               std::vector<stk::mesh::Entity *>::iterator & elements_iter ) const;

};


template <>
class TopologyType<shards::Hexahedron<8> > : public TopologyType<void>
{
  virtual int NumChildNodes() const { return 19; }

  virtual int NumChildElements() const { return 8; }

  virtual void FindChildNodes(RefineSet& hns,
                              const stk::mesh::Entity* e) const;

  virtual void HangingNodes( const RefineSet* rs,
                             Mesh* mesh,
                             const std::vector<stk::mesh::Entity*> & nodes,
                             std::vector<stk::mesh::Entity*> & hanging ) const;

  virtual void CreateElements( Mesh* mesh,
                               stk::mesh::EntityRank rank,
                               const std::pair<const unsigned *, const unsigned *> & part_ord,
                               std::vector<stk::mesh::Entity*> & nodes,
                               std::vector<stk::mesh::Entity*> & hanging,
                               std::vector<stk::mesh::Entity *>::iterator & elements_iter ) const;

};

}

#endif
#endif
