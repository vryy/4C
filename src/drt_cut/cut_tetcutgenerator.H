#ifndef CUT_TETCUTGENERATOR_H
#define CUT_TETCUTGENERATOR_H

#include "cut_mesh.H"
#include "cut_meshintersection.H"

namespace GEO
{
  namespace CUT
  {

class TetCutGenerator : public CellGenerator
{
public:

  TetCutGenerator( GEO::CUT::MeshIntersection & intersection,
                   CellGenerator * parent,
                   Mesh & cut_mesh,
                   Teuchos::RCP<PointPool> pp );

  virtual void NonCut( Element* element );

#ifdef QHULL
  virtual void Generate( Element* element, const tetgenio & out );
#endif

private:
  GEO::CUT::MeshIntersection & intersection_;
  CellGenerator * parent_;
  Mesh & cut_mesh_;
  Teuchos::RCP<PointPool> pp_;
};

class TetCutConverter : public CellGenerator
{
public:

  TetCutConverter( Mesh & mesh );

  virtual void NonCut( Element* element );

#ifdef QHULL
  virtual void Generate( Element* element, const tetgenio & out );
#endif

private:
  Mesh & mesh_;
};

  }
}

#endif
