#ifndef CUT_TEST_GENERATOR_H
#define CUT_TEST_GENERATOR_H

#include "../../src/drt_cut/cut_mesh.H"

class OutputGenerator : public GEO::CUT::CellGenerator
{
public:

  OutputGenerator() : counter_( 0 ) {}

#ifdef QHULL
  virtual void Generate( GEO::CUT::Element* element, const tetgenio & out );
#endif

private:
  int counter_;
};

#endif
