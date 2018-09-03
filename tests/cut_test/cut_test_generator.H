// Todo: Delete Me!!!

#ifndef CUT_TEST_GENERATOR_H
#define CUT_TEST_GENERATOR_H

#include "../../src/drt_cut/cut_mesh.H"

#if 0

class OutputGenerator : public GEO::CUT::CellGenerator
{
public:

  OutputGenerator() : counter_( 0 ) {}

  virtual void Generate( GEO::CUT::Element* element, const tetgenio & out );

private:
  int counter_;
};

#endif

#endif
