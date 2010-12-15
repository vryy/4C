#ifndef CELL_FACET_H
#define CELL_FACET_H

#include <vector>

#include "cell_line.H"

namespace GEO
{
namespace CELL
{

class Facet
{
public:

  Facet( const std::vector<Point*> & ps, const std::vector<Line*> & ls );

};

}
}

#endif
