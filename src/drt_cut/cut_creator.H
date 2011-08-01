#ifndef CUT_CREATOR_H
#define CUT_CREATOR_H

#include <vector>

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Side;
    class Element;
    class Mesh;

    /// delayed line object creation
    /*!
     * Unused.
     */
    class Creator
    {
    public:

      void NewLine( Point * p1, Point * p2, Side * side, Side * other, Element * element );

      void Execute( Mesh & mesh );

    private:

      struct LineCmd
      {
        void Assign( Point * p1, Point * p2, Side * side, Side * other, Element * element );

        void Execute( Mesh & mesh );

        Point * p1_;
        Point * p2_;
        Side * side_;
        Side * other_;
        Element * element_;
      };

      std::vector<LineCmd> lines_;
    };
  }
}

#endif
