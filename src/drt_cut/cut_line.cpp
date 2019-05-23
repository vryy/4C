/*---------------------------------------------------------------------*/
/*!

\brief cut line

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include "cut_line.H"
#include "cut_element.H"
#include "cut_side.H"

GEO::CUT::Line::Line(Point* p1, Point* p2, Side* cut_side1, Side* cut_side2, Element* cut_element)
    : p1_(p1), p2_(p2)
{
  if (cut_side1) AddSide(cut_side1);
  if (cut_side2) AddSide(cut_side2);
  if (cut_element) AddElement(cut_element);

  p1_->Register(this);
  p2_->Register(this);

  // self-register at all sides

  const plain_side_set& cut_sides1 = p1->CutSides();
  const plain_side_set& cut_sides2 = p2->CutSides();

  plain_side_set sides;

  std::set_intersection(cut_sides1.begin(), cut_sides1.end(), cut_sides2.begin(), cut_sides2.end(),
      std::inserter(sides, sides.begin()));

  for (plain_side_set::iterator i = sides.begin(); i != sides.end(); ++i)
  {
    Side* s = *i;
    AddSide(s);
  }

  // self-register at all elements

  const plain_element_set& cut_elements1 = p1->Elements();
  const plain_element_set& cut_elements2 = p2->Elements();

  plain_element_set elements;

  std::set_intersection(cut_elements1.begin(), cut_elements1.end(), cut_elements2.begin(),
      cut_elements2.end(), std::inserter(elements, elements.begin()));

  for (plain_element_set::iterator i = elements.begin(); i != elements.end(); ++i)
  {
    Element* s = *i;
    AddElement(s);
  }

#if 0
#ifdef DEBUGCUTLIBRARY
  double x1[] = { 1.0571400000000001906, 0.49999999999999994449, -0.024639335281227081609 };
  double x2[] = { 1.0571400000000001906, 0.49999999999999994449, -0.050000000000000009714 };

  LINALG::Matrix<3,1> px1( p1_->X() );
  LINALG::Matrix<3,1> px2( p2_->X() );

  LINALG::Matrix<3,1> mx1( x1 );
  LINALG::Matrix<3,1> mx2( x2 );

  LINALG::Matrix<3,1> d1;
  LINALG::Matrix<3,1> d2;

  d1.Update( 1, px1, -1, mx1, 0 );
  d2.Update( 1, px2, -1, mx2, 0 );

  if ( d1.Norm2() < 1e-12 and d2.Norm2() < 1e-12 )
  {
    std::cout << "offending line 1\n";
  }

  d1.Update( 1, px1, -1, mx2, 0 );
  d2.Update( 1, px2, -1, mx1, 0 );

  if ( d1.Norm2() < 1e-12 and d2.Norm2() < 1e-12 )
  {
    std::cout << "offending line 2\n";
  }
#endif
#endif
}

void GEO::CUT::Line::AddSide(Side* cut_side)
{
  p1_->AddSide(cut_side);
  p2_->AddSide(cut_side);
  cut_sides_.insert(cut_side);
  cut_side->AddLine(this);
}

void GEO::CUT::Line::AddElement(Element* cut_element)
{
  if (cut_element != NULL)
  {
#if 0
#ifdef DEBUGCUTLIBRARY
    Node * n1 = NULL;
    Node * n2 = NULL;
    const std::vector<Node*> & nodes = cut_element->Nodes();
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      if ( n->point()==p1_ )
      {
        n1 = n;
      }
      if ( n->point()==p2_ )
      {
        n2 = n;
      }
    }
    if ( n1!=NULL and n2!=NULL )
    {
      plain_edge_set edges;
      p1_->CommonEdge( p2_, edges );
      if ( edges.size()!=1 )
        throw std::runtime_error( "line does not belong to element" );
      Edge * e = *edges.begin();

      const plain_side_set & edge_sides = e->Sides();
      const std::vector<Side*> & element_sides = cut_element->Sides();

      bool found = false;
      for ( std::vector<Side*>::const_iterator i=element_sides.begin(); i!=element_sides.end(); ++i )
      {
        Side * s = *i;
        if ( edge_sides.count( s )>0 )
        {
          found = true;
        }
      }
      if ( not found )
      {
        throw std::runtime_error( "line does not belong to element" );
      }
    }
#endif
#endif
    cut_elements_.insert(cut_element);
#if 1
    if (not p1_->IsCut(cut_element) or not p2_->IsCut(cut_element))
    {
      throw std::runtime_error("cut line between non-cut points");
    }
#else
    p1_->AddElement(cut_element);
    p2_->AddElement(cut_element);
#endif
  }
}

bool GEO::CUT::Line::IsInternalCut(Side* side)
{
  return cut_sides_.count(side) > 0 and not side->OnEdge(this);
}
