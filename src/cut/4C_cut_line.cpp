/*---------------------------------------------------------------------*/
/*! \file

\brief cut line

\level 3


*----------------------------------------------------------------------*/

#include "4C_cut_line.hpp"

#include "4C_cut_element.hpp"
#include "4C_cut_side.hpp"

FOUR_C_NAMESPACE_OPEN

Core::Geo::Cut::Line::Line(
    Point* p1, Point* p2, Side* cut_side1, Side* cut_side2, Element* cut_element)
    : p1_(p1), p2_(p2)
{
  if (cut_side1) add_side(cut_side1);
  if (cut_side2) add_side(cut_side2);
  if (cut_element) add_element(cut_element);

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
    add_side(s);
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
    add_element(s);
  }
}

void Core::Geo::Cut::Line::add_side(Side* cut_side)
{
  p1_->add_side(cut_side);
  p2_->add_side(cut_side);
  cut_sides_.insert(cut_side);
  cut_side->AddLine(this);
}

void Core::Geo::Cut::Line::add_element(Element* cut_element)
{
  if (cut_element != nullptr)
  {
    if (not p1_->IsCut(cut_element) or not p2_->IsCut(cut_element))
    {
      FOUR_C_THROW("cut line between non-cut points");
    }

    cut_elements_.insert(cut_element);

    p1_->add_element(cut_element);
    p2_->add_element(cut_element);
  }
}

bool Core::Geo::Cut::Line::IsInternalCut(Side* side)
{
  return cut_sides_.count(side) > 0 and not side->OnEdge(this);
}

FOUR_C_NAMESPACE_CLOSE
