/*---------------------------------------------------------------------*/
/*! \file

\brief cut line

\level 3


*----------------------------------------------------------------------*/

#include "cut_line.H"
#include "cut_element.H"
#include "cut_side.H"

CORE::GEO::CUT::Line::Line(
    Point* p1, Point* p2, Side* cut_side1, Side* cut_side2, Element* cut_element)
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
}

void CORE::GEO::CUT::Line::AddSide(Side* cut_side)
{
  p1_->AddSide(cut_side);
  p2_->AddSide(cut_side);
  cut_sides_.insert(cut_side);
  cut_side->AddLine(this);
}

void CORE::GEO::CUT::Line::AddElement(Element* cut_element)
{
  if (cut_element != NULL)
  {
    if (not p1_->IsCut(cut_element) or not p2_->IsCut(cut_element))
    {
      throw std::runtime_error("cut line between non-cut points");
    }

    cut_elements_.insert(cut_element);

    p1_->AddElement(cut_element);
    p2_->AddElement(cut_element);
  }
}

bool CORE::GEO::CUT::Line::IsInternalCut(Side* side)
{
  return cut_sides_.count(side) > 0 and not side->OnEdge(this);
}
