/*!
\file
\brief Information about the elements.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is the central place where information about elements is
kept. When elements are added or changed this file needs to be
changed, too. The binary io module, the filters and maybe others
depend on it.

\author u.kue
\date 09/04

*/

#include "../headers/standardtypes.h"
#include "global_element_info.h"


/*----------------------------------------------------------------------*/
/*!
  \brief The table with all information about the elements.

  Major and minor number specify an element type. Here we have the
  table that contains everything we'd might need to know about the
  element.

  \warning This table is accessed by index. The entries must be in the
  same order as the major element number.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
ELEMENT_INFO element_info[el_count] = {
  { "none", {
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "shell1", {
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "shell8", {
      { 4, GiD_Hexahedra, 8, 0, 0 },
      { 8, GiD_Quadrilateral, 9, 0, 0 },
      { 9, GiD_Quadrilateral, 9, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "shell9", {
      { 4, GiD_Quadrilateral, 4, -1, 0 },
      { 4, GiD_Quadrilateral, 9, -1, 0 },
      { 8, GiD_Quadrilateral, 4, -1, 0 },
      { 8, GiD_Quadrilateral, 9, -1, 0 },
      { 9, GiD_Quadrilateral, 4, -1, 0 },
      { 9, GiD_Quadrilateral, 9, -1, 0 }
    } },
  { "brick1", {
      {  8, GiD_Hexahedra,  8, 0, 6 },
      { 20, GiD_Hexahedra, 27, 0, 6 },
      { 27, GiD_Hexahedra, 27, 0, 6 },
      {  0, GiD_NoElement,  0, 0, 0 },
      {  0, GiD_NoElement,  0, 0, 0 },
      {  0, GiD_NoElement,  0, 0, 0 }
    } },
  { "wall1", {
      { 3, GiD_Triangle,      1, 1*4, 0 },
      { 4, GiD_Quadrilateral, 4, 4*6, 0 },
      { 8, GiD_Quadrilateral, 9, 9*6, 0 },
      { 9, GiD_Quadrilateral, 9, 9*6, 0 },
      { 8, GiD_Quadrilateral, 4, 4*6, 0 },
      { 0, GiD_NoElement,     0,   0, 0 }
    } },
  { "beam3", {
      { 2, GiD_Linear,    1, 1*6, 0 },
      { 2, GiD_Linear,    2, 2*6, 0 },
      { 3, GiD_Linear,    2, 2*6, 0 },
      { 3, GiD_Linear,    3, 3*6, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "fluid2", {
      { 4, GiD_Quadrilateral, 4, 0, 0 },
      { 8, GiD_Quadrilateral, 9, 0, 0 },
      { 9, GiD_Quadrilateral, 9, 0, 0 },
      { 3, GiD_Triangle,      4, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "fluid2_pro", {
      { 4, GiD_NoElement, 4, 0, 0 },
      { 8, GiD_NoElement, 9, 0, 0 },
      { 9, GiD_NoElement, 9, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "fluid2_tu", {
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "fluid3", {
      {  8, GiD_Hexahedra,  8, 0, 6*2 },
      { 20, GiD_Hexahedra, 27, 0, 6*2 },
      { 27, GiD_Hexahedra, 27, 0, 6*2 },
      {  0, GiD_NoElement,  0, 0, 0 },
      {  0, GiD_NoElement,  0, 0, 0 },
      {  0, GiD_NoElement,  0, 0, 0 }
    } },
  { "ale2", {
      { 4, GiD_Quadrilateral, 1, 0, 0 },
      { 4, GiD_Quadrilateral, 4, 0, 0 },
      { 3, GiD_Triangle, 2, 0, 0 },
      { 3, GiD_Triangle, 3, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "ale3", {
      { 8, GiD_Hexahedra, 1, 0, 0 },
      { 8, GiD_Hexahedra, 8, 0, 0 },
      { 4, GiD_Tetrahedra, 1, 0, 0 },
      { 4, GiD_Tetrahedra, 4, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "axishell", {
      { 2, GiD_Linear,    1, 5, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "interf", {
      { 4, GiD_Quadrilateral, 4, 2*5, 0 },
      { 8, GiD_Quadrilateral, 9, 3*5, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } },
  { "wallge", {
      { 4, GiD_Quadrilateral, 4, 4*6, 0 },
      { 8, GiD_Quadrilateral, 9, 9*6, 0 },
      { 9, GiD_Quadrilateral, 9, 9*6, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 },
      { 0, GiD_NoElement, 0, 0, 0 }
    } }
};
