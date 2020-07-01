/*---------------------------------------------------------------------*/
/*! \file

\brief used in direct divergence line integration

\level 3


*----------------------------------------------------------------------*/

#ifndef BASE_H
#define BASE_H

#include <vector>
#include <cmath>
#include <iostream>

#include "cut_clnwrapper.H"
#include "../linalg/linalg_fixedsizematrix.H"

/*!
\brief Returns the function that is to be integrated along the bounding lines of the irregular
volume
*/
double base_func_line_int(LINALG::Matrix<2, 1> pt, int inte_num, std::vector<double> alfa)
{
  double basef_line = 0.0;
  if (inte_num == 1)  // f(x,y,z) = 1
  {
    basef_line =
        alfa[0] * pt(0, 0) + alfa[1] * pt(0, 0) * pt(0, 0) * 0.5 + alfa[2] * pt(0, 0) * pt(1, 0);
    return basef_line;
  }
  if (inte_num == 2)  // f(x,y,z) = x
  {
    if (fabs(alfa[1]) < 0.0000001)
      basef_line = 0.5 * pow((alfa[0] + alfa[2] * pt(1, 0)), 2) * pt(0, 0);
    else
      basef_line = std::pow((alfa[0] + alfa[1] * pt(0, 0) + alfa[2] * pt(1, 0)), 3) / 6.0 / alfa[1];
    return basef_line;
  }
  if (inte_num == 3)  // f(x,y,z) = y
  {
    basef_line =
        (0.5 * alfa[0] + alfa[1] * pt(0, 0) / 3.0 + 0.5 * alfa[2] * pt(1, 0)) * pt(0, 0) * pt(0, 0);
    return basef_line;
  }
  if (inte_num == 4)  // f(x,y,z) = z
  {
    basef_line = (alfa[0] + 0.5 * alfa[1] * pt(0, 0) + alfa[2] * pt(1, 0)) * pt(0, 0) * pt(1, 0);
    return basef_line;
  }

  if (inte_num == 5)  // f(x,y,z) = x^2
  {
    if (fabs(alfa[1]) < 0.0000001)
      basef_line = std::pow((alfa[0] + alfa[2] * pt(1, 0)), 3) * pt(0, 0) / 3.0;
    else
      basef_line =
          std::pow((alfa[0] + alfa[1] * pt(0, 0) + alfa[2] * pt(1, 0)), 4) / 12.0 / alfa[1];
    return basef_line;
  }
  if (inte_num == 6)  // f(x,y,z) = xy
  {
    basef_line =
        6 * pow((alfa[2] * pt(0, 0) * pt(1, 0)), 2) +
        (8 * alfa[1] * alfa[2] * pow(pt(0, 0), 3) + 12 * alfa[0] * alfa[2] * pow(pt(0, 0), 2)) *
            pt(1, 0) +
        3 * pow(alfa[1], 2) * pow(pt(0, 0), 4) + 8 * alfa[0] * alfa[1] * pow(pt(0, 0), 3) +
        6 * pow((alfa[0] * pt(0, 0)), 2);
    return basef_line / 24.0;
  }
  if (inte_num == 7)  // f(x,y,z) = xz
  {
    basef_line =
        pt(0, 0) * pow((alfa[2] * pt(1, 0)), 2) +
        2 * alfa[0] * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0)) +
        alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) * pt(1, 0) +
        pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) / 3.0 + alfa[0] * alfa[0] * pt(0, 0);
    return basef_line * pt(1, 0) * 0.5;
  }
  if (inte_num == 8)  // f(x,y,z) = y^2
  {
    basef_line = std::pow(pt(0, 0), 3) * (4 * alfa[2] * pt(1, 0) + 4 * alfa[0]) +
                 3 * alfa[1] * pow(pt(0, 0), 4);
    return basef_line / 12.0;
  }
  if (inte_num == 9)  // f(x,y,z) = yz
  {
    basef_line = (pt(0, 0) * pt(0, 0) * (3 * alfa[2] * pt(1, 0) + 3 * alfa[0]) +
                     2 * alfa[1] * pow(pt(0, 0), 3)) *
                 pt(1, 0);
    return basef_line / 6.0;
  }
  if (inte_num == 10)  // f(x,y,z) = z^2
  {
    basef_line =
        std::pow(pt(1, 0), 2) *
        (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0) + alfa[0] * pt(0, 0));
    return basef_line;
  }

  if (inte_num == 11)  // f(x,y,z) = x^3
  {
    basef_line = pt(0, 0) * pow((alfa[2] * pt(1, 0)), 4) +
                 4 * alfa[0] *
                     (pt(0, 0) * pow((alfa[2] * pt(1, 0)), 3) +
                         1.5 * alfa[1] * pow((alfa[2] * pt(0, 0) * pt(1, 0)), 2) +
                         alfa[1] * alfa[1] * alfa[2] * pow((pt(0, 0)), 3) * pt(1, 0) +
                         pow((alfa[1] * pt(0, 0)), 3) * 0.25 * pt(0, 0)) +
                 2 * alfa[1] * pt(0, 0) * pt(0, 0) * pow((alfa[2] * pt(1, 0)), 3) +
                 6 * alfa[0] * alfa[0] *
                     (pow((alfa[2] * pt(1, 0)), 2) * pt(0, 0) +
                         alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) * pt(1, 0) +
                         pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) / 3.0) +
                 2 * pow((alfa[1] * alfa[2] * pt(0, 0) * pt(1, 0)), 2) * pt(0, 0) +
                 4 * pow((alfa[0]), 3) *
                     (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0)) +
                 pow((alfa[1]), 3) * alfa[2] * pow((pt(0, 0)), 4) * pt(1, 0) +
                 pow((alfa[1] * pt(0, 0)), 4) * pt(0, 0) * 0.2 + pow(alfa[0], 4) * pt(0, 0);
    return basef_line * 0.25;
  }
  if (inte_num == 12)  // f(x,y,z) = x^2 y
  {
    basef_line =
        pt(0, 0) * pt(0, 0) *
            (10 * pow((alfa[2] * pt(1, 0)), 3) + 30 * alfa[0] * pow((alfa[2] * pt(1, 0)), 2) +
                30 * alfa[0] * alfa[0] * alfa[2] * pt(1, 0) + 10 * pow((alfa[0]), 3)) +
        pow((pt(0, 0)), 3) *
            (20 * alfa[1] * pow((alfa[2] * pt(1, 0)), 2) +
                40 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) + 20 * alfa[0] * alfa[0] * alfa[1]) +
        pow((pt(0, 0)), 4) *
            (15 * alfa[1] * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1] * alfa[1]) +
        4 * pow(alfa[1], 3) * pow(pt(0, 0), 5);
    return basef_line / 60.0;
  }
  if (inte_num == 13)  // f(x,y,z) = x^2 z
  {
    basef_line =
        pt(1, 0) * (pt(0, 0) * pow((alfa[2] * pt(1, 0)), 3) +
                       3 * alfa[0] *
                           (pt(0, 0) * pow((alfa[2] * pt(1, 0)), 2) +
                               alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) * pt(1, 0) +
                               pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) / 3.0) +
                       1.5 * alfa[1] * pow((alfa[2] * pt(0, 0) * pt(1, 0)), 2) +
                       3 * alfa[0] * alfa[0] *
                           (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0)) +
                       alfa[2] * pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) * pt(1, 0) +
                       0.25 * pow(alfa[1], 3) * pow(pt(0, 0), 4) + pow(alfa[0], 3) * pt(0, 0));
    return basef_line / 3.0;
  }
  if (inte_num == 14)  // f(x,y,z) = xy^2
  {
    basef_line =
        std::pow(pt(0, 0), 3) * (10 * pow((alfa[2] * pt(1, 0)), 2) +
                                    20 * alfa[0] * alfa[2] * pt(1, 0) + 10 * alfa[0] * alfa[0]) +
        pow(pt(0, 0), 4) * (15 * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1]) +
        6 * alfa[1] * alfa[1] * pow(pt(0, 0), 5);
    return basef_line / 60.0;
  }
  if (inte_num == 15)  // f(x,y,z) = xyz
  {
    basef_line =
        pt(1, 0) * (pt(0, 0) * pt(0, 0) *
                           (6 * pow((alfa[2] * pt(1, 0)), 2) + 12 * alfa[0] * alfa[2] * pt(1, 0) +
                               6 * alfa[0] * alfa[0]) +
                       pow(pt(0, 0), 3) * 8 * (alfa[1] * alfa[2] * pt(1, 0) + alfa[0] * alfa[1]) +
                       3 * alfa[1] * alfa[1] * pow(pt(0, 0), 4));
    return basef_line / 24.0;
  }
  if (inte_num == 16)  // f(x,y,z) = xz^2
  {
    basef_line =
        3 * alfa[2] * alfa[2] * pt(0, 0) * pow(pt(1, 0), 4) +
        pow(pt(1, 0), 3) *
            (3 * alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) + 6 * alfa[0] * alfa[2] * pt(0, 0)) +
        pt(1, 0) * pt(1, 0) *
            (alfa[1] * alfa[1] * pow(pt(0, 0), 3) + 3 * alfa[0] * alfa[1] * pt(0, 0) * pt(0, 0) +
                3 * alfa[0] * alfa[0] * pt(0, 0));
    return basef_line / 6.0;
  }
  if (inte_num == 17)  // f(x,y,z) = y^3
  {
    basef_line =
        std::pow(pt(0, 0), 4) * 5 * (alfa[2] * pt(1, 0) + alfa[0]) + 4 * alfa[1] * pow(pt(0, 0), 5);
    return basef_line * 0.05;
  }
  if (inte_num == 18)  // f(x,y,z) = y^2 z
  {
    basef_line = 4 * alfa[2] * pow((pt(0, 0) * pt(1, 0)), 2) * pt(0, 0) +
                 pt(1, 0) * (3 * alfa[1] * pow(pt(0, 0), 4) + 4 * alfa[0] * pow(pt(0, 0), 3));
    return basef_line / 12.0;
  }
  if (inte_num == 19)  // f(x,y,z) = yz^2
  {
    basef_line = pt(1, 0) * pt(1, 0) *
                 (pt(0, 0) * pt(0, 0) * (3 * alfa[2] * pt(1, 0) + 3 * alfa[0]) +
                     2 * alfa[1] * pow(pt(0, 0), 3));
    return basef_line / 6.0;
  }
  if (inte_num == 20)  // f(x,y,z) = z^3
  {
    basef_line =
        std::pow(pt(1, 0), 3) *
        (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0) + alfa[0] * pt(0, 0));
    return basef_line;
  }

  double a02 = pow(alfa[0], 2), a12 = pow(alfa[1], 2), a22 = pow(alfa[2], 2), y2 = pow(pt(0, 0), 2),
         z2 = pow(pt(1, 0), 2);
  double a03 = pow(alfa[0], 3), a13 = pow(alfa[1], 3), a23 = pow(alfa[2], 3), y3 = pow(pt(0, 0), 3),
         z3 = pow(pt(1, 0), 3);
  double a04 = pow(alfa[0], 4), a14 = pow(alfa[1], 4), a24 = pow(alfa[2], 4), y4 = pow(pt(0, 0), 4),
         z4 = pow(pt(1, 0), 4);
  double a05 = pow(alfa[0], 5), a15 = pow(alfa[1], 5), a25 = pow(alfa[2], 5), y5 = pow(pt(0, 0), 5),
         z5 = pow(pt(1, 0), 5);
  double y6 = pow(pt(0, 0), 6);
  if (inte_num == 21)  // f(x,y,z) = x^4
  {
    basef_line =
        6 * a25 * pt(0, 0) * z5 + (15 * alfa[1] * y2 + 30 * alfa[0] * pt(0, 0)) * a24 * z4 +
        (20 * a12 * y3 + 60 * alfa[0] * alfa[1] * y2 + 60 * a02 * pt(0, 0)) * a23 * z3 +
        (15 * a13 * y4 + 60 * alfa[0] * a12 * y3 + 90 * a02 * alfa[1] * y2 + 60 * a03 * pt(0, 0)) *
            a22 * z2 +
        (6 * a14 * y5 + 30 * alfa[0] * a13 * y4 + 60 * a02 * a12 * y3 + 60 * a03 * alfa[1] * y2 +
            30 * a04 * pt(0, 0)) *
            alfa[2] * pt(1, 0) +
        a15 * y6 + 6 * alfa[0] * a14 * y5 + 15 * a02 * a13 * y4 + 20 * a03 * a12 * y3 +
        15 * a04 * alfa[1] * y2 + 6 * a05 * pt(0, 0);
    return basef_line / 30.0;
  }
  if (inte_num == 22)  // f(x,y,z) = x^3*y
  {
    basef_line =
        y2 * (15 * a24 * z4 + 60 * alfa[0] * a23 * z3 + 90 * a02 * a22 * z2 +
                 60 * a03 * alfa[2] * pt(1, 0) + 15 * a04) +
        y3 * (40 * alfa[1] * a23 * z3 + 120 * alfa[0] * alfa[1] * a22 * z2 +
                 120 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 40 * a03 * alfa[1]) +
        y4 * (45 * a12 * a22 * z2 + 90 * alfa[0] * a12 * alfa[2] * pt(1, 0) + 45 * a02 * a12) +
        y5 * (24 * a13 * alfa[2] * pt(1, 0) + 24 * alfa[0] * a13) + 5 * a14 * y6;
    return basef_line / 120.0;
  }
  if (inte_num == 23)  // f(x,y,z) = x^3*z
  {
    basef_line =
        pt(1, 0) *
        (a24 * pt(0, 0) * z4 +
            4 * alfa[0] *
                (a23 * pt(0, 0) * z3 + 1.5 * alfa[1] * a22 * y2 * z2 +
                    a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4) +
            2 * alfa[1] * a23 * y2 * z3 +
            6 * a02 * (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
            2 * a12 * a22 * y3 * z2 +
            4 * a03 * (alfa[2] * pt(0, 0) * pt(1, 0) + alfa[1] * y2 * 0.5) +
            a13 * alfa[2] * y4 * pt(1, 0) + 0.2 * a14 * y5 + a04 * pt(0, 0));
    return 0.25 * basef_line;
  }
  if (inte_num == 24)  // f(x,y,z) = x^2*y^2
  {
    basef_line =
        y3 * (20 * a23 * z3 + 60 * alfa[0] * a22 * z2 + 60 * a02 * alfa[2] * pt(1, 0) + 20 * a03) +
        y4 * (45 * alfa[1] * a22 * z2 + 90 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                 45 * a02 * alfa[1]) +
        y5 * (36 * a12 * alfa[2] * pt(1, 0) + 36 * alfa[0] * a12) + 10 * a13 * y6;
    return basef_line / 180.0;
  }
  if (inte_num == 25)  // f(x,y,z) = x^2*yz
  {
    basef_line =
        pt(1, 0) *
        (y2 * (10 * a23 * z3 + 30 * alfa[0] * a22 * z2 + 30 * a02 * alfa[2] * pt(1, 0) + 10 * a03) +
            y3 * (20 * alfa[1] * a22 * z2 + 40 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                     20 * a02 * alfa[1]) +
            y4 * (15 * a12 * alfa[2] * pt(1, 0) + 15 * alfa[0] * a12) + 4 * a13 * y5);
    return basef_line / 60.0;
  }
  if (inte_num == 26)  // f(x,y,z) = x^2*z^2
  {
    basef_line =
        z2 * (a23 * pt(0, 0) * z3 +
                 3 * alfa[0] *
                     (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
                 1.5 * alfa[1] * a22 * y2 * z2 +
                 3 * a02 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                 a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4 + a03 * pt(0, 0));
    return basef_line / 3.0;
  }
  if (inte_num == 27)  // f(x,y,z) = x*y^3
  {
    basef_line = y4 * (15 * a22 * z2 + 30 * alfa[0] * alfa[2] * pt(1, 0) + 15 * a02) +
                 y5 * (24 * alfa[1] * alfa[2] * pt(1, 0) + 24 * alfa[0] * alfa[1]) + 10 * a12 * y6;
    return basef_line / 120.0;
  }
  if (inte_num == 28)  // f(x,y,z) = x*y^2*z
  {
    basef_line = pt(1, 0) * (y3 * (10 * a22 * z2 + 20 * alfa[0] * alfa[2] * pt(1, 0) + 10 * a02) +
                                y4 * (15 * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1]) +
                                6 * a12 * y5);
    return basef_line / 60.0;
  }
  if (inte_num == 29)  // f(x,y,z) = x*y*z^2
  {
    basef_line =
        z2 * (y2 * (6 * a22 * z2 + 12 * alfa[0] * alfa[2] * pt(1, 0) + 6 * a02) +
                 y3 * (8 * alfa[1] * alfa[2] * pt(1, 0) + 8 * alfa[0] * alfa[1]) + 3 * a12 * y4);
    return basef_line / 24.0;
  }
  if (inte_num == 30)  // f(x,y,z) = x*z^3
  {
    basef_line = z3 * (a22 * pt(0, 0) * z2 +
                          2 * alfa[0] * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                          alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0 + a02 * pt(0, 0));
    return basef_line * 0.5;
  }
  if (inte_num == 31)  // f(x,y,z) = y^4
  {
    basef_line = y5 * (6 * alfa[2] * pt(1, 0) + 6 * alfa[0]) + 5 * alfa[1] * y6;
    return basef_line / 30.0;
  }
  if (inte_num == 32)  // f(x,y,z) = y^3*z
  {
    basef_line = y4 * (5 * alfa[2] * pt(1, 0) + 5 * alfa[0]) + 4 * alfa[1] * y5;
    return basef_line * 0.05 * pt(1, 0);
  }
  if (inte_num == 33)  // f(x,y,z) = y^2*z^2
  {
    basef_line = 4 * y3 * (alfa[2] * pt(1, 0) + alfa[0]) + 3 * alfa[1] * y4;
    return basef_line / 12.0 * z2;
  }
  if (inte_num == 34)  // f(x,y,z) = y*z^3
  {
    basef_line = y2 * (3 * alfa[2] * pt(1, 0) + 3 * alfa[0]) + 2 * alfa[1] * y3;
    return basef_line * z3 / 6.0;
  }
  if (inte_num == 35)  // f(x,y,z) = z^4
  {
    basef_line = alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2 + alfa[0] * pt(0, 0);
    return basef_line * z4;
  }


  double a06 = pow(alfa[0], 6), a16 = pow(alfa[1], 6), a26 = pow(alfa[2], 6), z6 = pow(pt(1, 0), 6),
         y7 = pow(pt(0, 0), 7);
  if (inte_num == 36)  // f(x,y,z) = x^5
  {
    if (fabs(alfa[1]) < 0.0000001)
      basef_line = std::pow((alfa[0] + alfa[2] * pt(1, 0)), 6) * pt(0, 0) / 6.0;
    else
      basef_line =
          std::pow((alfa[0] + alfa[1] * pt(0, 0) + alfa[2] * pt(1, 0)), 7) / 42.0 / alfa[1];
    return basef_line;
  }
  if (inte_num == 37)  // f(x,y,z) = x^4*y
  {
    basef_line =
        pt(0, 0) * pt(0, 0) *
            (21 * a25 * z5 + 105 * alfa[0] * a24 * z4 + 210 * a02 * a23 * z3 +
                210 * a03 * a22 * z2 + 105 * a04 * alfa[2] * pt(1, 0) + 21 * a05) +
        pow(pt(0, 0), 3) * (70 * alfa[1] * a24 * z4 + 280 * alfa[0] * alfa[1] * a23 * z3 +
                               420 * a02 * alfa[1] * a22 * z2 +
                               280 * a03 * alfa[1] * alfa[2] * pt(1, 0) + 70 * a04 * alfa[1]) +
        pow(pt(0, 0), 4) * (105 * a12 * a23 * z3 + 315 * alfa[0] * a12 * a22 * z2 +
                               315 * a02 * a12 * alfa[2] * pt(1, 0) + 105 * a03 * a12) +
        pow(pt(0, 0), 5) *
            (84 * a13 * a22 * z2 + 168 * alfa[0] * a13 * alfa[2] * pt(1, 0) + 84 * a02 * a13) +
        pow(pt(0, 0), 6) * (35 * a14 * alfa[2] * pt(1, 0) + 35 * alfa[0] * a14) + 6 * a15 * y7;
    return basef_line / 210.0;
  }
  if (inte_num == 38)  // f(x,y,z) = x^4*z
  {
    basef_line =
        pt(1, 0) *
        (a25 * pt(0, 0) * z5 +
            5 * alfa[0] *
                (a24 * pt(0, 0) * z4 + 2 * alfa[1] * a23 * y2 * z3 + 2 * a12 * a22 * y3 * z2 +
                    a13 * alfa[2] * y4 * pt(1, 0) + 0.2 * a14 * y5) +
            2.5 * alfa[1] * a24 * y2 * z4 +
            10 * a02 *
                (a23 * pt(0, 0) * z3 + 1.5 * alfa[1] * a22 * y2 * z2 +
                    a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4) +
            10 * a12 * a23 * y3 * z3 / 3.0 +
            10 * a03 * (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
            2.5 * a13 * a22 * y4 * z2 +
            5 * a04 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
            a14 * alfa[2] * y5 * pt(1, 0) + a15 * y6 / 6.0 + a05 * pt(0, 0));
    return basef_line * 0.2;
  }
  if (inte_num == 39)  // f(x,y,z) = x^3*y^2
  {
    basef_line =
        std::pow(pt(0, 0), 3) * (35 * a24 * z4 + 140 * alfa[0] * a23 * z3 + 210 * a02 * a22 * z2 +
                                    140 * a03 * alfa[2] * pt(1, 0) + 35 * a04) +
        pow(pt(0, 0), 4) * (105 * alfa[1] * a23 * z3 + 315 * alfa[0] * alfa[1] * a22 * z2 +
                               315 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 105 * a03 * alfa[1]) +
        pow(pt(0, 0), 5) *
            (126 * a12 * a22 * z2 + 252 * alfa[0] * a12 * alfa[2] * pt(1, 0) + 126 * a02 * a12) +
        pow(pt(0, 0), 6) * (70 * a13 * alfa[2] * pt(1, 0) + 70 * alfa[0] * a13) + 15 * a14 * y7;
    return basef_line / 420.0;
  }
  if (inte_num == 40)  // f(x,y,z) = x^3*yz
  {
    basef_line =
        pt(1, 0) *
        (y2 * (15 * a24 * z4 + 60 * alfa[0] * a23 * z3 + 90 * a02 * a22 * z2 +
                  60 * a03 * alfa[2] * pt(1, 0) + 15 * a04) +
            y3 * (40 * alfa[1] * a23 * z3 + 120 * alfa[0] * alfa[1] * a22 * z2 +
                     120 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 40 * a03 * alfa[1]) +
            y4 * (45 * a12 * a22 * z2 + 90 * alfa[0] * a12 * alfa[2] * pt(1, 0) + 45 * a02 * a12) +
            y5 * (24 * a13 * alfa[2] * pt(1, 0) + 24 * alfa[0] * a13) + 5 * a14 * y6);
    return basef_line / 120.0;
  }
  if (inte_num == 41)  // f(x,y,z) = x^3*z^2
  {
    basef_line =
        z2 *
        (a24 * pt(0, 0) * z4 +
            4 * alfa[0] *
                (a23 * pt(0, 0) * z3 + 1.5 * alfa[1] * a22 * y2 * z2 +
                    a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4) +
            2 * alfa[1] * a23 * y2 * z3 +
            6 * a02 * (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
            2 * a12 * a22 * y3 * z2 +
            4 * a03 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
            a13 * alfa[2] * y4 * pt(1, 0) + 0.2 * a14 * y5 + a04 * pt(0, 0));
    return basef_line * 0.25;
  }
  if (inte_num == 42)  // f(x,y,z) = x^2*y^3
  {
    basef_line = y4 * (35 * a23 * z3 + 105 * alfa[0] * a22 * z2 + 105 * a02 * alfa[2] * pt(1, 0) +
                          35 * a03) +
                 y5 * (84 * alfa[1] * a22 * z2 + 168 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                          84 * a02 * alfa[1]) +
                 y6 * (70 * a12 * alfa[2] * pt(1, 0) + 70 * alfa[0] * a12) + 20 * a13 * y7;
    return basef_line / 420.0;
  }
  if (inte_num == 43)  // f(x,y,z) = x^2*y^2*z
  {
    basef_line =
        pt(1, 0) *
        (y3 * (20 * a23 * z3 + 60 * alfa[0] * a22 * z2 + 60 * a02 * alfa[2] * pt(1, 0) + 20 * a03) +
            y4 * (45 * alfa[1] * a22 * z2 + 90 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                     45 * a02 * alfa[1]) +
            y5 * (36 * a12 * alfa[2] * pt(1, 0) + 36 * alfa[0] * a12) + 10 * a13 * y6);
    return basef_line / 180.0;
  }
  if (inte_num == 44)  // f(x,y,z) = x^2*yz^2
  {
    basef_line =
        z2 *
        (y2 * (10 * a23 * z3 + 30 * alfa[0] * a22 * z2 + 30 * a02 * alfa[2] * pt(1, 0) + 10 * a03) +
            y3 * (20 * alfa[1] * a22 * z2 + 40 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                     20 * a02 * alfa[1]) +
            y4 * (15 * a12 * alfa[2] * pt(1, 0) + 15 * alfa[0] * a12) + 4 * a13 * y5);
    return basef_line / 60.0;
  }
  if (inte_num == 45)  // f(x,y,z) = x^2*z^3
  {
    basef_line =
        z3 * (a23 * pt(0, 0) * z3 +
                 3 * alfa[0] *
                     (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
                 1.5 * alfa[1] * a22 * y2 * z2 +
                 3 * a02 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                 a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4 + a03 * pt(0, 0));
    return basef_line / 3.0;
  }
  if (inte_num == 46)  // f(x,y,z) = xy^4
  {
    basef_line = y5 * (21 * a22 * z2 + 42 * alfa[0] * alfa[2] * pt(1, 0) + 21 * a02) +
                 y6 * (35 * alfa[1] * alfa[2] * pt(1, 0) + 35 * alfa[0] * alfa[1]) + 15 * a12 * y7;
    return basef_line / 210.0;
  }
  if (inte_num == 47)  // f(x,y,z) = xy^3*z
  {
    basef_line = pt(1, 0) * (y4 * (15 * a22 * z2 + 30 * alfa[0] * alfa[2] * pt(1, 0) + 15 * a02) +
                                y5 * (24 * alfa[1] * alfa[2] * pt(1, 0) + 24 * alfa[0] * alfa[1]) +
                                10 * a12 * y6);
    return basef_line / 120.0;
  }
  if (inte_num == 48)  // f(x,y,z) = xy^2*z^2
  {
    basef_line =
        z2 * (y3 * (10 * a22 * z2 + 20 * alfa[0] * alfa[2] * pt(1, 0) + 10 * a02) +
                 y4 * (15 * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1]) + 6 * a12 * y5);
    return basef_line / 60.0;
  }
  if (inte_num == 49)  // f(x,y,z) = xyz^3
  {
    basef_line =
        z3 * (y2 * (6 * a22 * z2 + 12 * alfa[0] * alfa[2] * pt(1, 0) + 6 * a02) +
                 y3 * (8 * alfa[1] * alfa[2] * pt(1, 0) + 8 * alfa[0] * alfa[1]) + 3 * a12 * y4);
    return basef_line / 24.0;
  }
  if (inte_num == 50)  // f(x,y,z) = xz^4
  {
    basef_line = z4 * (a22 * pt(0, 0) * z2 +
                          2 * alfa[0] * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                          alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0 + a02 * pt(0, 0));
    return basef_line * 0.5;
  }
  if (inte_num == 51)  // f(x,y,z) = y^5
  {
    basef_line = y6 * (7 * alfa[2] * pt(1, 0) + 7 * alfa[0]) + 6 * alfa[1] * y7;
    return basef_line / 42.0;
  }
  if (inte_num == 52)  // f(x,y,z) = y^4*z
  {
    basef_line = pt(1, 0) * (y5 * (6 * alfa[2] * pt(1, 0) + 6 * alfa[0]) + 5 * alfa[1] * y6);
    return basef_line / 30.0;
  }
  if (inte_num == 53)  // f(x,y,z) = y^3*z^2
  {
    basef_line = z2 * (y4 * (5 * alfa[2] * pt(1, 0) + 5 * alfa[0]) + 4 * alfa[1] * y5);
    return basef_line * 0.05;
  }
  if (inte_num == 54)  // f(x,y,z) = y^2*z^3
  {
    basef_line = z3 * (y3 * (4 * alfa[2] * pt(1, 0) + 4 * alfa[0]) + 3 * alfa[1] * y4);
    return basef_line / 12.0;
  }
  if (inte_num == 55)  // f(x,y,z) = yz^4
  {
    basef_line = z4 * (y2 * (3 * alfa[2] * pt(1, 0) + 3 * alfa[0]) + 2 * alfa[1] * y3);
    return basef_line / 6.0;
  }
  if (inte_num == 56)  // f(x,y,z) = z^5
  {
    basef_line = z5 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2 + alfa[0] * pt(0, 0));
    return basef_line;
  }

  double y8 = pow(pt(0, 0), 8);
  //      double
  //      a07=pow(alfa[0],7),a17=pow(alfa[1],7),a27=pow(alfa[2],7),z7=pow(pt(1,0),7),y8=pow(pt(0,0),8);
  if (inte_num == 57)  // f(x,y,z) = x^6
  {
    if (fabs(alfa[1]) < 0.0000001)
      basef_line = std::pow((alfa[0] + alfa[2] * pt(1, 0)), 7) * pt(0, 0) / 7.0;
    else
      basef_line =
          std::pow((alfa[0] + alfa[1] * pt(0, 0) + alfa[2] * pt(1, 0)), 8) / 56.0 / alfa[1];
    return basef_line;
  }
  if (inte_num == 58)  // f(x,y,z) = x^5*y
  {
    basef_line = y2 * (28.0 * a26 * z6 + 168.0 * alfa[0] * a25 * z5 + 420.0 * a02 * a24 * z4 +
                          560.0 * a03 * a23 * z3 + 420.0 * a04 * a22 * z2 +
                          168.0 * a05 * alfa[2] * pt(1, 0) + 28.0 * a06) +
                 y3 * (112.0 * alfa[1] * a25 * z5 + 560.0 * alfa[0] * alfa[1] * a24 * z4 +
                          1120.0 * a02 * alfa[1] * a23 * z3 + 1120.0 * a03 * alfa[1] * a22 * z2 +
                          560.0 * a04 * alfa[1] * alfa[2] * pt(1, 0) + 112.0 * a05 * alfa[1]) +
                 y4 * (210.0 * a12 * a24 * z4 + 840.0 * alfa[0] * a12 * a23 * z3 +
                          1260.0 * a02 * a12 * a22 * z2 + 840.0 * a03 * a12 * alfa[2] * pt(1, 0) +
                          210.0 * a04 * a12) +
                 y5 * (224.0 * a13 * a23 * z3 + 672.0 * alfa[0] * a13 * a22 * z2 +
                          672.0 * a02 * a13 * alfa[2] * pt(1, 0) + 224.0 * a03 * a13) +
                 y6 * (140.0 * a14 * a22 * z2 + 280.0 * alfa[0] * a14 * alfa[2] * pt(1, 0) +
                          140.0 * a02 * a14) +
                 y7 * (48.0 * a15 * alfa[2] * pt(1, 0) + 48.0 * alfa[0] * a15) + 7.0 * a16 * y8;

    return basef_line / 336.0;
  }
  if (inte_num == 59)  // f(x,y,z) = x^5*z
  {
    if (fabs(alfa[1]) < 0.0000001)
      basef_line = pt(0, 0) * pt(1, 0) * pow((alfa[2] * pt(1, 0) + alfa[0]), 6) / 6.0;
    else
      basef_line =
          pt(1, 0) * pow((alfa[2] * pt(1, 0) + alfa[1] * pt(0, 0) + alfa[0]), 7) / 42.0 / alfa[1];
    return basef_line;
  }
  if (inte_num == 60)  // f(x,y,z) = x^4*y^2
  {
    basef_line = y3 * (56.0 * a25 * z5 + 280.0 * alfa[0] * a24 * z4 + 560.0 * a02 * a23 * z3 +
                          560.0 * a03 * a22 * z2 + 280.0 * a04 * alfa[2] * pt(1, 0) + 56.0 * a05) +
                 y4 * (210.0 * alfa[1] * a24 * z4 + 840.0 * alfa[0] * alfa[1] * a23 * z3 +
                          1260.0 * a02 * alfa[1] * a22 * z2 +
                          840.0 * a03 * alfa[1] * alfa[2] * pt(1, 0) + 210.0 * a04 * alfa[1]) +
                 y5 * (336.0 * a12 * a23 * z3 + 1008.0 * alfa[0] * a12 * a22 * z2 +
                          1008.0 * a02 * a12 * alfa[2] * pt(1, 0) + 336.0 * a03 * a12) +
                 y6 * (280.0 * a13 * a22 * z2 + 560.0 * alfa[0] * a13 * alfa[2] * pt(1, 0) +
                          280.0 * a02 * a13) +
                 y7 * (120.0 * a14 * alfa[2] * pt(1, 0) + 120.0 * alfa[0] * a14) + 21.0 * a15 * y8;
    return basef_line / 840.0;
  }
  if (inte_num == 61)  // f(x,y,z) = x^4*yz
  {
    basef_line =
        pt(1, 0) *
        (y2 * (21.0 * a25 * z5 + 105.0 * alfa[0] * a24 * z4 + 210.0 * a02 * a23 * z3 +
                  210.0 * a03 * a22 * z2 + 105.0 * a04 * alfa[2] * pt(1, 0) + 21.0 * a05) +
            y3 * (70.0 * alfa[1] * a24 * z4 + 280.0 * alfa[0] * alfa[1] * a23 * z3 +
                     420.0 * a02 * alfa[1] * a22 * z2 + 280.0 * a03 * alfa[1] * alfa[2] * pt(1, 0) +
                     70.0 * a04 * alfa[1]) +
            y4 * (105.0 * a12 * a23 * z3 + 315.0 * alfa[0] * a12 * a22 * z2 +
                     315.0 * a02 * a12 * alfa[2] * pt(1, 0) + 105.0 * a03 * a12) +
            y5 * (84.0 * a13 * a22 * z2 + 168.0 * alfa[0] * a13 * alfa[2] * pt(1, 0) +
                     84.0 * a02 * a13) +
            y6 * (35.0 * a14 * alfa[2] * pt(1, 0) + 35.0 * alfa[0] * a14) + 6.0 * a15 * y7);
    return basef_line / 210.0;
  }
  if (inte_num == 62)  // f(x,y,z) = x^4*z^2
  {
    basef_line =
        z2 *
        (a25 * pt(0, 0) * z5 +
            5.0 * alfa[0] *
                (a24 * pt(0, 0) * z4 + 2.0 * alfa[1] * a23 * y2 * z3 + 2.0 * a12 * a22 * y3 * z2 +
                    a13 * alfa[2] * y4 * pt(1, 0) + a14 * y5 * 0.2) +
            2.5 * alfa[1] * a24 * y2 * z4 +
            10.0 * a02 *
                (a23 * pt(0, 0) * z3 + 1.5 * alfa[1] * a22 * y2 * z2 +
                    a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4) +
            10.0 / 3.0 * a12 * a23 * y3 * z3 +
            10.0 * a03 *
                (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
            2.5 * a13 * a22 * y4 * z2 +
            5.0 * a04 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
            a14 * alfa[2] * y5 * pt(1, 0) + a15 * y6 / 6.0 + a05 * pt(0, 0));
    return basef_line * 0.2;
  }
  if (inte_num == 63)  // f(x,y,z) = x^3*y^3
  {
    basef_line = y4 * (70.0 * a24 * z4 + 280.0 * alfa[0] * a23 * z3 + 420.0 * a02 * a22 * z2 +
                          280.0 * a03 * alfa[2] * pt(1, 0) + 70.0 * a04) +
                 y5 * (224.0 * alfa[1] * a23 * z3 + 672.0 * alfa[0] * alfa[1] * a22 * z2 +
                          672.0 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 224.0 * a03 * alfa[1]) +
                 y6 * (280.0 * a12 * a22 * z2 + 560.0 * alfa[0] * a12 * alfa[2] * pt(1, 0) +
                          280.0 * a02 * a12) +
                 y7 * (160.0 * a13 * alfa[2] * pt(1, 0) + 160.0 * alfa[0] * a13) + 35.0 * a14 * y8;
    return basef_line / 1120.0;
  }
  if (inte_num == 64)  // f(x,y,z) = x^3*y^2*z
  {
    basef_line =
        pt(1, 0) *
        (y3 * (35.0 * a24 * z4 + 140.0 * alfa[0] * a23 * z3 + 210.0 * a02 * a22 * z2 +
                  140.0 * a03 * alfa[2] * pt(1, 0) + 35.0 * a04) +
            y4 * (105.0 * alfa[1] * a23 * z3 + 315.0 * alfa[0] * alfa[1] * a22 * z2 +
                     315.0 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 105.0 * a03 * alfa[1]) +
            y5 * (126.0 * a12 * a22 * z2 + 252.0 * alfa[0] * a12 * alfa[2] * pt(1, 0) +
                     126.0 * a02 * a12) +
            y6 * (70.0 * a13 * alfa[2] * pt(1, 0) + 70.0 * alfa[0] * a13) + 15.0 * a14 * y7);
    return basef_line / 420.0;
  }
  if (inte_num == 65)  // f(x,y,z) = x^3*yz^2
  {
    basef_line =
        z2 * (y2 * (15.0 * a24 * z4 + 60.0 * alfa[0] * a23 * z3 + 90.0 * a02 * a22 * z2 +
                       60.0 * a03 * alfa[2] * pt(1, 0) + 15.0 * a04) +
                 y3 * (40.0 * alfa[1] * a23 * z3 + 120.0 * alfa[0] * alfa[1] * a22 * z2 +
                          120.0 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 40.0 * a03 * alfa[1]) +
                 y4 * (45.0 * a12 * a22 * z2 + 90.0 * alfa[0] * a12 * alfa[2] * pt(1, 0) +
                          45.0 * a02 * a12) +
                 y5 * (24.0 * a13 * alfa[2] * pt(1, 0) + 24.0 * alfa[0] * a13) + 5.0 * a14 * y6);
    return basef_line / 120.0;
  }
  if (inte_num == 66)  // f(x,y,z) = x^3*z^3
  {
    basef_line =
        z3 *
        (a24 * pt(0, 0) * z4 +
            4.0 * alfa[0] *
                (a23 * pt(0, 0) * z3 + 1.5 * alfa[1] * a22 * y2 * z2 +
                    a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4) +
            2.0 * alfa[1] * a23 * y2 * z3 +
            6.0 * a02 * (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
            2.0 * a12 * a22 * y3 * z2 +
            4.0 * a03 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
            a13 * alfa[2] * y4 * pt(1, 0) + 0.2 * a14 * y5 + a04 * pt(0, 0));
    return basef_line * 0.25;
  }
  if (inte_num == 67)  // f(x,y,z) = x^2*y^4
  {
    basef_line = y5 * (56.0 * a23 * z3 + 168.0 * alfa[0] * a22 * z2 +
                          168.0 * a02 * alfa[2] * pt(1, 0) + 56.0 * a03) +
                 y6 * (140.0 * alfa[1] * a22 * z2 + 280.0 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                          140.0 * a02 * alfa[1]) +
                 y7 * (120.0 * a12 * alfa[2] * pt(1, 0) + 120.0 * alfa[0] * a12) + 35.0 * a13 * y8;
    return basef_line / 840.0;
  }
  if (inte_num == 68)  // f(x,y,z) = x^2*y^3*z
  {
    basef_line =
        pt(1, 0) *
        (y4 * (35.0 * a23 * z3 + 105.0 * alfa[0] * a22 * z2 + 105.0 * a02 * alfa[2] * pt(1, 0) +
                  35.0 * a03) +
            y5 * (84.0 * alfa[1] * a22 * z2 + 168.0 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                     84.0 * a02 * alfa[1]) +
            y6 * (70.0 * a12 * alfa[2] * pt(1, 0) + 70.0 * alfa[0] * a12) + 20.0 * a13 * y7);
    return basef_line / 420.0;
  }
  if (inte_num == 69)  // f(x,y,z) = x^2*y^2*z^2
  {
    basef_line =
        z2 * (y3 * (20.0 * a23 * z3 + 60.0 * alfa[0] * a22 * z2 + 60.0 * a02 * alfa[2] * pt(1, 0) +
                       20.0 * a03) +
                 y4 * (45.0 * alfa[1] * a22 * z2 + 90.0 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                          45.0 * a02 * alfa[1]) +
                 y5 * (36.0 * a12 * alfa[2] * pt(1, 0) + 36.0 * alfa[0] * a12) + 10.0 * a13 * y6);
    return basef_line / 180.0;
  }
  if (inte_num == 70)  // f(x,y,z) = x^2*yz^3
  {
    basef_line =
        z3 * (y2 * (10.0 * a23 * z3 + 30.0 * alfa[0] * a22 * z2 + 30.0 * a02 * alfa[2] * pt(1, 0) +
                       10.0 * a03) +
                 y3 * (20.0 * alfa[1] * a22 * z2 + 40.0 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                          20.0 * a02 * alfa[1]) +
                 y4 * (15.0 * a12 * alfa[2] * pt(1, 0) + 15.0 * alfa[0] * a12) + 4.0 * a13 * y5);
    return basef_line / 60.0;
  }
  if (inte_num == 71)  // f(x,y,z) = x^2*z^4
  {
    basef_line =
        z4 * (a23 * pt(0, 0) * z3 +
                 3.0 * alfa[0] *
                     (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
                 1.5 * alfa[1] * a22 * y2 * z2 +
                 3.0 * a02 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                 a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4 + a03 * pt(0, 0));
    return basef_line / 3.0;
  }
  if (inte_num == 72)  // f(x,y,z) = x*y^5
  {
    basef_line = y6 * (28.0 * a22 * z2 + 56.0 * alfa[0] * alfa[2] * pt(1, 0) + 28.0 * a02) +
                 y7 * (48.0 * alfa[1] * alfa[2] * pt(1, 0) + 48.0 * alfa[0] * alfa[1]) +
                 21.0 * a12 * y8;
    return basef_line / 336.0;
  }
  if (inte_num == 73)  // f(x,y,z) = x*y^4*z
  {
    basef_line =
        pt(1, 0) * (y5 * (21.0 * a22 * z2 + 42.0 * alfa[0] * alfa[2] * pt(1, 0) + 21.0 * a02) +
                       y6 * (35.0 * alfa[1] * alfa[2] * pt(1, 0) + 35.0 * alfa[0] * alfa[1]) +
                       15.0 * a12 * y7);
    return basef_line / 210.0;
  }
  if (inte_num == 74)  // f(x,y,z) = x*y^3*z^2
  {
    basef_line = z2 * (y4 * (15.0 * a22 * z2 + 30.0 * alfa[0] * alfa[2] * pt(1, 0) + 15.0 * a02) +
                          y5 * (24.0 * alfa[1] * alfa[2] * pt(1, 0) + 24.0 * alfa[0] * alfa[1]) +
                          10.0 * a12 * y6);
    return basef_line / 120.0;
  }
  if (inte_num == 75)  // f(x,y,z) = x*y^2*z^3
  {
    basef_line = z3 * (y3 * (10.0 * a22 * z2 + 20.0 * alfa[0] * alfa[2] * pt(1, 0) + 10.0 * a02) +
                          y4 * (15.0 * alfa[1] * alfa[2] * pt(1, 0) + 15.0 * alfa[0] * alfa[1]) +
                          6.0 * a12 * y5);
    return basef_line / 60.0;
  }
  if (inte_num == 76)  // f(x,y,z) = xyz^4
  {
    basef_line = z4 * (y2 * (6.0 * a22 * z2 + 12.0 * alfa[0] * alfa[2] * pt(1, 0) + 6.0 * a02) +
                          y3 * (8.0 * alfa[1] * alfa[2] * pt(1, 0) + 8.0 * alfa[0] * alfa[1]) +
                          3.0 * a12 * y4);
    return basef_line / 24.0;
  }
  if (inte_num == 77)  // f(x,y,z) = xz^5
  {
    basef_line = z5 * (a22 * pt(0, 0) * z2 +
                          2.0 * alfa[0] * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                          alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0 + a02 * pt(0, 0));
    return basef_line * 0.5;
  }
  if (inte_num == 78)  // f(x,y) = y^6
  {
    basef_line = y7 * (8.0 * alfa[2] * pt(1, 0) + 8.0 * alfa[0]) + 7.0 * alfa[1] * y8;
    return basef_line / 56.0;
  }
  if (inte_num == 79)  // f(x,y,z) = y^5*z
  {
    basef_line = pt(1, 0) * (y6 * (7.0 * alfa[2] * pt(1, 0) + 7.0 * alfa[0]) + 6.0 * alfa[1] * y7);
    return basef_line / 42.0;
  }
  if (inte_num == 80)  // f(x,y,z) = y^4*z^2
  {
    basef_line = z2 * (y5 * (6.0 * alfa[2] * pt(1, 0) + 6.0 * alfa[0]) + 5.0 * alfa[1] * y6);
    return basef_line / 30.0;
  }
  if (inte_num == 81)  // f(x,y,z) = y^3*z^3
  {
    basef_line = z3 * (y4 * (5.0 * alfa[2] * pt(1, 0) + 5.0 * alfa[0]) + 4.0 * alfa[1] * y5);
    return basef_line / 20.0;
  }
  if (inte_num == 82)  // f(x,y,z) = y^2*z^4
  {
    basef_line = z4 * (y3 * (4.0 * alfa[2] * pt(1, 0) + 4.0 * alfa[0]) + 3.0 * alfa[1] * y4);
    return basef_line / 12.0;
  }
  if (inte_num == 83)  // f(x,y,z) = yz^5
  {
    basef_line = z5 * (y2 * (3.0 * alfa[2] * pt(1, 0) + 3.0 * alfa[0]) + 2.0 * alfa[1] * y3);
    return basef_line / 6.0;
  }
  if (inte_num == 84)  // f(x,y,z) = z^6
  {
    basef_line = z6 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2 + alfa[0] * pt(0, 0));
    return basef_line;
  }

  dserror("the base function required is not defined");
  return 0.0;
}
#endif
