/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over a line using one-dimensional Gauss quadrature
equations

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_line_integration.hpp"

#include "4C_cut_tolerance.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
  \brief Returns the function that is to be integrated along the bounding lines of the irregular
  volume
  */
  double base_func_line_int(Core::LinAlg::Matrix<2, 1> pt, int inte_num, std::vector<double> alfa)
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
        basef_line = 0.5 * std::pow((alfa[0] + alfa[2] * pt(1, 0)), 2) * pt(0, 0);
      else
        basef_line =
            std::pow((alfa[0] + alfa[1] * pt(0, 0) + alfa[2] * pt(1, 0)), 3) / 6.0 / alfa[1];
      return basef_line;
    }
    if (inte_num == 3)  // f(x,y,z) = y
    {
      basef_line = (0.5 * alfa[0] + alfa[1] * pt(0, 0) / 3.0 + 0.5 * alfa[2] * pt(1, 0)) *
                   pt(0, 0) * pt(0, 0);
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
      basef_line = 6 * std::pow((alfa[2] * pt(0, 0) * pt(1, 0)), 2) +
                   (8 * alfa[1] * alfa[2] * std::pow(pt(0, 0), 3) +
                       12 * alfa[0] * alfa[2] * std::pow(pt(0, 0), 2)) *
                       pt(1, 0) +
                   3 * std::pow(alfa[1], 2) * std::pow(pt(0, 0), 4) +
                   8 * alfa[0] * alfa[1] * std::pow(pt(0, 0), 3) +
                   6 * std::pow((alfa[0] * pt(0, 0)), 2);
      return basef_line / 24.0;
    }
    if (inte_num == 7)  // f(x,y,z) = xz
    {
      basef_line =
          pt(0, 0) * std::pow((alfa[2] * pt(1, 0)), 2) +
          2 * alfa[0] * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0)) +
          alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) * pt(1, 0) +
          std::pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) / 3.0 + alfa[0] * alfa[0] * pt(0, 0);
      return basef_line * pt(1, 0) * 0.5;
    }
    if (inte_num == 8)  // f(x,y,z) = y^2
    {
      basef_line = std::pow(pt(0, 0), 3) * (4 * alfa[2] * pt(1, 0) + 4 * alfa[0]) +
                   3 * alfa[1] * std::pow(pt(0, 0), 4);
      return basef_line / 12.0;
    }
    if (inte_num == 9)  // f(x,y,z) = yz
    {
      basef_line = (pt(0, 0) * pt(0, 0) * (3 * alfa[2] * pt(1, 0) + 3 * alfa[0]) +
                       2 * alfa[1] * std::pow(pt(0, 0), 3)) *
                   pt(1, 0);
      return basef_line / 6.0;
    }
    if (inte_num == 10)  // f(x,y,z) = z^2
    {
      basef_line =
          std::pow(pt(1, 0), 2) * (alfa[2] * pt(0, 0) * pt(1, 0) +
                                      0.5 * alfa[1] * pt(0, 0) * pt(0, 0) + alfa[0] * pt(0, 0));
      return basef_line;
    }

    if (inte_num == 11)  // f(x,y,z) = x^3
    {
      basef_line = pt(0, 0) * std::pow((alfa[2] * pt(1, 0)), 4) +
                   4 * alfa[0] *
                       (pt(0, 0) * std::pow((alfa[2] * pt(1, 0)), 3) +
                           1.5 * alfa[1] * std::pow((alfa[2] * pt(0, 0) * pt(1, 0)), 2) +
                           alfa[1] * alfa[1] * alfa[2] * std::pow((pt(0, 0)), 3) * pt(1, 0) +
                           std::pow((alfa[1] * pt(0, 0)), 3) * 0.25 * pt(0, 0)) +
                   2 * alfa[1] * pt(0, 0) * pt(0, 0) * std::pow((alfa[2] * pt(1, 0)), 3) +
                   6 * alfa[0] * alfa[0] *
                       (std::pow((alfa[2] * pt(1, 0)), 2) * pt(0, 0) +
                           alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) * pt(1, 0) +
                           std::pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) / 3.0) +
                   2 * std::pow((alfa[1] * alfa[2] * pt(0, 0) * pt(1, 0)), 2) * pt(0, 0) +
                   4 * std::pow((alfa[0]), 3) *
                       (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0)) +
                   std::pow((alfa[1]), 3) * alfa[2] * std::pow((pt(0, 0)), 4) * pt(1, 0) +
                   std::pow((alfa[1] * pt(0, 0)), 4) * pt(0, 0) * 0.2 +
                   std::pow(alfa[0], 4) * pt(0, 0);
      return basef_line * 0.25;
    }
    if (inte_num == 12)  // f(x,y,z) = x^2 y
    {
      basef_line =
          pt(0, 0) * pt(0, 0) *
              (10 * std::pow((alfa[2] * pt(1, 0)), 3) +
                  30 * alfa[0] * std::pow((alfa[2] * pt(1, 0)), 2) +
                  30 * alfa[0] * alfa[0] * alfa[2] * pt(1, 0) + 10 * std::pow((alfa[0]), 3)) +
          std::pow((pt(0, 0)), 3) *
              (20 * alfa[1] * std::pow((alfa[2] * pt(1, 0)), 2) +
                  40 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) + 20 * alfa[0] * alfa[0] * alfa[1]) +
          std::pow((pt(0, 0)), 4) *
              (15 * alfa[1] * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1] * alfa[1]) +
          4 * std::pow(alfa[1], 3) * std::pow(pt(0, 0), 5);
      return basef_line / 60.0;
    }
    if (inte_num == 13)  // f(x,y,z) = x^2 z
    {
      basef_line =
          pt(1, 0) * (pt(0, 0) * std::pow((alfa[2] * pt(1, 0)), 3) +
                         3 * alfa[0] *
                             (pt(0, 0) * std::pow((alfa[2] * pt(1, 0)), 2) +
                                 alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) * pt(1, 0) +
                                 std::pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) / 3.0) +
                         1.5 * alfa[1] * std::pow((alfa[2] * pt(0, 0) * pt(1, 0)), 2) +
                         3 * alfa[0] * alfa[0] *
                             (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * pt(0, 0) * pt(0, 0)) +
                         alfa[2] * std::pow((alfa[1] * pt(0, 0)), 2) * pt(0, 0) * pt(1, 0) +
                         0.25 * std::pow(alfa[1], 3) * std::pow(pt(0, 0), 4) +
                         std::pow(alfa[0], 3) * pt(0, 0));
      return basef_line / 3.0;
    }
    if (inte_num == 14)  // f(x,y,z) = xy^2
    {
      basef_line =
          std::pow(pt(0, 0), 3) * (10 * std::pow((alfa[2] * pt(1, 0)), 2) +
                                      20 * alfa[0] * alfa[2] * pt(1, 0) + 10 * alfa[0] * alfa[0]) +
          std::pow(pt(0, 0), 4) * (15 * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1]) +
          6 * alfa[1] * alfa[1] * std::pow(pt(0, 0), 5);
      return basef_line / 60.0;
    }
    if (inte_num == 15)  // f(x,y,z) = xyz
    {
      basef_line =
          pt(1, 0) *
          (pt(0, 0) * pt(0, 0) *
                  (6 * std::pow((alfa[2] * pt(1, 0)), 2) + 12 * alfa[0] * alfa[2] * pt(1, 0) +
                      6 * alfa[0] * alfa[0]) +
              std::pow(pt(0, 0), 3) * 8 * (alfa[1] * alfa[2] * pt(1, 0) + alfa[0] * alfa[1]) +
              3 * alfa[1] * alfa[1] * std::pow(pt(0, 0), 4));
      return basef_line / 24.0;
    }
    if (inte_num == 16)  // f(x,y,z) = xz^2
    {
      basef_line =
          3 * alfa[2] * alfa[2] * pt(0, 0) * std::pow(pt(1, 0), 4) +
          std::pow(pt(1, 0), 3) *
              (3 * alfa[1] * alfa[2] * pt(0, 0) * pt(0, 0) + 6 * alfa[0] * alfa[2] * pt(0, 0)) +
          pt(1, 0) * pt(1, 0) *
              (alfa[1] * alfa[1] * std::pow(pt(0, 0), 3) +
                  3 * alfa[0] * alfa[1] * pt(0, 0) * pt(0, 0) + 3 * alfa[0] * alfa[0] * pt(0, 0));
      return basef_line / 6.0;
    }
    if (inte_num == 17)  // f(x,y,z) = y^3
    {
      basef_line = std::pow(pt(0, 0), 4) * 5 * (alfa[2] * pt(1, 0) + alfa[0]) +
                   4 * alfa[1] * std::pow(pt(0, 0), 5);
      return basef_line * 0.05;
    }
    if (inte_num == 18)  // f(x,y,z) = y^2 z
    {
      basef_line =
          4 * alfa[2] * std::pow((pt(0, 0) * pt(1, 0)), 2) * pt(0, 0) +
          pt(1, 0) * (3 * alfa[1] * std::pow(pt(0, 0), 4) + 4 * alfa[0] * std::pow(pt(0, 0), 3));
      return basef_line / 12.0;
    }
    if (inte_num == 19)  // f(x,y,z) = yz^2
    {
      basef_line = pt(1, 0) * pt(1, 0) *
                   (pt(0, 0) * pt(0, 0) * (3 * alfa[2] * pt(1, 0) + 3 * alfa[0]) +
                       2 * alfa[1] * std::pow(pt(0, 0), 3));
      return basef_line / 6.0;
    }
    if (inte_num == 20)  // f(x,y,z) = z^3
    {
      basef_line =
          std::pow(pt(1, 0), 3) * (alfa[2] * pt(0, 0) * pt(1, 0) +
                                      0.5 * alfa[1] * pt(0, 0) * pt(0, 0) + alfa[0] * pt(0, 0));
      return basef_line;
    }

    double a02 = std::pow(alfa[0], 2), a12 = std::pow(alfa[1], 2), a22 = std::pow(alfa[2], 2),
           y2 = std::pow(pt(0, 0), 2), z2 = std::pow(pt(1, 0), 2);
    double a03 = std::pow(alfa[0], 3), a13 = std::pow(alfa[1], 3), a23 = std::pow(alfa[2], 3),
           y3 = std::pow(pt(0, 0), 3), z3 = std::pow(pt(1, 0), 3);
    double a04 = std::pow(alfa[0], 4), a14 = std::pow(alfa[1], 4), a24 = std::pow(alfa[2], 4),
           y4 = std::pow(pt(0, 0), 4), z4 = std::pow(pt(1, 0), 4);
    double a05 = std::pow(alfa[0], 5), a15 = std::pow(alfa[1], 5), a25 = std::pow(alfa[2], 5),
           y5 = std::pow(pt(0, 0), 5), z5 = std::pow(pt(1, 0), 5);
    double y6 = std::pow(pt(0, 0), 6);
    if (inte_num == 21)  // f(x,y,z) = x^4
    {
      basef_line = 6 * a25 * pt(0, 0) * z5 +
                   (15 * alfa[1] * y2 + 30 * alfa[0] * pt(0, 0)) * a24 * z4 +
                   (20 * a12 * y3 + 60 * alfa[0] * alfa[1] * y2 + 60 * a02 * pt(0, 0)) * a23 * z3 +
                   (15 * a13 * y4 + 60 * alfa[0] * a12 * y3 + 90 * a02 * alfa[1] * y2 +
                       60 * a03 * pt(0, 0)) *
                       a22 * z2 +
                   (6 * a14 * y5 + 30 * alfa[0] * a13 * y4 + 60 * a02 * a12 * y3 +
                       60 * a03 * alfa[1] * y2 + 30 * a04 * pt(0, 0)) *
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
      basef_line = y3 * (20 * a23 * z3 + 60 * alfa[0] * a22 * z2 + 60 * a02 * alfa[2] * pt(1, 0) +
                            20 * a03) +
                   y4 * (45 * alfa[1] * a22 * z2 + 90 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                            45 * a02 * alfa[1]) +
                   y5 * (36 * a12 * alfa[2] * pt(1, 0) + 36 * alfa[0] * a12) + 10 * a13 * y6;
      return basef_line / 180.0;
    }
    if (inte_num == 25)  // f(x,y,z) = x^2*yz
    {
      basef_line = pt(1, 0) *
                   (y2 * (10 * a23 * z3 + 30 * alfa[0] * a22 * z2 + 30 * a02 * alfa[2] * pt(1, 0) +
                             10 * a03) +
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
                   y5 * (24 * alfa[1] * alfa[2] * pt(1, 0) + 24 * alfa[0] * alfa[1]) +
                   10 * a12 * y6;
      return basef_line / 120.0;
    }
    if (inte_num == 28)  // f(x,y,z) = x*y^2*z
    {
      basef_line =
          pt(1, 0) *
          (y3 * (10 * a22 * z2 + 20 * alfa[0] * alfa[2] * pt(1, 0) + 10 * a02) +
              y4 * (15 * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1]) + 6 * a12 * y5);
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


    double a06 = std::pow(alfa[0], 6), a16 = std::pow(alfa[1], 6), a26 = std::pow(alfa[2], 6),
           z6 = std::pow(pt(1, 0), 6), y7 = std::pow(pt(0, 0), 7);
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
          std::pow(pt(0, 0), 3) *
              (70 * alfa[1] * a24 * z4 + 280 * alfa[0] * alfa[1] * a23 * z3 +
                  420 * a02 * alfa[1] * a22 * z2 + 280 * a03 * alfa[1] * alfa[2] * pt(1, 0) +
                  70 * a04 * alfa[1]) +
          std::pow(pt(0, 0), 4) * (105 * a12 * a23 * z3 + 315 * alfa[0] * a12 * a22 * z2 +
                                      315 * a02 * a12 * alfa[2] * pt(1, 0) + 105 * a03 * a12) +
          std::pow(pt(0, 0), 5) *
              (84 * a13 * a22 * z2 + 168 * alfa[0] * a13 * alfa[2] * pt(1, 0) + 84 * a02 * a13) +
          std::pow(pt(0, 0), 6) * (35 * a14 * alfa[2] * pt(1, 0) + 35 * alfa[0] * a14) +
          6 * a15 * y7;
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
              10 * a03 *
                  (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
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
          std::pow(pt(0, 0), 4) *
              (105 * alfa[1] * a23 * z3 + 315 * alfa[0] * alfa[1] * a22 * z2 +
                  315 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 105 * a03 * alfa[1]) +
          std::pow(pt(0, 0), 5) *
              (126 * a12 * a22 * z2 + 252 * alfa[0] * a12 * alfa[2] * pt(1, 0) + 126 * a02 * a12) +
          std::pow(pt(0, 0), 6) * (70 * a13 * alfa[2] * pt(1, 0) + 70 * alfa[0] * a13) +
          15 * a14 * y7;
      return basef_line / 420.0;
    }
    if (inte_num == 40)  // f(x,y,z) = x^3*yz
    {
      basef_line =
          pt(1, 0) * (y2 * (15 * a24 * z4 + 60 * alfa[0] * a23 * z3 + 90 * a02 * a22 * z2 +
                               60 * a03 * alfa[2] * pt(1, 0) + 15 * a04) +
                         y3 * (40 * alfa[1] * a23 * z3 + 120 * alfa[0] * alfa[1] * a22 * z2 +
                                  120 * a02 * alfa[1] * alfa[2] * pt(1, 0) + 40 * a03 * alfa[1]) +
                         y4 * (45 * a12 * a22 * z2 + 90 * alfa[0] * a12 * alfa[2] * pt(1, 0) +
                                  45 * a02 * a12) +
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
      basef_line = pt(1, 0) *
                   (y3 * (20 * a23 * z3 + 60 * alfa[0] * a22 * z2 + 60 * a02 * alfa[2] * pt(1, 0) +
                             20 * a03) +
                       y4 * (45 * alfa[1] * a22 * z2 + 90 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                                45 * a02 * alfa[1]) +
                       y5 * (36 * a12 * alfa[2] * pt(1, 0) + 36 * alfa[0] * a12) + 10 * a13 * y6);
      return basef_line / 180.0;
    }
    if (inte_num == 44)  // f(x,y,z) = x^2*yz^2
    {
      basef_line =
          z2 * (y2 * (10 * a23 * z3 + 30 * alfa[0] * a22 * z2 + 30 * a02 * alfa[2] * pt(1, 0) +
                         10 * a03) +
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
                   y6 * (35 * alfa[1] * alfa[2] * pt(1, 0) + 35 * alfa[0] * alfa[1]) +
                   15 * a12 * y7;
      return basef_line / 210.0;
    }
    if (inte_num == 47)  // f(x,y,z) = xy^3*z
    {
      basef_line =
          pt(1, 0) *
          (y4 * (15 * a22 * z2 + 30 * alfa[0] * alfa[2] * pt(1, 0) + 15 * a02) +
              y5 * (24 * alfa[1] * alfa[2] * pt(1, 0) + 24 * alfa[0] * alfa[1]) + 10 * a12 * y6);
      return basef_line / 120.0;
    }
    if (inte_num == 48)  // f(x,y,z) = xy^2*z^2
    {
      basef_line = z2 * (y3 * (10 * a22 * z2 + 20 * alfa[0] * alfa[2] * pt(1, 0) + 10 * a02) +
                            y4 * (15 * alfa[1] * alfa[2] * pt(1, 0) + 15 * alfa[0] * alfa[1]) +
                            6 * a12 * y5);
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

    double y8 = std::pow(pt(0, 0), 8);
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
        basef_line = pt(0, 0) * pt(1, 0) * std::pow((alfa[2] * pt(1, 0) + alfa[0]), 6) / 6.0;
      else
        basef_line = pt(1, 0) * std::pow((alfa[2] * pt(1, 0) + alfa[1] * pt(0, 0) + alfa[0]), 7) /
                     42.0 / alfa[1];
      return basef_line;
    }
    if (inte_num == 60)  // f(x,y,z) = x^4*y^2
    {
      basef_line =
          y3 * (56.0 * a25 * z5 + 280.0 * alfa[0] * a24 * z4 + 560.0 * a02 * a23 * z3 +
                   560.0 * a03 * a22 * z2 + 280.0 * a04 * alfa[2] * pt(1, 0) + 56.0 * a05) +
          y4 * (210.0 * alfa[1] * a24 * z4 + 840.0 * alfa[0] * alfa[1] * a23 * z3 +
                   1260.0 * a02 * alfa[1] * a22 * z2 + 840.0 * a03 * alfa[1] * alfa[2] * pt(1, 0) +
                   210.0 * a04 * alfa[1]) +
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
                       420.0 * a02 * alfa[1] * a22 * z2 +
                       280.0 * a03 * alfa[1] * alfa[2] * pt(1, 0) + 70.0 * a04 * alfa[1]) +
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
                   y7 * (160.0 * a13 * alfa[2] * pt(1, 0) + 160.0 * alfa[0] * a13) +
                   35.0 * a14 * y8;
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
          z3 * (a24 * pt(0, 0) * z4 +
                   4.0 * alfa[0] *
                       (a23 * pt(0, 0) * z3 + 1.5 * alfa[1] * a22 * y2 * z2 +
                           a12 * alfa[2] * y3 * pt(1, 0) + 0.25 * a13 * y4) +
                   2.0 * alfa[1] * a23 * y2 * z3 +
                   6.0 * a02 *
                       (a22 * pt(0, 0) * z2 + alfa[1] * alfa[2] * y2 * pt(1, 0) + a12 * y3 / 3.0) +
                   2.0 * a12 * a22 * y3 * z2 +
                   4.0 * a03 * (alfa[2] * pt(0, 0) * pt(1, 0) + 0.5 * alfa[1] * y2) +
                   a13 * alfa[2] * y4 * pt(1, 0) + 0.2 * a14 * y5 + a04 * pt(0, 0));
      return basef_line * 0.25;
    }
    if (inte_num == 67)  // f(x,y,z) = x^2*y^4
    {
      basef_line =
          y5 * (56.0 * a23 * z3 + 168.0 * alfa[0] * a22 * z2 + 168.0 * a02 * alfa[2] * pt(1, 0) +
                   56.0 * a03) +
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
          z2 * (y3 * (20.0 * a23 * z3 + 60.0 * alfa[0] * a22 * z2 +
                         60.0 * a02 * alfa[2] * pt(1, 0) + 20.0 * a03) +
                   y4 * (45.0 * alfa[1] * a22 * z2 + 90.0 * alfa[0] * alfa[1] * alfa[2] * pt(1, 0) +
                            45.0 * a02 * alfa[1]) +
                   y5 * (36.0 * a12 * alfa[2] * pt(1, 0) + 36.0 * alfa[0] * a12) + 10.0 * a13 * y6);
      return basef_line / 180.0;
    }
    if (inte_num == 70)  // f(x,y,z) = x^2*yz^3
    {
      basef_line =
          z3 * (y2 * (10.0 * a23 * z3 + 30.0 * alfa[0] * a22 * z2 +
                         30.0 * a02 * alfa[2] * pt(1, 0) + 10.0 * a03) +
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
      basef_line =
          pt(1, 0) * (y6 * (7.0 * alfa[2] * pt(1, 0) + 7.0 * alfa[0]) + 6.0 * alfa[1] * y7);
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

    FOUR_C_THROW("the base function required is not defined");
    return 0.0;
  }

  /*!
  \brief Returns the base function when the boundaycell is projected in y-z plane
  */
  // pt(0,0) --> y
  // pt(1,0) --> z
  double base_func_surf_x(
      const Core::LinAlg::Matrix<2, 1> &pt, int inte_num, std::vector<double> alfa)
  {
    double basef_surf = 0.0;
    if (inte_num == 1)  // f(x,y,z) = 1
    {
      basef_surf = pt(0, 0);
      return basef_surf;
    }
    double y = pt(0, 0), y2 = std::pow(pt(0, 0), 2), z = pt(1, 0);
    double a0 = alfa[0], a1 = alfa[1], a2 = alfa[2];
    if (inte_num == 2)  // f(x,y,z) = x
    {
      basef_surf = a2 * y * z + 0.5 * a1 * y2 + a0 * y;
      return basef_surf;
    }
    if (inte_num == 3)  // f(x,y,z) = y
    {
      basef_surf = 0.5 * y * y;
      return basef_surf;
    }
    if (inte_num == 4)  // f(x,y,z) = z
    {
      basef_surf = z * y;
      return basef_surf;
    }

    double a22 = std::pow(a2, 2), a02 = std::pow(a0, 2), a12 = std::pow(a1, 2);
    double z2 = std::pow(z, 2), y3 = std::pow(y, 3);
    if (inte_num == 5)  // f(x,y,z) = x^2
    {
      basef_surf = a22 * y * z2 + 2 * a0 * (a2 * y * z + 0.5 * a1 * y2) + a1 * a2 * y2 * z +
                   a12 * y3 / 3.0 + a02 * y;
      return basef_surf;
    }
    if (inte_num == 6)  // f(x,y,z) = xy
    {
      basef_surf = y2 * (3 * a2 * z + 3 * a0) + 2 * a1 * y3;
      return basef_surf / 6.0;
    }
    if (inte_num == 7)  // f(x,y,z) = xz
    {
      basef_surf = (a2 * y * z + 0.5 * a1 * y2 + a0 * y) * z;
      return basef_surf;
    }
    if (inte_num == 8)  // f(x,y,z) = y^2
    {
      basef_surf = y3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 9)  // f(x,y,z) = yz
    {
      basef_surf = 0.5 * y2 * z;
      return basef_surf;
    }
    if (inte_num == 10)  // f(x,y,z) = z^2
    {
      basef_surf = z2 * y;
      return basef_surf;
    }

    double a03 = std::pow(a0, 3), a13 = std::pow(a1, 3), a23 = std::pow(a2, 3);
    double y4 = std::pow(y, 4), z3 = std::pow(z, 3);
    if (inte_num == 11)  // f(x,y,z) = x^3
    {
      basef_surf = a23 * y * z3 + 3.0 * a0 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
                   1.5 * a1 * a22 * y2 * z2 + 3.0 * a02 * (a2 * y * z + 0.5 * a1 * y2) +
                   a12 * a2 * y3 * z + 0.25 * a13 * y4 + a03 * y;
      return basef_surf;
    }
    if (inte_num == 12)  // f(x,y,z) = x^2*y
    {
      basef_surf = y2 * (6.0 * a22 * z2 + 12.0 * a0 * a2 * z + 6.0 * a02) +
                   y3 * (8.0 * a1 * a2 * z + 8.0 * a0 * a1) + 3.0 * a12 * y4;
      return basef_surf / 12.0;
    }
    if (inte_num == 13)  // f(x,y,z) = x^2*z
    {
      basef_surf = z * (a22 * y * z2 + 2.0 * a0 * (a2 * y * z + 0.5 * a1 * y2) + a1 * a2 * y2 * z +
                           a12 * y3 / 3.0 + a02 * y);
      return basef_surf;
    }
    if (inte_num == 14)  // f(x,y,z) = xy^2
    {
      basef_surf = y3 * (4.0 * a2 * z + 4.0 * a0) + 3.0 * a1 * y4;
      return basef_surf / 12.0;
    }
    if (inte_num == 15)  // f(x,y,z) = xyz
    {
      basef_surf = z * (y2 * (3.0 * a2 * z + 3.0 * a0) + 2.0 * a1 * y3);
      return basef_surf / 6.0;
    }
    if (inte_num == 16)  // f(x,y,z) = xz^2
    {
      basef_surf = z2 * (a2 * y * z + 0.5 * a1 * y2 + a0 * y);
      return basef_surf;
    }
    if (inte_num == 17)  // f(x,y,z) = y^3
    {
      basef_surf = 0.25 * y4;
      return basef_surf;
    }
    if (inte_num == 18)  // f(x,y,z) = y^2*z
    {
      basef_surf = y3 * z / 3.0;
      return basef_surf;
    }
    if (inte_num == 19)  // f(x,y,z) = yz^2
    {
      basef_surf = y2 * z2 * 0.5;
      return basef_surf;
    }
    if (inte_num == 20)  // f(x,y,z) = z^3
    {
      basef_surf = y * z3;
      return basef_surf;
    }

    double a04 = std::pow(a0, 4), a14 = std::pow(a1, 4), a24 = std::pow(a2, 4);
    double y5 = std::pow(y, 5), z4 = std::pow(z, 4);
    if (inte_num == 21)  // f(x,y,z) = x^4
    {
      basef_surf =
          a24 * y * z4 +
          4.0 * a0 *
              (a23 * y * z3 + 1.5 * a1 * a22 * y2 * z2 + a12 * a2 * y3 * z + 0.25 * a13 * y4) +
          2.0 * a1 * a23 * y2 * z3 +
          6.0 * a02 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
          2.0 * a12 * a22 * y3 * z2 + 4.0 * a03 * (a2 * y * z + 0.5 * a1 * y2) + a13 * a2 * y4 * z +
          0.2 * a14 * y5 + a04 * y;
      return basef_surf;
    }
    if (inte_num == 22)  // f(x,y,z) = x^3*y
    {
      basef_surf =
          y2 * (10.0 * a23 * z3 + 30.0 * a0 * a22 * z2 + 30.0 * a02 * a2 * z + 10.0 * a03) +
          y3 * (20.0 * a1 * a22 * z2 + 40.0 * a0 * a1 * a2 * z + 20.0 * a02 * a1) +
          y4 * (15.0 * a12 * a2 * z + 15.0 * a0 * a12) + 4.0 * a13 * y5;
      return basef_surf * 0.05;
    }
    if (inte_num == 23)  // f(x,y,z) = x^3*z
    {
      basef_surf =
          z * (a23 * y * z3 + 3.0 * a0 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
                  1.5 * a1 * a22 * y2 * z2 + 3.0 * a02 * (a2 * y * z + 0.5 * a1 * y2) +
                  a12 * a2 * y3 * z + 0.25 * a13 * y4 + a03 * y);
      return basef_surf;
    }
    if (inte_num == 24)  // f(x,y,z) = x^2*y^2
    {
      basef_surf = y3 * (10.0 * a22 * z2 + 20.0 * a0 * a2 * z + 10.0 * a02) +
                   y4 * (15.0 * a1 * a2 * z + 15.0 * a0 * a1) + 6.0 * a12 * y5;
      return basef_surf / 30.0;
    }
    if (inte_num == 25)  // f(x,y,z) = x^2*yz
    {
      basef_surf = z * (y2 * (6.0 * a22 * z2 + 12.0 * a0 * a2 * z + 6.0 * a02) +
                           y3 * (8.0 * a1 * a2 * z + 8.0 * a0 * a1) + 3.0 * a12 * y4);
      return basef_surf / 12.0;
    }
    if (inte_num == 26)  // f(x,y,z) = x^2*z^2
    {
      basef_surf = z2 * (a22 * y * z2 + 2.0 * a0 * (a2 * y * z + 0.5 * a1 * y2) + a1 * a2 * y2 * z +
                            a12 * y3 / 3.0 + a02 * y);
      return basef_surf;
    }
    if (inte_num == 27)  // f(x,y,z) = x*y^3
    {
      basef_surf = y4 * (5.0 * a2 * z + 5.0 * a0) + 4.0 * a1 * y5;
      return basef_surf * 0.05;
    }
    if (inte_num == 28)  // f(x,y,z) = x*y^2*z
    {
      basef_surf = z * (y3 * (4.0 * a2 * z + 4.0 * a0) + 3.0 * a1 * y4);
      return basef_surf / 12.0;
    }
    if (inte_num == 29)  // f(x,y,z) = xyz^2
    {
      basef_surf = z2 * (y2 * (3.0 * a2 * z + 3.0 * a0) + 2.0 * a1 * y3);
      return basef_surf / 6.0;
    }
    if (inte_num == 30)  // f(x,y,z) = xz^3
    {
      basef_surf = z3 * (a2 * y * z + 0.5 * a1 * y2 + a0 * y);
      return basef_surf;
    }
    if (inte_num == 31)  // f(x,y,z) = y^4
    {
      basef_surf = 0.2 * y5;
      return basef_surf;
    }
    if (inte_num == 32)  // f(x,y,z) = y^3*z
    {
      basef_surf = 0.25 * y4 * z;
      return basef_surf;
    }
    if (inte_num == 33)  // f(x,y,z) = y^2*z^2
    {
      basef_surf = y3 * z2 / 3.0;
      return basef_surf;
    }
    if (inte_num == 34)  // f(x,y,z) = yz^3
    {
      basef_surf = 0.5 * y2 * z3;
      return basef_surf;
    }
    if (inte_num == 35)  // f(x,y,z) = z^4
    {
      basef_surf = y * z4;
      return basef_surf;
    }

    double a05 = std::pow(a0, 5), a15 = std::pow(a1, 5), a25 = std::pow(a2, 5);
    double y6 = std::pow(y, 6), z5 = std::pow(z, 5);
    if (inte_num == 36)  // f(x,y,z) = x^5
    {
      basef_surf =
          a25 * y * z5 +
          5.0 * a0 *
              (a24 * y * z4 + 2.0 * a1 * a23 * y2 * z3 + 2.0 * a12 * a22 * y3 * z2 +
                  a13 * a2 * y4 * z + 0.2 * a14 * y5) +
          2.5 * a1 * a24 * y2 * z4 +
          10.0 * a02 *
              (a23 * y * z3 + 1.5 * a1 * a22 * y2 * z2 + a12 * a2 * y3 * z + 0.25 * a13 * y4) +
          10.0 * a12 * a23 * y3 * z3 / 3.0 +
          10.0 * a03 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
          2.5 * a13 * a22 * y4 * z2 + 5.0 * a04 * (a2 * y * z + 0.5 * a1 * y2) + a14 * a2 * y5 * z +
          a15 * y6 / 6.0 + a05 * y;
      return basef_surf;
    }
    if (inte_num == 37)  // f(x,y,z) = x^4*y
    {
      basef_surf = y2 * (15.0 * a24 * z4 + 60.0 * a0 * a23 * z3 + 90.0 * a02 * a22 * z2 +
                            60.0 * a03 * a2 * z + 15.0 * a04) +
                   y3 * (40.0 * a1 * a23 * z3 + 120.0 * a0 * a1 * a22 * z2 +
                            120.0 * a02 * a1 * a2 * z + 40.0 * a03 * a1) +
                   y4 * (45.0 * a12 * a22 * z2 + 90.0 * a0 * a12 * a2 * z + 45.0 * a02 * a12) +
                   y5 * (24.0 * a13 * a2 * z + 24.0 * a0 * a13) + 5.0 * a14 * y6;
      return basef_surf / 30.0;
    }
    if (inte_num == 38)  // f(x,y,z) = x^4*z
    {
      basef_surf =
          z *
          (a24 * y * z4 +
              4.0 * a0 *
                  (a23 * y * z3 + 1.5 * a1 * a22 * y2 * z2 + a12 * a2 * y3 * z + 0.25 * a13 * y4) +
              2.0 * a1 * a23 * y2 * z3 +
              6.0 * a02 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
              2.0 * a12 * a22 * y3 * z2 + 4.0 * a03 * (a2 * y * z + 0.5 * a1 * y2) +
              a13 * a2 * y4 * z + 0.2 * a14 * y5 + a04 * y);
      return basef_surf;
    }
    if (inte_num == 39)  // f(x,y,z) = x^3*y^2
    {
      basef_surf =
          y3 * (20.0 * a23 * z3 + 60.0 * a0 * a22 * z2 + 60.0 * a02 * a2 * z + 20.0 * a03) +
          y4 * (45.0 * a1 * a22 * z2 + 90.0 * a0 * a1 * a2 * z + 45.0 * a02 * a1) +
          y5 * (36.0 * a12 * a2 * z + 36.0 * a0 * a12) + 10.0 * a13 * y6;
      return basef_surf / 60.0;
    }
    if (inte_num == 40)  // f(x,y,z) = x^3*yz
    {
      basef_surf =
          z * (y2 * (10.0 * a23 * z3 + 30.0 * a0 * a22 * z2 + 30.0 * a02 * a2 * z + 10.0 * a03) +
                  y3 * (20.0 * a1 * a22 * z2 + 40.0 * a0 * a1 * a2 * z + 20.0 * a02 * a1) +
                  y4 * (15.0 * a12 * a2 * z + 15.0 * a0 * a12) + 4.0 * a13 * y5);
      return basef_surf * 0.05;
    }
    if (inte_num == 41)  // f(x,y,z) = x^3*z^2
    {
      basef_surf =
          z2 * (a23 * y * z3 + 3.0 * a0 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
                   1.5 * a1 * a22 * y2 * z2 + 3.0 * a02 * (a2 * y * z + 0.5 * a1 * y2) +
                   a12 * a2 * y3 * z + 0.25 * a13 * y4 + a03 * y);
      return basef_surf;
    }
    if (inte_num == 42)  // f(x,y,z) = x^2*y^3
    {
      basef_surf = y4 * (15.0 * a22 * z2 + 30.0 * a0 * a2 * z + 15.0 * a02) +
                   y5 * (24.0 * a1 * a2 * z + 24.0 * a0 * a1) + 10.0 * a12 * y6;
      return basef_surf / 60.0;
    }
    if (inte_num == 43)  // f(x,y,z) = x^2*y^2*z
    {
      basef_surf = z * (y3 * (10.0 * a22 * z2 + 20.0 * a0 * a2 * z + 10.0 * a02) +
                           y4 * (15.0 * a1 * a2 * z + 15.0 * a0 * a1) + 6.0 * a12 * y5);
      return basef_surf / 30.0;
    }
    if (inte_num == 44)  // f(x,y,z) = x^2*yz^2
    {
      basef_surf = z2 * (y2 * (6.0 * a22 * z2 + 12.0 * a0 * a2 * z + 6.0 * a02) +
                            y3 * (8.0 * a1 * a2 * z + 8.0 * a0 * a1) + 3.0 * a12 * y4);
      return basef_surf / 12.0;
    }
    if (inte_num == 45)  // f(x,y,z) = x^2*z^3
    {
      basef_surf = z3 * (a22 * y * z2 + 2.0 * a0 * (a2 * y * z + 0.5 * a1 * y2) + a1 * a2 * y2 * z +
                            a12 * y3 / 3.0 + a02 * y);
      return basef_surf;
    }
    if (inte_num == 46)  // f(x,y,z) = xy^4
    {
      basef_surf = y5 * (6.0 * a2 * z + 6.0 * a0) + 5.0 * a1 * y6;
      return basef_surf / 30.0;
    }
    if (inte_num == 47)  // f(x,y,z) = xy^3*z
    {
      basef_surf = z * (y4 * (5.0 * a2 * z + 5.0 * a0) + 4.0 * a1 * y5);
      return basef_surf * 0.05;
    }
    if (inte_num == 48)  // f(x,y,z) = xy^2*z^2
    {
      basef_surf = z2 * (y3 * (4.0 * a2 * z + 4.0 * a0) + 3.0 * a1 * y4);
      return basef_surf / 12.0;
    }
    if (inte_num == 49)  // f(x,y,z) = xyz^3
    {
      basef_surf = z3 * (y2 * (3.0 * a2 * z + 3.0 * a0) + 2.0 * a1 * y3);
      return basef_surf / 6.0;
    }
    if (inte_num == 50)  // f(x,y,z) = xz^4
    {
      basef_surf = z4 * (a2 * y * z + 0.5 * a1 * y2 + a0 * y);
      return basef_surf;
    }
    if (inte_num == 51)  // f(x,y,z) = y^5
    {
      basef_surf = y6 / 6.0;
      return basef_surf;
    }
    if (inte_num == 52)  // f(x,y,z) = y^4*z
    {
      basef_surf = 0.2 * y5 * z;
      return basef_surf;
    }
    if (inte_num == 53)  // f(x,y,z) = y^3*z^2
    {
      basef_surf = 0.25 * y4 * z2;
      return basef_surf;
    }
    if (inte_num == 54)  // f(x,y,z) = y^2*z^3
    {
      basef_surf = y3 * z3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 55)  // f(x,y,z) = yz^4
    {
      basef_surf = y2 * z4 * 0.5;
      return basef_surf;
    }
    if (inte_num == 56)  // f(x,y,z) = z^5
    {
      basef_surf = y * z5;
      return basef_surf;
    }

    double y7 = std::pow(y, 7), z6 = std::pow(z, 6);
    if (inte_num == 57)  // f(x,y,z) = x^6
    {
      if (fabs(a1) > 0.00000001)
        basef_surf = std::pow((a2 * z + a1 * y + a0), 7) / 7.0 / a1;
      else
        basef_surf = std::pow((a0 + a2 * z), 6) * y;
      return basef_surf;
    }
    if (inte_num == 58)  // f(x,y,z) = x^5*y
    {
      basef_surf =
          y2 * (21.0 * a25 * z5 + 105.0 * a0 * a24 * z4 + 210.0 * a02 * a23 * z3 +
                   210.0 * a03 * a22 * z2 + 105.0 * a04 * a2 * z + 21.0 * a05) +
          y3 * (70.0 * a1 * a24 * z4 + 280.0 * a0 * a1 * a23 * z3 + 420.0 * a02 * a1 * a22 * z2 +
                   280.0 * a03 * a1 * a2 * z + 70.0 * a04 * a1) +
          y4 * (105.0 * a12 * a23 * z3 + 315.0 * a0 * a12 * a22 * z2 + 315.0 * a02 * a12 * a2 * z +
                   105.0 * a03 * a12) +
          y5 * (84.0 * a13 * a22 * z2 + 168.0 * a0 * a13 * a2 * z + 84.0 * a02 * a13) +
          y6 * (35.0 * a14 * a2 * z + 35.0 * a0 * a14) + 6.0 * a15 * y7;
      return basef_surf / 42.0;
    }
    if (inte_num == 59)  // f(x,y,z) = x^5*z
    {
      basef_surf =
          z *
          (a25 * y * z5 +
              5.0 * a0 *
                  (a24 * y * z4 + 2.0 * a1 * a23 * y2 * z3 + 2.0 * a12 * a22 * y3 * z2 +
                      a13 * a2 * y4 * z + 0.2 * a14 * y5) +
              2.5 * a1 * a24 * y2 * z4 +
              10.0 * a02 *
                  (a23 * y * z3 + 1.5 * a1 * a22 * y2 * z2 + a12 * a2 * y3 * z + 0.25 * a13 * y4) +
              10.0 / 3.0 * a12 * a23 * y3 * z3 +
              10.0 * a03 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
              2.5 * a13 * a22 * y4 * z2 + 5.0 * a04 * (a2 * y * z + 0.5 * a1 * y2) +
              a14 * a2 * y5 * z + a15 * y6 / 6.0 + a05 * y);
      return basef_surf;
    }
    if (inte_num == 60)  // f(x,y,z) = x^4*y^2
    {
      basef_surf = y3 * (35.0 * a24 * z4 + 140.0 * a0 * a23 * z3 + 210.0 * a02 * a22 * z2 +
                            140.0 * a03 * a2 * z + 35.0 * a04) +
                   y4 * (105.0 * a1 * a23 * z3 + 315.0 * a0 * a1 * a22 * z2 +
                            315.0 * a02 * a1 * a2 * z + 105.0 * a03 * a1) +
                   y5 * (126.0 * a12 * a22 * z2 + 252.0 * a0 * a12 * a2 * z + 126.0 * a02 * a12) +
                   y6 * (70.0 * a13 * a2 * z + 70.0 * a0 * a13) + 15.0 * a14 * y7;
      return basef_surf / 105.0;
    }
    if (inte_num == 61)  // f(x,y,z) = x^4*yz
    {
      basef_surf =
          z * (y2 * (15.0 * a24 * z4 + 60.0 * a0 * a23 * z3 + 90.0 * a02 * a22 * z2 +
                        60.0 * a03 * a2 * z + 15.0 * a04) +
                  y3 * (40.0 * a1 * a23 * z3 + 120.0 * a0 * a1 * a22 * z2 +
                           120.0 * a02 * a1 * a2 * z + 40.0 * a03 * a1) +
                  y4 * (45.0 * a12 * a22 * z2 + 90.0 * a0 * a12 * a2 * z + 45.0 * a02 * a12) +
                  y5 * (24.0 * a13 * a2 * z + 24.0 * a0 * a13) + 5.0 * a14 * y6);
      return basef_surf / 30.0;
    }
    if (inte_num == 62)  // f(x,y,z) = x^4*z^2
    {
      basef_surf =
          z2 *
          (a24 * y * z4 +
              4.0 * a0 *
                  (a23 * y * z3 + 1.5 * a1 * a22 * y2 * z2 + a12 * a2 * y3 * z + 0.25 * a13 * y4) +
              2.0 * a1 * a23 * y2 * z3 +
              6.0 * a02 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
              2.0 * a12 * a22 * y3 * z2 + 4.0 * a03 * (a2 * y * z + 0.5 * a1 * y2) +
              a13 * a2 * y4 * z + 0.2 * a14 * y5 + a04 * y);
      return basef_surf;
    }
    if (inte_num == 63)  // f(x,y,z) = x^3*y^3
    {
      basef_surf =
          y4 * (35.0 * a23 * z3 + 105.0 * a0 * a22 * z2 + 105.0 * a02 * a2 * z + 35.0 * a03) +
          y5 * (84.0 * a1 * a22 * z2 + 168.0 * a0 * a1 * a2 * z + 84.0 * a02 * a1) +
          y6 * (70.0 * a12 * a2 * z + 70.0 * a0 * a12) + 20.0 * a13 * y7;
      return basef_surf / 140.0;
    }
    if (inte_num == 64)  // f(x,y,z) = x^3*y^2*z
    {
      basef_surf =
          z * (y3 * (20.0 * a23 * z3 + 60.0 * a0 * a22 * z2 + 60.0 * a02 * a2 * z + 20.0 * a03) +
                  y4 * (45.0 * a1 * a22 * z2 + 90.0 * a0 * a1 * a2 * z + 45.0 * a02 * a1) +
                  y5 * (36.0 * a12 * a2 * z + 36.0 * a0 * a12) + 10.0 * a13 * y6);
      return basef_surf / 60.0;
    }
    if (inte_num == 65)  // f(x,y,z) = x^3*y*z^2
    {
      basef_surf =
          z2 * (y2 * (10.0 * a23 * z3 + 30.0 * a0 * a22 * z2 + 30.0 * a02 * a2 * z + 10.0 * a03) +
                   y3 * (20.0 * a1 * a22 * z2 + 40.0 * a0 * a1 * a2 * z + 20.0 * a02 * a1) +
                   y4 * (15.0 * a12 * a2 * z + 15.0 * a0 * a12) + 4.0 * a13 * y5);
      return basef_surf * 0.05;
    }
    if (inte_num == 66)  // f(x,y,z) = x^3*z^3
    {
      basef_surf =
          z3 * (a23 * y * z3 + 3.0 * a0 * (a22 * y * z2 + a1 * a2 * y2 * z + a12 * y3 / 3.0) +
                   1.5 * a1 * a22 * y2 * z2 + 3.0 * a02 * (a2 * y * z + 0.5 * a1 * y2) +
                   a12 * a2 * y3 * z + 0.25 * a13 * y4 + a03 * y);
      return basef_surf;
    }
    if (inte_num == 67)  // f(x,y,z) = x^2*y^4
    {
      basef_surf = y5 * (21.0 * a22 * z2 + 42.0 * a0 * a2 * z + 21.0 * a02) +
                   y6 * (35.0 * a1 * a2 * z + 35.0 * a0 * a1) + 15.0 * a12 * y7;
      return basef_surf / 105.0;
    }
    if (inte_num == 68)  // f(x,y,z) = x^2*y^3*z
    {
      basef_surf = z * (y4 * (15.0 * a22 * z2 + 30.0 * a0 * a2 * z + 15.0 * a02) +
                           y5 * (24.0 * a1 * a2 * z + 24.0 * a0 * a1) + 10.0 * a12 * y6);
      return basef_surf / 60.0;
    }
    if (inte_num == 69)  // f(x,y,z) = x^2*y^2*z^2
    {
      basef_surf = z2 * (y3 * (10.0 * a22 * z2 + 20.0 * a0 * a2 * z + 10.0 * a02) +
                            y4 * (15.0 * a1 * a2 * z + 15.0 * a0 * a1) + 6.0 * a12 * y5);
      return basef_surf / 30.0;
    }
    if (inte_num == 70)  // f(x,y,z) = x^2*yz^3
    {
      basef_surf = z3 * (y2 * (6.0 * a22 * z2 + 12.0 * a0 * a2 * z + 6.0 * a02) +
                            y3 * (8.0 * a1 * a2 * z + 8.0 * a0 * a1) + 3.0 * a12 * y4);
      return basef_surf / 12.0;
    }
    if (inte_num == 71)  // f(x,y,z) = x^2*z^4
    {
      basef_surf = z4 * (a22 * y * z2 + 2.0 * a0 * (a2 * y * z + 0.5 * a1 * y2) + a1 * a2 * y2 * z +
                            a12 * y3 / 3.0 + a02 * y);
      return basef_surf;
    }
    if (inte_num == 72)  // f(x,y,z) = xy^5
    {
      basef_surf = y6 * (7.0 * a2 * z + 7.0 * a0) + 6.0 * a1 * y7;
      return basef_surf / 42.0;
    }
    if (inte_num == 73)  // f(x,y,z) = xy^4*z
    {
      basef_surf = z * (y5 * (6.0 * a2 * z + 6.0 * a0) + 5.0 * a1 * y6);
      return basef_surf / 30.0;
    }
    if (inte_num == 74)  // f(x,y,z) = xy^3*z^2
    {
      basef_surf = z2 * (y4 * (5.0 * a2 * z + 5.0 * a0) + 4.0 * a1 * y5);
      return basef_surf * 0.05;
    }
    if (inte_num == 75)  // f(x,y,z) = xy^2*z^3
    {
      basef_surf = z3 * (y3 * (4.0 * a2 * z + 4.0 * a0) + 3.0 * a1 * y4);
      return basef_surf / 12.0;
    }
    if (inte_num == 76)  // f(x,y,z) = xyz^4
    {
      basef_surf = z4 * (y2 * (3.0 * a2 * z + 3.0 * a0) + 2.0 * a1 * y3);
      return basef_surf / 6.0;
    }
    if (inte_num == 77)  // f(x,y,z) = xz^5
    {
      basef_surf = z5 * (a2 * y * z + 0.5 * a1 * y2 + a0 * y);
      return basef_surf;
    }
    if (inte_num == 78)  // f(x,y,z) = y^6
    {
      basef_surf = y7 / 7.0;
      return basef_surf;
    }
    if (inte_num == 79)  // f(x,y,z) = y^5*z
    {
      basef_surf = y6 * z / 6.0;
      return basef_surf;
    }
    if (inte_num == 80)  // f(x,y,z) = y^4*z^2
    {
      basef_surf = 0.2 * y5 * z2;
      return basef_surf;
    }
    if (inte_num == 81)  // f(x,y,z) = y^3*z^3
    {
      basef_surf = 0.25 * y4 * z3;
      return basef_surf;
    }
    if (inte_num == 82)  // f(x,y,z) = y^2*z^4
    {
      basef_surf = y3 * z4 / 3.0;
      return basef_surf;
    }
    if (inte_num == 83)  // f(x,y,z) = yz^5
    {
      basef_surf = 0.5 * y2 * z5;
      return basef_surf;
    }
    if (inte_num == 84)  // f(x,y,z) = z^6
    {
      basef_surf = y * z6;
      return basef_surf;
    }

    FOUR_C_THROW("The base function for boundarycell integration undefined");
    exit(1);
    return 0.0;
  }

  /*!
  \brief Returns the base function when the boundaycell is projected in z-x plane
  */
  // pt(0,0) --> z
  // pt(1,0) --> x
  double base_func_surf_y(
      const Core::LinAlg::Matrix<2, 1> &pt, int inte_num, std::vector<double> alfa)
  {
    double basef_surf = 0.0;

    if (inte_num == 1)  // f(x,y,z) = 1
    {
      basef_surf = pt(0, 0);
      return basef_surf;
    }

    double a0 = alfa[0], a1 = alfa[1], a2 = alfa[2];
    double z = pt(0, 0), x = pt(1, 0), z2 = std::pow(pt(0, 0), 2);
    if (inte_num == 2)  // f(x,y,z) = x
    {
      basef_surf = x * z;
      return basef_surf;
    }
    if (inte_num == 3)  // f(x,y,z) = y
    {
      basef_surf = 0.5 * a1 * z2 + a2 * x * z + a0 * z;
      return basef_surf;
    }
    if (inte_num == 4)  // f(x,y,z) = z
    {
      basef_surf = 0.5 * z * z;
      return basef_surf;
    }

    double a22 = std::pow(a2, 2), a02 = std::pow(a0, 2), a12 = std::pow(a1, 2);
    double x2 = std::pow(x, 2), z3 = std::pow(z, 3);
    if (inte_num == 5)  // f(x,y,z) = x^2
    {
      basef_surf = x2 * z;
      return basef_surf;
    }
    if (inte_num == 6)  // f(x,y,z) = xy
    {
      basef_surf = x * (a1 * z2 * 0.5 + a2 * x * z + a0 * z);
      return basef_surf;
    }
    if (inte_num == 7)  // f(x,y,z) = xz
    {
      basef_surf = 0.5 * x * z2;
      return basef_surf;
    }
    if (inte_num == 8)  // f(x,y,z) = y^2
    {
      basef_surf = a12 * z3 / 3.0 + 2.0 * a0 * (0.5 * a1 * z2 + a2 * x * z) + a1 * a2 * x * z2 +
                   a22 * x2 * z + a02 * z;
      return basef_surf;
    }
    if (inte_num == 9)  // f(x,y,z) = yz
    {
      basef_surf = 2.0 * a1 * z3 + (3.0 * a2 * x + 3.0 * a0) * z2;
      return basef_surf / 6.0;
    }
    if (inte_num == 10)  // f(x,y,z) = z^2
    {
      basef_surf = z3 / 3.0;
      return basef_surf;
    }

    double a03 = std::pow(a0, 3), a13 = std::pow(a1, 3), a23 = std::pow(a2, 3);
    double x3 = std::pow(x, 3), z4 = std::pow(z, 4);
    if (inte_num == 11)  // f(x,y,z) = x^3
    {
      basef_surf = x3 * z;
      return basef_surf;
    }
    if (inte_num == 12)  // f(x,y,z) = x^2*y
    {
      basef_surf = x2 * (0.5 * a1 * z2 + a2 * x * z + a0 * z);
      return basef_surf;
    }
    if (inte_num == 13)  // f(x,y,z) = x^2*z
    {
      basef_surf = x2 * z2 * 0.5;
      return basef_surf;
    }
    if (inte_num == 14)  // f(x,y,z) = xy^2
    {
      basef_surf = x * (a12 * z3 / 3.0 + 2.0 * a0 * (a1 * z2 * 0.5 + a2 * x * z) +
                           a1 * a2 * x * z2 + a22 * x2 * z + a02 * z);
      return basef_surf;
    }
    if (inte_num == 15)  // f(x,y,z) = xyz
    {
      basef_surf = x * (2.0 * a1 * z3 + (3.0 * a2 * x + 3.0 * a0) * z2);
      return basef_surf / 6.0;
    }
    if (inte_num == 16)  // f(x,y,z) = xz^2
    {
      basef_surf = x * z3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 17)  // f(x,y,z) = y^3
    {
      basef_surf = 0.25 * a13 * z4 + 3.0 * a0 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
                   a12 * a2 * x * z3 + 3.0 * a02 * (0.5 * a1 * z2 + a2 * x * z) +
                   1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z + a03 * z;
      return basef_surf;
    }
    if (inte_num == 18)  // f(x,y,z) = y^2*z
    {
      basef_surf = 3.0 * a12 * z4 + (8.0 * a1 * a2 * x + 8.0 * a0 * a1) * z3 +
                   (6.0 * a22 * x2 + 12.0 * a0 * a2 * x + 6.0 * a02) * z2;
      return basef_surf / 12.0;
    }
    if (inte_num == 19)  // f(x,y,z) = yz^2
    {
      basef_surf = 3.0 * a1 * z4 + (4.0 * a2 * x + 4.0 * a0) * z3;
      return basef_surf / 12.0;
    }
    if (inte_num == 20)  // f(x,y,z) = z^3
    {
      basef_surf = 0.25 * z4;
      return basef_surf;
    }

    double a04 = std::pow(a0, 4), a14 = std::pow(a1, 4), a24 = std::pow(a2, 4);
    double x4 = std::pow(x, 4), z5 = std::pow(z, 5);
    if (inte_num == 21)  // f(x,y,z) = x^4
    {
      basef_surf = x4 * z;
      return basef_surf;
    }
    if (inte_num == 22)  // f(x,y,z) = x^3*y
    {
      basef_surf = x3 * (0.5 * a1 * z2 + a2 * x * z + a0 * z);
      return basef_surf;
    }
    if (inte_num == 23)  // f(x,y,z) = x^3*z
    {
      basef_surf = x3 * z2 * 0.5;
      return basef_surf;
    }
    if (inte_num == 24)  // f(x,y,z) x^2*y^2
    {
      basef_surf = x2 * (a12 * z3 / 3.0 + 2.0 * a0 * (0.5 * a1 * z2 + a2 * x * z) +
                            a1 * a2 * x * z2 + a22 * x2 * z + a02 * z);
      return basef_surf;
    }
    if (inte_num == 25)  // f(x,y,z) = x^2*yz
    {
      basef_surf = x2 * (2.0 * a1 * z3 + (3.0 * a2 * x + 3.0 * a0) * z2);
      return basef_surf / 6.0;
    }
    if (inte_num == 26)  // f(x,y,z) = x^2*z^2
    {
      basef_surf = x2 * z3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 27)  // f(x,y,z) = xy^3
    {
      basef_surf =
          x * (a13 * z4 * 0.25 + 3.0 * a0 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
                  a12 * a2 * x * z3 + 3.0 * a02 * (0.5 * a1 * z2 + a2 * x * z) +
                  1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z + a03 * z);
      return basef_surf;
    }
    if (inte_num == 28)  // f(x,y,z) = xy^2*z
    {
      basef_surf = x * (3.0 * a12 * z4 + (8.0 * a1 * a2 * x + 8.0 * a0 * a1) * z3 +
                           (6.0 * a22 * x2 + 12.0 * a0 * a2 * x + 6.0 * a02) * z2);
      return basef_surf / 12.0;
    }
    if (inte_num == 29)  // f(x,y,z) = xyz^2
    {
      basef_surf = x * (3.0 * a1 * z4 + (4.0 * a2 * x + 4.0 * a0) * z3);
      return basef_surf / 12.0;
    }
    if (inte_num == 30)  // f(x,y,z) = xz^3
    {
      basef_surf = 0.25 * x * z4;
      return basef_surf;
    }
    if (inte_num == 31)  // f(x,y,z) = y^4
    {
      basef_surf =
          0.2 * a14 * z5 +
          4.0 * a0 *
              (0.25 * a13 * z4 + a12 * a2 * x * z3 + 1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z) +
          a13 * a2 * x * z4 + 6.0 * a02 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
          2.0 * a12 * a22 * x2 * z3 + 4.0 * a03 * (0.5 * a1 * z2 + a2 * x * z) +
          2.0 * a1 * a23 * x3 * z2 + a24 * x4 * z + a04 * z;
      return basef_surf;
    }
    if (inte_num == 32)  // f(x,y,z) = y^3*z
    {
      basef_surf = 4.0 * a13 * z5 + (15.0 * a12 * a2 * x + 15.0 * a0 * a12) * z4 +
                   (20.0 * a1 * a22 * x2 + 40.0 * a0 * a1 * a2 * x + 20.0 * a02 * a1) * z3 +
                   (10.0 * a23 * x3 + 30.0 * a0 * a22 * x2 + 30.0 * a02 * a2 * x + 10.0 * a03) * z2;
      return basef_surf * 0.05;
    }
    if (inte_num == 33)  // f(x,y,z) = y^2*z^2
    {
      basef_surf = 6.0 * a12 * z5 + (15.0 * a1 * a2 * x + 15.0 * a0 * a1) * z4 +
                   (10.0 * a22 * x2 + 20.0 * a0 * a2 * x + 10.0 * a02) * z3;
      return basef_surf / 30.0;
    }
    if (inte_num == 34)  // f(x,y,z) = yz^3
    {
      basef_surf = 4.0 * a1 * z5 + (5.0 * a2 * x + 5.0 * a0) * z4;
      return basef_surf * 0.05;
    }
    if (inte_num == 35)  // f(x,y,z) = z^4
    {
      basef_surf = 0.2 * z5;
      return basef_surf;
    }

    double a05 = std::pow(a0, 5), a15 = std::pow(a1, 5), a25 = std::pow(a2, 5);
    double x5 = std::pow(x, 5), z6 = std::pow(z, 6);
    if (inte_num == 36)  // f(x,y,z) = x^5
    {
      basef_surf = x5 * z;
      return basef_surf;
    }
    if (inte_num == 37)  // f(x,y,z) = x^4*y
    {
      basef_surf = x4 * (0.5 * a1 * z2 + a2 * x * z + a0 * z);
      return basef_surf;
    }
    if (inte_num == 38)  // f(x,y,z) = x^4*z
    {
      basef_surf = x4 * z2 * 0.5;
      return basef_surf;
    }
    if (inte_num == 39)  // f(x,y,z) = x^3*y^2
    {
      basef_surf = x3 * (a12 * z3 / 3.0 + 2.0 * a0 * (0.5 * a1 * z2 + a2 * x * z) +
                            a1 * a2 * x * z2 + a22 * x2 * z + a02 * z);
      return basef_surf;
    }
    if (inte_num == 40)  // f(x,y,z) = x^3*yz
    {
      basef_surf = x3 * (2.0 * a1 * z3 + (3.0 * a2 * x + 3.0 * a0) * z2);
      return basef_surf / 6.0;
    }
    if (inte_num == 41)  // f(x,y,z) = x^3*z^2
    {
      basef_surf = x3 * z3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 42)  // f(x,y,z) = x^2*y^3
    {
      basef_surf =
          x2 * (0.25 * a13 * z4 + 3.0 * a0 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
                   a12 * a2 * x * z3 + 3.0 * a02 * (0.5 * a1 * z2 + a2 * x * z) +
                   1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z + a03 * z);
      return basef_surf;
    }
    if (inte_num == 43)  // f(x,y,z) = x^2*y^2*z
    {
      basef_surf = x2 * (3.0 * a12 * z4 + (8.0 * a1 * a2 * x + 8.0 * a0 * a1) * z3 +
                            (6.0 * a22 * x2 + 12.0 * a0 * a2 * x + 6.0 * a02) * z2);
      return basef_surf / 12.0;
    }
    if (inte_num == 44)  // f(x,y,z) = x^2*yz^2
    {
      basef_surf = x2 * (3.0 * a1 * z4 + (4.0 * a2 * x + 4.0 * a0) * z3);
      return basef_surf / 12.0;
    }
    if (inte_num == 45)  // f(x,y,z) = x^2*z^3
    {
      basef_surf = x2 * z4 * 0.25;
      return basef_surf;
    }
    if (inte_num == 46)  // f(x,y,z) = xy^4
    {
      basef_surf =
          x *
          (0.2 * a14 * z5 +
              4.0 * a0 *
                  (0.25 * a13 * z4 + a12 * a2 * x * z3 + 1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z) +
              a13 * a2 * x * z4 + 6.0 * a02 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
              2.0 * a12 * a22 * x2 * z3 + 4.0 * a03 * (0.5 * a1 * z2 + a2 * x * z) +
              2.0 * a1 * a23 * x3 * z2 + a24 * x4 * z + a04 * z);
      return basef_surf;
    }
    if (inte_num == 47)  // f(x,y,z) = xy^3*z
    {
      basef_surf =
          x * (4.0 * a13 * z5 + (15.0 * a12 * a2 * x + 15.0 * a0 * a12) * z4 +
                  (20.0 * a1 * a22 * x2 + 40.0 * a0 * a1 * a2 * x + 20.0 * a02 * a1) * z3 +
                  (10.0 * a23 * x3 + 30.0 * a0 * a22 * x2 + 30.0 * a02 * a2 * x + 10.0 * a03) * z2);
      return basef_surf * 0.05;
    }
    if (inte_num == 48)  // f(x,y,z) = xy^2*z^2
    {
      basef_surf = x * (6.0 * a12 * z5 + (15.0 * a1 * a2 * x + 15.0 * a0 * a1) * z4 +
                           (10.0 * a22 * x2 + 20.0 * a0 * a2 * x + 10.0 * a02) * z3);
      return basef_surf / 30.0;
    }
    if (inte_num == 49)  // f(x,y,z) = xyz^3
    {
      basef_surf = x * (4.0 * a1 * z5 + (5.0 * a2 * x + 5.0 * a0) * z4);
      return basef_surf * 0.05;
    }
    if (inte_num == 50)  // f(x,y,z) = xz^4
    {
      basef_surf = x * z5 * 0.2;
      return basef_surf;
    }
    if (inte_num == 51)  // f(x,y,z) = y^5
    {
      basef_surf =
          a15 * z6 / 6.0 +
          5.0 * a0 *
              (0.2 * a14 * z5 + a13 * a2 * x * z4 + 2.0 * a12 * a22 * x2 * z3 +
                  2.0 * a1 * a23 * x3 * z2 + a24 * x4 * z) +
          a14 * a2 * x * z5 +
          10.0 * a02 *
              (0.25 * a13 * z4 + a12 * a2 * x * z3 + 1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z) +
          2.5 * a13 * a22 * x2 * z4 +
          10.0 * a03 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
          10.0 / 3.0 * a12 * a23 * x3 * z3 + 5.0 * a04 * (0.5 * a1 * z2 + a2 * x * z) +
          2.5 * a1 * a24 * x4 * z2 + a25 * x5 * z + a05 * z;
      return basef_surf;
    }
    if (inte_num == 52)  // f(x,y,z) = y^4*z
    {
      basef_surf = 5.0 * a14 * z6 + (24.0 * a13 * a2 * x + 24.0 * a0 * a13) * z5 +
                   (45.0 * a12 * a22 * x2 + 90.0 * a0 * a12 * a2 * x + 45.0 * a02 * a12) * z4 +
                   (40.0 * a1 * a23 * x3 + 120.0 * a0 * a1 * a22 * x2 + 120.0 * a02 * a1 * a2 * x +
                       40.0 * a03 * a1) *
                       z3 +
                   (15.0 * a24 * x4 + 60.0 * a0 * a23 * x3 + 90.0 * a02 * a22 * x2 +
                       60.0 * a03 * a2 * x + 15.0 * a04) *
                       z2;
      return basef_surf / 30.0;
    }
    if (inte_num == 53)  // f(x,y,z) = y^3*z^2
    {
      basef_surf = 10.0 * a13 * z6 + (36.0 * a12 * a2 * x + 36.0 * a0 * a12) * z5 +
                   (45.0 * a1 * a22 * x2 + 90.0 * a0 * a1 * a2 * x + 45.0 * a02 * a1) * z4 +
                   (20.0 * a23 * x3 + 60.0 * a0 * a22 * x2 + 60.0 * a02 * a2 * x + 20.0 * a03) * z3;
      return basef_surf / 60.0;
    }
    if (inte_num == 54)  // f(x,y,z) = y^2*z^3
    {
      basef_surf = 10.0 * a12 * z6 + (24.0 * a1 * a2 * x + 24.0 * a0 * a1) * z5 +
                   (15.0 * a22 * x2 + 30.0 * a0 * a2 * x + 15.0 * a02) * z4;
      return basef_surf / 60.0;
    }
    if (inte_num == 55)  // f(x,y,z) = yz^4
    {
      basef_surf = 5.0 * a1 * z6 + (6.0 * a2 * x + 6.0 * a0) * z5;
      return basef_surf / 30.0;
    }
    if (inte_num == 56)  // f(x,y,z) = z^5
    {
      basef_surf = z6 / 6.0;
      return basef_surf;
    }

    double x6 = std::pow(x, 6), z7 = std::pow(z, 7);
    if (inte_num == 57)  // f(x,y,z) = x^6
    {
      basef_surf = x6 * z;
      return basef_surf;
    }
    if (inte_num == 58)  // f(x,y,z) = x^5*y
    {
      basef_surf = x5 * (0.5 * a1 * z2 + a2 * x * z + a0 * z);
      return basef_surf;
    }
    if (inte_num == 59)  // f(x,y,z) = x^5*z
    {
      basef_surf = x5 * z2 * 0.5;
      return basef_surf;
    }
    if (inte_num == 60)  // f(x,y,z) = x^4*y^2
    {
      basef_surf = x4 * (a12 * z3 / 3.0 + 2.0 * a0 * (0.5 * a1 * z2 + a2 * x * z) +
                            a1 * a2 * x * z2 + a22 * x2 * z + a02 * z);
      return basef_surf;
    }
    if (inte_num == 61)  // f(x,y,z) = x^4*yz
    {
      basef_surf = x4 * (2.0 * a1 * z3 + (3.0 * a2 * x + 3.0 * a0) * z2);
      return basef_surf / 6.0;
    }
    if (inte_num == 62)  // f(x,y,z) = x^4*z^2
    {
      basef_surf = x4 * z3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 63)  // f(x,y,z) = x^3*y^3
    {
      basef_surf =
          x3 * (0.25 * a13 * z4 + 3.0 * a0 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
                   a12 * a2 * x * z3 + 3.0 * a02 * (0.5 * a1 * z2 + a2 * x * z) +
                   1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z + a03 * z);
      return basef_surf;
    }
    if (inte_num == 64)  // f(x,y,z) = x^3*y^2*z
    {
      basef_surf = x3 * (3.0 * a12 * z4 + (8.0 * a1 * a2 * x + 8.0 * a0 * a1) * z3 +
                            (6.0 * a22 * x2 + 12.0 * a0 * a2 * x + 6.0 * a02) * z2);
      return basef_surf / 12.0;
    }
    if (inte_num == 65)  // f(x,y,z) = x^3*y*z^2
    {
      basef_surf = x3 * (3.0 * a1 * z4 + (4.0 * a2 * x + 4.0 * a0) * z3);
      return basef_surf / 12.0;
    }
    if (inte_num == 66)  // f(x,y,z) = x^3*z^3
    {
      basef_surf = x3 * z4 * 0.25;
      return basef_surf;
    }
    if (inte_num == 67)  // f(x,y,z) = x^2*y^4
    {
      basef_surf =
          x2 *
          (0.2 * a14 * z5 +
              4.0 * a0 *
                  (0.25 * a13 * z4 + a12 * a2 * x * z3 + 1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z) +
              a13 * a2 * x * z4 + 6.0 * a02 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
              2.0 * a12 * a22 * x2 * z3 + 4.0 * a03 * (0.5 * a1 * z2 + a2 * x * z) +
              2.0 * a1 * a23 * x3 * z2 + a24 * x4 * z + a04 * z);
      return basef_surf;
    }
    if (inte_num == 68)  // f(x,y,z) = x^2*y^3*z
    {
      basef_surf =
          x2 *
          (4.0 * a13 * z5 + (15.0 * a12 * a2 * x + 15.0 * a0 * a12) * z4 +
              (20.0 * a1 * a22 * x2 + 40.0 * a0 * a1 * a2 * x + 20.0 * a02 * a1) * z3 +
              (10.0 * a23 * x3 + 30.0 * a0 * a22 * x2 + 30.0 * a02 * a2 * x + 10.0 * a03) * z2);
      return basef_surf * 0.05;
    }
    if (inte_num == 69)  // f(x,y,z) = x^2*y^2*z^2
    {
      basef_surf = x2 * (6.0 * a12 * z5 + (15.0 * a1 * a2 * x + 15.0 * a0 * a1) * z4 +
                            (10.0 * a22 * x2 + 20.0 * a0 * a2 * x + 10.0 * a02) * z3);
      return basef_surf / 30.0;
    }
    if (inte_num == 70)  // f(x,y,z) = x^2*yz^3
    {
      basef_surf = x2 * (4.0 * a1 * z5 + (5.0 * a2 * x + 5.0 * a0) * z4);
      return basef_surf * 0.05;
    }
    if (inte_num == 71)  // f(x,y,z) = x^2*z^4
    {
      basef_surf = x2 * z5 * 0.2;
      return basef_surf;
    }
    if (inte_num == 72)  // f(x,y,z) = xy^5
    {
      basef_surf =
          x *
          (a15 * z6 / 6.0 +
              5.0 * a0 *
                  (a14 * z5 * 0.2 + a13 * a2 * x * z4 + 2.0 * a12 * a22 * x2 * z3 +
                      2.0 * a1 * a23 * x3 * z2 + a24 * x4 * z) +
              a14 * a2 * x * z5 +
              10.0 * a02 *
                  (0.25 * a13 * z4 + a12 * a2 * x * z3 + 1.5 * a1 * a22 * x2 * z2 + a23 * x3 * z) +
              2.5 * a13 * a22 * x2 * z4 +
              10.0 * a03 * (a12 * z3 / 3.0 + a1 * a2 * x * z2 + a22 * x2 * z) +
              10.0 / 3.0 * a12 * a23 * x3 * z3 + 5.0 * a04 * (0.5 * a1 * z2 + a2 * x * z) +
              2.5 * a1 * a24 * x4 * z2 + a25 * x5 * z + a05 * z);
      return basef_surf;
    }
    if (inte_num == 73)  // f(x,y,z) = xy^4*z
    {
      basef_surf =
          x * (5.0 * a14 * z6 + (24.0 * a13 * a2 * x + 24.0 * a0 * a13) * z5 +
                  (45.0 * a12 * a22 * x2 + 90.0 * a0 * a12 * a2 * x + 45.0 * a02 * a12) * z4 +
                  (40.0 * a1 * a23 * x3 + 120.0 * a0 * a1 * a22 * x2 + 120.0 * a02 * a1 * a2 * x +
                      40.0 * a03 * a1) *
                      z3 +
                  (15.0 * a24 * x4 + 60.0 * a0 * a23 * x3 + 90.0 * a02 * a22 * x2 +
                      60.0 * a03 * a2 * x + 15.0 * a04) *
                      z2);
      return basef_surf / 30.0;
    }
    if (inte_num == 74)  // f(x,y,z) = xy^3*z^2
    {
      basef_surf =
          x * (10.0 * a13 * z6 + (36.0 * a12 * a2 * x + 36.0 * a0 * a12) * z5 +
                  (45.0 * a1 * a22 * x2 + 90.0 * a0 * a1 * a2 * x + 45.0 * a02 * a1) * z4 +
                  (20.0 * a23 * x3 + 60.0 * a0 * a22 * x2 + 60.0 * a02 * a2 * x + 20.0 * a03) * z3);
      return basef_surf / 60.0;
    }
    if (inte_num == 75)  // f(x,y,z) = xy^2*z^3
    {
      basef_surf = x * (10.0 * a12 * z6 + (24.0 * a1 * a2 * x + 24.0 * a0 * a1) * z5 +
                           (15.0 * a22 * x2 + 30.0 * a0 * a2 * x + 15.0 * a02) * z4);
      return basef_surf / 60.0;
    }
    if (inte_num == 76)  // f(x,y,z) = xyz^4
    {
      basef_surf = x * (5.0 * a1 * z6 + (6.0 * a2 * x + 6.0 * a0) * z5);
      return basef_surf / 30.0;
    }
    if (inte_num == 77)  // f(x,y,z) = xz^5
    {
      basef_surf = x * z6 / 6.0;
      return basef_surf;
    }
    if (inte_num == 78)  // f(x,y,z) = y^6
    {
      if (fabs(a1) > 0.00000001)
        basef_surf = std::pow((a0 + a2 * x + a1 * z), 7) / 7.0 / a1;
      else
        basef_surf = std::pow((a0 + a1 * x), 6) * z;
      return basef_surf;
    }
    if (inte_num == 79)  // f(x,y,z) = y^5*z
    {
      basef_surf = 6.0 * a15 * z7 + (35.0 * a14 * a2 * x + 35.0 * a0 * a14) * z6 +
                   (84.0 * a13 * a22 * x2 + 168.0 * a0 * a13 * a2 * x + 84.0 * a02 * a13) * z5 +
                   (105.0 * a12 * a23 * x3 + 315.0 * a0 * a12 * a22 * x2 +
                       315.0 * a02 * a12 * a2 * x + 105.0 * a03 * a12) *
                       z4 +
                   (70.0 * a1 * a24 * x4 + 280.0 * a0 * a1 * a23 * x3 +
                       420.0 * a02 * a1 * a22 * x2 + 280.0 * a03 * a1 * a2 * x + 70.0 * a04 * a1) *
                       z3 +
                   (21.0 * a25 * x5 + 105.0 * a0 * a24 * x4 + 210.0 * a02 * a23 * x3 +
                       210.0 * a03 * a22 * x2 + 105.0 * a04 * a2 * x + 21.0 * a05) *
                       z2;
      return basef_surf / 42.0;
    }
    if (inte_num == 80)  // f(x,y,z) = y^4*z^2
    {
      basef_surf = 15.0 * a14 * z7 + (70.0 * a13 * a2 * x + 70.0 * a0 * a13) * z6 +
                   (126.0 * a12 * a22 * x2 + 252.0 * a0 * a12 * a2 * x + 126.0 * a02 * a12) * z5 +
                   (105.0 * a1 * a23 * x3 + 315.0 * a0 * a1 * a22 * x2 + 315.0 * a02 * a1 * a2 * x +
                       105.0 * a03 * a1) *
                       z4 +
                   (35.0 * a24 * x4 + 140.0 * a0 * a23 * x3 + 210.0 * a02 * a22 * x2 +
                       140.0 * a03 * a2 * x + 35.0 * a04) *
                       z3;
      return basef_surf / 105.0;
    }
    if (inte_num == 81)  // f(x,y,z) = y^3*z^3
    {
      basef_surf =
          20.0 * a13 * z7 + (70.0 * a12 * a2 * x + 70.0 * a0 * a12) * z6 +
          (84.0 * a1 * a22 * x2 + 168.0 * a0 * a1 * a2 * x + 84.0 * a02 * a1) * z5 +
          (35.0 * a23 * x3 + 105.0 * a0 * a22 * x2 + 105.0 * a02 * a2 * x + 35.0 * a03) * z4;
      return basef_surf / 140.0;
    }
    if (inte_num == 82)  // f(x,y,z) = y^2*z^4
    {
      basef_surf = 15.0 * a12 * z7 + (35.0 * a1 * a2 * x + 35.0 * a0 * a1) * z6 +
                   (21.0 * a22 * x2 + 42.0 * a0 * a2 * x + 21.0 * a02) * z5;
      return basef_surf / 105.0;
    }
    if (inte_num == 83)  // f(x,y,z) = yz^5
    {
      basef_surf = 6.0 * a1 * z7 + (7.0 * a2 * x + 7.0 * a0) * z6;
      return basef_surf / 42.0;
    }
    if (inte_num == 84)  // f(x,y,z) = z^6
    {
      basef_surf = z7 / 7.0;
      return basef_surf;
    }
    FOUR_C_THROW("The base function for boundarycell integration undefined");
    exit(1);
    return 0.0;
  }

  /*!
  \brief Returns the base function when the boundaycell is projected in x-y plane
  */
  // pt(0,0) --> x
  // pt(1,0) --> y
  double base_func_surf_z(
      const Core::LinAlg::Matrix<2, 1> &pt, int inte_num, std::vector<double> alfa)
  {
    double basef_surf = 0.0;
    if (inte_num == 1)  // f(x,y,z) = 1
    {
      basef_surf = pt(0, 0);
      return basef_surf;
    }

    double a0 = alfa[0], a1 = alfa[1], a2 = alfa[2];
    double x = pt(0, 0), y = pt(1, 0), x2 = std::pow(pt(0, 0), 2);
    if (inte_num == 2)  // f(x,y,z) = x
    {
      basef_surf = 0.5 * x * x;
      return basef_surf;
    }
    if (inte_num == 3)  // f(x,y,z) = y
    {
      basef_surf = y * x;
      return basef_surf;
    }
    if (inte_num == 4)  // f(x,y,z) = z
    {
      basef_surf = a2 * x * y + a1 * x2 * 0.5 + a0 * x;
      return basef_surf;
    }

    double a22 = std::pow(a2, 2), a02 = std::pow(a0, 2), a12 = std::pow(a1, 2);
    double x3 = std::pow(x, 3), y2 = std::pow(y, 2);
    if (inte_num == 5)  // f(x,y,z) = x^2
    {
      basef_surf = x3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 6)  // f(x,y,z) = xy
    {
      basef_surf = 0.5 * x2 * y;
      return basef_surf;
    }
    if (inte_num == 7)  // f(x,y,z) = xz
    {
      basef_surf = x2 * (3 * a2 * y + 3 * a0) + 2 * a1 * x3;
      return basef_surf / 6.0;
    }
    if (inte_num == 8)  // f(x,y,z) = y^2
    {
      basef_surf = x * y2;
      return basef_surf;
    }
    if (inte_num == 9)  // f(x,y,z) = yz
    {
      basef_surf = y * (a2 * x * y + 0.5 * a1 * x2 + a0 * x);
      return basef_surf;
    }
    if (inte_num == 10)  // f(x,y,z) = z^2
    {
      basef_surf = a22 * x * y2 + 2.0 * a0 * (a2 * x * y + 0.5 * a1 * x2) + a1 * a2 * x2 * y +
                   a12 * x3 / 3.0 + a02 * x;
      return basef_surf;
    }

    double x4 = std::pow(x, 4), y3 = std::pow(y, 3);
    double a03 = std::pow(a0, 3), a13 = std::pow(a1, 3), a23 = std::pow(a2, 3);
    if (inte_num == 11)  // f(x,y,z) = x^3
    {
      basef_surf = 0.25 * x4;
      return basef_surf;
    }
    if (inte_num == 12)  // f(x,y,z) = x^2*y
    {
      basef_surf = x3 * y / 3.0;
      return basef_surf;
    }
    if (inte_num == 13)  // f(x,y,z) = x^2*z
    {
      basef_surf = x3 * (4.0 * a2 * y + 4.0 * a0) + 3.0 * a1 * x4;
      return basef_surf / 12.0;
    }
    if (inte_num == 14)  // f(x,y,z) = xy^2
    {
      basef_surf = 0.5 * x2 * y2;
      return basef_surf;
    }
    if (inte_num == 15)  // f(x,y,z) = xyz
    {
      basef_surf = y * (x2 * (3.0 * a2 * y + 3.0 * a0) + 2.0 * a1 * x3);
      return basef_surf / 6.0;
    }
    if (inte_num == 16)  // f(x,y,z) = xz^2
    {
      basef_surf = x2 * (6.0 * a22 * y2 + 12.0 * a0 * a2 * y + 6.0 * a02) +
                   x3 * (8.0 * a1 * a2 * y + 8.0 * a0 * a1) + 3.0 * a12 * x4;
      return basef_surf / 12.0;
    }
    if (inte_num == 17)  // f(x,y,z) = y^3
    {
      basef_surf = y3 * x;
      return basef_surf;
    }
    if (inte_num == 18)  // f(x,y,z) = y^2*z
    {
      basef_surf = y2 * (a2 * x * y + 0.5 * a1 * x2 + a0 * x);
      return basef_surf;
    }
    if (inte_num == 19)  // f(x,y,z) = yz^2
    {
      basef_surf = y * (a22 * x * y2 + 2.0 * a0 * (a2 * x * y + 0.5 * a1 * x2) + a1 * a2 * x2 * y +
                           a12 * x3 / 3.0 + a02 * x);
      return basef_surf;
    }
    if (inte_num == 20)  // f(x,y,z) = z^3
    {
      basef_surf = a23 * x * y3 + 3.0 * a0 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
                   1.5 * a1 * a22 * x2 * y2 + 3.0 * a02 * (a2 * x * y + 0.5 * a1 * x2) +
                   a12 * a2 * x3 * y + 0.25 * a13 * x4 + a03 * x;
      return basef_surf;
    }

    double x5 = std::pow(x, 5), y4 = std::pow(y, 4);
    double a04 = std::pow(a0, 4), a14 = std::pow(a1, 4), a24 = std::pow(a2, 4);
    if (inte_num == 21)  // f(x,y,z) = x^4
    {
      basef_surf = 0.2 * x5;
      return basef_surf;
    }
    if (inte_num == 22)  // f(x,y,z) = x^3*y
    {
      basef_surf = 0.25 * x4 * y;
      return basef_surf;
    }
    if (inte_num == 23)  // f(x,y,z) = x^3*z
    {
      basef_surf = x4 * (5.0 * a2 * y + 5.0 * a0) + 4.0 * a1 * x5;
      return basef_surf * 0.05;
    }
    if (inte_num == 24)  // f(x,y,z) = x^2*y^2
    {
      basef_surf = x3 * y2 / 3.0;
      return basef_surf;
    }
    if (inte_num == 25)  // f(x,y,z) = x^2*yz
    {
      basef_surf = y * (x3 * (4.0 * a2 * y + 4.0 * a0) + 3.0 * a1 * x4);
      return basef_surf / 12.0;
    }
    if (inte_num == 26)  // f(x,y,z) = x^2*z^2
    {
      basef_surf = x3 * (10.0 * a22 * y2 + 20.0 * a0 * a2 * y + 10.0 * a02) +
                   x4 * (15.0 * a1 * a2 * y + 15.0 * a0 * a1) + 6.0 * a12 * x5;
      return basef_surf / 30.0;
    }
    if (inte_num == 27)  // f(x,y,z) = x*y^3
    {
      basef_surf = x2 * y3 * 0.5;
      return basef_surf;
    }
    if (inte_num == 28)  // f(x,y,z) = xy^2*z
    {
      basef_surf = y2 * (x2 * (3.0 * a2 * y + 3.0 * a0) + 2.0 * a1 * x3);
      return basef_surf / 6.0;
    }
    if (inte_num == 29)  // f(x,y,z) = xyz^2
    {
      basef_surf = y * (x2 * (6.0 * a22 * y2 + 12.0 * a0 * a2 * y + 6.0 * a02) +
                           x3 * (8.0 * a1 * a2 * y + 8.0 * a0 * a1) + 3.0 * a12 * x4);
      return basef_surf / 12.0;
    }
    if (inte_num == 30)  // f(x,y,z) = xz^3
    {
      basef_surf =
          x2 * (10.0 * a23 * y3 + 30.0 * a0 * a22 * y2 + 30.0 * a02 * a2 * y + 10.0 * a03) +
          x3 * (20.0 * a1 * a22 * y2 + 40.0 * a0 * a1 * a2 * y + 20.0 * a02 * a1) +
          x4 * (15.0 * a12 * a2 * y + 15.0 * a0 * a12) + 4.0 * a13 * x5;
      return basef_surf * 0.05;
    }
    if (inte_num == 31)  // f(x,y,z) = y^4
    {
      basef_surf = x * y4;
      return basef_surf;
    }
    if (inte_num == 32)  // f(x,y,z) = y^3*z
    {
      basef_surf = y3 * (a2 * x * y + 0.5 * a1 * x2 + a0 * x);
      return basef_surf;
    }
    if (inte_num == 33)  // f(x,y,z) = y^2*z^2
    {
      basef_surf = y2 * (a22 * x * y2 + 2.0 * a0 * (a2 * x * y + 0.5 * a1 * x2) + a1 * a2 * x2 * y +
                            a12 * x3 / 3.0 + a02 * x);
      return basef_surf;
    }
    if (inte_num == 34)  // f(x,y,z) = yz^3
    {
      basef_surf =
          y * (a23 * x * y3 + 3.0 * a0 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
                  1.5 * a1 * a22 * x2 * y2 + 3.0 * a02 * (a2 * x * y + 0.5 * a1 * x2) +
                  a12 * a2 * x3 * y + 0.25 * a13 * x4 + a03 * x);
      return basef_surf;
    }
    if (inte_num == 35)  // f(x,y,z) = z^4
    {
      basef_surf =
          a24 * x * y4 +
          4.0 * a0 *
              (a23 * x * y3 + 1.5 * a1 * a22 * x2 * y2 + a12 * a2 * x3 * y + 0.25 * a13 * x4) +
          2.0 * a1 * a23 * x2 * y3 +
          6.0 * a02 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
          2.0 * a12 * a22 * x3 * y2 + 4.0 * a03 * (a2 * x * y + 0.5 * a1 * x2) + a13 * a2 * x4 * y +
          0.2 * a14 * x5 + a04 * x;
      return basef_surf;
    }

    double x6 = std::pow(x, 6), y5 = std::pow(y, 5);
    double a05 = std::pow(a0, 5), a15 = std::pow(a1, 5), a25 = std::pow(a2, 5);
    if (inte_num == 36)  // f(x,y,z) = x^5
    {
      basef_surf = x6 / 6.0;
      return basef_surf;
    }
    if (inte_num == 37)  // f(x,y,z) = x^4*y
    {
      basef_surf = 0.2 * x5 * y;
      return basef_surf;
    }
    if (inte_num == 38)  // f(x,y,z) = x^4*z
    {
      basef_surf = x5 * (6.0 * a2 * y + 6.0 * a0) + 5.0 * a1 * x6;
      return basef_surf / 30.0;
    }
    if (inte_num == 39)  // f(x,y,z) = x^3*y^2
    {
      basef_surf = 0.25 * x4 * y2;
      return basef_surf;
    }
    if (inte_num == 40)  // f(x,y,z) = x^3*yz
    {
      basef_surf = y * (x4 * (5.0 * a2 * y + 5.0 * a0) + 4.0 * a1 * x5);
      return basef_surf * 0.05;
    }
    if (inte_num == 41)  // f(x,y,z) = x^3*z^2
    {
      basef_surf = x4 * (15.0 * a22 * y2 + 30.0 * a0 * a2 * y + 15.0 * a02) +
                   x5 * (24.0 * a1 * a2 * y + 24.0 * a0 * a1) + 10.0 * a12 * x6;
      return basef_surf / 60.0;
    }
    if (inte_num == 42)  // f(x,y,z) = x^2*y^3
    {
      basef_surf = x3 * y3 / 3.0;
      return basef_surf;
    }
    if (inte_num == 43)  // f(x,y,z) = x^2*y^2*z
    {
      basef_surf = y2 * (x3 * (4.0 * a2 * y + 4.0 * a0) + 3.0 * a1 * x4);
      return basef_surf / 12.0;
    }
    if (inte_num == 44)  // f(x,y,z) = x^2*yz^2
    {
      basef_surf = y * (x3 * (10.0 * a22 * y2 + 20.0 * a0 * a2 * y + 10.0 * a02) +
                           x4 * (15.0 * a1 * a2 * y + 15.0 * a0 * a1) + 6.0 * a12 * x5);
      return basef_surf / 30.0;
    }
    if (inte_num == 45)  // f(x,y,z) = x^2*z^3
    {
      basef_surf =
          x3 * (20.0 * a23 * y3 + 60.0 * a0 * a22 * y2 + 60.0 * a02 * a2 * y + 20.0 * a03) +
          x4 * (45.0 * a1 * a22 * y2 + 90.0 * a0 * a1 * a2 * y + 45.0 * a02 * a1) +
          x5 * (36.0 * a12 * a2 * y + 36.0 * a0 * a12) + 10.0 * a13 * x6;
      return basef_surf / 60.0;
    }
    if (inte_num == 46)  // f(x,y,z) = xy^4
    {
      basef_surf = x2 * y4 * 0.5;
      return basef_surf;
    }
    if (inte_num == 47)  // f(x,y,z) = xy^3*z
    {
      basef_surf = y3 * (x2 * (3.0 * a2 * y + 3.0 * a0) + 2.0 * a1 * x3);
      return basef_surf / 6.0;
    }
    if (inte_num == 48)  // f(x,y,z) = xy^2*z^2
    {
      basef_surf = y2 * (x2 * (6.0 * a22 * y2 + 12.0 * a0 * a2 * y + 6.0 * a02) +
                            x3 * (8.0 * a1 * a2 * y + 8.0 * a0 * a1) + 3.0 * a12 * x4);
      return basef_surf / 12.0;
    }
    if (inte_num == 49)  // f(x,y,z) = xyz^3
    {
      basef_surf =
          y * (x2 * (10.0 * a23 * y3 + 30.0 * a0 * a22 * y2 + 30.0 * a02 * a2 * y + 10.0 * a03) +
                  x3 * (20.0 * a1 * a22 * y2 + 40.0 * a0 * a1 * a2 * y + 20.0 * a02 * a1) +
                  x4 * (15.0 * a12 * a2 * y + 15.0 * a0 * a12) + 4.0 * a13 * x5);
      return basef_surf * 0.05;
    }
    if (inte_num == 50)  // f(x,y,z) = xz^4
    {
      basef_surf = x2 * (15.0 * a24 * y4 + 60.0 * a0 * a23 * y3 + 90.0 * a02 * a22 * y2 +
                            60.0 * a03 * a2 * y + 15.0 * a04) +
                   x3 * (40.0 * a1 * a23 * y3 + 120.0 * a0 * a1 * a22 * y2 +
                            120.0 * a02 * a1 * a2 * y + 40.0 * a03 * a1) +
                   x4 * (45.0 * a12 * a22 * y2 + 90.0 * a0 * a12 * a2 * y + 45.0 * a02 * a12) +
                   x5 * (24.0 * a13 * a2 * y + 24.0 * a0 * a13) + 5.0 * a14 * x6;
      return basef_surf / 30.0;
    }
    if (inte_num == 51)  // f(x,y,z) = y^5
    {
      basef_surf = y5 * x;
      return basef_surf;
    }
    if (inte_num == 52)  // f(x,y,z) = y^4*z
    {
      basef_surf = y4 * (a2 * x * y + 0.5 * a1 * x2 + a0 * x);
      return basef_surf;
    }
    if (inte_num == 53)  // f(x,y,z) = y^3*z^2
    {
      basef_surf = y3 * (a22 * x * y2 + 2.0 * a0 * (a2 * x * y + 0.5 * a1 * x2) + a1 * a2 * x2 * y +
                            a12 * x3 / 3.0 + a02 * x);
      return basef_surf;
    }
    if (inte_num == 54)  // f(x,y,z) = y^2*z^3
    {
      basef_surf =
          y2 * (a23 * x * y3 + 3.0 * a0 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
                   1.5 * a1 * a22 * x2 * y2 + 3.0 * a02 * (a2 * x * y + 0.5 * a1 * x2) +
                   a12 * a2 * x3 * y + 0.25 * a13 * x4 + a03 * x);
      return basef_surf;
    }
    if (inte_num == 55)  // f(x,y,z) = yz^4
    {
      basef_surf =
          y *
          (a24 * x * y4 +
              4.0 * a0 *
                  (a23 * x * y3 + 1.5 * a1 * a22 * x2 * y2 + a12 * a2 * x3 * y + 0.25 * a13 * x4) +
              2.0 * a1 * a23 * x2 * y3 +
              6.0 * a02 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
              2.0 * a12 * a22 * x3 * y2 + 4.0 * a03 * (a2 * x * y + 0.5 * a1 * x2) +
              a13 * a2 * x4 * y + 0.2 * a14 * x5 + a04 * x);
      return basef_surf;
    }
    if (inte_num == 56)  // f(x,y,z) = z^5
    {
      basef_surf =
          a25 * x * y5 +
          5.0 * a0 *
              (a24 * x * y4 + 2.0 * a1 * a23 * x2 * y3 + 2.0 * a12 * a22 * x3 * y2 +
                  a13 * a2 * x4 * y + a14 * x5 * 0.2) +
          2.5 * a1 * a24 * x2 * y4 +
          10.0 * a02 *
              (a23 * x * y3 + 1.5 * a1 * a22 * x2 * y2 + a12 * a2 * x3 * y + 0.25 * a13 * x4) +
          10.0 / 3.0 * a12 * a23 * x3 * y3 +
          10.0 * a03 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
          2.5 * a13 * a22 * x4 * y2 + 5.0 * a04 * (a2 * x * y + 0.5 * a1 * x2) + a14 * a2 * x5 * y +
          a15 * x6 / 6.0 + a05 * x;
      return basef_surf;
    }

    double x7 = std::pow(x, 7), y6 = std::pow(y, 6);
    if (inte_num == 57)  // f(x,y,z) = x^6
    {
      basef_surf = x7 / 7.0;
      return basef_surf;
    }
    if (inte_num == 58)  // f(x,y,z) = x^5*y
    {
      basef_surf = x6 * y / 6.0;
      return basef_surf;
    }
    if (inte_num == 59)  // f(x,y,z) = x^5*z
    {
      basef_surf = x6 * (7.0 * a2 * y + 7.0 * a0) + 6.0 * a1 * x7;
      return basef_surf / 42.0;
    }
    if (inte_num == 60)  // f(x,y,z) = x^4*y^2
    {
      basef_surf = 0.2 * x5 * y2;
      return basef_surf;
    }
    if (inte_num == 61)  // f(x,y,z) = x^4*yz
    {
      basef_surf = y * (x5 * (6.0 * a2 * y + 6.0 * a0) + 5.0 * a1 * x6);
      return basef_surf / 30.0;
    }
    if (inte_num == 62)  // f(x,y,z) = x^4*z^2
    {
      basef_surf = x5 * (21.0 * a22 * y2 + 42.0 * a0 * a2 * y + 21.0 * a02) +
                   x6 * (35.0 * a1 * a2 * y + 35.0 * a0 * a1) + 15.0 * a12 * x7;
      return basef_surf / 105.0;
    }
    if (inte_num == 63)  // f(x,y,z) = x^3*y^3
    {
      basef_surf = 0.25 * x4 * y3;
      return basef_surf;
    }
    if (inte_num == 64)  // f(x,y,z) = x^3*y^2*z
    {
      basef_surf = y2 * (x4 * (5.0 * a2 * y + 5.0 * a0) + 4.0 * a1 * x5);
      return basef_surf * 0.05;
    }
    if (inte_num == 65)  // f(x,y,z) = x^3*y*z^2
    {
      basef_surf = y * (x4 * (15.0 * a22 * y2 + 30.0 * a0 * a2 * y + 15.0 * a02) +
                           x5 * (24.0 * a1 * a2 * y + 24.0 * a0 * a1) + 10.0 * a12 * x6);
      return basef_surf / 60.0;
    }
    if (inte_num == 66)  // f(x,y,z) = x^3*z^3
    {
      basef_surf =
          x4 * (35.0 * a23 * y3 + 105.0 * a0 * a22 * y2 + 105.0 * a02 * a2 * y + 35.0 * a03) +
          x5 * (84.0 * a1 * a22 * y2 + 168.0 * a0 * a1 * a2 * y + 84.0 * a02 * a1) +
          x6 * (70.0 * a12 * a2 * y + 70.0 * a0 * a12) + 20.0 * a13 * x7;
      return basef_surf / 140.0;
    }
    if (inte_num == 67)  // f(x,y,z) = x^2*y^4
    {
      basef_surf = x3 * y4 / 3.0;
      return basef_surf;
    }
    if (inte_num == 68)  // f(x,y,z) = x^2*y^3*z
    {
      basef_surf = 4.0 * a2 * x3 * y4 + (3.0 * a1 * x4 + 4.0 * a0 * x3) * y3;
      return basef_surf / 12.0;
    }
    if (inte_num == 69)  // f(x,y,z) = x^2*y^2*z^2
    {
      basef_surf = y2 * (x3 * (10.0 * a22 * y2 + 20.0 * a0 * a2 * y + 10.0 * a02) +
                            x4 * (15.0 * a1 * a2 * y + 15.0 * a0 * a1) + 6.0 * a12 * x5);
      return basef_surf / 30.0;
    }
    if (inte_num == 70)  // f(x,y,z) = x^2*yz^3
    {
      basef_surf =
          y * (x3 * (20.0 * a23 * y3 + 60.0 * a0 * a22 * y2 + 60.0 * a02 * a2 * y + 20.0 * a03) +
                  x4 * (45.0 * a1 * a22 * y2 + 90.0 * a0 * a1 * a2 * y + 45.0 * a02 * a1) +
                  x5 * (36.0 * a12 * a2 * y + 36.0 * a0 * a12) + 10.0 * a13 * x6);
      return basef_surf / 60.0;
    }
    if (inte_num == 71)  // f(x,y,z) = x^2*z^4
    {
      basef_surf = x3 * (35.0 * a24 * y4 + 140.0 * a0 * a23 * y3 + 210.0 * a02 * a22 * y2 +
                            140.0 * a03 * a2 * y + 35.0 * a04) +
                   x4 * (105.0 * a1 * a23 * y3 + 315.0 * a0 * a1 * a22 * y2 +
                            315.0 * a02 * a1 * a2 * y + 105.0 * a03 * a1) +
                   x5 * (126.0 * a12 * a22 * y2 + 252.0 * a0 * a12 * a2 * y + 126.0 * a02 * a12) +
                   x6 * (70.0 * a13 * a2 * y + 70.0 * a0 * a13) + 15.0 * a14 * x7;
      return basef_surf / 105.0;
    }
    if (inte_num == 72)  // f(x,y,z) = xy^5
    {
      basef_surf = 0.5 * x2 * y5;
      return basef_surf;
    }
    if (inte_num == 73)  // f(x,y,z) = xy^4*z
    {
      basef_surf = y4 * (x2 * (3.0 * a2 * y + 3.0 * a0) + 2.0 * a1 * x3);
      return basef_surf / 6.0;
    }
    if (inte_num == 74)  // f(x,y,z) = xy^3*z^2
    {
      basef_surf = y3 * (x2 * (6.0 * a22 * y2 + 12.0 * a0 * a2 * y + 6.0 * a02) +
                            x3 * (8.0 * a1 * a2 * y + 8.0 * a0 * a1) + 3.0 * a12 * x4);
      return basef_surf / 12.0;
    }
    if (inte_num == 75)  // f(x,y,z) = xy^2*z^3
    {
      basef_surf =
          y2 * (x2 * (10.0 * a23 * y3 + 30.0 * a0 * a22 * y2 + 30.0 * a02 * a2 * y + 10.0 * a03) +
                   x3 * (20.0 * a1 * a22 * y2 + 40.0 * a0 * a1 * a2 * y + 20.0 * a02 * a1) +
                   x4 * (15.0 * a12 * a2 * y + 15.0 * a0 * a12) + 4.0 * a13 * x5);
      return basef_surf * 0.05;
    }
    if (inte_num == 76)  // f(x,y,z) = xyz^4
    {
      basef_surf =
          y * (x2 * (15.0 * a24 * y4 + 60.0 * a0 * a23 * y3 + 90.0 * a02 * a22 * y2 +
                        60.0 * a03 * a2 * y + 15.0 * a04) +
                  x3 * (40.0 * a1 * a23 * y3 + 120.0 * a0 * a1 * a22 * y2 +
                           120.0 * a02 * a1 * a2 * y + 40.0 * a03 * a1) +
                  x4 * (45.0 * a12 * a22 * y2 + 90.0 * a0 * a12 * a2 * y + 45.0 * a02 * a12) +
                  x5 * (24.0 * a13 * a2 * y + 24.0 * a0 * a13) + 5.0 * a14 * x6);
      return basef_surf / 30.0;
    }
    if (inte_num == 77)  // f(x,y,z) = xz^5
    {
      basef_surf =
          x2 * (21.0 * a25 * y5 + 105.0 * a0 * a24 * y4 + 210.0 * a02 * a23 * y3 +
                   210.0 * a03 * a22 * y2 + 105.0 * a04 * a2 * y + 21.0 * a05) +
          x3 * (70.0 * a1 * a24 * y4 + 280.0 * a0 * a1 * a23 * y3 + 420.0 * a02 * a1 * a22 * y2 +
                   280.0 * a03 * a1 * a2 * y + 70.0 * a04 * a1) +
          x4 * (105.0 * a12 * a23 * y3 + 315.0 * a0 * a12 * a22 * y2 + 315.0 * a02 * a12 * a2 * y +
                   105.0 * a03 * a12) +
          x5 * (84.0 * a13 * a22 * y2 + 168.0 * a0 * a13 * a2 * y + 84.0 * a02 * a13) +
          x6 * (35.0 * a14 * a2 * y + 35.0 * a0 * a14) + 6.0 * a15 * x7;
      return basef_surf / 42.0;
    }
    if (inte_num == 78)  // f(x,y,z) = y^6
    {
      basef_surf = y6 * x;
      return basef_surf;
    }
    if (inte_num == 79)  // f(x,y,z) = y^5*z
    {
      basef_surf = y5 * (a2 * x * y + 0.5 * a1 * x2 + a0 * x);
      return basef_surf;
    }
    if (inte_num == 80)  // f(x,y,z) = y^4*z^2
    {
      basef_surf = y4 * (a22 * x * y2 + 2.0 * a0 * (a2 * x * y + 0.5 * a1 * x2) + a1 * a2 * x2 * y +
                            a12 * x3 / 3.0 + a02 * x);
      return basef_surf;
    }
    if (inte_num == 81)  // f(x,y,z) = y^3*z^3
    {
      basef_surf =
          y3 * (a23 * x * y3 + 3.0 * a0 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
                   1.5 * a1 * a22 * x2 * y2 + 3.0 * a02 * (a2 * x * y + 0.5 * a1 * x2) +
                   a12 * a2 * x3 * y + 0.25 * a13 * x4 + a03 * x);
      return basef_surf;
    }
    if (inte_num == 82)  // f(x,y,z) = y^2*z^4
    {
      basef_surf =
          y2 *
          (a24 * x * y4 +
              4.0 * a0 *
                  (a23 * x * y3 + 1.5 * a1 * a22 * x2 * y2 + a12 * a2 * x3 * y + 0.25 * a13 * x4) +
              2.0 * a1 * a23 * x2 * y3 +
              6.0 * a02 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
              2.0 * a12 * a22 * x3 * y2 + 4.0 * a03 * (a2 * x * y + 0.5 * a1 * x2) +
              a13 * a2 * x4 * y + 0.2 * a14 * x5 + a04 * x);
      return basef_surf;
    }
    if (inte_num == 83)  // f(x,y,z) = yz^5
    {
      basef_surf =
          y *
          (a25 * x * y5 +
              5.0 * a0 *
                  (a24 * x * y4 + 2.0 * a1 * a23 * x2 * y3 + 2.0 * a12 * a22 * x3 * y2 +
                      a13 * a2 * x4 * y + 0.2 * a14 * x5) +
              2.5 * a1 * a24 * x2 * y4 +
              10.0 * a02 *
                  (a23 * x * y3 + 1.5 * a1 * a22 * x2 * y2 + a12 * a2 * x3 * y + 0.25 * a13 * x4) +
              10.0 / 3.0 * a12 * a23 * x3 * y3 +
              10.0 * a03 * (a22 * x * y2 + a1 * a2 * x2 * y + a12 * x3 / 3.0) +
              2.5 * a13 * a22 * x4 * y2 + 5.0 * a04 * (a2 * x * y + 0.5 * a1 * x2) +
              a14 * a2 * x5 * y + a15 * x6 / 6.0 + a05 * x);
      return basef_surf;
    }
    if (inte_num == 84)  // f(x,y,z) = z^6
    {
      if (fabs(a1) > 0.00000001)
        basef_surf = std::pow((a0 + a1 * x + a2 * y), 7) / 7.0 / a1;
      else
        basef_surf = std::pow((a0 + a2 * y), 6) * x;
      return basef_surf;
    }
    FOUR_C_THROW("The base function for boundarycell integration undefined");
    exit(1);
    return 0.0;
  }
}  // namespace

/*---------------------------------------------------------------------------------------------------------------------*
 *      Compute normal vector for the line. if normal(0)==0, the line need not be integrated
 *(Divergence theorem)       *
 *----------------------------------------------------------------------------------------------------------------------*/
Core::LinAlg::Matrix<2, 1> LineIntegration::compute_normal()
{
  Core::LinAlg::Matrix<2, 1> normal;
  double dy = end_pts_(1, 1) - end_pts_(1, 0);
  double dx = -end_pts_(0, 1) + end_pts_(0, 0);
  double modd = sqrt(dx * dx + dy * dy);

  normal(0, 0) = dy / modd;
  normal(1, 0) = dx / modd;

  return normal;
}

/*-----------------------------------------------------------------------------*
 *                performs integration over the given line                      *
 *------------------------------------------------------------------------------*/
double LineIntegration::integrate_line()
{
  Core::LinAlg::Matrix<2, 1> normal;
  normal = compute_normal();

  if (fabs(normal(0, 0)) < TOL_LINE_NORMAL) return 0.0;

  double inte = 0.0;

  // 8 is the order of Gauss integration used in the line integration
  // since we integrate 6th order polynomial in volume, 8th order must be used for line
  Core::FE::GaussIntegration gi(Core::FE::CellType::line2, (DIRECTDIV_GAUSSRULE + 1));

  for (Core::FE::GaussIntegration::iterator iquad = gi.begin(); iquad != gi.end(); ++iquad)
  {
    const Core::LinAlg::Matrix<1, 1> eta(iquad.point());
    double weight = iquad.weight();
    Core::LinAlg::Matrix<2, 1> normaltemp, actCoord;
    double drs = 0.0;
    transform(end_pts_, eta(0, 0), actCoord, normaltemp, drs);

    if (bcell_int_ == false)  // integration over volumecell
    {
      double linein = base_func_line_int(actCoord, inte_num_, alpha_);
      inte = inte + weight * linein * drs;
    }
    else  // integration over boundarycell
    {
      double linein = 0.0;
      if (int_type_ == Cut::proj_x)
        linein = base_func_surf_x(actCoord, inte_num_, alpha_);
      else if (int_type_ == Cut::proj_y)
        linein = base_func_surf_y(actCoord, inte_num_, alpha_);
      else if (int_type_ == Cut::proj_z)
        linein = base_func_surf_z(actCoord, inte_num_, alpha_);
      else
        FOUR_C_THROW("Integration type unspecified");
      inte = inte + weight * linein * drs;
    }
  }
  inte = inte * normal(0, 0);

  return inte;
}

/*---------------------------------------------------------------------------------------------------------------------*
 *     Transform the Gaussian point coordinates and weight from (-1,1) interval to actual coordinate
 *of the lines       *
 *----------------------------------------------------------------------------------------------------------------------*/
void LineIntegration::transform(const Core::LinAlg::Matrix<2, 2> &xyze, const double &eta,
    Core::LinAlg::Matrix<2, 1> &x_gp_lin, Core::LinAlg::Matrix<2, 1> &normal, double &drs)
{
  const int numnodes = Core::FE::num_nodes<Core::FE::CellType::line2>;
  Core::LinAlg::Matrix<numnodes, 1> funct;
  Core::LinAlg::Matrix<1, numnodes> deriv;
  Core::LinAlg::Matrix<1, 1> metrictensor;

  Core::FE::shape_function_1d(funct, eta, Core::FE::CellType::line2);
  Core::FE::shape_function_1d_deriv1(deriv, eta, Core::FE::CellType::line2);
  Core::FE::compute_metric_tensor_for_boundary_ele<Core::FE::CellType::line2>(
      xyze, deriv, metrictensor, drs, &normal);

  x_gp_lin.multiply(xyze, funct);

  return;
}

FOUR_C_NAMESPACE_CLOSE
