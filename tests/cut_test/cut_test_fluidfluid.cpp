/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_fluidfluid.cpp

\level 1

\maintainer Christoph Ager
*----------------------------------------------------------------------*/

#include "cut_test_loader.H"

void test_fluidfluid()
{
  MeshLoader loader;

  loader.GetCutNode(96, -0.0883883, -0.276777, 0.0883883, 0);
  loader.GetCutNode(99, -0.0258883, -0.188388, 0.150888, 0);
  loader.GetCutNode(98, -0.0883883, -0.1, 0.0883883, 0);
  loader.GetCutNode(97, -0.150888, -0.188388, 0.0258883, 0);
  loader.CreateSide(0, 96, 99, 98, 97, DRT::Element::quad4);
  loader.GetCutNode(-1, -0.0883883, -0.188388, 0.0883883, 0);
  loader.GetCutNode(101, -0.0625, -0.188388, -0.0625, 0);
  loader.GetCutNode(100, 0, -0.276777, 0, 0);
  loader.CreateSide(1, 96, 97, 101, 100, DRT::Element::quad4);
  loader.GetCutNode(-2, -0.0754442, -0.232583, 0.0129442, 0);
  loader.GetCutNode(103, 0.0625, -0.188388, 0.0625, 0);
  loader.CreateSide(2, 96, 100, 103, 99, DRT::Element::quad4);
  loader.GetCutNode(-3, -0.0129442, -0.232583, 0.0754442, 0);
  loader.GetCutNode(105, -0.150888, -0.0116117, 0.0258883, 0);
  loader.GetCutNode(104, -0.213388, -0.1, -0.0366117, 0);
  loader.CreateSide(3, 97, 98, 105, 104, DRT::Element::quad4);
  loader.GetCutNode(-4, -0.150888, -0.1, 0.0258883, 0);
  loader.GetCutNode(106, -0.125, -0.1, -0.125, 0);
  loader.CreateSide(4, 97, 104, 106, 101, DRT::Element::quad4);
  loader.GetCutNode(-5, -0.137944, -0.144194, -0.0495558, 0);
  loader.GetCutNode(109, 0.0366117, -0.1, 0.213388, 0);
  loader.GetCutNode(108, -0.0258883, -0.0116117, 0.150888, 0);
  loader.CreateSide(5, 99, 109, 108, 98, DRT::Element::quad4);
  loader.GetCutNode(-6, -0.0258883, -0.1, 0.150888, 0);
  loader.GetCutNode(112, -0.0883883, 0.0767767, 0.0883883, 0);
  loader.CreateSide(6, 98, 108, 112, 105, DRT::Element::quad4);
  loader.GetCutNode(-7, -0.0883883, -0.0116117, 0.0883883, 0);
  loader.GetCutNode(111, 0.125, -0.1, 0.125, 0);
  loader.CreateSide(7, 99, 103, 111, 109, DRT::Element::quad4);
  loader.GetCutNode(-8, 0.0495558, -0.144194, 0.137944, 0);
  loader.GetCutNode(115, 0.0258883, -0.188388, -0.150888, 0);
  loader.GetCutNode(114, 0.0883883, -0.276777, -0.0883883, 0);
  loader.CreateSide(8, 100, 101, 115, 114, DRT::Element::quad4);
  loader.GetCutNode(-9, 0.0129442, -0.232583, -0.0754442, 0);
  loader.GetCutNode(117, 0.150888, -0.188388, -0.0258883, 0);
  loader.CreateSide(9, 100, 114, 117, 103, DRT::Element::quad4);
  loader.GetCutNode(-10, 0.0754442, -0.232583, -0.0129442, 0);
  loader.GetCutNode(118, -0.0366117, -0.1, -0.213388, 0);
  loader.CreateSide(10, 101, 106, 118, 115, DRT::Element::quad4);
  loader.GetCutNode(-11, -0.0495558, -0.144194, -0.137944, 0);
  loader.GetCutNode(121, 0.213388, -0.1, 0.0366117, 0);
  loader.CreateSide(11, 103, 117, 121, 111, DRT::Element::quad4);
  loader.GetCutNode(-12, 0.137944, -0.144194, 0.0495558, 0);
  loader.GetCutNode(107, -0.0625, -0.0116117, -0.0625, 0);
  loader.CreateSide(12, 104, 105, 107, 106, DRT::Element::quad4);
  loader.GetCutNode(-13, -0.137944, -0.0558058, -0.0495558, 0);
  loader.GetCutNode(113, 0, 0.0767767, 0, 0);
  loader.CreateSide(13, 105, 112, 113, 107, DRT::Element::quad4);
  loader.GetCutNode(-14, -0.0754442, 0.0325825, 0.0129442, 0);
  loader.GetCutNode(119, 0.0258883, -0.0116117, -0.150888, 0);
  loader.CreateSide(14, 106, 107, 119, 118, DRT::Element::quad4);
  loader.GetCutNode(-15, -0.0495558, -0.0558058, -0.137944, 0);
  loader.GetCutNode(122, 0.0883883, 0.0767767, -0.0883883, 0);
  loader.CreateSide(15, 107, 113, 122, 119, DRT::Element::quad4);
  loader.GetCutNode(-16, 0.0129442, 0.0325825, -0.0754442, 0);
  loader.GetCutNode(110, 0.0625, -0.0116117, 0.0625, 0);
  loader.CreateSide(16, 108, 109, 111, 110, DRT::Element::quad4);
  loader.GetCutNode(-17, 0.0495558, -0.0558058, 0.137944, 0);
  loader.CreateSide(17, 112, 108, 110, 113, DRT::Element::quad4);
  loader.GetCutNode(-18, -0.0129442, 0.0325825, 0.0754442, 0);
  loader.GetCutNode(120, 0.150888, -0.0116117, -0.0258883, 0);
  loader.CreateSide(18, 110, 111, 121, 120, DRT::Element::quad4);
  loader.GetCutNode(-19, 0.137944, -0.0558058, 0.0495558, 0);
  loader.CreateSide(19, 113, 110, 120, 122, DRT::Element::quad4);
  loader.GetCutNode(-20, 0.0754442, 0.0325825, -0.0129442, 0);
  loader.GetCutNode(116, 0.0883883, -0.1, -0.0883883, 0);
  loader.CreateSide(20, 114, 115, 116, 117, DRT::Element::quad4);
  loader.GetCutNode(-21, 0.0883883, -0.188388, -0.0883883, 0);
  loader.CreateSide(21, 115, 118, 119, 116, DRT::Element::quad4);
  loader.GetCutNode(-22, 0.0258883, -0.1, -0.150888, 0);
  loader.CreateSide(22, 117, 116, 120, 121, DRT::Element::quad4);
  loader.GetCutNode(-23, 0.150888, -0.1, -0.0258883, 0);
  loader.CreateSide(23, 116, 119, 122, 120, DRT::Element::quad4);
  loader.GetCutNode(-24, 0.0883883, -0.0116117, -0.0883883, 0);
  loader.GetNode(3, -0.25, -0.3, 0.25, 0);
  loader.GetNode(2, -0.25, -0.3, 0.0833333, 0);
  loader.GetNode(16, -0.25, -0.1, 0.0833333, 0);
  loader.GetNode(17, -0.25, -0.1, 0.25, 0);
  loader.GetNode(7, -0.0833333, -0.3, 0.25, 0);
  loader.GetNode(6, -0.0833333, -0.3, 0.0833333, 0);
  loader.GetNode(18, -0.0833333, -0.1, 0.0833333, 0);
  loader.GetNode(19, -0.0833333, -0.1, 0.25, 0);
  loader.CreateElement(11, 3, 2, 16, 17, 7, 6, 18, 19, DRT::Element::hex8);
  loader.GetNode(9, -0.25, -0.3, -0.0833333, 0);
  loader.GetNode(20, -0.25, -0.1, -0.0833333, 0);
  loader.GetNode(11, -0.0833333, -0.3, -0.0833333, 0);
  loader.GetNode(21, -0.0833333, -0.1, -0.0833333, 0);
  loader.CreateElement(12, 2, 9, 20, 16, 6, 11, 21, 18, DRT::Element::hex8);
  loader.GetNode(13, -0.25, -0.3, -0.25, 0);
  loader.GetNode(22, -0.25, -0.1, -0.25, 0);
  loader.GetNode(15, -0.0833333, -0.3, -0.25, 0);
  loader.GetNode(23, -0.0833333, -0.1, -0.25, 0);
  loader.CreateElement(13, 9, 13, 22, 20, 11, 15, 23, 21, DRT::Element::hex8);
  loader.GetNode(24, -0.25, 0.1, 0.0833333, 0);
  loader.GetNode(25, -0.25, 0.1, 0.25, 0);
  loader.GetNode(26, -0.0833333, 0.1, 0.0833333, 0);
  loader.GetNode(27, -0.0833333, 0.1, 0.25, 0);
  loader.CreateElement(14, 17, 16, 24, 25, 19, 18, 26, 27, DRT::Element::hex8);
  loader.GetNode(28, -0.25, 0.1, -0.0833333, 0);
  loader.GetNode(29, -0.0833333, 0.1, -0.0833333, 0);
  loader.CreateElement(15, 16, 20, 28, 24, 18, 21, 29, 26, DRT::Element::hex8);
  loader.GetNode(30, -0.25, 0.1, -0.25, 0);
  loader.GetNode(31, -0.0833333, 0.1, -0.25, 0);
  loader.CreateElement(16, 20, 22, 30, 28, 21, 23, 31, 29, DRT::Element::hex8);
  loader.GetNode(51, 0.0833333, -0.3, 0.25, 0);
  loader.GetNode(50, 0.0833333, -0.3, 0.0833333, 0);
  loader.GetNode(56, 0.0833333, -0.1, 0.0833333, 0);
  loader.GetNode(57, 0.0833333, -0.1, 0.25, 0);
  loader.CreateElement(26, 7, 6, 18, 19, 51, 50, 56, 57, DRT::Element::hex8);
  loader.GetNode(53, 0.0833333, -0.3, -0.0833333, 0);
  loader.GetNode(58, 0.0833333, -0.1, -0.0833333, 0);
  loader.CreateElement(27, 6, 11, 21, 18, 50, 53, 58, 56, DRT::Element::hex8);
  loader.GetNode(55, 0.0833333, -0.3, -0.25, 0);
  loader.GetNode(59, 0.0833333, -0.1, -0.25, 0);
  loader.CreateElement(28, 11, 15, 23, 21, 53, 55, 59, 58, DRT::Element::hex8);
  loader.GetNode(60, 0.0833333, 0.1, 0.0833333, 0);
  loader.GetNode(61, 0.0833333, 0.1, 0.25, 0);
  loader.CreateElement(29, 19, 18, 26, 27, 57, 56, 60, 61, DRT::Element::hex8);
  loader.GetNode(62, 0.0833333, 0.1, -0.0833333, 0);
  loader.CreateElement(30, 18, 21, 29, 26, 56, 58, 62, 60, DRT::Element::hex8);
  loader.GetNode(63, 0.0833333, 0.1, -0.25, 0);
  loader.CreateElement(31, 21, 23, 31, 29, 58, 59, 63, 62, DRT::Element::hex8);
  loader.GetNode(75, 0.25, -0.3, 0.25, 0);
  loader.GetNode(74, 0.25, -0.3, 0.0833333, 0);
  loader.GetNode(80, 0.25, -0.1, 0.0833333, 0);
  loader.GetNode(81, 0.25, -0.1, 0.25, 0);
  loader.CreateElement(41, 51, 50, 56, 57, 75, 74, 80, 81, DRT::Element::hex8);
  loader.GetNode(77, 0.25, -0.3, -0.0833333, 0);
  loader.GetNode(82, 0.25, -0.1, -0.0833333, 0);
  loader.CreateElement(42, 50, 53, 58, 56, 74, 77, 82, 80, DRT::Element::hex8);
  loader.GetNode(79, 0.25, -0.3, -0.25, 0);
  loader.GetNode(83, 0.25, -0.1, -0.25, 0);
  loader.CreateElement(43, 53, 55, 59, 58, 77, 79, 83, 82, DRT::Element::hex8);
  loader.GetNode(84, 0.25, 0.1, 0.0833333, 0);
  loader.GetNode(85, 0.25, 0.1, 0.25, 0);
  loader.CreateElement(44, 57, 56, 60, 61, 81, 80, 84, 85, DRT::Element::hex8);
  loader.GetNode(86, 0.25, 0.1, -0.0833333, 0);
  loader.CreateElement(45, 56, 58, 62, 60, 80, 82, 86, 84, DRT::Element::hex8);
  loader.GetNode(87, 0.25, 0.1, -0.25, 0);
  loader.CreateElement(46, 58, 59, 63, 62, 82, 83, 87, 86, DRT::Element::hex8);
  loader.GetNode(0, -0.0883883, -0.276777, 0.0883883, 0);
  loader.GetNode(1, -0.0883883, -0.1, 0.0883883, 0);
  loader.GetNode(2, -0.0883883, -0.188388, 0.0883883, 0);
  loader.GetNode(3, -0.25, -0.3, 0.25, 0);
  loader.GetNode(4, -0.25, -0.3, 0.0833333, 0);
  loader.GetNode(5, -0.25, -0.1, 0.0833333, 0);
  loader.GetNode(6, -0.25, -0.1, 0.25, 0);
  loader.GetNode(7, -0.0833333, -0.3, 0.25, 0);
  loader.GetNode(8, -0.0833333, -0.3, 0.0833333, 0);
  loader.GetNode(9, -0.0833333, -0.1, 0.0833333, 0);
  loader.GetNode(10, -0.0833333, -0.1, 0.25, 0);
  loader.GetNode(11, -0.0833333, -0.1, 0.0934434, 0);
  loader.GetNode(12, -0.0833333, -0.107149, 0.0934434, 0);
  loader.GetNode(13, -0.0934434, -0.1, 0.0833333, 0);
  loader.GetNode(14, -0.0934434, -0.107149, 0.0833333, 0);
  loader.GetNode(15, -0.0833333, -0.269628, 0.0934434, 0);
  loader.GetNode(16, -0.0833333, -0.273816, 0.087521, 0);
  loader.GetNode(17, -0.0833333, -0.276777, 0.0833333, 0);
  loader.GetNode(18, -0.0934434, -0.269628, 0.0833333, 0);
  loader.GetNode(19, -0.087521, -0.273816, 0.0833333, 0);
  loader.GetNode(20, -0.0934434, -0.188388, 0.0833333, 0);
  loader.GetNode(21, -0.0833333, -0.188388, 0.0934434, 0);
  loader.GetCutNode(96, -0.0883883, -0.276777, 0.0883883, 0);
  loader.GetCutNode(99, -0.0258883, -0.188388, 0.150888, 0);
  loader.GetCutNode(-1, -0.0883883, -0.188388, 0.0883883, 0);
  loader.GetCutNode(98, -0.0883883, -0.1, 0.0883883, 0);
  loader.GetCutNode(97, -0.150888, -0.188388, 0.0258883, 0);
  loader.GetCutNode(-2, -0.0754442, -0.232583, 0.0129442, 0);
  loader.GetCutNode(100, 0, -0.276777, 0, 0);
  loader.GetCutNode(-3, -0.0129442, -0.232583, 0.0754442, 0);
  loader.GetCutNode(-4, -0.150888, -0.1, 0.0258883, 0);
  loader.GetCutNode(-6, -0.0258883, -0.1, 0.150888, 0);
  loader.GetNode(0, -0.125, -0.1, -0.125, 0);
  loader.GetNode(1, -0.25, -0.3, -0.0833333, 0);
  loader.GetNode(2, -0.25, -0.1, -0.0833333, 0);
  loader.GetNode(3, -0.0833333, -0.3, -0.0833333, 0);
  loader.GetNode(4, -0.0833333, -0.1, -0.0833333, 0);
  loader.GetNode(5, -0.25, -0.3, -0.25, 0);
  loader.GetNode(6, -0.25, -0.1, -0.25, 0);
  loader.GetNode(7, -0.0833333, -0.3, -0.25, 0);
  loader.GetNode(8, -0.0833333, -0.1, -0.25, 0);
  loader.GetNode(9, -0.0833333, -0.1, -0.166667, 0);
  loader.GetNode(10, -0.166667, -0.1, -0.0833333, 0);
  loader.GetNode(11, -0.0833333, -0.158926, -0.0833333, 0);
  loader.GetNode(12, -0.0833333, -0.124408, -0.132149, 0);
  loader.GetNode(13, -0.132149, -0.124408, -0.0833333, 0);
  loader.GetCutNode(104, -0.213388, -0.1, -0.0366117, 0);
  loader.GetCutNode(106, -0.125, -0.1, -0.125, 0);
  loader.GetCutNode(-5, -0.137944, -0.144194, -0.0495558, 0);
  loader.GetCutNode(101, -0.0625, -0.188388, -0.0625, 0);
  loader.GetCutNode(-11, -0.0495558, -0.144194, -0.137944, 0);
  loader.GetCutNode(118, -0.0366117, -0.1, -0.213388, 0);
  loader.GetNode(0, -0.0883883, -0.1, 0.0883883, 0);
  loader.GetNode(1, -0.0883883, 0.0767767, 0.0883883, 0);
  loader.GetNode(2, -0.0883883, -0.0116117, 0.0883883, 0);
  loader.GetNode(3, -0.25, -0.1, 0.0833333, 0);
  loader.GetNode(4, -0.25, -0.1, 0.25, 0);
  loader.GetNode(5, -0.0833333, -0.1, 0.0833333, 0);
  loader.GetNode(6, -0.0833333, -0.1, 0.25, 0);
  loader.GetNode(7, -0.25, 0.1, 0.0833333, 0);
  loader.GetNode(8, -0.25, 0.1, 0.25, 0);
  loader.GetNode(9, -0.0833333, 0.1, 0.0833333, 0);
  loader.GetNode(10, -0.0833333, 0.1, 0.25, 0);
  loader.GetNode(11, -0.0833333, 0.0696278, 0.0934434, 0);
  loader.GetNode(12, -0.0833333, 0.0738155, 0.087521, 0);
  loader.GetNode(13, -0.0833333, 0.0767767, 0.0833333, 0);
  loader.GetNode(14, -0.0934434, 0.0696278, 0.0833333, 0);
  loader.GetNode(15, -0.087521, 0.0738155, 0.0833333, 0);
  loader.GetNode(16, -0.0934434, -0.0928511, 0.0833333, 0);
  loader.GetNode(17, -0.0934434, -0.0116117, 0.0833333, 0);
  loader.GetNode(18, -0.0833333, -0.0928511, 0.0934434, 0);
  loader.GetNode(19, -0.0833333, -0.0116117, 0.0934434, 0);
  loader.GetNode(20, -0.0833333, -0.1, 0.0934434, 0);
  loader.GetNode(21, -0.0934434, -0.1, 0.0833333, 0);
  loader.GetCutNode(98, -0.0883883, -0.1, 0.0883883, 0);
  loader.GetCutNode(105, -0.150888, -0.0116117, 0.0258883, 0);
  loader.GetCutNode(-4, -0.150888, -0.1, 0.0258883, 0);
  loader.GetCutNode(108, -0.0258883, -0.0116117, 0.150888, 0);
  loader.GetCutNode(-6, -0.0258883, -0.1, 0.150888, 0);
  loader.GetCutNode(-7, -0.0883883, -0.0116117, 0.0883883, 0);
  loader.GetCutNode(112, -0.0883883, 0.0767767, 0.0883883, 0);
  loader.GetCutNode(-14, -0.0754442, 0.0325825, 0.0129442, 0);
  loader.GetCutNode(113, 0, 0.0767767, 0, 0);
  loader.GetCutNode(-18, -0.0129442, 0.0325825, 0.0754442, 0);
  loader.GetNode(0, -0.125, -0.1, -0.125, 0);
  loader.GetNode(1, -0.25, -0.1, -0.0833333, 0);
  loader.GetNode(2, -0.0833333, -0.1, -0.0833333, 0);
  loader.GetNode(3, -0.25, -0.1, -0.25, 0);
  loader.GetNode(4, -0.0833333, -0.1, -0.25, 0);
  loader.GetNode(5, -0.25, 0.1, -0.0833333, 0);
  loader.GetNode(6, -0.0833333, 0.1, -0.0833333, 0);
  loader.GetNode(7, -0.25, 0.1, -0.25, 0);
  loader.GetNode(8, -0.0833333, 0.1, -0.25, 0);
  loader.GetNode(9, -0.0833333, -0.0410744, -0.0833333, 0);
  loader.GetNode(10, -0.0833333, -0.0755922, -0.132149, 0);
  loader.GetNode(11, -0.0833333, -0.1, -0.166667, 0);
  loader.GetNode(12, -0.166667, -0.1, -0.0833333, 0);
  loader.GetNode(13, -0.132149, -0.0755922, -0.0833333, 0);
  loader.GetCutNode(107, -0.0625, -0.0116117, -0.0625, 0);
  loader.GetCutNode(106, -0.125, -0.1, -0.125, 0);
  loader.GetCutNode(-13, -0.137944, -0.0558058, -0.0495558, 0);
  loader.GetCutNode(104, -0.213388, -0.1, -0.0366117, 0);
  loader.GetCutNode(-15, -0.0495558, -0.0558058, -0.137944, 0);
  loader.GetCutNode(118, -0.0366117, -0.1, -0.213388, 0);
  loader.GetNode(0, 0.125, -0.1, 0.125, 0);
  loader.GetNode(1, 0.0833333, -0.3, 0.25, 0);
  loader.GetNode(2, 0.0833333, -0.3, 0.0833333, 0);
  loader.GetNode(3, 0.0833333, -0.1, 0.0833333, 0);
  loader.GetNode(4, 0.0833333, -0.1, 0.25, 0);
  loader.GetNode(5, 0.25, -0.3, 0.25, 0);
  loader.GetNode(6, 0.25, -0.3, 0.0833333, 0);
  loader.GetNode(7, 0.25, -0.1, 0.0833333, 0);
  loader.GetNode(8, 0.25, -0.1, 0.25, 0);
  loader.GetNode(9, 0.166667, -0.1, 0.0833333, 0);
  loader.GetNode(10, 0.0833333, -0.1, 0.166667, 0);
  loader.GetNode(11, 0.0833333, -0.158926, 0.0833333, 0);
  loader.GetNode(12, 0.132149, -0.124408, 0.0833333, 0);
  loader.GetNode(13, 0.0833333, -0.124408, 0.132149, 0);
  loader.GetCutNode(103, 0.0625, -0.188388, 0.0625, 0);
  loader.GetCutNode(111, 0.125, -0.1, 0.125, 0);
  loader.GetCutNode(-8, 0.0495558, -0.144194, 0.137944, 0);
  loader.GetCutNode(109, 0.0366117, -0.1, 0.213388, 0);
  loader.GetCutNode(121, 0.213388, -0.1, 0.0366117, 0);
  loader.GetCutNode(-12, 0.137944, -0.144194, 0.0495558, 0);
  loader.GetNode(0, 0.0883883, -0.276777, -0.0883883, 0);
  loader.GetNode(1, 0.0883883, -0.1, -0.0883883, 0);
  loader.GetNode(2, 0.0883883, -0.188388, -0.0883883, 0);
  loader.GetNode(3, 0.0833333, -0.3, -0.0833333, 0);
  loader.GetNode(4, 0.0833333, -0.1, -0.0833333, 0);
  loader.GetNode(5, 0.0833333, -0.3, -0.25, 0);
  loader.GetNode(6, 0.0833333, -0.1, -0.25, 0);
  loader.GetNode(7, 0.25, -0.3, -0.0833333, 0);
  loader.GetNode(8, 0.25, -0.1, -0.0833333, 0);
  loader.GetNode(9, 0.25, -0.3, -0.25, 0);
  loader.GetNode(10, 0.25, -0.1, -0.25, 0);
  loader.GetNode(11, 0.0934434, -0.1, -0.0833333, 0);
  loader.GetNode(12, 0.0934434, -0.107149, -0.0833333, 0);
  loader.GetNode(13, 0.0833333, -0.1, -0.0934434, 0);
  loader.GetNode(14, 0.0833333, -0.107149, -0.0934434, 0);
  loader.GetNode(15, 0.0833333, -0.269628, -0.0934434, 0);
  loader.GetNode(16, 0.0833333, -0.188388, -0.0934434, 0);
  loader.GetNode(17, 0.0934434, -0.269628, -0.0833333, 0);
  loader.GetNode(18, 0.0934434, -0.188388, -0.0833333, 0);
  loader.GetNode(19, 0.0833333, -0.276777, -0.0833333, 0);
  loader.GetNode(20, 0.087521, -0.273816, -0.0833333, 0);
  loader.GetNode(21, 0.0833333, -0.273816, -0.087521, 0);
  loader.GetCutNode(115, 0.0258883, -0.188388, -0.150888, 0);
  loader.GetCutNode(114, 0.0883883, -0.276777, -0.0883883, 0);
  loader.GetCutNode(-9, 0.0129442, -0.232583, -0.0754442, 0);
  loader.GetCutNode(100, 0, -0.276777, 0, 0);
  loader.GetCutNode(-10, 0.0754442, -0.232583, -0.0129442, 0);
  loader.GetCutNode(117, 0.150888, -0.188388, -0.0258883, 0);
  loader.GetCutNode(-21, 0.0883883, -0.188388, -0.0883883, 0);
  loader.GetCutNode(116, 0.0883883, -0.1, -0.0883883, 0);
  loader.GetCutNode(-22, 0.0258883, -0.1, -0.150888, 0);
  loader.GetCutNode(-23, 0.150888, -0.1, -0.0258883, 0);
  loader.GetNode(0, 0.125, -0.1, 0.125, 0);
  loader.GetNode(1, 0.0833333, -0.1, 0.0833333, 0);
  loader.GetNode(2, 0.0833333, -0.1, 0.25, 0);
  loader.GetNode(3, 0.0833333, 0.1, 0.0833333, 0);
  loader.GetNode(4, 0.0833333, 0.1, 0.25, 0);
  loader.GetNode(5, 0.25, -0.1, 0.0833333, 0);
  loader.GetNode(6, 0.25, -0.1, 0.25, 0);
  loader.GetNode(7, 0.25, 0.1, 0.0833333, 0);
  loader.GetNode(8, 0.25, 0.1, 0.25, 0);
  loader.GetNode(9, 0.0833333, -0.0410744, 0.0833333, 0);
  loader.GetNode(10, 0.132149, -0.0755922, 0.0833333, 0);
  loader.GetNode(11, 0.166667, -0.1, 0.0833333, 0);
  loader.GetNode(12, 0.0833333, -0.1, 0.166667, 0);
  loader.GetNode(13, 0.0833333, -0.0755922, 0.132149, 0);
  loader.GetCutNode(109, 0.0366117, -0.1, 0.213388, 0);
  loader.GetCutNode(111, 0.125, -0.1, 0.125, 0);
  loader.GetCutNode(-17, 0.0495558, -0.0558058, 0.137944, 0);
  loader.GetCutNode(110, 0.0625, -0.0116117, 0.0625, 0);
  loader.GetCutNode(-19, 0.137944, -0.0558058, 0.0495558, 0);
  loader.GetCutNode(121, 0.213388, -0.1, 0.0366117, 0);
  loader.GetNode(0, 0.0883883, 0.0767767, -0.0883883, 0);
  loader.GetNode(1, 0.0883883, -0.1, -0.0883883, 0);
  loader.GetNode(2, 0.0883883, -0.0116117, -0.0883883, 0);
  loader.GetNode(3, 0.0833333, -0.1, -0.0833333, 0);
  loader.GetNode(4, 0.0833333, -0.1, -0.25, 0);
  loader.GetNode(5, 0.0833333, 0.1, -0.0833333, 0);
  loader.GetNode(6, 0.0833333, 0.1, -0.25, 0);
  loader.GetNode(7, 0.25, -0.1, -0.0833333, 0);
  loader.GetNode(8, 0.25, -0.1, -0.25, 0);
  loader.GetNode(9, 0.25, 0.1, -0.0833333, 0);
  loader.GetNode(10, 0.25, 0.1, -0.25, 0);
  loader.GetNode(11, 0.0833333, -0.0928511, -0.0934434, 0);
  loader.GetNode(12, 0.0833333, -0.0116117, -0.0934434, 0);
  loader.GetNode(13, 0.0934434, -0.0928511, -0.0833333, 0);
  loader.GetNode(14, 0.0934434, -0.0116117, -0.0833333, 0);
  loader.GetNode(15, 0.0833333, 0.0696278, -0.0934434, 0);
  loader.GetNode(16, 0.0934434, 0.0696278, -0.0833333, 0);
  loader.GetNode(17, 0.0934434, -0.1, -0.0833333, 0);
  loader.GetNode(18, 0.0833333, -0.1, -0.0934434, 0);
  loader.GetNode(19, 0.0833333, 0.0767767, -0.0833333, 0);
  loader.GetNode(20, 0.087521, 0.0738155, -0.0833333, 0);
  loader.GetNode(21, 0.0833333, 0.0738155, -0.087521, 0);
  loader.GetCutNode(113, 0, 0.0767767, 0, 0);
  loader.GetCutNode(122, 0.0883883, 0.0767767, -0.0883883, 0);
  loader.GetCutNode(-16, 0.0129442, 0.0325825, -0.0754442, 0);
  loader.GetCutNode(119, 0.0258883, -0.0116117, -0.150888, 0);
  loader.GetCutNode(120, 0.150888, -0.0116117, -0.0258883, 0);
  loader.GetCutNode(-20, 0.0754442, 0.0325825, -0.0129442, 0);
  loader.GetCutNode(116, 0.0883883, -0.1, -0.0883883, 0);
  loader.GetCutNode(-22, 0.0258883, -0.1, -0.150888, 0);
  loader.GetCutNode(-23, 0.150888, -0.1, -0.0258883, 0);
  loader.GetCutNode(-24, 0.0883883, -0.0116117, -0.0883883, 0);

  loader.CutTest_Cut(true, true);
}


void test_fluidfluid2()
{
  MeshLoader loader;

  loader.GetCutNode(32, -0.25, -0.25, 0.025, 0);
  loader.GetCutNode(35, -0.25, 0.25, 0.025, 0);
  loader.GetCutNode(34, -0.25, 0.25, -0.025, 0);
  loader.GetCutNode(33, -0.25, -0.25, -0.025, 0);
  loader.CreateSide(0, 32, 35, 34, 33, DRT::Element::quad4);
  loader.GetCutNode(-1, -0.25, 0, 0, 0);
  loader.GetCutNode(37, 0.25, -0.25, -0.025, 0);
  loader.GetCutNode(36, 0.25, -0.25, 0.025, 0);
  loader.CreateSide(1, 32, 33, 37, 36, DRT::Element::quad4);
  loader.GetCutNode(-2, 0, -0.25, 0, 0);
  loader.GetCutNode(39, 0.25, 0.25, 0.025, 0);
  loader.CreateSide(2, 32, 36, 39, 35, DRT::Element::quad4);
  loader.GetCutNode(-3, 0, 0, 0.025, 0);
  loader.GetCutNode(38, 0.25, 0.25, -0.025, 0);
  loader.CreateSide(3, 33, 34, 38, 37, DRT::Element::quad4);
  loader.GetCutNode(-4, 0, 0, -0.025, 0);
  loader.CreateSide(4, 34, 35, 39, 38, DRT::Element::quad4);
  loader.GetCutNode(-5, 0, 0.25, 0, 0);
  loader.CreateSide(5, 36, 37, 38, 39, DRT::Element::quad4);
  loader.GetCutNode(-6, 0.25, 0, 0, 0);
  loader.CreateSide(6, 32, 35, 34, 33, DRT::Element::quad4);
  loader.CreateSide(7, 32, 33, 37, 36, DRT::Element::quad4);
  loader.CreateSide(8, 32, 36, 39, 35, DRT::Element::quad4);
  loader.CreateSide(9, 33, 34, 38, 37, DRT::Element::quad4);
  loader.CreateSide(10, 34, 35, 39, 38, DRT::Element::quad4);
  loader.CreateSide(11, 36, 37, 38, 39, DRT::Element::quad4);
  loader.GetNode(0, -0.75, -0.75, 0.025, 0);
  loader.GetNode(1, -0.75, -0.75, -0.025, 0);
  loader.GetNode(2, -0.75, -0.25, -0.025, 0);
  loader.GetNode(3, -0.75, -0.25, 0.025, 0);
  loader.GetNode(4, -0.25, -0.75, 0.025, 0);
  loader.GetNode(5, -0.25, -0.75, -0.025, 0);
  loader.GetNode(6, -0.25, -0.25, -0.025, 0);
  loader.GetNode(7, -0.25, -0.25, 0.025, 0);
  loader.CreateElement(0, 0, 1, 2, 3, 4, 5, 6, 7, DRT::Element::hex8);
  loader.GetNode(8, -0.75, 0.25, -0.025, 0);
  loader.GetNode(9, -0.75, 0.25, 0.025, 0);
  loader.GetNode(10, -0.25, 0.25, -0.025, 0);
  loader.GetNode(11, -0.25, 0.25, 0.025, 0);
  loader.CreateElement(1, 3, 2, 8, 9, 7, 6, 10, 11, DRT::Element::hex8);
  loader.GetNode(12, -0.75, 0.75, -0.025, 0);
  loader.GetNode(13, -0.75, 0.75, 0.025, 0);
  loader.GetNode(14, -0.25, 0.75, -0.025, 0);
  loader.GetNode(15, -0.25, 0.75, 0.025, 0);
  loader.CreateElement(2, 9, 8, 12, 13, 11, 10, 14, 15, DRT::Element::hex8);
  loader.GetNode(16, 0.25, -0.75, 0.025, 0);
  loader.GetNode(17, 0.25, -0.75, -0.025, 0);
  loader.GetNode(18, 0.25, -0.25, -0.025, 0);
  loader.GetNode(19, 0.25, -0.25, 0.025, 0);
  loader.CreateElement(3, 4, 5, 6, 7, 16, 17, 18, 19, DRT::Element::hex8);
  loader.GetNode(20, 0.25, 0.25, -0.025, 0);
  loader.GetNode(21, 0.25, 0.25, 0.025, 0);
  loader.CreateElement(4, 7, 6, 10, 11, 19, 18, 20, 21, DRT::Element::hex8);
  loader.GetNode(22, 0.25, 0.75, -0.025, 0);
  loader.GetNode(23, 0.25, 0.75, 0.025, 0);
  loader.CreateElement(5, 11, 10, 14, 15, 21, 20, 22, 23, DRT::Element::hex8);
  loader.GetNode(24, 0.75, -0.75, 0.025, 0);
  loader.GetNode(25, 0.75, -0.75, -0.025, 0);
  loader.GetNode(26, 0.75, -0.25, -0.025, 0);
  loader.GetNode(27, 0.75, -0.25, 0.025, 0);
  loader.CreateElement(6, 16, 17, 18, 19, 24, 25, 26, 27, DRT::Element::hex8);
  loader.GetNode(28, 0.75, 0.25, -0.025, 0);
  loader.GetNode(29, 0.75, 0.25, 0.025, 0);
  loader.CreateElement(7, 19, 18, 20, 21, 27, 26, 28, 29, DRT::Element::hex8);
  loader.GetNode(30, 0.75, 0.75, -0.025, 0);
  loader.GetNode(31, 0.75, 0.75, 0.025, 0);
  loader.CreateElement(8, 21, 20, 22, 23, 29, 28, 30, 31, DRT::Element::hex8);

  loader.CutTest_Cut(true, true);
}
