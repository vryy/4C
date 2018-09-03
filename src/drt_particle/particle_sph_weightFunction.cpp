/*----------------------------------------------------------------------*/
/*!
\file particle_sph_weightFunction.cpp

\brief weight functions for SPH methods

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
/* headers */
#include "particle_sph_weightFunction.H"

/*----------------------------------------------------------------------*
 | gradient, r_ij = r_i - r_j                         cattabiani 08/16  |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 1> PARTICLE::WeightFunction_Base::GradW(
    const LINALG::Matrix<3, 1> &rRelVersor, const double &dw)
{
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Attention: If we want to split
  // the gradient according to grad(W) = dW/dr_{ij}=:F_{ij}*r_{ij} (see e.g. Monaghan2005), we get
  // the identity:F_{ij}=DW/||r_{ij}||!!! Attention: In Espanol2003, F_{ij} ist defined with the
  // opposite sign convention as applied here or in Monaghan2005!!!
  LINALG::Matrix<3, 1> gradW(rRelVersor);
  gradW.Scale(dw);

  return gradW;
}

/*----------------------------------------------------------------------*
 | gradient, r_ij = r_i - r_j (FAD version)                meier 03/17  |
 *----------------------------------------------------------------------*/
LINALG::TMatrix<FAD, 3, 1> PARTICLE::WeightFunction_Base::GradW(
    const LINALG::TMatrix<FAD, 3, 1> &rRelVersor, const FAD &dw)
{
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Attention: If we want to split
  // the gradient according to grad(W) = dW/dr_{ij}=:F_{ij}*r_{ij} (see e.g. Monaghan2005), we get
  // the identity:F_{ij}=DW/||r_{ij}||!!! Attention: In Espanol2003, F_{ij} ist defined with the
  // opposite sign convention as applied here or in Monaghan2005!!!
  LINALG::TMatrix<FAD, 3, 1> gradW(rRelVersor);
  gradW.Scale(dw);

  return gradW;
}


/*----------------------------------------------------------------------*
 | hessian, r_ij = r_i - r_j                          cattabiani 08/16  |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 3> PARTICLE::WeightFunction_Base::HessW(const LINALG::Matrix<3, 1> &rRelVersor,
    const double &dw, const double &rNorm2, const double &ddw)
{
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Following the product rule,
  // the second derivative reads (x represents a dyadic product and I_3 a proper identify matrix):
  // d^2W/dr_{ij}^2 = (d^2W/dq^2* 1/h * e_{ij}) * 1/h * e_{ij} + dW/dq * 1/h * de_{ij}/ddr_{ij}
  //                =1/h^2*d^2W/dq^2*(e_{ij} x e_{ij}^T) + 1/h*dW/dq*(I_3-e_{ij} x e_{ij}^T)
  //                =DDW*(e_{ij} x e_{ij}^T) + DW*(I_3-e_{ij} x e_{ij}^T)
  // Here, the definition DDW=1/h^2*d^2W/dq^2 as applied in the funciton DDW() has been introduced.

  // compute the resized first derivative
  const double dw_r = dw / rNorm2;
  // diadic product
  LINALG::Matrix<3, 3> rRelVersor2;
  rRelVersor2.MultiplyNT(rRelVersor, rRelVersor);
  // create the hessian and add the various terms
  LINALG::Matrix<3, 3> hessW;
  for (int ii = 0; ii < 3; ++ii)
  {
    hessW(ii, ii) = dw_r;
  }
  hessW.Update(ddw - dw_r, rRelVersor2, 1.0);

  return hessW;
}

/*----------------------------------------------------------------------*
 | gradient check                                     cattabiani 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::WeightFunction_Base::DBG_GradW()
{
  LINALG::Matrix<3, 1> r;
  const double radius = 1;
  const double sampl = 100.0;  // double because of the divisions
  const double incr = 1e-8;
  const double tol = 1e-6;

  for (int ri = 1; ri < sampl; ++ri)
  {
    r(0) = 2 * radius * ri / sampl;
    for (int rj = 1; rj < sampl; ++rj)
    {
      r(1) = 2 * radius * rj / sampl;
      for (int rk = 1; rk < sampl; ++rk)
      {
        r(2) = 2 * radius * rk / sampl;

        LINALG::Matrix<3, 1> rVersor(r);
        const double rNorm2 = rVersor.Norm2();
        rVersor.Scale(1 / rNorm2);

        // compute the gradient in the standard method
        LINALG::Matrix<3, 1> gradW1;
        gradW1 = GradW(rVersor, DW(rNorm2, radius));

        // compute the gradient with the finite differences
        LINALG::Matrix<3, 1> gradW2;
        double w0 = W(rNorm2, radius);
        for (int dim = 0; dim < 3; ++dim)
        {
          LINALG::Matrix<3, 1> rfd(r);

          rfd(dim) += incr;
          double w1 = W(rfd.Norm2(), radius);
          gradW2(dim) = (w1 - w0) / incr;
        }

        LINALG::Matrix<3, 1> gradWdiff(gradW1);
        gradWdiff.Update(1.0, gradW2, -1.0);

        if (gradWdiff.Norm2() / 3 > tol)
        {
          std::cout << "The check of the gradient failed\n";
          std::cout << "Standard:\n";
          std::cout << gradW1 << std::endl;
          std::cout << "Finite differences:\n";
          std::cout << gradW2 << std::endl;
          std::cout << "Difference:\n";
          std::cout << gradWdiff << std::endl;
          std::cout << "Difference norm:\n";
          std::cout << gradWdiff.Norm2() << std::endl;
          std::cin.get();
        }
      }
    }
  }
  std::cout << "Gradient check passed!\n";
  std::cin.get();
}

/*----------------------------------------------------------------------*
 | gradient check                                     cattabiani 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::WeightFunction_Base::DBG_HessW()
{
  LINALG::Matrix<3, 1> r;
  const double radius = 1;
  const double sampl = 100.0;  // double because of the divisions
  const double incr = 1e-8;
  const double tol = 1e-6;

  for (int ri = 1; ri < sampl; ++ri)
  {
    r(0) = radius * ri / sampl;
    for (int rj = 1; rj < sampl; ++rj)
    {
      r(1) = radius * rj / sampl;
      for (int rk = 1; rk < sampl; ++rk)
      {
        r(2) = radius * rk / sampl;

        LINALG::Matrix<3, 1> rVersor(r);
        const double rNorm2 = rVersor.Norm2();
        rVersor.Scale(1 / rNorm2);

        // compute the gradient in the standard method
        LINALG::Matrix<3, 3> hessW1;
        hessW1 = HessW(rVersor, DW(rNorm2, radius), rNorm2, DDW(rNorm2, radius));

        // compute the gradient with the finite differences
        LINALG::Matrix<3, 3> hessW2;
        LINALG::Matrix<3, 1> gradW0 = GradW(rVersor, DW(rNorm2, radius));

        for (int dimj = 0; dimj < 3; ++dimj)
        {
          LINALG::Matrix<3, 1> rfd(r);
          rfd(dimj) += incr;
          LINALG::Matrix<3, 1> rfdVersor(rfd);
          const double rfdNorm2 = rfd.Norm2();
          rfdVersor.Scale(1 / rfdNorm2);

          for (int dimi = 0; dimi < 3; ++dimi)
          {
            LINALG::Matrix<3, 1> gradW1 = GradW(rfdVersor, DW(rfdNorm2, radius));
            hessW2(dimi, dimj) = (gradW1(dimi) - gradW0(dimi)) / incr;
          }
        }


        LINALG::Matrix<3, 3> hessWdiff(hessW1);
        hessWdiff.Update(1.0, hessW2, -1.0);

        if (hessWdiff.Norm2() / 9 > tol)
        {
          std::cout << "The check of the hessian failed\n";
          std::cout << "Standard:\n";
          std::cout << hessW1 << std::endl;
          std::cout << "Finite differences:\n";
          std::cout << hessW2 << std::endl;
          std::cout << "Difference:\n";
          std::cout << hessWdiff << std::endl;
          std::cout << "Difference norm:\n";
          std::cout << hessWdiff.Norm2() << std::endl;
          std::cin.get();
        }
      }
    }
  }
  std::cout << "Hessian check passed!\n";
  std::cin.get();
}


/*----------------------------------------------------------------------*
 | compute the cubicBspline w function                cattabiani 08/16  |
 *----------------------------------------------------------------------*/

// The 3D variant can be found in Monaghan2005, Eq. (2.6)
// [see also Eq. (9) in Antoci2007]

double PARTICLE::WeightFunction_CubicBspline::W(const double &disRel, const double &radius)
{
  // Attention: The support of our kernel functions is (in 3D) defined by a sphere with radius
  // measured by our variable "radius". In the SLM literature, typically the smoothing length
  // h:=radius/2 is used for defining the kernel functions!!!
  // rszDisRel=q:=||r_{ij}||/h=2*||r_{ij}||/radius
  const double rszDisRel = RszDisRel(disRel, radius);

  double w = 0;
  if (rszDisRel < 1)
  {
    w = (2.0 / 3.0) - std::pow(rszDisRel, 2) + 0.5 * std::pow(rszDisRel, 3);
  }
  else if (rszDisRel < 2)
  {
    w = std::pow(2 - rszDisRel, 3) / 6;
  }

  // resizing to have an integral = 1
  w *= RszDim(radius);

  return w;
}

/*----------------------------------------------------------------------*
 | compute the cubicBspline w function (FAD version)       meier 03/17  |
 *----------------------------------------------------------------------*/

// The 3D variant can be found in Monaghan2005, Eq. (2.6)
// [see also Eq. (9) in Antoci2007]

FAD PARTICLE::WeightFunction_CubicBspline::W(const FAD &disRel, const double &radius)
{
  // Attention: The support of our kernel functions is (in 3D) defined by a sphere with radius
  // measured by our variable "radius". In the SLM literature, typically the smoothing length
  // h:=radius/2 is used for defining the kernel functions!!!
  // rszDisRel=q:=||r_{ij}||/h=2*||r_{ij}||/radius
  const FAD rszDisRel = RszDisRel(disRel, radius);

  FAD w = 0;
  if (rszDisRel.val() < 1)
  {
    w = (2.0 / 3.0) - std::pow(rszDisRel, 2) + 0.5 * std::pow(rszDisRel, 3);
  }
  else if (rszDisRel.val() < 2)
  {
    w = std::pow(2 - rszDisRel, 3) / 6;
  }

  // resizing to have an integral = 1
  w *= RszDim(radius);

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the cubicBspline w function derivative       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D
// See again Monaghan2005, Eq. (2.6)

double PARTICLE::WeightFunction_CubicBspline::DW(const double &disRel, const double &radius)
{
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Attention: If we want to split
  // the gradient according to grad(W) = dW/dr_{ij}=:F_{ij}*r_{ij} (see e.g. Monaghan2005), we get
  // the identity:F_{ij}=DW/||r_{ij}||!!! Attention: In Espanol2003, F_{ij} ist defined with the
  // opposite sign convention as applied here or in Monaghan2005!!!
  const double rszDisRel = RszDisRel(disRel, radius);

  double dw = 0;
  if (rszDisRel < 1)
  {
    dw = (-4 * rszDisRel + 3 * std::pow(rszDisRel, 2)) / radius;
  }
  else if (rszDisRel < 2)
  {
    dw = -std::pow(2 - rszDisRel, 2) / radius;
  }

  // resizing to have an integral = 1
  dw *= RszDim(radius);

  return dw;
}

/*--------------------------------------------------------------------------------*
 | compute the cubicBspline w function derivative (FAD version)      meier 03/17  |
 *-------------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D
// See again Monaghan2005, Eq. (2.6)

FAD PARTICLE::WeightFunction_CubicBspline::DW(const FAD &disRel, const double &radius)
{
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Attention: If we want to split
  // the gradient according to grad(W) = dW/dr_{ij}=:F_{ij}*r_{ij} (see e.g. Monaghan2005), we get
  // the identity:F_{ij}=DW/||r_{ij}||!!! Attention: In Espanol2003, F_{ij} ist defined with the
  // opposite sign convention as applied here or in Monaghan2005!!!
  const FAD rszDisRel = RszDisRel(disRel, radius);

  FAD dw = 0;
  if (rszDisRel.val() < 1)
  {
    dw = (-4 * rszDisRel + 3 * std::pow(rszDisRel, 2)) / radius;
  }
  else if (rszDisRel.val() < 2)
  {
    dw = -std::pow(2 - rszDisRel, 2) / radius;
  }

  // resizing to have an integral = 1
  dw *= RszDim(radius);

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the cubicBspline w function derivative       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D
// See again Monaghan2005, Eq. (2.6)

double PARTICLE::WeightFunction_CubicBspline::DDW(const double &disRel, const double &radius)
{
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Following the product rule,
  // the second derivative reads (x represents a dyadic product and I_3 a proper identify matrix):
  // d^2W/dr_{ij}^2 = (d^2W/dq^2* 1/h * e_{ij}) * 1/h * e_{ij} + dW/dq * 1/h * de_{ij}/ddr_{ij}
  //                =1/h^2*d^2W/dq^2*(e_{ij} x e_{ij}^T) + 1/h*dW/dq*(I_3-e_{ij} x e_{ij}^T)
  //                =DDW*(e_{ij} x e_{ij}^T) + DW*(I_3-e_{ij} x e_{ij}^T)
  // Here, the definition DDW=1/h^2*d^2W/dq^2 as applied in the funciton DDW() has been introduced.

  const double rszDisRel = RszDisRel(disRel, radius);

  double ddw = 0;
  if (rszDisRel < 1)
  {
    ddw = (-8 + 12 * rszDisRel) / std::pow(radius, 2);
  }
  else if (rszDisRel < 2)
  {
    ddw = 4 * (2 - rszDisRel) / std::pow(radius, 2);
  }

  // resizing to have an integral = 1
  ddw *= RszDim(radius);

  return ddw;
}


/*-----------------------------------------------------------------------------*
 | rsz in case of different dimensions                       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_CubicBspline::RszDim(const double &radius)
{
  double rszdim = 0.0;
  switch (WF_DIM_)
  {
    case INPAR::PARTICLE::WF_3D:
      rszdim = 12.0 * M_1_PI / std::pow(radius, 3);
      break;
    case INPAR::PARTICLE::WF_2D:
      rszdim = 60.0 * M_1_PI / (7.0 * std::pow(radius, 2));
      break;
    case INPAR::PARTICLE::WF_1D:
      rszdim = 2.0 / radius;
      break;
    default:
      dserror("Only the problem dimensions 1, 2 and 3 are possible!");
  }

  return rszdim;
}



//*************************************************************************************************************************
//*************************************************************************************************************************

/*----------------------------------------------------------------------*
 | compute the QuinticBspline w function                   meier 04/17  |
 *----------------------------------------------------------------------*/

double PARTICLE::WeightFunction_QuinticBspline::W(const double &disRel, const double &radius)
{
  // See Liu2010, Eq. (38) and Morris1997, Eq. (24) for the definition of the quintic spline in
  // different space dimensions! Attention: The support of our kernel functions is (in 3D) defined by
  // a sphere with radius measured by our variable "radius". In the SLM literature, typically the
  // smoothing length !!!h:=radius/3!!! is used for defining the quintic kernel functions!!!
  // q:=||r_{ij}||/h=3*||r_{ij}||/radius
  const double h = radius / 3.0;
  const double q = disRel / h;

  double w = 0;
  if (q < 1)
  {
    w = std::pow(3 - q, 5) - 6 * std::pow(2 - q, 5) + 15 * std::pow(1 - q, 5);
    // w(q=0)=RszDim(h)*(3^5-6*2^5+15*1^5)=66*RszDim(h) is evaluated in the method W0()
  }
  else if (q < 2)
  {
    w = std::pow(3 - q, 5) - 6 * std::pow(2 - q, 5);
  }
  else if (q < 3)
  {
    w = std::pow(3 - q, 5);
  }

  // resizing to have an integral = 1
  w *= RszDim(h);

  return w;
}

/*----------------------------------------------------------------------*
 | compute the QuinticBspline w function (FAD version)     meier 04/17  |
 *----------------------------------------------------------------------*/

FAD PARTICLE::WeightFunction_QuinticBspline::W(const FAD &disRel, const double &radius)
{
  // See Liu2010, Eq. (38) and Morris1997, Eq. (24) for the definition of the quintic spline in
  // different space dimensions! Attention: The support of our kernel functions is (in 3D) defined by
  // a sphere with radius measured by our variable "radius". In the SLM literature, typically the
  // smoothing length !!!h:=radius/3!!! is used for defining the quintic kernel functions!!!
  // q:=||r_{ij}||/h=3*||r_{ij}||/radius
  const double h = radius / 3.0;
  const FAD q = disRel / h;

  FAD w = 0;
  if (q < 1)
  {
    w = std::pow(3 - q, 5) - 6 * std::pow(2 - q, 5) + 15 * std::pow(1 - q, 5);
    // w(q=0)=RszDim(h)*(3^5-6*2^5+15*1^5)=66*RszDim(h) is evaluated in the method W0()
  }
  else if (q < 2)
  {
    w = std::pow(3 - q, 5) - 6 * std::pow(2 - q, 5);
  }
  else if (q < 3)
  {
    w = std::pow(3 - q, 5);
  }

  // resizing to have an integral = 1
  w *= RszDim(h);

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the QuinticBspline w function derivative               meier 04/17  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_QuinticBspline::DW(const double &disRel, const double &radius)
{
  // See Liu2010, Eq. (38) and Morris1997, Eq. (24) for the definition of the quintic spline in
  // different space dimensions!
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Attention: If we want to split
  // the gradient according to grad(W) = dW/dr_{ij}=:F_{ij}*r_{ij} (see e.g. Monaghan2005), we get
  // the identity:F_{ij}=DW/||r_{ij}||!!! Attention: In Espanol2003, F_{ij} ist defined with the
  // opposite sign convention as applied here or in Monaghan2005!!!
  const double h = radius / 3.0;
  const double q = disRel / h;

  double dw = 0;
  if (q < 1)
  {
    dw = -5.0 * std::pow(3 - q, 4) + 30 * std::pow(2 - q, 4) - 75 * std::pow(1 - q, 4);
    // dw(q=0)=0 for symmetry reasons
  }
  else if (q < 2)
  {
    dw = -5.0 * std::pow(3 - q, 4) + 30.0 * std::pow(2 - q, 4);
  }
  else if (q < 3)
  {
    dw = -5.0 * std::pow(3 - q, 4);
  }

  // resizing to have an integral = 1 and consideration of factor 1/h in DW = dW/dq * 1/h (see
  // comment above)
  dw *= RszDim(h) / h;

  return dw;
}

/*--------------------------------------------------------------------------------*
 | compute the QuinticBspline w function derivative (FAD version)    meier 04/17  |
 *-------------------------------------------------------------------------------*/

FAD PARTICLE::WeightFunction_QuinticBspline::DW(const FAD &disRel, const double &radius)
{
  // See Liu2010, Eq. (38) and Morris1997, Eq. (24) for the definition of the quintic spline in
  // different space dimensions!
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Attention: If we want to split
  // the gradient according to grad(W) = dW/dr_{ij}=:F_{ij}*r_{ij} (see e.g. Monaghan2005), we get
  // the identity:F_{ij}=DW/||r_{ij}||!!! Attention: In Espanol2003, F_{ij} ist defined with the
  // opposite sign convention as applied here or in Monaghan2005!!!
  const double h = radius / 3.0;
  const FAD q = disRel / h;

  FAD dw = 0;
  if (q < 1)
  {
    dw = -5.0 * std::pow(3 - q, 4) + 30 * std::pow(2 - q, 4) - 75 * std::pow(1 - q, 4);
    // dw(q=0)=0 for symmetry reasons
  }
  else if (q < 2)
  {
    dw = -5.0 * std::pow(3 - q, 4) + 30.0 * std::pow(2 - q, 4);
  }
  else if (q < 3)
  {
    dw = -5.0 * std::pow(3 - q, 4);
  }

  // resizing to have an integral = 1 and consideration of factor 1/h in DW = dW/dq * 1/h (see
  // comment above)
  dw *= RszDim(h) / h;

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the QuinticBspline w function derivative               meier 04/17  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_QuinticBspline::DDW(const double &disRel, const double &radius)
{
  // See Liu2010, Eq. (38) and Morris1997, Eq. (24) for the definition of the quintic spline in
  // different space dimensions!
  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} =
  // dW/dq * 1/h * e_{ij} The part DW = dW/dq * 1/h is determined in the function DW() while the
  // part e_{ij}:=r_{ij}/||r_{ij}|| has to be multiplied in addition! Following the product rule,
  // the second derivative reads (x represents a dyadic product and I_3 a proper identify matrix):
  // d^2W/dr_{ij}^2 = (d^2W/dq^2* 1/h * e_{ij}) * 1/h * e_{ij} + dW/dq * 1/h * de_{ij}/ddr_{ij}
  //                =1/h^2*d^2W/dq^2*(e_{ij} x e_{ij}^T) + 1/h*dW/dq*(I_3-e_{ij} x e_{ij}^T)
  //                =DDW*(e_{ij} x e_{ij}^T) + DW*(I_3-e_{ij} x e_{ij}^T)
  // Here, the definition DDW=1/h^2*d^2W/dq^2 as applied in the funciton DDW() has been introduced.

  const double h = radius / 3.0;
  const double q = disRel / h;

  double ddw = 0;
  if (q < 1)
  {
    ddw = 20.0 * std::pow(3 - q, 3) - 120 * std::pow(2 - q, 3) + 300 * std::pow(1 - q, 3);
    // ddw(q=0)=RszDim(h)/(h*h)*(20*3^3-120*2^3+300*1^3)=-120*RszDim(h)/(h*h) is evaluated in the
    // method DDW0()
  }
  else if (q < 2)
  {
    ddw = 20.0 * std::pow(3 - q, 3) - 120.0 * std::pow(2 - q, 3);
  }
  else if (q < 3)
  {
    ddw = 20.0 * std::pow(3 - q, 3);
  }

  // resizing to have an integral = 1 and consideration of factor 1/h in DW = dW/dq * 1/h (see
  // comment above)
  ddw *= RszDim(h) / (h * h);

  return ddw;
}


/*-----------------------------------------------------------------------------*
 | rsz in case of different dimensions                            meier 04/17  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_QuinticBspline::RszDim(const double &h)
{
  // See Liu2010, Eq. (38) and Morris1997, Eq. (24) for the definition of the quintic spline in
  // different space dimensions!
  double rszdim = 0.0;
  switch (WF_DIM_)
  {
    case INPAR::PARTICLE::WF_3D:
      rszdim = 3.0 * M_1_PI / (359 * std::pow(h, 3));
      break;
    case INPAR::PARTICLE::WF_2D:
      rszdim = 7.0 * M_1_PI / (478 * std::pow(h, 2));
      break;
    case INPAR::PARTICLE::WF_1D:
      rszdim = 1.0 / (120 * h);
      break;
    default:
      dserror("Only the problem dimensions 1, 2 and 3 are possible!");
  }

  return rszdim;
}

//*************************************************************************************************************************
//*************************************************************************************************************************



/*----------------------------------------------------------------------*
 | compute the SqrtHyperbola w function          cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::W(const double &disRel, const double &radius)
{
  double w = 0;
  if (disRel < radius)
  {
    w = (std::pow(radius / disRel, 0.5) - 1) * RszDim(radius);
  }

  return w;
}

/*----------------------------------------------------------------------*
 | compute the SqrtHyperbola w function (FAD version)      meier 03/17  |
 *----------------------------------------------------------------------*/
FAD PARTICLE::WeightFunction_SqrtHyperbola::W(const FAD &disRel, const double &radius)
{
  FAD w = 0;
  if (disRel.val() < radius)
  {
    w = (std::pow(radius / disRel, 0.5) - 1) * RszDim(radius);
  }

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola 1 function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::DW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel < radius)
  {
    dw = (-radius / (2 * std::pow(radius / disRel, 0.5) * disRel)) * RszDim(radius);
  }

  return dw;
}

/*------------------------------------------------------------------------------*
 | compute the SqrtHyperbola 1 function derivative (FAD version)   meier 03/17  |
 *------------------------------------------------------------------------------*/
FAD PARTICLE::WeightFunction_SqrtHyperbola::DW(const FAD &disRel, const double &radius)
{
  FAD dw = 0;
  if (disRel.val() < radius)
  {
    dw = (-radius / (2 * std::pow(radius / disRel, 0.5) * disRel)) * RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola 2 function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::DDW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel < radius)
  {
    dw = (-std::pow(radius, 2) / (4 * std::pow(radius / disRel, 3.0 / 2.0) * std::pow(disRel, 4)) +
             radius / (std::pow(radius / disRel, 0.5) * std::pow(disRel, 3))) *
         RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | rsz in case of different dimensions                       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_SqrtHyperbola::RszDim(const double &radius)
{
  double rszDim = 0.0;
  switch (WF_DIM_)
  {
    case INPAR::PARTICLE::WF_3D:
    {
      rszDim = 15.0 * M_1_PI / (4.0 * std::pow(radius, 3));
      break;
    }
    case INPAR::PARTICLE::WF_2D:
    {
      rszDim = 3.0 * M_1_PI / (std::pow(radius, 2));
      break;
    }
    case INPAR::PARTICLE::WF_1D:
    {
      rszDim = 0.5 / radius;
      break;
    }
    default:
    {
      dserror("Only the problem dimensions 1, 2 and 3 are possible!");
      break;
    }
  }

  return rszDim;
}



/*----------------------------------------------------------------------*
 | compute the HyperbolaNoRsz w function         cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::W(const double &disRel, const double &radius)
{
  double w = 0;
  if (disRel < radius)
  {
    w = (radius / disRel - 1) * RszDim(radius);
  }

  return w;
}

/*----------------------------------------------------------------------*
 | compute the HyperbolaNoRsz w function (FAD version)     meier 03/17  |
 *----------------------------------------------------------------------*/
FAD PARTICLE::WeightFunction_HyperbolaNoRsz::W(const FAD &disRel, const double &radius)
{
  FAD w = 0;
  if (disRel.val() < radius)
  {
    w = (radius / disRel - 1) * RszDim(radius);
  }

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the HyperbolaNoRsz w function derivative     cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::DW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel < radius)
  {
    dw = -radius / (disRel * disRel) * RszDim(radius);
  }

  return dw;
}

/*---------------------------------------------------------------------------------*
 | compute the HyperbolaNoRsz w function derivative (FAD version)     meier 03/17  |
 *---------------------------------------------------------------------------------*/
FAD PARTICLE::WeightFunction_HyperbolaNoRsz::DW(const FAD &disRel, const double &radius)
{
  FAD dw = 0;
  if (disRel.val() < radius)
  {
    dw = -radius / (disRel * disRel) * RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola 2 function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::DDW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel < radius)
  {
    dw = 2 * radius / std::pow(disRel, 3) * RszDim(radius);
  }

  return dw;
}
