/*----------------------------------------------------------------------*/
/*!
\file particle_sph_weightFunction.H

\brief weight functions for SPH methods

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLE_SPH_WEIGHTFUNCTION_H_
#define SRC_DRT_PARTICLE_PARTICLE_SPH_WEIGHTFUNCTION_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_particle.H"

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> FAD;

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

class WeightFunction_Base
{
public :

  // destructor
  virtual ~WeightFunction_Base() { ; };

  // get smoothing length h from radius r
  virtual double SmoothingLength(const double &radius) = 0;

  // weight function
  virtual double W(const double &disRel, const double &radius) = 0;

  // weight function (FAD)
  virtual FAD W(const FAD &disRel, const double &radius) = 0;

  // weight function 1 derivative
  virtual double DW(const double &disRel, const double &radius) = 0;

  // weight function 1 derivative (FAD)
  virtual FAD DW(const FAD &disRel, const double &radius) = 0;

  // weight function 2 derivative
  virtual double DDW(const double &disRel, const double &radius) = 0;

  // gradient, r = r_i - r_j
  virtual LINALG::Matrix<3,1> GradW(const LINALG::Matrix<3,1> &rRelVersor, const double& dw);

  // gradient, r = r_i - r_j (FAD)
  virtual LINALG::TMatrix<FAD,3,1> GradW(const LINALG::TMatrix<FAD,3,1> &rRelVersor, const FAD& dw);

  // hessian, r = r_i - r_j
  virtual LINALG::Matrix<3,3> HessW(const LINALG::Matrix<3,1> &rRelVersor, const double& dw, const double& rNorm2, const double& ddw);

  // weight function in 0 (auto-interactions)
  virtual double W0(const double& radius) = 0;

  // hessian of the weight function in 0 (auto-interactions)
  virtual double DDW0(const double& radius) = 0;

  // weight function type
  virtual INPAR::PARTICLE::WeightFunction Name() = 0;

  // check if the gradient is correct
  virtual void DBG_GradW();

  // check if the gradient is correct
  virtual void DBG_HessW();

private :

  // resizer due to the dimension integral
  virtual double RszDim(const double& radius) = 0;

};

class WeightFunction_CubicBspline : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

   //! constructor
   WeightFunction_CubicBspline(INPAR::PARTICLE::WeightFunctionDim dim):
    WF_DIM_(dim){}

   // get smoothing length h from radius r
   virtual inline double SmoothingLength(const double &radius){return radius/2.0; }

   //! compute the cubicBspline weight function
   virtual double W(const double &disRel, const double &radius);

   //! compute the cubicBspline weight function (FAD)
   virtual FAD W(const FAD &disRel, const double &radius);

   //! compute the weight function 1 derivative of the cubicBspline weight function
   virtual double DW(const double &disRel, const double &radius);

   //! compute the weight function 1 derivative of the cubicBspline weight function (FAD)
   virtual FAD DW(const FAD &disRel, const double &radius);

   //! compute the weight function 2 derivative of the cubicBspline weight function
   virtual double DDW(const double &disRel, const double &radius);

   // weight function in 0 (auto-interactions)
   virtual inline double W0(const double& radius) { return 2.0/3.0 * RszDim(radius); }

   // hessian of the weight function in 0 (auto-interactions)
   virtual inline double DDW0(const double& radius) { return -8 * RszDim(radius) / std::pow(radius,2); }

   // weight function type
   virtual inline INPAR::PARTICLE::WeightFunction Name() { return INPAR::PARTICLE::CubicBspline; }

   //@}

private :

   virtual double RszDisRel(const double& disRel, const double& radius) {return 2 * disRel / radius; }

   virtual FAD RszDisRel(const FAD& disRel, const double& radius) {return 2 * disRel / radius; }

   virtual double RszDim(const double& radius);

protected:

   // store problem dimension required by weight functions
   INPAR::PARTICLE::WeightFunctionDim WF_DIM_;

};

class WeightFunction_QuinticBspline : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

   //! constructor
   WeightFunction_QuinticBspline(INPAR::PARTICLE::WeightFunctionDim dim):
    WF_DIM_(dim){}

   // get smoothing length h from radius r
   virtual inline double SmoothingLength(const double &radius){return radius/3.0; }

   //! compute the QuinticBspline weight function
   virtual double W(const double &disRel, const double &radius);

   //! compute the QuinticBspline weight function (FAD)
   virtual FAD W(const FAD &disRel, const double &radius);

   //! compute the weight function 1 derivative of the QuinticBspline weight function
   virtual double DW(const double &disRel, const double &radius);

   //! compute the weight function 1 derivative of the QuinticBspline weight function (FAD)
   virtual FAD DW(const FAD &disRel, const double &radius);

   //! compute the weight function 2 derivative of the QuinticBspline weight function
   virtual double DDW(const double &disRel, const double &radius);

   // weight function in 0 (auto-interactions)
   virtual inline double W0(const double& radius)
   {
     double h=radius/3.0;
     return 66*RszDim(h);
   }

   // hessian of the weight function in 0 (auto-interactions)
   virtual inline double DDW0(const double& radius)
   {
     double h=radius/3.0;
     return -120*RszDim(h)/(h*h);
   }

   // weight function type
   virtual inline INPAR::PARTICLE::WeightFunction Name() { return INPAR::PARTICLE::QuinticBspline; }

   //@}

private :

   virtual double RszDim(const double& h);

protected:

   // store problem dimension required by weight functions
   INPAR::PARTICLE::WeightFunctionDim WF_DIM_;

};

class WeightFunction_SqrtHyperbola : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

   //! constructor
   WeightFunction_SqrtHyperbola(INPAR::PARTICLE::WeightFunctionDim dim):
    WF_DIM_(dim){}

   // get smoothing length h from radius r
   virtual inline double SmoothingLength(const double &radius){return radius; }

   //! compute the SqrtHyperbola weight function
   virtual double W(const double &disRel, const double &radius);

   //! compute the SqrtHyperbola weight function
   virtual FAD W(const FAD &disRel, const double &radius);

   //! compute the weight function 1 derivative of the SqrtHyperbola weight function
   virtual double DW(const double &disRel, const double &radius);

   //! compute the weight function 1 derivative of the SqrtHyperbola weight function (FAD)
   virtual FAD DW(const FAD &disRel, const double &radius);

   //! compute the weight function 2 derivative of the SqrtHyperbola weight function
   virtual double DDW(const double &disRel, const double &radius);

   // weight function in 0 (auto-interactions) - if it is not allowed, it is 0
   virtual inline double W0(const double& radius) { return RszDim(radius); }

   // hessian of the weight function in 0 (auto-interactions)
   virtual inline double DDW0(const double& radius) { return RszDim(radius); }

   // weight function type
   virtual inline INPAR::PARTICLE::WeightFunction Name() { return INPAR::PARTICLE::SqrtHyperbola; }

   //@}

private :

   virtual double RszDim(const double& radius);

protected:

   // store problem dimension required by weight functions
   INPAR::PARTICLE::WeightFunctionDim WF_DIM_;

};


class WeightFunction_HyperbolaNoRsz : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

  // get smoothing length h from radius r
  virtual inline double SmoothingLength(const double &radius){return radius; }

   //! compute the HyperbolaNoRsz weight function
   virtual double W(const double &disRel, const double &radius);

   //! compute the HyperbolaNoRsz weight function
   virtual FAD W(const FAD &disRel, const double &radius);

   //! compute the weight function 1 derivative of the HyperbolaNoRsz weight function
   virtual double DW(const double &disRel, const double &radius);

   //! compute the weight function 1 derivative of the HyperbolaNoRsz weight function (FAD)
   virtual FAD DW(const FAD &disRel, const double &radius);

   //! compute the weight function 1 derivative of the HyperbolaNoRsz weight function
   virtual double DDW(const double &disRel, const double &radius);

   // weight function in 0 (auto-interactions) - if it is not allowed, it is 0
   virtual inline double W0(const double& radius) { return RszDim(radius); }

   // hessian of the weight function in 0 (auto-interactions)
   virtual inline double DDW0(const double& radius) { return RszDim(radius); }

   // weight function type
   virtual inline INPAR::PARTICLE::WeightFunction Name() { return INPAR::PARTICLE::HyperbolaNoRsz; }

   //@}

private :

   virtual double RszDim(const double& radius) {return 1; };

};


}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLE_SPH_WEIGHTFUNCTION_H_ */
