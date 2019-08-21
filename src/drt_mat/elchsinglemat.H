/*----------------------------------------------------------------------*/
/*! \file
\brief abstract interface for electrode and electrolyte materials carrying concentration and
electric potential as degrees of freedom

\maintainer Christoph Schmidt

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef ELCHSINGLEMAT_H
#define ELCHSINGLEMAT_H

#include "material.H"
#include "matpar_parameter.H"

namespace MAT
{
  namespace PAR
  {
    //! parameters for abstract battery material
    class ElchSingleMat : public Parameter
    {
     public:
      //! create instance of electrode or electrolyte material
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() = 0;

      //! @name parameters for abstract battery material
      //@{
      //! curve number for diffusion coefficient
      const int diffcurve_;

      //! number of parameters for diffusion coefficient
      const int diffparanum_;

      //! parameters for diffusion coefficient
      const std::vector<double> diffpara_;

      //! curve number for conductivity
      const int condcurve_;

      //! number of parameters for conductivity
      const int condparanum_;

      //! parameters for conductivity
      const std::vector<double> condpara_;
      //@}

     protected:
      //! constructor
      ElchSingleMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~ElchSingleMat() { return; };

      //! check whether number of parameters is consistent with curve number
      void CheckProvidedParams(const int functnr, const std::vector<double>& functparams);
    };  // class MAT::PAR::ElchSingleMat
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  //! wrapper for abstract battery material
  class ElchSingleMat : public Material
  {
   public:
    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const = 0;

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const = 0;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    virtual void Unpack(const std::vector<char>& data) = 0;
    //@}

    //! compute diffusion coefficient according to curve number
    double ComputeDiffusionCoefficient(const double cint) const;

    //! compute first derivative of diffusion coefficient according to curve number
    double ComputeFirstDerivDiffCoeff(const double cint) const;

    //! compute conductivity according to curve number
    double ComputeConductivity(const double cint) const;

    //! compute first derivative of conductivity according to curve number
    double ComputeFirstDerivCond(const double cint) const;

   protected:
    //! return curve number for diffusion coefficient
    int DiffCurve() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())->diffcurve_;
    };

    //! return curve number for conductivity
    int CondCurve() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())->condcurve_;
    };

    //! return parameters for diffusion coefficient
    const std::vector<double>& DiffParams() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())->diffpara_;
    };

    //! return parameters for conductivity
    const std::vector<double>& CondParams() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())->condpara_;
    };

    //! evaluate value of predefined function
    double EvalFunctValue(
        const int functnr, const double cint, const std::vector<double>& functparams) const;

    //! evaluate first derivative of predefined function
    double EvalFirstDerivFunctValue(
        const int functnr, const double cint, const std::vector<double>& functparams) const;
  };
}  // namespace MAT
#endif
