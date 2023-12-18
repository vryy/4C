/*----------------------------------------------------------------------*/
/*! \file

 \brief main routine of the Ensight filter

 \level 1

 */


#ifndef BACI_POST_COMMON_SINGLE_FIELD_WRITERS_H
#define BACI_POST_COMMON_SINGLE_FIELD_WRITERS_H

#include "baci_post_common_filter_base.H"
#include "baci_post_common_writer_base.H"

#include <Epetra_MultiVector.h>

BACI_NAMESPACE_OPEN

class PostWriterBase;

/*!
 \brief Writer for structural problems
 */
class StructureFilter : public PostFilterBase
{
 public:
  StructureFilter(PostField* field, std::string name, std::string stresstype = "none",
      std::string straintype = "none", std::string optquantitytype = "none")
      : PostFilterBase(field, name),
        stresstype_(stresstype),
        straintype_(straintype),
        optquantitytype_(optquantitytype)
  {
  }

 protected:
  void WriteAllResults(PostField* field) override;

  void WriteAllResultsOneTimeStep(PostResult& result, bool firststep, bool laststep) override;

  /*!
  \brief postprocess gauss point stresses and write results
  \author lw
  \date 02/08
  */
  void PostStress(const std::string groupname, const std::string stresstype);
  void WriteStress(const std::string groupname, PostResult& result, const ResultType stresskind);
  void WriteEigenStress(
      const std::string groupname, PostResult& result, const ResultType stresskind);

  std::string stresstype_;
  std::string straintype_;
  std::string optquantitytype_;
};

/*!
 \brief Writer for mortar interface problems

 Each mortar interface is written as its own discretization. The MortarFilter will process only one
 of these interfaces, i.e. there will be as many MortarFilers as there are mortar interfaces.

 \author mayr.mt \date 08/2018
 */
class MortarFilter : public PostFilterBase
{
 public:
  /*!
  \brief Constructor

  @param field Field to be processed
  @param name ??
  */
  MortarFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  /*!
  \brief Write all results of a field

  So far, we don't need special control over quantities to be filtered. We just filter every result
  field.

  @param[in] field Field to be processed

  \sa WriteDofResults(), WriteNodeResults(), WriteElementResults()
  */
  void WriteAllResults(PostField* field) final;
};

/*!
 \brief Writer for fluid problems
 */
class FluidFilter : public PostFilterBase
{
 public:
  FluidFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};

/*!
 \brief Writer for xfluid problems
 */
class XFluidFilter : public PostFilterBase
{
 public:
  XFluidFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};

/*!
 \brief Writer for ale problems
 */
class AleFilter : public PostFilterBase
{
 public:
  AleFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};


/*!
 \brief Writer for interface fields in XFEM
 */
class InterfaceFilter : public PostFilterBase
{
 public:
  InterfaceFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};


/*!
 \brief Writer for lubrication problems

 \author wirtz
 \date 11/15
*/
class LubricationFilter : public PostFilterBase
{
 public:
  /// constructor
  LubricationFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};

/*!
 \brief Writer for multiphase porous fluid  problems

 \author vuong
 \date 08/16
*/
class PoroFluidMultiPhaseFilter : public PostFilterBase
{
 public:
  /// constructor
  PoroFluidMultiPhaseFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};

/*!
 \brief Writer for scalar transport problems

 \author gjb
 \date 12/07
*/
class ScaTraFilter : public PostFilterBase
{
 public:
  /// constructor
  ScaTraFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};


/*!
 \brief Writer for electrochemistry problems

 \author gjb
 \date 09/08
*/
class ElchFilter : public PostFilterBase
{
 public:
  /// constructor
  ElchFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};


/*!
 \brief Writer for (in)stationary heat conduction

 \author bborn
 \date 09/09
*/
class ThermoFilter : public PostFilterBase
{
 public:
  /// constructor
  ThermoFilter(PostField* field, std::string name, std::string heatfluxtype = "none",
      std::string tempgradtype = "none")
      : PostFilterBase(field, name), heatfluxtype_(heatfluxtype), tempgradtype_(tempgradtype)
  {
  }

 protected:
  void WriteAllResults(PostField* field) override;

  /*!
  \brief postprocess gauss point heatfluxes and write results
  \author originally by lw
  \date 11/09
  */
  void PostHeatflux(const std::string groupname, const std::string heatfluxtype);
  void WriteHeatflux(const std::string groupname, PostResult& result, const ResultType kind);

  std::string heatfluxtype_;  ///< type of heat flux output
  std::string tempgradtype_;  ///< type of spatial temperature gradient output
};

/*!
 \brief Writer for electromagnetic problems

 \author gravemeier
 \date 06/17
 */
class ElemagFilter : public PostFilterBase
{
 public:
  ElemagFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};

/// Writer for undefined problem types
/*
  Just write all the vectors we have.
 */
class AnyFilter : public PostFilterBase
{
 public:
  AnyFilter(PostField* field, std::string name) : PostFilterBase(field, name) {}

 protected:
  void WriteAllResults(PostField* field) override;
};

BACI_NAMESPACE_CLOSE

#endif  // POST_COMMON_SINGLE_FIELD_WRITERS_H
