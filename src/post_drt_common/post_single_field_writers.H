/*!
 \file post_single_field_writers.H

 \brief main routine of the Ensight filter

 \maintainer Martin Kronbichler
 */


#ifndef POST_SINGLE_FIELD_WRITER_H_
#define POST_SINGLE_FIELD_WRITER_H_

#include "post_writer_base.H"
#include "post_filter_base.H"

#include <Epetra_MultiVector.h>

class PostWriterBase;

/*!
 \brief Writer for structural problems
 */
class StructureFilter : public PostFilterBase
{
public:
  StructureFilter(PostField *field,
                  std::string name,
                  std::string stresstype="none",
                  std::string straintype="none")
  :
    PostFilterBase(field, name),
    stresstype_(stresstype),
    straintype_(straintype)
  {
  }

protected:
  virtual void WriteAllResults(PostField* field);

  virtual void WriteAllResultsOneTimeStep(PostResult& result, bool firststep, bool laststep );

  /*!
  \brief postprocess gauss point stresses and write results
  \author lw
  \date 02/08
  */
  void PostStress(const std::string groupname, const std::string stresstype);
  void WriteStress(const std::string groupname,
                   PostResult& result,
                   const ResultType stresskind);
  void WriteEigenStress(const std::string groupname,
                        PostResult& result,
                        const ResultType stresskind);

  std::string stresstype_;
  std::string straintype_;
};

/*!
 \brief Writer for fluid problems
 */
class FluidFilter : public PostFilterBase
{
public:
    FluidFilter(PostField *field,
                std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*!
 \brief Writer for xfluid problems
 */
class XFluidFilter : public PostFilterBase
{
public:
    XFluidFilter(PostField *field,
                 std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*!
 \brief Writer for ale problems
 */
class AleFilter : public PostFilterBase
{
public:
    AleFilter(PostField *field,
              std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};


/*!
 \brief Writer for interface fields in XFEM
 */
class InterfaceFilter : public PostFilterBase
{
public:
    InterfaceFilter(PostField *field,
                    std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
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
    LubricationFilter(PostField *field,
                      std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
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
    ScaTraFilter(PostField *field,
                 std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
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
    ElchFilter(PostField *field,
               std::string name)
    :
      PostFilterBase (field, name)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
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
    ThermoFilter(PostField *field,
                 std::string name,
                 std::string heatfluxtype="none",
                 std::string tempgradtype="none")
    :
      PostFilterBase (field, name),
      heatfluxtype_(heatfluxtype),
      tempgradtype_(tempgradtype)
    {
    }

protected:

    virtual void WriteAllResults(PostField* field);

    /*!
    \brief postprocess gauss point heatfluxes and write results
    \author originally by lw
    \date 11/09
    */
    void PostHeatflux(const std::string groupname, const std::string heatfluxtype);
    void WriteHeatflux(const std::string groupname,
                       PostResult& result,
                       const ResultType kind);

    std::string heatfluxtype_;  ///< type of heat flux output
    std::string tempgradtype_;  ///< type of spatial temperature gradient output
};


/*!
 \brief Writer for particle problems

 \author ghamm
 \date 02/13
 */
class ParticleFilter : public PostFilterBase
{
public:
  ParticleFilter(
            PostField* field,
            std::string name);

protected:

  virtual void WriteAllResults(PostField* field);
};

/*!
 \brief Writer for acoustical problems

 \author schoeder
 \date 01/14
 */
class AcouFilter : public PostFilterBase
{
public:
  AcouFilter(
            PostField* field,
            std::string name) :
        PostFilterBase (field, name)
  {
  }

protected:

  virtual void WriteAllResults(PostField* field);
};

/// Writer for undefined problem types
/*
  Just write all the vectors we have.
 */
class AnyFilter : public PostFilterBase
{
public:
  AnyFilter(PostField* field,
            std::string name)
    :
    PostFilterBase (field, name)
  {}

protected:

  virtual void WriteAllResults(PostField* field);
};


#endif /*POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_*/
