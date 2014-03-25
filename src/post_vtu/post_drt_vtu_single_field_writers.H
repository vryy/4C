/*!
 \file post_drt_vtu_single_field_writers.H

 \brief field writers for VTU filter

 <pre>
 Maintainer: Martin Kronbichler
 kronbichler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15235
 </pre>

 */


#ifndef POST_DRT_VTU_SINGLE_FIELD_WRITER_H_
#define POST_DRT_VTU_SINGLE_FIELD_WRITER_H_

#include "post_drt_vtu_writer.H"
#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"

/*!
 \brief Writer for structural problems
 */
class StructureVtuWriter : public VtuWriter
{
public:
  StructureVtuWriter(PostField* field,
      std::string filename,
      std::string stresstype="none",
      std::string straintype="none") :
  VtuWriter(field, filename), stresstype_(stresstype), straintype_(straintype)
  {
  }

protected:

  virtual void WriteAllResults(PostField* field);

  /*!
  \brief postprocess gauss point stresses and write results
  \author lw
  \date 02/08
  */
  void PostStress(const std::string groupname, const std::string stresstype)
  {
    CurrentlyNotImplemented();
  }

  std::string stresstype_;
  std::string straintype_;
};

/*!
 \brief Writer for fluid problems
 */
class FluidVtuWriter : public VtuWriter
{
public:
    FluidVtuWriter(
            PostField* field,
            std::string filename) :
        VtuWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*!
 \brief Writer for xfluid problems
 */
class XFluidVtuWriter : public VtuWriter
{
public:
    XFluidVtuWriter(
            PostField* field,
            std::string filename) :
        VtuWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*!
 \brief Writer for ale problems
 */
class AleVtuWriter : public VtuWriter
{
public:
    AleVtuWriter(
            PostField* field,
            std::string filename) :
        VtuWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};


/*!
 \brief Writer for interface fields in XFEM
 */
class InterfaceVtuWriter : public VtuWriter
{
public:
    InterfaceVtuWriter(
            PostField* field,
            std::string filename) :
        VtuWriter(field, filename)
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
class ScaTraVtuWriter : public VtuWriter
{
public:
    /// constructor
    ScaTraVtuWriter(
            PostField* field,
            std::string filename) :
        VtuWriter(field, filename)
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
class ElchVtuWriter : public VtuWriter
{
public:
    /// constructor
    ElchVtuWriter(
            PostField* field,
            std::string filename) :
        VtuWriter(field, filename)
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
class ThermoVtuWriter : public VtuWriter
{
public:
    /// constructor
    ThermoVtuWriter(PostField* field,
        std::string filename,
        std::string heatfluxtype="none",
        std::string tempgradtype="none") :
    VtuWriter(field, filename),
    heatfluxtype_(heatfluxtype),
    tempgradtype_(tempgradtype)
    {
      CurrentlyNotImplemented();
    }

    virtual void WriteAllResults(PostField* field);

protected:
    /*!
    \brief postprocess gauss point heatfluxes and write results
    \author originally by lw
    \date 11/09
    */
    void PostHeatflux(const std::string groupname, const std::string heatfluxtype)
    {
      CurrentlyNotImplemented();
    }

    std::string heatfluxtype_;  ///< type of heat flux output
    std::string tempgradtype_;  ///< type of spatial temperature gradient output
};


/*!
 \brief Writer for particle problems

 \author ghamm
 \date 02/13
 */
class ParticleVtuWriter : public VtuWriter
{
public:
  ParticleVtuWriter(
            PostField* field,
            std::string filename) :
  VtuWriter(field, filename)
  {
  }

protected:

  virtual void WriteAllResults(PostField* field);

  virtual void WriteCells(
    std::ofstream& geofile,
    const Teuchos::RCP<DRT::Discretization> dis,
    const Teuchos::RCP<Epetra_Map>& proc0map
    ) const;
};

/*!
 \brief Writer for acoustical problems

 \author schoeder
 \date 01/14
 */
class AcouVtuWriter : public VtuWriter
{
public:
  AcouVtuWriter(
            PostField* field,
            std::string filename) :
  VtuWriter(field, filename)
  {
  }

protected:

  virtual void WriteAllResults(PostField* field);

};


#endif /*POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_*/
