/*!
 \file post_drt_ensight_single_field_writers.H

 \brief main routine of the Ensight filter

 <pre>
 Maintainer: Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/gerstenberger
 089 - 289-15236
 </pre>

 */


#ifndef POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_
#define POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_

#include "post_drt_ensight_writer.H"

/*!
 \brief Writer for structural problems
 */
class StructureEnsightWriter : public EnsightWriter
{
public:
  StructureEnsightWriter(PostField* field,
      std::string filename,
      std::string stresstype="none",
      std::string straintype="none") :
  EnsightWriter(field, filename), stresstype_(stresstype), straintype_(straintype)
  {
  }

protected:

  virtual void WriteAllResults(PostField* field);

  /*!
  \brief postprocess gauss point stresses and write results
  \author lw
  \date 02/08
  */
  void PostStress(const std::string groupname, const std::string stresstype);
  void WriteNodalStress(const std::string groupname,
                        PostResult& result);
  void WriteNodalStressStep(std::ofstream& file,
                            PostResult& result,
                            std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                            const std::string groupname,
                            const std::string name) const;
  void WriteElementCenterStress(const std::string groupname,
                                PostResult& result);
  void WriteElementCenterStressStep(std::ofstream& file,
                                    PostResult& result,
                                    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                    const std::string groupname,
                                    const std::string name) const;
  void WriteNodalEigenStress(const std::string groupname,
                             PostResult& result);
  void WriteNodalEigenStressStep(std::vector<Teuchos::RCP<std::ofstream> > files,
                                 PostResult& result,
                                 std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                 const std::string groupname,
                                 std::vector<std::string> name);
  void WriteElementCenterEigenStress(const std::string groupname,
                                     PostResult& result);
  void WriteElementCenterEigenStressStep(std::vector<Teuchos::RCP<std::ofstream> > files,
                                         PostResult& result,
                                         std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                         const std::string groupname,
                                         std::vector<std::string> name);

  std::string stresstype_;
  std::string straintype_;
};

/*!
 \brief Writer for fluid problems
 */
class FluidEnsightWriter : public EnsightWriter
{
public:
    FluidEnsightWriter(
            PostField* field,
            std::string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*!
 \brief Writer for xfluid problems
 */
class XFluidEnsightWriter : public EnsightWriter
{
public:
    XFluidEnsightWriter(
            PostField* field,
            std::string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*!
 \brief Writer for ale problems
 */
class AleEnsightWriter : public EnsightWriter
{
public:
    AleEnsightWriter(
            PostField* field,
            std::string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};


/*!
 \brief Writer for interface fields in XFEM
 */
class InterfaceEnsightWriter : public EnsightWriter
{
public:
    InterfaceEnsightWriter(
            PostField* field,
            std::string filename) :
        EnsightWriter(field, filename)
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
class ScaTraEnsightWriter : public EnsightWriter
{
public:
    /// constructor
    ScaTraEnsightWriter(
            PostField* field,
            std::string filename) :
        EnsightWriter(field, filename)
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
class ElchEnsightWriter : public EnsightWriter
{
public:
    /// constructor
    ElchEnsightWriter(
            PostField* field,
            std::string filename) :
        EnsightWriter(field, filename)
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
class ThermoEnsightWriter : public EnsightWriter
{
public:
    /// constructor
    ThermoEnsightWriter(PostField* field,
        std::string filename,
        std::string heatfluxtype="none",
        std::string tempgradtype="none") :
    EnsightWriter(field, filename), heatfluxtype_(heatfluxtype), tempgradtype_(tempgradtype)
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
    void WriteNodalHeatflux(const std::string groupname,
                            PostResult& result);
    void WriteNodalHeatfluxStep(std::ofstream& file,
                                PostResult& result,
                                std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                const std::string groupname,
                                const std::string name,
                                const int numdf) const;
    void WriteElementCenterHeatflux(const std::string groupname,
                                    PostResult& result);
    void WriteElementCenterHeatfluxStep(std::ofstream& file,
                                        PostResult& result,
                                        std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                        const std::string groupname,
                                        const std::string name,
                                        const int numdf) const;

    std::string heatfluxtype_;  ///< type of heat flux output
    std::string tempgradtype_;  ///< type of spatial temperature gradient output
};


/*!
 \brief Writer for particle problems

 \author ghamm
 \date 02/13
 */
class ParticleEnsightWriter : public EnsightWriter
{
public:
  ParticleEnsightWriter(
            PostField* field,
            std::string filename) :
  EnsightWriter(field, filename)
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


/// Writer for undefined problem types
/*
  Just write all the vectors we have.
 */
class AnyEnsightWriter : public EnsightWriter
{
public:
  AnyEnsightWriter(PostField* field, std::string filename)
    : EnsightWriter(field, filename) {}

protected:

  virtual void WriteAllResults(PostField* field);
};


#endif /*POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_*/
