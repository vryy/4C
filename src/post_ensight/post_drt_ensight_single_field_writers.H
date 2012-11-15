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
                         string filename,
                         string stresstype="none",
                         string straintype="none") :
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
  void PostStress(const string groupname, const string stresstype);
  void WriteNodalStress(const string groupname,
                        PostResult& result);
  void WriteNodalStressStep(std::ofstream& file,
                            PostResult& result,
                            map<string, vector<ofstream::pos_type> >& resultfilepos,
                            const string groupname,
                            const string name) const;
  void WriteElementCenterStress(const string groupname,
                                PostResult& result);
  void WriteElementCenterStressStep(std::ofstream& file,
                                    PostResult& result,
                                    map<string, vector<ofstream::pos_type> >& resultfilepos,
                                    const string groupname,
                                    const string name) const;
  void WriteNodalEigenStress(const string groupname,
                             PostResult& result);
  void WriteNodalEigenStressStep(std::vector<Teuchos::RCP<ofstream> > files,
                                 PostResult& result,
                                 map<string, vector<ofstream::pos_type> >& resultfilepos,
                                 const string groupname,
                                 vector<string> name);
  void WriteElementCenterEigenStress(const string groupname,
                                     PostResult& result);
  void WriteElementCenterEigenStressStep(std::vector<Teuchos::RCP<ofstream> > files,
                                         PostResult& result,
                                         map<string, vector<ofstream::pos_type> >& resultfilepos,
                                         const string groupname,
                                         vector<string> name);

  string stresstype_;
  string straintype_;
};

/*!
 \brief Writer for fluid problems
 */
class FluidEnsightWriter : public EnsightWriter
{
public:
    FluidEnsightWriter(
            PostField* field,
            string filename) :
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
            string filename) :
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
            string filename) :
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
            string filename) :
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
            string filename) :
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
            string filename) :
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
                        string filename,
                        string heatfluxtype="none",
                        string tempgradtype="none") :
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
    void PostHeatflux(const string groupname, const string heatfluxtype);
    void WriteNodalHeatflux(const string groupname,
                            PostResult& result);
    void WriteNodalHeatfluxStep(std::ofstream& file,
                                PostResult& result,
                                map<string, vector<ofstream::pos_type> >& resultfilepos,
                                const string groupname,
                                const string name,
                                const int numdf) const;
    void WriteElementCenterHeatflux(const string groupname,
                                    PostResult& result);
    void WriteElementCenterHeatfluxStep(std::ofstream& file,
                                        PostResult& result,
                                        map<string, vector<ofstream::pos_type> >& resultfilepos,
                                        const string groupname,
                                        const string name,
                                        const int numdf) const;

    string heatfluxtype_;  ///< type of heat flux output
    string tempgradtype_;  ///< type of spatial temperature gradient output
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
