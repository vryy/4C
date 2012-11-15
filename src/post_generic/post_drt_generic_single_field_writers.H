/*!
 \file post_drt_generic_single_field_writers.H

 \brief a generic filter

 <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

 */


#ifndef POST_DRT_GENERIC_SINGLE_FIELD_WRITER_H_
#define POST_DRT_GENERIC_SINGLE_FIELD_WRITER_H_

#include "post_drt_generic_writer.H"
#include "../drt_lib/drt_condition_utils.H"

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
  void WriteNodalEigenStressStep(std::vector<RCP<ofstream> > files,
                                 PostResult& result,
                                 map<string, vector<ofstream::pos_type> >& resultfilepos,
                                 const string groupname,
                                 vector<string> name);
  void WriteElementCenterEigenStress(const string groupname,
                                     PostResult& result);
  void WriteElementCenterEigenStressStep(std::vector<RCP<ofstream> > files,
                                         PostResult& result,
                                         map<string, vector<ofstream::pos_type> >& resultfilepos,
                                         const string groupname,
                                         vector<string> name);

  string stresstype_;
  string straintype_;
};



#endif /*POST_DRT_GENERIC_SINGLE_FIELD_WRITER_H_*/
