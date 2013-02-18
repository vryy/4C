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
                            std::map<std::string, std::vector< std::ofstream::pos_type> >& resultfilepos,
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



#endif /*POST_DRT_GENERIC_SINGLE_FIELD_WRITER_H_*/
