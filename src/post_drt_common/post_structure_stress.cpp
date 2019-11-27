/*----------------------------------------------------------------------*/
/*! \file
\brief postprocessing of structural stresses

\level 1

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------*/



#include "post_writer_base.H"
#include "post_drt_common.H"
#include <string>
#include "post_single_field_writers.H"
#include "../linalg/linalg_utils_densematrix_eigen.H"
#include "../pss_full/pss_cpp.h"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureFilter::PostStress(const std::string groupname, const std::string stresstype)
{
  PostField* field = writer_->GetField();
  PostResult result = PostResult(field);
  result.next_result();

  if (!map_has_map(result.group(), groupname.c_str())) return;

  //--------------------------------------------------------------------
  // calculation and output of nodal stresses in xyz-reference frame
  //--------------------------------------------------------------------

  if (stresstype == "ndxyz")
  {
    WriteStress(groupname, result, nodebased);
  }

  //-------------------------------------------------------------------------
  // calculation and output of element center stresses in xyz-reference frame
  //-------------------------------------------------------------------------

  else if (stresstype == "cxyz")
  {
    WriteStress(groupname, result, elementbased);
  }

  //-----------------------------------------------------------------------------------
  // calculation and output of nodal and element center stresses in xyz-reference frame
  //-----------------------------------------------------------------------------------

  else if (stresstype == "cxyz_ndxyz")
  {
    WriteStress(groupname, result, nodebased);

    // reset result for postprocessing and output of element center stresses
    PostResult resultelestress = PostResult(field);
    resultelestress.next_result();
    WriteStress(groupname, resultelestress, elementbased);
  }

  else if (stresstype == "nd123")
  {
    WriteEigenStress(groupname, result, nodebased);
  }

  else if (stresstype == "c123")
  {
    WriteEigenStress(groupname, result, elementbased);
  }

  else if (stresstype == "c123_nd123")
  {
    WriteEigenStress(groupname, result, nodebased);

    // reset result for postprocessing and output of element center stresses
    PostResult resultelestress = PostResult(field);
    resultelestress.next_result();
    WriteEigenStress(groupname, resultelestress, elementbased);
  }

  else
  {
    dserror("Unknown stress/strain type");
  }

  return;
}



//--------------------------------------------------------------------
// calculate nodal stresses from gauss point stresses
//--------------------------------------------------------------------
struct WriteNodalStressStep : public SpecialFieldInterface
{
  WriteNodalStressStep(StructureFilter& filter) : filter_(filter) {}

  virtual std::vector<int> NumDfMap() { return std::vector<int>(1, 6); }

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name)
  {
    dsassert(name.size() == 1, "Unexpected number of names");

    const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);

    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();
    const Epetra_Map* noderowmap = dis->NodeRowMap();

    Teuchos::ParameterList p;
    p.set("action", "postprocess_stress");
    p.set("stresstype", "ndxyz");
    p.set("gpstressmap", data);
    Epetra_MultiVector* tmp = new Epetra_MultiVector(*noderowmap, 6, true);
    Teuchos::RCP<Epetra_MultiVector> nodal_stress = Teuchos::rcp(tmp);
    p.set("poststress", nodal_stress);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (nodal_stress == Teuchos::null)
    {
      dserror("vector containing element center stresses/strains not available");
    }

    filter_.GetWriter().WriteNodalResultStep(
        *files[0], nodal_stress, resultfilepos, groupname, name[0], 6);
  }

  StructureFilter& filter_;
};



//--------------------------------------------------------------------
// calculate element center stresses from gauss point stresses
//--------------------------------------------------------------------
struct WriteElementCenterStressStep : public SpecialFieldInterface
{
  WriteElementCenterStressStep(StructureFilter& filter) : filter_(filter) {}

  virtual std::vector<int> NumDfMap() { return std::vector<int>(1, 6); }

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name)
  {
    dsassert(name.size() == 1, "Unexpected number of names");
    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();
    const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);
    Teuchos::ParameterList p;
    p.set("action", "postprocess_stress");
    p.set("stresstype", "cxyz");
    p.set("gpstressmap", data);
    Teuchos::RCP<Epetra_MultiVector> elestress =
        Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()), 6));
    p.set("poststress", elestress);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (elestress == Teuchos::null)
    {
      dserror("vector containing element center stresses/strains not available");
    }
    filter_.GetWriter().WriteElementResultStep(
        *files[0], elestress, resultfilepos, groupname, name[0], 6, 0);
  }

  StructureFilter& filter_;
};



//--------------------------------------------------------------------
// Get structural rotation tensor R at element center     pfaller may17
//--------------------------------------------------------------------
struct WriteElementCenterRotation : public SpecialFieldInterface
{
  WriteElementCenterRotation(StructureFilter& filter) : filter_(filter) {}

  virtual std::vector<int> NumDfMap() { return std::vector<int>(1, 9); }

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name)
  {
    dsassert(name.size() == 1, "Unexpected number of names");
    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();
    const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);
    Teuchos::ParameterList p;
    p.set("action", "postprocess_stress");
    p.set("stresstype", "rotation");
    // this is not really stress at gp! we are misusing "postprocess_stress" here
    // we rather store the rotation tensor R at the element center
    p.set("gpstressmap", data);
    Teuchos::RCP<Epetra_MultiVector> elestress =
        Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()), 9));
    p.set("poststress", elestress);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (elestress == Teuchos::null)
    {
      dserror("vector containing element center stresses/strains not available");
    }
    filter_.GetWriter().WriteElementResultStep(
        *files[0], elestress, resultfilepos, groupname, name[0], 9, 0);
  }

  StructureFilter& filter_;
};



//--------------------------------------------------------------------
// calculate nodal membrane thickness from gauss point membrane thickness
//--------------------------------------------------------------------
struct WriteNodalMembraneThicknessStep : public SpecialFieldInterface
{
  WriteNodalMembraneThicknessStep(StructureFilter& filter) : filter_(filter) {}

  virtual std::vector<int> NumDfMap() { return std::vector<int>(1, 1); }

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name)
  {
    dsassert(name.size() == 1, "Unexpected number of names");

    const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);

    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();
    const Epetra_Map* noderowmap = dis->NodeRowMap();

    Teuchos::ParameterList p;
    p.set("action", "postprocess_thickness");
    p.set("optquantitytype", "ndxyz");
    p.set("gpthickmap", data);
    Epetra_MultiVector* tmp = new Epetra_MultiVector(*noderowmap, 1, true);
    Teuchos::RCP<Epetra_MultiVector> nodal_thickness = Teuchos::rcp(tmp);
    p.set("postthick", nodal_thickness);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (nodal_thickness == Teuchos::null)
    {
      dserror("vector containing nodal thickness not available");
    }

    filter_.GetWriter().WriteNodalResultStep(
        *files[0], nodal_thickness, resultfilepos, groupname, name[0], 1);
  }

  StructureFilter& filter_;
};



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureFilter::WriteStress(
    const std::string groupname, PostResult& result, const ResultType stresskind)
{
  std::string name;
  std::string out;

  if (groupname == "gauss_2PK_stresses_xyz")
  {
    name = "2PK_stresses_xyz";
    out = "2nd Piola-Kirchhoff stresses";
  }
  else if (groupname == "gauss_cauchy_stresses_xyz")
  {
    name = "cauchy_stresses_xyz";
    out = "Cauchy stresses";
  }
  else if (groupname == "gauss_2PK_coupling_stresses_xyz")
  {
    name = "2PK_coupling_stresses_xyz";
    out = "2nd Piola-Kirchhoff coupling stresses";
  }
  else if (groupname == "gauss_cauchy_coupling_stresses_xyz")
  {
    name = "cauchy_coupling_stresses_xyz";
    out = "Cauchy coupling stresses";
  }
  else if (groupname == "gauss_GL_strains_xyz")
  {
    name = "GL_strains_xyz";
    out = "Green-Lagrange strains";
  }
  else if (groupname == "gauss_EA_strains_xyz")
  {
    name = "EA_strains_xyz";
    out = "Euler-Almansi strains";
  }
  else if (groupname == "gauss_LOG_strains_xyz")
  {
    name = "LOG_strains_xyz";
    out = "Logarithmic strains";
  }
  else if (groupname == "gauss_pl_GL_strains_xyz")
  {
    name = "pl_GL_strains_xyz";
    out = "Plastic Green-Lagrange strains";
  }
  else if (groupname == "gauss_pl_EA_strains_xyz")
  {
    name = "pl_EA_strains_xyz";
    out = "Plastic Euler-Almansi strains";
  }
  else if (groupname == "gauss_membrane_thickness")
  {
    name = "membrane_thickness";
    out = "membrane thickness";
  }
  else if (groupname == "rotation")
  {
    name = "rotation";
    out = "structural rotation tensor";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  if (groupname == "rotation")
  {
    name = "element_" + name;
    WriteElementCenterRotation stresses(*this);
    writer_->WriteSpecialField(
        stresses, result, elementbased, groupname, std::vector<std::string>(1, name), out);
  }
  if (groupname == "gauss_membrane_thickness")
  {
    if (stresskind == nodebased)
    {
      name = "nodal_" + name;
      WriteNodalMembraneThicknessStep thickness(*this);
      writer_->WriteSpecialField(
          thickness, result, nodebased, groupname, std::vector<std::string>(1, name), out);
    }
    else if (stresskind == elementbased)
    {
      dserror("element based membrane thickness postprocessed anyway!");
    }
    else
      dserror("Unknown stress type");
  }
  else
  {
    if (stresskind == nodebased)
    {
      name = "nodal_" + name;
      WriteNodalStressStep stresses(*this);
      writer_->WriteSpecialField(
          stresses, result, nodebased, groupname, std::vector<std::string>(1, name), out);
    }
    else if (stresskind == elementbased)
    {
      name = "element_" + name;
      WriteElementCenterStressStep stresses(*this);
      writer_->WriteSpecialField(
          stresses, result, elementbased, groupname, std::vector<std::string>(1, name), out);
    }
    else
      dserror("Unknown stress type");
  }
}



//--------------------------------------------------------------------
// calculate nodal eigen stresses from gauss point stresses
//--------------------------------------------------------------------
struct WriteNodalEigenStressStep : public SpecialFieldInterface
{
  WriteNodalEigenStressStep(StructureFilter& filter) : filter_(filter) {}

  virtual std::vector<int> NumDfMap()
  {
    std::vector<int> map(3, 1);
    for (int i = 0; i < 3; ++i) map.push_back(3);
    return map;
  }

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name)
  {
    dsassert(name.size() == 6, "Unexpected number of names");

    const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);

    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();
    const Epetra_Map* noderowmap = dis->NodeRowMap();

    Teuchos::ParameterList p;
    p.set("action", "postprocess_stress");
    p.set("stresstype", "ndxyz");
    p.set("gpstressmap", data);
    Teuchos::RCP<Epetra_MultiVector> nodal_stress =
        Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 6, true));
    p.set("poststress", nodal_stress);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (nodal_stress == Teuchos::null)
    {
      dserror("vector containing nodal stresses/strains not available");
    }

    // Epetra_MultiVector with eigenvalues (3) and eigenvectors (9 components) in each row (=node)
    std::vector<Teuchos::RCP<Epetra_MultiVector>> nodal_eigen_val_vec(6);
    for (int i = 0; i < 3; ++i)
      nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 1));
    for (int i = 3; i < 6; ++i)
      nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3));

    const int numnodes = dis->NumMyRowNodes();
    bool threedim = true;
    if (result.field()->problem()->num_dim() == 2) threedim = false;

    // the three-dimensional case
    if (threedim)
    {
      for (int i = 0; i < numnodes; ++i)
      {
        Epetra_SerialDenseMatrix eigenvec(3, 3);
        Epetra_SerialDenseVector eigenval(3);

        eigenvec(0, 0) = (*((*nodal_stress)(0)))[i];
        eigenvec(0, 1) = (*((*nodal_stress)(3)))[i];
        eigenvec(0, 2) = (*((*nodal_stress)(5)))[i];
        eigenvec(1, 0) = eigenvec(0, 1);
        eigenvec(1, 1) = (*((*nodal_stress)(1)))[i];
        eigenvec(1, 2) = (*((*nodal_stress)(4)))[i];
        eigenvec(2, 0) = eigenvec(0, 2);
        eigenvec(2, 1) = eigenvec(1, 2);
        eigenvec(2, 2) = (*((*nodal_stress)(2)))[i];

        LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

        for (int d = 0; d < 3; ++d)
        {
          (*((*nodal_eigen_val_vec[d])(0)))[i] = eigenval(d);
          for (int e = 0; e < 3; ++e) (*((*nodal_eigen_val_vec[d + 3])(e)))[i] = eigenvec(e, d);
        }
      }
    }
    // the two-dimensional case
    else
    {
      for (int i = 0; i < numnodes; ++i)
      {
        Epetra_SerialDenseMatrix eigenvec(2, 2);
        Epetra_SerialDenseVector eigenval(2);

        eigenvec(0, 0) = (*((*nodal_stress)(0)))[i];
        eigenvec(0, 1) = (*((*nodal_stress)(3)))[i];
        eigenvec(1, 0) = eigenvec(0, 1);
        eigenvec(1, 1) = (*((*nodal_stress)(1)))[i];

        LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

        (*((*nodal_eigen_val_vec[0])(0)))[i] = eigenval(0);
        (*((*nodal_eigen_val_vec[1])(0)))[i] = eigenval(1);
        (*((*nodal_eigen_val_vec[2])(0)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[3])(0)))[i] = eigenvec(0, 0);
        (*((*nodal_eigen_val_vec[3])(1)))[i] = eigenvec(1, 0);
        (*((*nodal_eigen_val_vec[3])(2)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[4])(0)))[i] = eigenvec(0, 1);
        (*((*nodal_eigen_val_vec[4])(1)))[i] = eigenvec(1, 1);
        (*((*nodal_eigen_val_vec[4])(2)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[5])(0)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[5])(1)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[5])(2)))[i] = 0.0;
      }
    }

    for (int i = 0; i < 3; ++i)
      filter_.GetWriter().WriteNodalResultStep(
          *files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 1);
    for (int i = 3; i < 6; ++i)
      filter_.GetWriter().WriteNodalResultStep(
          *files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 3);
  }

  StructureFilter& filter_;
};



//--------------------------------------------------------------------
// calculate element center eigen stresses from gauss point stresses
//--------------------------------------------------------------------
struct WriteElementCenterEigenStressStep : public SpecialFieldInterface
{
  WriteElementCenterEigenStressStep(StructureFilter& filter) : filter_(filter) {}

  virtual std::vector<int> NumDfMap()
  {
    std::vector<int> map(3, 1);
    for (int i = 0; i < 3; ++i) map.push_back(3);
    return map;
  }

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name)
  {
    const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);

    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();

    Teuchos::ParameterList p;
    p.set("action", "postprocess_stress");
    p.set("stresstype", "cxyz");
    p.set("gpstressmap", data);
    Teuchos::RCP<Epetra_MultiVector> elestress =
        Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()), 6));
    p.set("poststress", elestress);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (elestress == Teuchos::null)
    {
      dserror("vector containing element center stresses/strains not available");
    }

    std::vector<Teuchos::RCP<Epetra_MultiVector>> nodal_eigen_val_vec(6);
    for (int i = 0; i < 3; ++i)
      nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()), 1));
    for (int i = 3; i < 6; ++i)
      nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()), 3));

    const int numnodes = dis->NumMyRowNodes();
    bool threedim = true;
    if (result.field()->problem()->num_dim() == 2) threedim = false;

    // the three-dimensional case
    if (threedim)
    {
      for (int i = 0; i < dis->NumMyRowElements(); ++i)
      {
        Epetra_SerialDenseMatrix eigenvec(3, 3);
        Epetra_SerialDenseVector eigenval(3);

        eigenvec(0, 0) = (*((*elestress)(0)))[i];
        eigenvec(0, 1) = (*((*elestress)(3)))[i];
        eigenvec(0, 2) = (*((*elestress)(5)))[i];
        eigenvec(1, 0) = eigenvec(0, 1);
        eigenvec(1, 1) = (*((*elestress)(1)))[i];
        eigenvec(1, 2) = (*((*elestress)(4)))[i];
        eigenvec(2, 0) = eigenvec(0, 2);
        eigenvec(2, 1) = eigenvec(1, 2);
        eigenvec(2, 2) = (*((*elestress)(2)))[i];

        LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

        for (int d = 0; d < 3; ++d)
        {
          (*((*nodal_eigen_val_vec[d])(0)))[i] = eigenval(d);
          for (int e = 0; e < 3; ++e) (*((*nodal_eigen_val_vec[d + 3])(e)))[i] = eigenvec(e, d);
        }
      }
    }
    // the two-dimensional case
    else
    {
      for (int i = 0; i < numnodes; ++i)
      {
        Epetra_SerialDenseMatrix eigenvec(2, 2);
        Epetra_SerialDenseVector eigenval(2);

        eigenvec(0, 0) = (*((*elestress)(0)))[i];
        eigenvec(0, 1) = (*((*elestress)(3)))[i];
        eigenvec(1, 0) = eigenvec(0, 1);
        eigenvec(1, 1) = (*((*elestress)(1)))[i];

        LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

        (*((*nodal_eigen_val_vec[0])(0)))[i] = eigenval(0);
        (*((*nodal_eigen_val_vec[1])(0)))[i] = eigenval(1);
        (*((*nodal_eigen_val_vec[2])(0)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[3])(0)))[i] = eigenvec(0, 0);
        (*((*nodal_eigen_val_vec[3])(1)))[i] = eigenvec(1, 0);
        (*((*nodal_eigen_val_vec[3])(2)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[4])(0)))[i] = eigenvec(0, 1);
        (*((*nodal_eigen_val_vec[4])(1)))[i] = eigenvec(1, 1);
        (*((*nodal_eigen_val_vec[4])(2)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[5])(0)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[5])(1)))[i] = 0.0;
        (*((*nodal_eigen_val_vec[5])(2)))[i] = 0.0;
      }
    }

    for (int i = 0; i < 3; ++i)
      filter_.GetWriter().WriteElementResultStep(
          *files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 1, 0);
    for (int i = 3; i < 6; ++i)
      filter_.GetWriter().WriteElementResultStep(
          *files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 3, 0);
  }

  StructureFilter& filter_;
};



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureFilter::WriteEigenStress(
    const std::string groupname, PostResult& result, const ResultType stresskind)
{
  std::vector<std::string> name(6);
  std::string out;

  if (groupname == "gauss_2PK_stresses_xyz")
  {
    name[0] = "2PK_stresses_eigenval1";
    name[1] = "2PK_stresses_eigenval2";
    name[2] = "2PK_stresses_eigenval3";
    name[3] = "2PK_stresses_eigenvec1";
    name[4] = "2PK_stresses_eigenvec2";
    name[5] = "2PK_stresses_eigenvec3";
    out = "principal 2nd Piola-Kirchhoff stresses";
  }
  else if (groupname == "gauss_cauchy_stresses_xyz")
  {
    name[0] = "cauchy_stresses_eigenval1";
    name[1] = "cauchy_stresses_eigenval2";
    name[2] = "cauchy_stresses_eigenval3";
    name[3] = "cauchy_stresses_eigenvec1";
    name[4] = "cauchy_stresses_eigenvec2";
    name[5] = "cauchy_stresses_eigenvec3";
    out = "principal Cauchy stresses";
  }
  else if (groupname == "gauss_2PK_coupling_stresses_xyz")
  {
    name[0] = "2PK_coupling_stresses_eigenval1";
    name[1] = "2PK_coupling_stresses_eigenval2";
    name[2] = "2PK_coupling_stresses_eigenval3";
    name[3] = "2PK_coupling_stresses_eigenvec1";
    name[4] = "2PK_coupling_stresses_eigenvec2";
    name[5] = "2PK_coupling_stresses_eigenvec3";
    out = "principal 2nd Piola-Kirchhoff coupling stresses";
  }
  else if (groupname == "gauss_cauchy_coupling_stresses_xyz")
  {
    name[0] = "cauchy_coupling_stresses_eigenval1";
    name[1] = "cauchy_coupling_stresses_eigenval2";
    name[2] = "cauchy_coupling_stresses_eigenval3";
    name[3] = "cauchy_coupling_stresses_eigenvec1";
    name[4] = "cauchy_coupling_stresses_eigenvec2";
    name[5] = "cauchy_coupling_stresses_eigenvec3";
    out = "principal Cauchy coupling stresses";
  }
  else if (groupname == "gauss_GL_strains_xyz")
  {
    name[0] = "GL_strains_eigenval1";
    name[1] = "GL_strains_eigenval2";
    name[2] = "GL_strains_eigenval3";
    name[3] = "GL_strains_eigenvec1";
    name[4] = "GL_strains_eigenvec2";
    name[5] = "GL_strains_eigenvec3";
    out = "principal Green-Lagrange strains";
  }
  else if (groupname == "gauss_EA_strains_xyz")
  {
    name[0] = "EA_strains_eigenval1";
    name[1] = "EA_strains_eigenval2";
    name[2] = "EA_strains_eigenval3";
    name[3] = "EA_strains_eigenvec1";
    name[4] = "EA_strains_eigenvec2";
    name[5] = "EA_strains_eigenvec3";
    out = "principal Euler-Almansi strains";
  }
  else if (groupname == "gauss_LOG_strains_xyz")
  {
    name[0] = "LOG_strains_eigenval1";
    name[1] = "LOG_strains_eigenval2";
    name[2] = "LOG_strains_eigenval3";
    name[3] = "LOG_strains_eigenvec1";
    name[4] = "LOG_strains_eigenvec2";
    name[5] = "LOG_strains_eigenvec3";
    out = "principal Logarithmic strains";
  }
  else if (groupname == "gauss_pl_GL_strains_xyz")
  {
    name[0] = "pl_GL_strains_eigenval1";
    name[1] = "pl_GL_strains_eigenval2";
    name[2] = "pl_GL_strains_eigenval3";
    name[3] = "pl_GL_strains_eigenvec1";
    name[4] = "pl_GL_strains_eigenvec2";
    name[5] = "pl_GL_strains_eigenvec3";
    out = "principal plastic Green-Lagrange strains";
  }
  else if (groupname == "gauss_pl_EA_strains_xyz")
  {
    name[0] = "pl_EA_strains_eigenval1";
    name[1] = "pl_EA_strains_eigenval2";
    name[2] = "pl_EA_strains_eigenval3";
    name[3] = "pl_EA_strains_eigenvec1";
    name[4] = "pl_EA_strains_eigenvec2";
    name[5] = "pl_EA_strains_eigenvec3";
    out = "principal plastic Euler-Almansi strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }


  if (stresskind == nodebased)
  {
    for (int i = 0; i < 6; ++i) name[i] = "nodal_" + name[i];
    WriteNodalEigenStressStep stresses(*this);
    writer_->WriteSpecialField(stresses, result, nodebased, groupname, name, out);
  }
  else if (stresskind == elementbased)
  {
    for (int i = 0; i < 6; ++i) name[i] = "element_" + name[i];
    WriteElementCenterEigenStressStep stresses(*this);
    writer_->WriteSpecialField(stresses, result, elementbased, groupname, name, out);
  }
  else
    dserror("Unknown heatflux type");
}
