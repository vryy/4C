/*----------------------------------------------------------------------------*/
/*! \file
\brief convert results for shell 8 elements

\level 3

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/

#include "../post_drt_common/post_drt_common.H"
#include "../drt_s8/shell8.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include <Teuchos_RCP.hpp>
#include <vector>

#include "../pss_full/pss_cpp.h"

class Converter
{
 public:
  Converter(PostField* field);

  void ConvertResults();

  void write_vector_result(std::string result_name, PostField* field, PostResult* result);
  void write_element_result(std::string result_name, PostField* field, PostResult* result);

 private:
  PostField* field_;
  Teuchos::RCP<DRT::Discretization> h8dis_;
  Teuchos::RCP<DRT::Discretization> s8dis_;
  Teuchos::RCP<IO::OutputControl> control_;
  Teuchos::RCP<IO::DiscretizationWriter> writer_;
  std::map<int, int> nodeids_;
  int maxgid_;
};


Converter::Converter(PostField* field) : field_(field)
{
  s8dis_ = field->discretization();
  h8dis_ =
      Teuchos::rcp(new DRT::Discretization("structure", Teuchos::rcp(new Epetra_SerialComm())));

  int numnodes = s8dis_->NumGlobalNodes();

  Epetra_SerialDenseMatrix coords(2 * numnodes, 3);
  for (int i = 0; i < numnodes; ++i)
  {
    DRT::Node* actnode = s8dis_->lRowNode(i);

    nodeids_[actnode->Id()] = i;

    coords(i, 0) = actnode->X()[0];
    coords(i, 1) = actnode->X()[1];
    coords(i, 2) = actnode->X()[2];
    coords(numnodes + i, 0) = actnode->X()[0];
    coords(numnodes + i, 1) = actnode->X()[1];
    coords(numnodes + i, 2) = actnode->X()[2];
  }

  std::set<int> donenodes;

  int numelements = s8dis_->NumGlobalElements();
  for (int i = 0; i < numelements; ++i)
  {
    DRT::Element* actele = s8dis_->lRowElement(i);
    DRT::ELEMENTS::Shell8* actshell = dynamic_cast<DRT::ELEMENTS::Shell8*>(actele);
    if (actshell == NULL) dserror("not a shell8 element");

    const Epetra_SerialDenseMatrix* dirs = actshell->GetDirectors();
    const std::vector<double>* thick = actshell->GetThickness();
    int enumnode = actshell->NumNode();
    for (int n = 0; n < enumnode; ++n)
    {
      DRT::Node* node = actshell->Nodes()[n];
      if (donenodes.find(node->Id()) == donenodes.end())
      {
        donenodes.insert(node->Id());
        int nlid = nodeids_[node->Id()];
        coords(nlid, 0) -= (*dirs)(0, n) * (*thick)[n];
        coords(nlid, 1) -= (*dirs)(1, n) * (*thick)[n];
        coords(nlid, 2) -= (*dirs)(2, n) * (*thick)[n];
        coords(numnodes + nlid, 0) += (*dirs)(0, n) * (*thick)[n];
        coords(numnodes + nlid, 1) += (*dirs)(1, n) * (*thick)[n];
        coords(numnodes + nlid, 2) += (*dirs)(2, n) * (*thick)[n];
      }
    }
  }

  // create new nodes
  maxgid_ = 0;
  for (std::map<int, int>::iterator i = nodeids_.begin(); i != nodeids_.end(); ++i)
  {
    int lid = i->first;
    int gid = i->second;
    maxgid_ = std::max(maxgid_, gid);
    double x[3];
    x[0] = coords(lid, 0);
    x[1] = coords(lid, 1);
    x[2] = coords(lid, 2);
    h8dis_->AddNode(Teuchos::rcp(new DRT::Node(gid, x, 0)));
  }
  for (std::map<int, int>::iterator i = nodeids_.begin(); i != nodeids_.end(); ++i)
  {
    int lid = i->first;
    int gid = i->second;
    double x[3];
    x[0] = coords(numnodes + lid, 0);
    x[1] = coords(numnodes + lid, 1);
    x[2] = coords(numnodes + lid, 2);
    h8dis_->AddNode(Teuchos::rcp(new DRT::Node(maxgid_ + gid + 1, x, 0)));
  }

  // create elements
  for (int i = 0; i < numelements; ++i)
  {
    DRT::Element* actele = s8dis_->lRowElement(i);
    DRT::ELEMENTS::Shell8* actshell = dynamic_cast<DRT::ELEMENTS::Shell8*>(actele);
    if (actshell == NULL) dserror("not a shell8 element");

    std::vector<int> h8nodeids;
    copy(actele->NodeIds(), actele->NodeIds() + 4, back_inserter(h8nodeids));
    for (int j = 0; j < 4; ++j)
    {
      h8nodeids.push_back(h8nodeids[j] + numnodes);
    }
    Teuchos::RCP<DRT::ELEMENTS::So_hex8> h8 =
        Teuchos::rcp(new DRT::ELEMENTS::So_hex8(actshell->Id(), 0));
    h8->SetNodeIds(8, &h8nodeids[0]);
    h8dis_->AddElement(h8);
  }

  h8dis_->FillComplete();

  control_ = Teuchos::rcp(new IO::OutputControl(h8dis_->Comm(), "Structure",
      ShapeFunctionType::shapefunction_polynomial, "generated", "s8convert", 3, 0, 1000, 1));
  writer_ = h8dis_->Writer();
  writer_->SetOutput(control_);
  writer_->WriteMesh(0, 0);
}


void Converter::ConvertResults()
{
  PostResult result = PostResult(field_);
  while (result.next_result())
  {
    writer_->NewStep(result.step(), result.time());

    if (map_has_map(result.group(), "displacement"))
    {
      write_vector_result("displacement", field_, &result);
    }
    if (map_has_map(result.group(), "velocity"))
    {
      write_vector_result("velocity", field_, &result);
    }
    if (map_has_map(result.group(), "acceleration"))
    {
      write_vector_result("acceleration", field_, &result);
    }
    if (map_has_map(result.group(), "gauss_EA_strains_xyz"))
    {
      write_element_result("gauss_EA_strains_xyz", field_, &result);
    }
  }
}


void Converter::write_vector_result(std::string result_name, PostField* field, PostResult* result)
{
  Teuchos::RCP<Epetra_Vector> s8data = result->read_result(result_name);
  const Epetra_BlockMap& s8map = s8data->Map();

  Teuchos::RCP<Epetra_Vector> h8data = Teuchos::rcp(new Epetra_Vector(*h8dis_->DofRowMap()));
  const Epetra_BlockMap& h8map = h8data->Map();

  int numnodes = s8dis_->NumGlobalNodes();
  for (int i = 0; i < numnodes; ++i)
  {
    DRT::Node* s8node = s8dis_->lRowNode(i);
    int llid = nodeids_[s8node->Id()];
    int ulid = llid + numnodes;
    std::vector<int> s8dof = s8dis_->Dof(s8node);

    DRT::Node* lh8node = h8dis_->lRowNode(llid);
    DRT::Node* uh8node = h8dis_->lRowNode(ulid);
    std::vector<int> lh8dof = h8dis_->Dof(lh8node);
    std::vector<int> uh8dof = h8dis_->Dof(uh8node);
    for (int i = 0; i < field->problem()->num_dim(); ++i)
    {
      (*h8data)[h8map.LID(lh8dof[i])] =
          (*s8data)[s8map.LID(s8dof[i])] - (*s8data)[s8map.LID(s8dof[i + 3])];
      (*h8data)[h8map.LID(uh8dof[i])] =
          (*s8data)[s8map.LID(s8dof[i])] + (*s8data)[s8map.LID(s8dof[i + 3])];
    }
  }

  writer_->WriteVector(result_name, h8data);
}


void Converter::write_element_result(std::string result_name, PostField* field, PostResult* result)
{
  Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> s8data =
      result->read_result_serialdensematrix(result_name);

  const Epetra_Map* h8map = h8dis_->ElementRowMap();

  DRT::PackBuffer h8data;
  for (int i = 0; i < h8map->NumMyElements(); ++i)
  {
    int gid = h8map->GID(i);
    if (s8data->find(gid) == s8data->end()) dserror("gid %d not in s8data", gid);
    if ((*s8data)[gid] == Teuchos::null) dserror("null matrix in s8data");
    DRT::ParObject::AddtoPack(h8data, *(*s8data)[gid]);
  }

  h8data.StartPacking();

  for (int i = 0; i < h8map->NumMyElements(); ++i)
  {
    int gid = h8map->GID(i);
    if (s8data->find(gid) == s8data->end()) dserror("gid %d not in s8data", gid);
    if ((*s8data)[gid] == Teuchos::null) dserror("null matrix in s8data");
    DRT::ParObject::AddtoPack(h8data, *(*s8data)[gid]);
  }

  writer_->WriteVector(result_name, h8data(), *h8map);
}

int main(int argc, char** argv)
{
  Teuchos::CommandLineProcessor My_CLP;
  My_CLP.setDocString("Post DRT shell8 conversion Filter\n");

  PostProblem problem(My_CLP, argc, argv);

  switch (problem.Problemtype())
  {
    case prb_structure:
    {
      PostField* field = problem.get_discretization(0);
      Converter conv(field);
      conv.ConvertResults();
      break;
    }
    default:
      dserror("problem type %d not yet supported", problem.Problemtype());
  }

  return 0;
}
