/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of methods to modify patient specific geometries

\maintainer Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/


#include "patspec.H"
#include "../drt_mat/material.H"
#include "../drt_mat/elasthyper.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <iostream>
#include "../drt_io/io_pstream.H"  // has to go before io.H
#include "../drt_io/io.H"
#include <Epetra_Time.h>

//#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_matelast/elast_summand.H"
#include "../drt_matelast/elast_isovolaaagasser.H"


/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::PatientSpecificGeometry(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Teuchos::ParameterList> params)
{
  if (!dis->Comm().MyPID())
  {
    IO::cout << "------------------------------------------------------------\n";
    IO::cout << "Entering patient specific structural preprocessing (PATSPEC)\n";
    IO::cout << "\n";
  }

  //------------- test discretization for presence of the Gasser ILT material
  int lfoundit = 0;
  for (int i = 0; i < dis->ElementRowMap()->NumMyElements(); ++i)
  {
    INPAR::MAT::MaterialType type = dis->lRowElement(i)->Material()->MaterialType();
    if (type == INPAR::MAT::m_aaagasser)
    {
      lfoundit = 1;
      break;
    }

    if (type == INPAR::MAT::m_aaa_mixedeffects)
    {
      lfoundit = 1;
      break;
    }

    if (type == INPAR::MAT::m_elasthyper)
    {
      lfoundit = 1;
      break;
    }
  }

  int gfoundit = 0;
  dis->Comm().SumAll(&lfoundit, &gfoundit, 1);

  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();

  int maxhulumen = pslist.get<int>("MAXHULUMEN");
  params->set("max hu lumen", maxhulumen);
  if (!dis->Comm().MyPID())
  {
    if (maxhulumen == 0)
      IO::cout << "max HU in lumen not set " << IO::endl;
    else
      IO::cout << "max HU in lumen " << maxhulumen << IO::endl;
  }

  int calc_strength = DRT::INPUT::IntegralValue<int>(pslist, "CALCSTRENGTH");
  bool calc_shear = DRT::INPUT::IntegralValue<bool>(pslist, "CALCINNERRADIUS");

  if (gfoundit or calc_strength)
  {
    if (!dis->Comm().MyPID()) IO::cout << "Computing distance functions...\n";
    PATSPEC::ComputeEleNormalizedLumenDistance(dis, params);
    PATSPEC::ComputeEleLocalRadius(dis);
  }

  if (calc_strength)
  {
    if (!dis->Comm().MyPID()) IO::cout << "Computing strength model...\n";
    PATSPEC::ComputeEleStrength(dis, params);
  }

  if (calc_shear)
  {
    if (!dis->Comm().MyPID()) IO::cout << "Computing inner radius...\n";
    PATSPEC::InitializeMappingInnerSurface(dis);
  }


  if (!dis->Comm().MyPID())
  {
    IO::cout << "\n";
    IO::cout << "Leaving patient specific structural preprocessing (PATSPEC)\n";
    IO::cout << "-----------------------------------------------------------\n";
  }


  return;
}


/*----------------------------------------------------------------------*
 |                                                          amaier 07/11|
 *----------------------------------------------------------------------*/
void PATSPEC::ComputeEleStrength(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Teuchos::ParameterList> params)
{
  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  double subrendia = pslist.get<double>("AAA_SUBRENDIA");
  int is_male = DRT::INPUT::IntegralValue<int>(pslist, "MALE_PATIENT");
  int has_familyhist = DRT::INPUT::IntegralValue<int>(pslist, "FAMILYHIST");
  double spatialconst = 922000.0;  // spatial constant strength contrib
                                   // acc. Vande Geest [Pa]

  double maxiltthick = params->get<double>("max ilt thick");


  if (!dis->Comm().MyPID())
  {
    if (subrendia == 22.01)
      IO::cout << "Subrenal diameter not specified, taking default value (22mm).\n";
    else
      IO::cout << "Subrenal diameter " << subrendia << " mm" << IO::endl;

    if (is_male)
      IO::cout << "Male patient " << IO::endl;
    else
      IO::cout << "Female patient." << IO::endl;

    if (!has_familyhist)
      IO::cout << "No AAA familiy history." << IO::endl;
    else
      IO::cout << "Patient has AAA family history!!" << IO::endl;
  }


  if (has_familyhist) spatialconst -= 213000;
  if (!is_male) spatialconst -= 193000;

  Teuchos::RCP<Epetra_Vector> elestrength = LINALG::CreateVector(*(dis->ElementRowMap()), true);

  std::vector<DRT::Condition*> mypatspeccond;
  dis->GetCondition("PatientSpecificData", mypatspeccond);

  if (!mypatspeccond.size()) dserror("Cannot find the Patient Specific Data Conditions :-(");


  for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
  {
    const Epetra_Vector* ilt = mypatspeccond[i]->Get<Epetra_Vector>("normalized ilt thickness");
    if (ilt)
    {
      for (int j = 0; j < elestrength->MyLength(); ++j)
      {
        // contribution of local ilt thickness to strength; lower and
        // upper bounds for the ilt thickness are 0 and 36
        // Careful! ilt thickness is still normalized! Multiply with max
        // ilt thickness.
        (*elestrength)[j] =
            spatialconst -
            379000 *
                (std::pow((std::max(0.,
                               std::min(36., (*ilt)[ilt->Map().LID(dis->ElementRowMap()->GID(j))] *
                                                 maxiltthick)) /
                              10),
                     0.5) -
                    0.81);  // from Vande Geest strength formula
      }
    }
  }

  for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
  {
    const Epetra_Vector* locrad = mypatspeccond[i]->Get<Epetra_Vector>("local radius");
    if (locrad)
    {
      for (int j = 0; j < elestrength->MyLength(); ++j)
      {
        // contribution of local diameter to strength; lower and upper bounds for the normalized
        // diameter are 1.0 and 3.9
        (*elestrength)[j] -=
            156000 *
            (std::max(1.,
                 std::min(3.9,
                     2 * (*locrad)[locrad->Map().LID(dis->ElementRowMap()->GID(j))] / subrendia)) -
                2.46);  // from Vande Geeststrength formula
      }
    }
  }

  Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis->ElementColMap()), true);
  LINALG::Export(*elestrength, *tmp);
  elestrength = tmp;

  // put this column vector of element strength in a condition and store in
  // discretization
  Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
      new DRT::Condition(0, DRT::Condition::PatientSpecificData, false, DRT::Condition::Volume));
  cond->Add("elestrength", *elestrength);

  // check whether discretization has been filled before putting condition to dis
  bool filled = dis->Filled();
  dis->SetCondition("PatientSpecificData", cond);
  if (filled && !dis->Filled()) dis->FillComplete();

  if (!dis->Comm().MyPID()) IO::cout << "Strength calculation completed." << IO::endl;

  return;
}

/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::ComputeEleNormalizedLumenDistance(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Teuchos::ParameterList> params)
{
  // find out whether we have a orthopressure or FSI condition
  std::vector<DRT::Condition*> conds;
  std::vector<DRT::Condition*> ortho(0);
  std::vector<DRT::Condition*> fsi(0);
  dis->GetCondition("SurfaceNeumann", ortho);
  for (int i = 0; i < (int)ortho.size(); ++i)
  {
    if (ortho[i]->GType() != DRT::Condition::Surface) continue;
    const std::string* type = ortho[i]->Get<std::string>("type");
    if (*type == "neum_orthopressure" || *type == "neum_pseudo_orthopressure")
      conds.push_back(ortho[i]);
  }

  dis->GetCondition("FSICoupling", fsi);
  for (int i = 0; i < (int)fsi.size(); ++i)
  {
    if (fsi[i]->GType() != DRT::Condition::Surface) continue;
    conds.push_back(fsi[i]);
  }

  if (!conds.size()) dserror("There is no orthopressure nor FSI condition in this discretization");

  // measure time as there is a brute force search in here
  Epetra_Time timer(dis->Comm());
  timer.ResetStartTime();

  // start building the distance function
  std::set<int> allnodes;
  for (int i = 0; i < (int)conds.size(); ++i)
  {
    const std::vector<int>* nodes = conds[i]->Nodes();
    if (!nodes) dserror("Cannot find node ids in condition");
    for (int j = 0; j < (int)nodes->size(); ++j) allnodes.insert((*nodes)[j]);
  }

  // create coordinates for all these nodes
  const int nnodes = (int)allnodes.size();
  std::vector<double> lcoords(nnodes * 3, 0.0);
  std::vector<double> gcoords(nnodes * 3, 0.0);
  std::set<int>::iterator fool;
  int count = 0;
  for (fool = allnodes.begin(); fool != allnodes.end(); ++fool)
  {
    if (!dis->NodeRowMap()->MyGID(*fool))
    {
      count++;
      continue;
    }
    DRT::Node* node = dis->gNode(*fool);
    lcoords[count * 3] = node->X()[0];
    lcoords[count * 3 + 1] = node->X()[1];
    lcoords[count * 3 + 2] = node->X()[2];
    count++;
  }
  dis->Comm().SumAll(&lcoords[0], &gcoords[0], nnodes * 3);
  lcoords.clear();
  allnodes.clear();

  // compute distance of all of my nodes to these nodes
  // create vector for nodal values of ilt thickness
  // WARNING: This is a brute force expensive minimum distance search!
  const Epetra_Map* nrowmap = dis->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> iltthick = LINALG::CreateVector(*nrowmap, true);
  for (int i = 0; i < nrowmap->NumMyElements(); ++i)
  {
    const double* x = dis->gNode(nrowmap->GID(i))->X();
    // loop nodes from the condition and find minimum distance
    double mindist = 1.0e12;
    for (int j = 0; j < nnodes; ++j)
    {
      double* xorth = &gcoords[j * 3];
      double dist =
          sqrt((x[0] - xorth[0]) * (x[0] - xorth[0]) + (x[1] - xorth[1]) * (x[1] - xorth[1]) +
               (x[2] - xorth[2]) * (x[2] - xorth[2]));
      if (dist < mindist) mindist = dist;
    }
    (*iltthick)[i] = mindist;
  }
  gcoords.clear();

  // now we need to compute the maximum ILT thickness
  // historically the maximum ilt thickness was computed based on distance to orthopressure/fsi
  // surface on luminal side of the ilt. From the maximum an approximate wall thickness
  // of 1.0 mm is subtrated (hardcoded). This obviously can cause problems
  // when the wall thickness is not constant e.g. during UQ analysis.
  // Therefore, a new more accurate method was added which needs the luminal
  // and the outer surface of the ILT as AAA surface condition.
  // This can be used by setting the CALC_ACCURATE_MAX_ILT_THICK flag to yes in
  // the input file.
  // The old way is kept here only to allow evaluation of the AAA database.
  // As a safety measure the accurate method must be used if the problem type is UQ

  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  bool accurate_ilt_thick_calc =
      DRT::INPUT::IntegralValue<int>(pslist, "CALC_ACCURATE_MAX_ILT_THICK");

  double maxiltthick = 0.0;
  if (accurate_ilt_thick_calc)
  {
    maxiltthick = ComputeMaxILTThickness(dis);
  }
  else
  {
    if (DRT::Problem::Instance()->GetProblemType() == prb_uq)
      dserror(
          "UQ requires accurate computation of max ILT thickness, Set CALC_ACCURATE_MAX_ILT_THICK "
          "to yes");
    else
      iltthick->MaxValue(&maxiltthick);
    maxiltthick -= 1.0;  // subtract an approximate arterial wall thickness
    if (!dis->Comm().MyPID())
    {
      IO::cout << "WARNING: WALL THICKNESS IS CURRENTLY HARDCODED TO 1.0 mm FOR ILT THICK CALC"
               << IO::endl;
      IO::cout << "WARNING: THINK ABOUT SWITCHING TO ACCURATE ILT THICK CALCULATION" << IO::endl;
    }
  }

  iltthick->Scale(1.0 / maxiltthick);
  if (!dis->Comm().MyPID()) IO::cout << "Max ILT thickness " << maxiltthick << IO::endl;
  params->set("max ilt thick", maxiltthick);

  // export nodal distances to column map
  Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis->NodeColMap()), true);
  LINALG::Export(*iltthick, *tmp);
  iltthick = tmp;

  // compute element-wise mean distance from nodal distance
  Teuchos::RCP<Epetra_Vector> iltele = LINALG::CreateVector(*(dis->ElementRowMap()), true);
  for (int i = 0; i < dis->ElementRowMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = dis->gElement(dis->ElementRowMap()->GID(i));
    double mean = 0.0;
    for (int j = 0; j < actele->NumNode(); ++j)
      mean += (*iltthick)[iltthick->Map().LID(actele->Nodes()[j]->Id())];
    mean /= actele->NumNode();
    (*iltele)[i] = mean;
  }

  tmp = LINALG::CreateVector(*(dis->ElementColMap()), true);
  LINALG::Export(*iltele, *tmp);
  iltele = tmp;

  // put this column vector of element ilt thickness in a condition and store in
  // discretization
  Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
      new DRT::Condition(0, DRT::Condition::PatientSpecificData, false, DRT::Condition::Volume));
  cond->Add("normalized ilt thickness", *iltele);

  // check materials if they need normalized ilt thickness

  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();
  std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
  for (curr = mats.begin(); curr != mats.end(); ++curr)
  {
    if (curr->second->Type() == INPAR::MAT::mes_isovolaaagasser)
    {
      // set thrombus thickness
      MAT::ELASTIC::PAR::IsoVolAAAGasser* params =
          dynamic_cast<MAT::ELASTIC::PAR::IsoVolAAAGasser*>(curr->second->Parameter());
      curr->second->Parameter()->SetParameter(params->normdist, iltele);
    }
  }
  // check whether discretization has been filled before putting condition to dis
  bool filled = dis->Filled();
  dis->SetCondition("PatientSpecificData", cond);
  if (filled && !dis->Filled()) dis->FillComplete();

  if (!dis->Comm().MyPID())
    IO::cout << "Normalized ILT thickness computed in " << timer.ElapsedTime() << " sec"
             << IO::endl;


  return;
}

/*----------------------------------------------------------------------*
 |                                                           maier 11/10|
 *----------------------------------------------------------------------*/
void PATSPEC::ComputeEleLocalRadius(Teuchos::RCP<DRT::Discretization> dis)
{
  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  std::string filename = pslist.get<std::string>("CENTERLINEFILE");

  if (filename == "name.txt")
  {
    if (!dis->Comm().MyPID()) IO::cout << "No centerline file provided" << IO::endl;
    // set element-wise mean distance to zero
    Teuchos::RCP<Epetra_Vector> locradele = LINALG::CreateVector(*(dis->ElementRowMap()), true);
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis->ElementColMap()), true);
    LINALG::Export(*locradele, *tmp);
    locradele = tmp;

    // put this column vector of element local radius in a condition and store in
    // discretization
    Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
        new DRT::Condition(0, DRT::Condition::PatientSpecificData, false, DRT::Condition::Volume));
    cond->Add("local radius", *locradele);

    // check whether discretization has been filled before putting condition to dis
    bool filled = dis->Filled();
    dis->SetCondition("PatientSpecificData", cond);
    if (filled && !dis->Filled()) dis->FillComplete();

    if (!dis->Comm().MyPID()) IO::cout << "No local radii computed" << IO::endl;

    return;
  }

  std::vector<double> clcoords = PATSPEC::GetCenterline(filename);

  // compute distance of all of my nodes to the clpoints
  // create vector for nodal values of ilt thickness
  // WARNING: This is a brute force expensive minimum distance search!
  const Epetra_Map* nrowmap = dis->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> localrad = LINALG::CreateVector(*nrowmap, true);
  for (int i = 0; i < nrowmap->NumMyElements(); ++i)
  {
    const double* x = dis->gNode(nrowmap->GID(i))->X();
    // loop nodes from the condition and find minimum distance
    double mindist = 1.0e12;
    for (int j = 0; j < (int)(clcoords.size() / 3); ++j)
    {
      double* xcl = &clcoords[j * 3];
      double dist = sqrt((x[0] - xcl[0]) * (x[0] - xcl[0]) + (x[1] - xcl[1]) * (x[1] - xcl[1]) +
                         (x[2] - xcl[2]) * (x[2] - xcl[2]));
      if (dist < mindist) mindist = dist;
    }
    (*localrad)[i] = mindist;
  }
  clcoords.clear();
  // max local rad just for information purpose:
  double maxlocalrad;
  localrad->MaxValue(&maxlocalrad);
  if (!dis->Comm().MyPID()) IO::cout << "Max local radius:  " << maxlocalrad << IO::endl;

  // export nodal distances to column map
  Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis->NodeColMap()), true);
  LINALG::Export(*localrad, *tmp);
  localrad = tmp;

  // compute element-wise mean distance from nodal distance
  Teuchos::RCP<Epetra_Vector> locradele = LINALG::CreateVector(*(dis->ElementRowMap()), true);
  for (int i = 0; i < dis->ElementRowMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = dis->gElement(dis->ElementRowMap()->GID(i));
    double mean = 0.0;
    for (int j = 0; j < actele->NumNode(); ++j)
      mean += (*localrad)[localrad->Map().LID(actele->Nodes()[j]->Id())];
    mean /= actele->NumNode();
    (*locradele)[i] = mean;
  }

  tmp = LINALG::CreateVector(*(dis->ElementColMap()), true);
  LINALG::Export(*locradele, *tmp);
  locradele = tmp;

  // put this column vector of element local radius in a condition and store in
  // discretization
  Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
      new DRT::Condition(0, DRT::Condition::PatientSpecificData, false, DRT::Condition::Volume));
  cond->Add("local radius", *locradele);

  // check whether discretization has been filled before putting condition to dis
  bool filled = dis->Filled();
  dis->SetCondition("PatientSpecificData", cond);
  if (filled && !dis->Filled()) dis->FillComplete();

  if (!dis->Comm().MyPID()) IO::cout << "Local radii computed. " << IO::endl;

  return;
}

/*----------------------------------------------------------------------*
 |                                                           tinkl 11/13|
 *----------------------------------------------------------------------*/
void PATSPEC::InitializeMappingInnerSurface(Teuchos::RCP<DRT::Discretization> dis)
{
  // find out whether we have a orthopressure or FSI condition
  std::vector<DRT::Condition*> conds;
  std::vector<DRT::Condition*> ortho(0);
  std::vector<DRT::Condition*> fsi(0);
  dis->GetCondition("SurfaceNeumann", ortho);
  for (int i = 0; i < (int)ortho.size(); ++i)
  {
    if (ortho[i]->GType() != DRT::Condition::Surface) continue;
    const std::string* type = ortho[i]->Get<std::string>("type");
    if (*type == "neum_orthopressure" || *type == "neum_pseudo_orthopressure")
      conds.push_back(ortho[i]);
  }

  dis->GetCondition("FSICoupling", fsi);
  for (int i = 0; i < (int)fsi.size(); ++i)
  {
    if (fsi[i]->GType() != DRT::Condition::Surface) continue;
    conds.push_back(fsi[i]);
  }

  if (!conds.size()) dserror("There is no orthopressure nor FSI condition in this discretization");

  // collect all nodes
  std::set<int> allnodes;
  for (int i = 0; i < (int)conds.size(); ++i)
  {
    const std::vector<int>* nodes = conds[i]->Nodes();
    if (!nodes) dserror("Cannot find node ids in condition");
    for (int j = 0; j < (int)nodes->size(); ++j) allnodes.insert((*nodes)[j]);
  }

  // create coordinates for all these nodes
  const int nnodes = (int)allnodes.size();
  std::vector<double> lcoords(nnodes * 3, 0.0);
  std::vector<double> gcoords(nnodes * 3, 0.0);
  std::set<int>::iterator fool;
  int count = 0;
  for (fool = allnodes.begin(); fool != allnodes.end(); ++fool)
  {
    if (!dis->NodeRowMap()->MyGID(*fool))
    {
      count++;
      continue;
    }
    DRT::Node* node = dis->gNode(*fool);
    lcoords[count * 3] = node->X()[0];
    lcoords[count * 3 + 1] = node->X()[1];
    lcoords[count * 3 + 2] = node->X()[2];
    count++;
  }
  dis->Comm().SumAll(&lcoords[0], &gcoords[0], nnodes * 3);
  lcoords.clear();
  allnodes.clear();

  // compute distance of all of my elements to these nodes
  // WARNING: This is a brute force expensive minimum distance search!
  Teuchos::RCP<Epetra_Vector> innerid = LINALG::CreateVector(*(dis->NodeRowMap()), true);
  for (int i = 0; i < dis->NodeRowMap()->NumMyElements(); ++i)
  {
    const double* x = dis->gNode(dis->NodeRowMap()->GID(i))->X();
    // loop nodes from the condition and find minimum distance
    double mindist = 1.0e12;
    int inneridloc = 0;
    for (int j = 0; j < nnodes; ++j)
    {
      double* xorth = &gcoords[j * 3];
      double dist =
          sqrt((x[0] - xorth[0]) * (x[0] - xorth[0]) + (x[1] - xorth[1]) * (x[1] - xorth[1]) +
               (x[2] - xorth[2]) * (x[2] - xorth[2]));
      if (dist < mindist)
      {
        mindist = dist;
        inneridloc = j;
      }
    }
    (*innerid)[i] = inneridloc;
  }
  gcoords.clear();

  // export nodal distances to column map
  Teuchos::RCP<Epetra_Vector> tmpnode = LINALG::CreateVector(*(dis->NodeColMap()), true);
  LINALG::Export(*innerid, *tmpnode);
  innerid = tmpnode;

  // discretization
  Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
      new DRT::Condition(0, DRT::Condition::PatientSpecificData, false, DRT::Condition::Volume));
  cond->Add("inner surface node id", *innerid);

  // check whether discretization has been filled before putting condition to dis
  bool filled = dis->Filled();
  dis->SetCondition("PatientSpecificData", cond);
  if (filled && !dis->Filled()) dis->FillComplete();

  return;
}

/*----------------------------------------------------------------------*
 |                                                           tinkl 11/13|
 *----------------------------------------------------------------------*/
void PATSPEC::ComputeEleInnerRadius(Teuchos::RCP<DRT::Discretization> dis)
{
  // calculation only needed when PatientSpecificData available
  std::vector<DRT::Condition*> mypatspeccond;
  dis->GetCondition("PatientSpecificData", mypatspeccond);
  if (!mypatspeccond.size()) return;
  // check if inner id has been computed otherwise leave
  for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
  {
    const Epetra_Vector* fool = mypatspeccond[i]->Get<Epetra_Vector>("inner surface node id");
    if (!fool) return;
  }

  // get centerline
  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  bool lin_centerline = DRT::INPUT::IntegralValue<bool>(pslist, "LINEARCENTERLINE");
  std::vector<double> clcoords;
  std::vector<double> direction;
  double normdirection = 0.0;
  std::vector<double> support;
  if (lin_centerline)
  {
    // centerline is straight line
    double word;
    std::istringstream directionstream(
        Teuchos::getNumericStringParameter(pslist, "CENTERLINEDIRECTION"));
    while (directionstream >> word) direction.push_back(word);
    if ((int)direction.size() != 3)
      dserror("Length of direction is wrong %f!\n", (int)direction.size());

    normdirection = sqrt(
        direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);

    std::istringstream supportstream(Teuchos::getNumericStringParameter(pslist, "CENTERLINEPOINT"));
    while (supportstream >> word) support.push_back(word);
    if ((int)support.size() != 3) dserror("Length of direction is wrong!\n");
  }
  else
  {
    std::string filename = pslist.get<std::string>("CENTERLINEFILE");
    if (filename == "name.txt")
    {
      dserror("no centerline file provided");
      // this function makes no sense without centerline
    }
    clcoords = PATSPEC::GetCenterline(filename);
  }

  // find out whether we have a orthopressure or FSI condition
  std::vector<DRT::Condition*> conds;
  std::vector<DRT::Condition*> ortho(0);
  std::vector<DRT::Condition*> fsi(0);
  dis->GetCondition("SurfaceNeumann", ortho);
  for (int i = 0; i < (int)ortho.size(); ++i)
  {
    if (ortho[i]->GType() != DRT::Condition::Surface) continue;
    const std::string* type = ortho[i]->Get<std::string>("type");
    if (*type == "neum_orthopressure" || *type == "neum_pseudo_orthopressure")
      conds.push_back(ortho[i]);
  }

  dis->GetCondition("FSICoupling", fsi);
  for (int i = 0; i < (int)fsi.size(); ++i)
  {
    if (fsi[i]->GType() != DRT::Condition::Surface) continue;
    conds.push_back(fsi[i]);
  }

  if (!conds.size()) dserror("There is no orthopressure nor FSI condition in this discretization");

  // collect all nodes
  std::set<int> allnodes;
  for (int i = 0; i < (int)conds.size(); ++i)
  {
    const std::vector<int>* nodes = conds[i]->Nodes();
    if (!nodes) dserror("Cannot find node ids in condition");
    for (int j = 0; j < (int)nodes->size(); ++j) allnodes.insert((*nodes)[j]);
  }

  // compute current coordinates and inner radius for all these nodes
  Teuchos::RCP<const Epetra_Vector> disp = dis->GetState("displacement");
  if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
  // compute inner radius for all these nodes
  const int nnodes = (int)allnodes.size();
  std::vector<double> linnerradius(nnodes, 0.0);
  std::vector<double> ginnerradius(nnodes, 0.0);
  std::set<int>::iterator fool;
  int count = 0;
  for (fool = allnodes.begin(); fool != allnodes.end(); ++fool)
  {
    if (!dis->NodeRowMap()->MyGID(*fool))
    {
      count++;
      continue;
    }
    DRT::Node* node = dis->gNode(*fool);
    std::vector<double> mydisp(1);
    std::vector<int> lm;
    dis->Dof(node, lm);
    DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
    const double* x = node->X();
    std::vector<double> xcurr(3);
    xcurr[0] = x[0] + mydisp[0];
    xcurr[1] = x[1] + mydisp[1];
    xcurr[2] = x[2] + mydisp[2];
    // inner radius
    double mindist = 1.0e12;
    if (normdirection == 0.0)
    {
      // centerline from file
      for (int j = 0; j < (int)(clcoords.size() / 3); ++j)
      {
        double* xcl = &clcoords[j * 3];
        double dist = sqrt((xcurr[0] - xcl[0]) * (xcurr[0] - xcl[0]) +
                           (xcurr[1] - xcl[1]) * (xcurr[1] - xcl[1]) +
                           (xcurr[2] - xcl[2]) * (xcurr[2] - xcl[2]));
        if (dist < mindist) mindist = dist;
      }
    }
    else
    {
      // centerline is a straight line
      std::vector<double> diff(3);
      diff[0] = xcurr[0] - support[0];
      diff[1] = xcurr[1] - support[1];
      diff[2] = xcurr[2] - support[2];
      std::vector<double> product(3);
      product[0] = direction[1] * diff[2] - direction[2] * diff[1];
      product[1] = direction[2] * diff[0] - direction[0] * diff[2];
      product[2] = direction[0] * diff[1] - direction[1] * diff[0];
      mindist = sqrt(product[0] * product[0] + product[1] * product[1] + product[2] * product[2]) /
                normdirection;
    }

    linnerradius[count] = mindist;
    count++;
  }
  dis->Comm().SumAll(&linnerradius[0], &ginnerradius[0], nnodes);

  // add inner radius to condition
  for (unsigned int j = 0; j < mypatspeccond.size(); ++j)
  {
    const Epetra_Vector* foolnode = mypatspeccond[j]->Get<Epetra_Vector>("inner surface node id");
    if (foolnode)
    {
      Teuchos::RCP<Epetra_Vector> innerradius = LINALG::CreateVector(*(dis->NodeRowMap()), true);
      for (int i = 0; i < dis->NodeRowMap()->NumMyElements(); ++i)
      {
        if (!foolnode->Map().MyGID(dis->NodeRowMap()->GID(i))) dserror("I do not have this node");

        int innerid = (*foolnode)[foolnode->Map().LID(dis->NodeRowMap()->GID(i))];
        (*innerradius)[i] = ginnerradius[innerid];
      }
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis->NodeColMap()), true);
      LINALG::Export(*innerradius, *tmp);
      innerradius = tmp;

      // compute element-wise mean distance from nodal distance
      Teuchos::RCP<Epetra_Vector> innerradiusele =
          LINALG::CreateVector(*(dis->ElementRowMap()), true);
      for (int i = 0; i < dis->ElementRowMap()->NumMyElements(); ++i)
      {
        DRT::Element* actele = dis->gElement(dis->ElementRowMap()->GID(i));
        double mean = 0.0;
        for (int k = 0; k < actele->NumNode(); ++k)
          mean += (*innerradius)[innerradius->Map().LID(actele->Nodes()[k]->Id())];
        mean /= actele->NumNode();
        (*innerradiusele)[i] = mean;
      }

      tmp = LINALG::CreateVector(*(dis->ElementColMap()), true);
      LINALG::Export(*innerradiusele, *tmp);
      innerradiusele = tmp;
      // add this column vector of inner radius to the condition and store in
      // discretization
      mypatspeccond[j]->Add("inner radius", *innerradiusele);
    }
  }
  ginnerradius.clear();

  return;
}

/*----------------------------------------------------------------------*
 |                                                           maier 06/11|
 *----------------------------------------------------------------------*/
std::vector<double> PATSPEC::GetCenterline(std::string filename)
{
  std::ifstream file(filename.c_str());
  if (file.fail()) dserror("Error opening centerline file");
  std::string sLine;
  std::vector<double> clcoords;
  while (std::getline(file, sLine))
  {
    if (sLine.size() != 0)
    {
      std::stringstream ss(sLine);
      double myDouble = 0;
      int count = 0;
      while (ss >> myDouble && count < 3)
      {
        clcoords.push_back(myDouble);
        count++;
      }
    }
  }
  return clcoords;
}

/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::GetILTDistance(
    const int eleid, Teuchos::ParameterList& params, DRT::Discretization& dis)
{
  std::vector<DRT::Condition*> mypatspeccond;
  dis.GetCondition("PatientSpecificData", mypatspeccond);
  if (!mypatspeccond.size()) return;

  for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
  {
    const Epetra_Vector* fool = mypatspeccond[i]->Get<Epetra_Vector>("normalized ilt thickness");
    if (fool)
    {
      if (!fool->Map().MyGID(eleid)) dserror("I do not have this element");

      double meaniltthick = (*fool)[fool->Map().LID(eleid)];
      params.set("iltthick meanvalue", meaniltthick);

      DRT::Element* actele = dis.gElement(eleid);
      Teuchos::RCP<MAT::Material> actmat = actele->Material();

      if (actmat->MaterialType() == INPAR::MAT::m_elasthyper)
      {
        MAT::ElastHyper* hyper = static_cast<MAT::ElastHyper*>(actmat.get());
        hyper->SetupAAA(params, eleid);
      }

      return;
    }
  }

  // if ilt thickness not found in condition just return
  return;
}


/*----------------------------------------------------------------------*
 |                                                           maier 11/10|
 *----------------------------------------------------------------------*/
void PATSPEC::GetLocalRadius(
    const int eleid, Teuchos::ParameterList& params, DRT::Discretization& dis)
{
  std::vector<DRT::Condition*> mypatspeccond;
  dis.GetCondition("PatientSpecificData", mypatspeccond);
  if (!mypatspeccond.size()) return;

  for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
  {
    const Epetra_Vector* fool = mypatspeccond[i]->Get<Epetra_Vector>("local radius");
    if (fool)
    {
      if (!fool->Map().MyGID(eleid)) dserror("I do not have this element");

      double meanlocalrad = (*fool)[fool->Map().LID(eleid)];
      params.set("localrad meanvalue", meanlocalrad);
      return;
    }
  }

  // if local radius not found in condition just return
  return;
}

/*----------------------------------------------------------------------*
 |                                                           tinkl 11/13|
 *----------------------------------------------------------------------*/
void PATSPEC::GetInnerRadius(
    const int eleid, Teuchos::ParameterList& params, DRT::Discretization& dis)
{
  std::vector<DRT::Condition*> mypatspeccond;
  dis.GetCondition("PatientSpecificData", mypatspeccond);
  if (!mypatspeccond.size()) return;

  for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
  {
    const Epetra_Vector* fool = mypatspeccond[i]->Get<Epetra_Vector>("inner radius");
    if (fool)
    {
      if (!fool->Map().MyGID(eleid)) dserror("I do not have this element");

      double innerradius = (*fool)[fool->Map().LID(eleid)];
      params.set("inner radius", innerradius);

      return;
    }
  }

  // if inner radius not found in condition just return
  return;
}


/*----------------------------------------------------------------------*
 |                                                           kehl 03/12|
 *----------------------------------------------------------------------*/
void PATSPEC::PatspecOutput(Teuchos::RCP<IO::DiscretizationWriter> output_,
    Teuchos::RCP<DRT::Discretization> discret_, Teuchos::RCP<Teuchos::ParameterList> params)
{
  std::vector<DRT::Condition*> mypatspeccond;
  discret_->GetCondition("PatientSpecificData", mypatspeccond);
  IO::VectorType vt = IO::elementvector;

  // also check if parameter is initialised. For monte carlo we want to setup patspec only once and
  // do not want output hence the parameterlist is only initialized in the first run and then kept
  // teuchos:null
  if (mypatspeccond.size() && params != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> patspecstuff =
        LINALG::CreateVector(*(discret_->ElementRowMap()), true);
    for (unsigned int i = 0; i < mypatspeccond.size(); ++i)
    {
      const Epetra_Vector* actcond =
          mypatspeccond[i]->Get<Epetra_Vector>("normalized ilt thickness");
      if (actcond)
      {
        double maxiltthick = params->get<double>("max ilt thick");
        for (int j = 0; j < patspecstuff->MyLength(); ++j)
          (*patspecstuff)[j] = (*actcond)[actcond->Map().LID(discret_->ElementRowMap()->GID(j))];
        patspecstuff->Scale(maxiltthick);
        output_->WriteVector("thrombus_thickness", patspecstuff, vt);
      }

      actcond = mypatspeccond[i]->Get<Epetra_Vector>("local radius");
      if (actcond)
      {
        for (int j = 0; j < patspecstuff->MyLength(); ++j)
          (*patspecstuff)[j] = (*actcond)[actcond->Map().LID(discret_->ElementRowMap()->GID(j))];
        output_->WriteVector("local_radius", patspecstuff, vt);
      }

      actcond = mypatspeccond[i]->Get<Epetra_Vector>("elestrength");
      if (actcond)
      {
        for (int j = 0; j < patspecstuff->MyLength(); ++j)
          (*patspecstuff)[j] = (*actcond)[actcond->Map().LID(discret_->ElementRowMap()->GID(j))];
        output_->WriteVector("strength", patspecstuff, vt);
      }

      actcond = mypatspeccond[i]->Get<Epetra_Vector>("inner radius");
      if (actcond)
      {
        for (int j = 0; j < patspecstuff->MyLength(); ++j)
          (*patspecstuff)[j] = (*actcond)[actcond->Map().LID(discret_->ElementRowMap()->GID(j))];
        output_->WriteVector("inner_radius", patspecstuff, vt);
      }
    }
    Teuchos::RCP<Epetra_Vector> eleID = LINALG::CreateVector(*(discret_->ElementRowMap()), true);
    for (int i = 0; i < eleID->MyLength(); ++i)
      (*eleID)[i] = (discret_->ElementRowMap()->GID(i)) + 1;
    output_->WriteVector("eleID", eleID, vt);
  }
}
/*----------------------------------------------------------------------*
 |                                                          biehler 09/14|
 *----------------------------------------------------------------------*/
double PATSPEC::ComputeMaxILTThickness(Teuchos::RCP<DRT::Discretization> dis)
{
  // get the AAA surface conditions
  std::vector<DRT::Condition*> aaa_surface;
  dis->GetCondition("AAASurface", aaa_surface);
  // check whether length of condition is two first condition luminal ILT surface second interface
  // ilt/wall
  if (aaa_surface.size() == 0)
    dserror(
        "AAA Surface condition condition not set! It must contain two surfaces. One surface of "
        "luminal ILT, one interface ILT/wall");
  // check whether length of condition is two first condition luminal ILT surface second interface
  // ilt/wall
  if (aaa_surface.size() != 2)
    dserror(
        "AAA Surface condition must contain two surfaces. One surface of luminal ILT, one "
        "interface ILT/wall");
  const std::vector<int>* luminal_nodes = aaa_surface[0]->Nodes();
  if (!luminal_nodes) dserror("Cannot find node ids in condition");
  const std::vector<int>* wall_nodes = aaa_surface[1]->Nodes();
  if (!wall_nodes) dserror("Cannot find node ids in condition");

  // now do a nice brute force search for maximum distance
  const int n_lum_itl_nodes = (int)luminal_nodes->size();
  const int n_lum_wall_nodes = (int)wall_nodes->size();
  // get coordinates for all nodes on luminal surface of ILT
  // and communicate them to all procs
  std::vector<double> lcoords_lum_ilt(n_lum_itl_nodes * 3, 0.0);
  std::vector<double> gcoords_lum_ilt(n_lum_itl_nodes * 3, 0.0);
  std::set<int>::iterator fool_lum_ilt;
  int count_lum_ilt = 0;
  for (int i = 0; i < (int)luminal_nodes->size(); i++)
  {
    if (!dis->NodeRowMap()->MyGID(luminal_nodes->at(i)))
    {
      count_lum_ilt++;
      continue;
    }
    DRT::Node* node = dis->gNode(luminal_nodes->at(i));
    lcoords_lum_ilt[count_lum_ilt * 3] = node->X()[0];
    lcoords_lum_ilt[count_lum_ilt * 3 + 1] = node->X()[1];
    lcoords_lum_ilt[count_lum_ilt * 3 + 2] = node->X()[2];
    count_lum_ilt++;
  }
  dis->Comm().Barrier();
  dis->Comm().SumAll(&lcoords_lum_ilt[0], &gcoords_lum_ilt[0], n_lum_itl_nodes * 3);
  lcoords_lum_ilt.clear();
  // because we have all nodes of luminal side stored redundant on all procs,
  // we only need to loop over local nodes of the wall surface
  // create vector for nodal values of potential maximum ilt thickness
  // WARNING: This is a brute force expensive minimum distance search!
  const Epetra_Map* nrowmap = dis->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> iltthick = LINALG::CreateVector(*nrowmap, true);
  iltthick->PutScalar(-10.0e12);

  for (int i = 0; i < n_lum_wall_nodes; i++)
  {
    if (!dis->NodeRowMap()->MyGID(wall_nodes->at(i)))
    {
      continue;
    }
    const double* x = dis->gNode(wall_nodes->at(i))->X();
    // loop nodes from the condition and find minimum distance
    double mindist = 10e12;

    for (int j = 0; j < n_lum_itl_nodes; ++j)
    {
      double* xorth = &gcoords_lum_ilt[j * 3];
      double dist =
          sqrt((x[0] - xorth[0]) * (x[0] - xorth[0]) + (x[1] - xorth[1]) * (x[1] - xorth[1]) +
               (x[2] - xorth[2]) * (x[2] - xorth[2]));
      if (dist < mindist) mindist = dist;
    }
    (*iltthick)[dis->gNode(wall_nodes->at(i))->LID()] = mindist;
  }
  dis->Comm().Barrier();
  gcoords_lum_ilt.clear();
  double maxiltthick;
  iltthick->MaxValue(&maxiltthick);
  return maxiltthick;
}
