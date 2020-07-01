/*----------------------------------------------------------------------*/
/*! \file
\brief manage material parameters during optimization

\level 3


!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/growth.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include "../drt_mat/elasthyper.H"  //to fit elasthyper-materials
#include "../drt_mat/growth_law.H"

/*----------------------------------------------------------------------*/
INVANA::MatParManager::MatParManager(Teuchos::RCP<DRT::Discretization> discret)
    : optparams_(Teuchos::null),
      optparams_initial_(Teuchos::null),
      paramlayoutmap_(Teuchos::null),
      paramlayoutmapunique_(Teuchos::null),
      paramapextractor_(Teuchos::null),
      metaparams_(),
      parapos_(),
      paramap_(),
      paraposGIDtoLID_(),
      numparams_(0),
      params_(Teuchos::null),
      discret_(discret)
{
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::Init(
    const Teuchos::ParameterList& invp, Teuchos::RCP<INVANA::ObjectiveFunct> objfunct)
{
  inpar_invana_ = invp;
  objfunct_ = objfunct;

  // want metaparametrization
  metaparams_.SetParametrization(
      DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMetaParamsType>(inpar_invana_, "METAPARAMS"));

  // set up maps to link materials their respective parameters
  // and the set of final optimization parameters
  SetupMatOptMap(inpar_invana_);

  // parameter initialization
  params_ = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementColMap()), numparams_, true));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::InitParams()
{
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();

  std::map<int, std::vector<int>>::const_iterator it;
  for (it = paramap_.begin(); it != paramap_.end(); it++)
  {
    Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(it->first);

    switch (actmat->Parameter()->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      case INPAR::MAT::m_scatra:
      case INPAR::MAT::m_growth_const:
        // case INPAR::MAT::m_growth_linsimple:
        {
          std::vector<int>::const_iterator jt;
          for (jt = it->second.begin(); jt != it->second.end(); jt++)
          {
            double val = metaparams_.Material2Meta(actmat->Parameter()->GetParameter(*jt, 0));
            InitParameters(parapos_.at(it->first).at(jt - it->second.begin()), val);
          }
          break;
        }
      // elasthyper materials
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i = 0; i < nummat; ++i)
        {
          const int id = (*matids)[i];
          const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;

          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_coupneohooke:
            case INPAR::MAT::mes_coup1pow:
              // add additional summands here
              {
                std::vector<int>::const_iterator jt;
                for (jt = it->second.begin(); jt != it->second.end(); jt++)
                {
                  double val =
                      metaparams_.Material2Meta(actelastmat->Parameter()->GetParameter(*jt, 0));
                  InitParameters(parapos_.at(it->first).at(jt - it->second.begin()), val);
                }
                break;
              }
            default:
              dserror("Elasthyper-Material not provided by the Material Manager for Optimization");
          }
        }
        break;
      }  // end elasthyper
      default:
      {
        std::cout << "This material " << *actmat << std::endl;
        dserror(".. is not provided by the Material Manager for Optimization");
      }
      break;
    }
  }

  // keep the inital set of optimization parameters
  optparams_initial_->Scale(1.0, *optparams_);
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::SetupMatOptMap(const Teuchos::ParameterList& invp)
{
  // the materials of the problem
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();

  // parameters to be optimized
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(invp, "PARAMLIST"));
  int matid;
  int actmatid = 0;
  char* pEnd;
  while (pstream >> word2)
  {
    matid = std::strtol(&word2[0], &pEnd, 10);
    if (*pEnd == '\0')  // if (matid != 0)
    {
      actmatid = matid;
      continue;
    }

    int localcount = 0;
    if (word2 != "none" && actmatid != 0)
    {
      // check whether this material exists in the problem
      if (mats.find(actmatid) == mats.end()) dserror("material %d not found in matset", actmatid);

      // check if this material has parameters to be optimized:
      std::map<std::string, int> optparams;
      mats.at(actmatid)->Parameter()->OptParams(&optparams);
      if (optparams.find(word2) == optparams.end())
        dserror("parameter %s is not prepared to be optimized for mat %s", word2.c_str(),
            mats.at(actmatid)->Name().c_str());

      paramap_[actmatid].push_back(optparams.at(word2));
      parapos_[actmatid].push_back(numparams_);
      paraposGIDtoLID_[numparams_] = localcount;
      localcount += 1;
      numparams_ += 1;
    }
    else
      dserror("Give the parameters for the respective materials");
  }

  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "MatParManager Initialization:" << std::endl;
    std::cout << "  Optimizing materials: ";
    for (std::map<int, std::vector<int>>::iterator it = paramap_.begin(); it != paramap_.end();
         it++)
      std::cout << it->first << " ";
    std::cout << "" << std::endl;
    std::cout << "  Total number of material parameters: " << numparams_ << std::endl;
  }
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::SetParams()
{
  // get the actual set of elementwise material parameters from the derived classes
  Teuchos::RCP<Epetra_MultiVector> getparams =
      Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementRowMap()), numparams_, false));
  FillParameters(getparams);

  // export to column layout to be able to run column elements
  LINALG::Export(*getparams, *params_);

  // set parameters to the elements
  PushParamsToElements();
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::PushParamsToElements()
{
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();
  Teuchos::RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*params_));

  metaparams_.Meta2Material(params_, tmp);

  // loop materials to be optimized
  std::map<int, std::vector<int>>::const_iterator curr;
  for (curr = paramap_.begin(); curr != paramap_.end(); curr++)
  {
    Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(curr->first);

    // in case of elasthyper-materials --> actmat is summand
    if (actmat->Parameter()->Type() == INPAR::MAT::m_elasthyper)
    {
      MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
      if (!params) dserror("Cannot cast material parameters");
      const int nummat = params->nummat_;
      const std::vector<int>* matids = params->matids_;
      for (int i = 0; i < nummat; ++i)
      {
        const int id = (*matids)[i];
        const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;

        // loop the parameters to be optimized
        std::vector<int> actparams = paramap_.at(curr->first);
        std::vector<int>::const_iterator it;
        for (it = actparams.begin(); it != actparams.end(); it++)
        {
          actelastmat->Parameter()->SetParameter(*it,
              Teuchos::rcp((*tmp)(parapos_.at(curr->first).at(it - actparams.begin())), false));
        }
      }
    }
    // for all other materials
    else
    {
      // loop the parameters to be optimized
      std::vector<int> actparams = paramap_.at(curr->first);
      std::vector<int>::const_iterator it;
      for (it = actparams.begin(); it != actparams.end(); it++)
      {
        actmat->Parameter()->SetParameter(
            *it, Teuchos::rcp((*tmp)(parapos_.at(curr->first).at(it - actparams.begin())), false));
      }
    }  // loop optimized materials
  }
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> INVANA::MatParManager::GetRawParams()
{
  Teuchos::RCP<Epetra_MultiVector> getparams =
      Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementRowMap()), numparams_, false));

  // get the actual set of elementwise parameters
  // from the optimization parameter layout
  FillParameters(getparams);

  return getparams;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> INVANA::MatParManager::GetMatParams()
{
  // get parameters mapped to the element row layout
  Teuchos::RCP<Epetra_MultiVector> getparams = GetRawParams();

  // export to column layout to be able to run column elements
  LINALG::Export(*getparams, *params_);

  Teuchos::RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*params_));
  metaparams_.Meta2Material(params_, tmp);

  return tmp;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> INVANA::MatParManager::GetMatrixDiagonal(DcsMatrix& matrix)
{
  Teuchos::RCP<Epetra_MultiVector> diagonals =
      Teuchos::rcp(new Epetra_MultiVector(*Discret()->ElementRowMap(), numparams_, true));

  ApplyParametrization(matrix, diagonals);

  return diagonals;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::ReplaceParams(const Epetra_MultiVector& toreplace)
{
  optparams_->Update(1.0, toreplace, 0.0);

  // bring updated parameters to the elements
  SetParams();
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::ReplaceParamsShallow(const Epetra_MultiVector& toreplace)
{
  optparams_->Update(1.0, toreplace, 0.0);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  Teuchos::RCP<const Epetra_Vector> disdual = discret_->GetState("dual displacement");

  // get the current elementwise material parameters
  SetParams();

  // use an intermediate overlapping vector
  // to assemble the gradient
  Teuchos::RCP<Epetra_MultiVector> dfint_tmp =
      Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_, 1, true));

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  double dt = sdyn.get<double>("TIMESTEP");

  // the reason not to do this loop via a discretizations evaluate call is that if done as is,
  // all elements have to be looped only once and evaluation is done only in case when an
  // element really has materials with parameters to be optimized. And the chain-rule application
  // with respect to the parameters to be optimized can be done without setting up the whole
  // gradient dR/dp_m and postmultiply it with dp_m\dp_o
  // with R: Residual forces; p_m: material params; p_o parametrization of p_m
  // TODO: Think of a more standard baci design way!
  for (int i = 0; i < discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    int elematid = ElementOptMat(actele);

    if (elematid == -1) continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set("total time", time);
    p.set("delta time", dt);

    std::vector<int> actparams = paramap_.at(elematid);
    std::vector<int>::const_iterator it;
    for (it = actparams.begin(); it != actparams.end(); it++)
    {
      p.set("action", "calc_struct_nlnstiff");
      p.set("matparderiv", *it);

      // initialize element vectors
      DRT::Element::LocationArray la(discret_->NumDofSets());
      actele->LocationVector(*discret_, la, false);
      int ndof = la[0].lm_.size();
      Epetra_SerialDenseMatrix elematrix1(ndof, ndof, false);
      Epetra_SerialDenseMatrix elematrix2(ndof, ndof, false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      actele->Evaluate(
          p, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);

      // dont forget product rule in case of parametrized material parameters!
      double metaval =
          (*(*params_)(parapos_.at(elematid).at(it - actparams.begin())))[actele->LID()];
      double val1 = metaparams_.DMaterialDMeta(metaval);
      elevector1.Scale(val1);

      // dualstate^T*(dR/dp_m)
      for (int l = 0; l < ndof; l++)
      {
        int lid = disdual->Map().LID(la[0].lm_.at(l));
        if (lid == -1) dserror("not found on this processor");
        elevector2[l] = (*disdual)[lid];
      }

      // functional differentiated wrt to this elements material parameters
      double val2 = elevector2.Dot(elevector1);

      // Assemble the final gradient; this is parametrization class business
      // (i.e contraction to (optimization)-parameter space:
      ContractGradient(dfint_tmp, val2, actele->Id(),
          parapos_.at(elematid).at(it - actparams.begin()), it - actparams.begin());

    }  // loop this elements material parameters (only the ones to be optimized)

  }  // loop elements

  // Finalize assembly, i.e. account for off-processor components
  Finalize(dfint_tmp, dfint);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManager::AddEvaluateFD(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  if (Comm().NumProc() > 1) dserror("this does probably not run in parallel");

  Teuchos::RCP<const Epetra_Vector> disdual = discret_->GetState("dual displacement");

  // get the actual set of elementwise material parameters from the derived classes
  Teuchos::RCP<Epetra_MultiVector> getparams =
      Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementRowMap()), numparams_, false));
  FillParameters(getparams);

  // export to column layout to be able to run column elements
  Comm().Barrier();
  LINALG::Export(*getparams, *params_);

  // a backup copy
  Epetra_MultiVector paramsbak(*params_);

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  double dt = sdyn.get<double>("TIMESTEP");

  for (int i = 0; i < discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    int elematid = ElementOptMat(actele);

    if (elematid == -1) continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set("total time", time);
    p.set("delta time", dt);

    std::vector<int> actparams = paramap_.at(elematid);
    std::vector<int>::const_iterator it;
    for (it = actparams.begin(); it != actparams.end(); it++)
    {
      p.set("action", "calc_struct_nlnstiff");

      double pa = 1.0e-6;
      double pb = 1.0e-12;

      // initialize element vectors
      DRT::Element::LocationArray la(discret_->NumDofSets());
      actele->LocationVector(*discret_, la, false);
      int ndof = la[0].lm_.size();
      Epetra_SerialDenseMatrix elematrix1(ndof, ndof, false);
      Epetra_SerialDenseMatrix elematrix2(ndof, ndof, false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);


      double actp = (*(*params_)(parapos_.at(elematid).at(it - actparams.begin())))[actele->LID()];
      double perturb = actp + pb + actp * pa;

      for (int j = 0; j < 2; j++)
      {
        if (j == 1)
        {
          params_->ReplaceMyValue(
              actele->LID(), parapos_.at(elematid).at(it - actparams.begin()), perturb);
          PushParamsToElements();
        }

        if (j == 0)
        {
          actele->Evaluate(
              p, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector2);
          // elevector1.Print(std::cout);
        }
        else if (j == 1)
        {
          actele->Evaluate(
              p, *discret_, la, elematrix1, elematrix2, elevector3, elevector2, elevector2);
          // elevector3.Print(std::cout);
        }

        // reset params
        if (j == 1)
        {
          params_->Update(1.0, paramsbak, 0.0);
          PushParamsToElements();
        }
      }

      // build fd approx
      elevector1.Scale(-1.0);
      elevector1 += elevector3;
      elevector1.Scale(1.0 / (pb + actp * pa));

      // reuse elevector2
      for (int l = 0; l < (int)la[0].lm_.size(); l++)
      {
        int lid = disdual->Map().LID(la[0].lm_.at(l));
        if (lid == -1) dserror("not found on this processor");
        elevector2[l] = (*disdual)[lid];
      }
      double val2 = elevector2.Dot(elevector1);

      // Assemble the final gradient; this is parametrization class business
      // (i.e contraction to (optimization)-parameter space:
      ContractGradient(dfint, val2, actele->Id(), parapos_.at(elematid).at(it - actparams.begin()),
          it - actparams.begin());

    }  // loop this elements material parameters (only the ones to be optimized)

  }  // loop elements
}

/*----------------------------------------------------------------------*/
int INVANA::MatParManager::GetParameterLocation(int eleid, std::string name)
{
  int loc = -1;

  if (!discret_->HaveGlobalElement(eleid))
    dserror("provide only ids of elements on this processor");

  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();

  DRT::Element* actele = discret_->gElement(eleid);
  int matid = actele->Material()->Parameter()->Id();

  std::map<std::string, int> optparams;
  mats.at(matid)->Parameter()->OptParams(&optparams);
  if (optparams.find(name) == optparams.end())
    dserror("parameter %s is not prepared to be optimized for mat %s", name.c_str(),
        mats.at(matid)->Name().c_str());

  if (paramap_.find(matid) == paramap_.end())
    dserror("Material with matid %d is not given for optimization in datfile", matid);
  else
  {
    // this is the vector index
    std::vector<int> actparams = paramap_.at(matid);
    std::vector<int>::const_iterator it;
    for (it = actparams.begin(); it != actparams.end(); it++)
    {
      if (actparams.at(it - actparams.begin()) == optparams.at(name))
        loc = parapos_.at(matid).at(it - actparams.begin());
    }
  }

  return loc;
}

/*----------------------------------------------------------------------*/
int INVANA::MatParManager::ElementOptMat(DRT::Element* ele)
{
  int elematid = ele->Material()->Parameter()->Id();

  if (paramap_.find(elematid) == paramap_.end())
  {
    MAT::PAR::Growth* mgrowth = dynamic_cast<MAT::PAR::Growth*>(ele->Material()->Parameter());
    if (mgrowth == NULL) return -1;

    int mgrowthid = mgrowth->growthlaw_->Parameter()->Id();
    if (paramap_.find(mgrowthid) == paramap_.end())
      return -1;
    else
      elematid = mgrowthid;
  }

  return elematid;
}

/*----------------------------------------------------------------------*/
const Epetra_Comm& INVANA::MatParManager::Comm() { return discret_->Comm(); }

/*----------------------------------------------------------------------*/
const Epetra_Comm& INVANA::MatParManager::GComm()
{
  Epetra_Comm& gcomm = *(DRT::Problem::Instance()->GetNPGroup()->GlobalComm());
  return gcomm;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<INVANA::ConnectivityData> INVANA::MatParManager::GetConnectivityData()
{
  Teuchos::RCP<Epetra_CrsMatrix> graph = Graph();
  Teuchos::RCP<Epetra_CrsMatrix> projector = Projector();

  // store graph in a container
  Teuchos::RCP<ConnectivityData> connectivity =
      Teuchos::rcp(new ConnectivityData(graph, projector));

  return connectivity;
}

/*----------------------------------------------------------------------*/
double INVANA::MetaParametrization::Material2Meta(double matval)
{
  double val = 0.0;

  switch (metatype_)
  {
    case INPAR::INVANA::stat_inv_meta_none:
    {
      val = matval;
      break;
    }
    case INPAR::INVANA::stat_inv_meta_quad:
    {
      val = sqrt(2 * (matval - 0.1));
      break;
    }
    case INPAR::INVANA::stat_inv_meta_exp:
    {
      val = log(matval);
      break;
    }
    case INPAR::INVANA::stat_inv_meta_arctan:
    {
      val = tan(PI * (matval - 0.5));
      break;
    }
    default:
      dserror("metaparams only implemented for none/quad/arctan");
      break;
  }

  return val;
}

/*----------------------------------------------------------------------*/
double INVANA::MetaParametrization::DMaterialDMeta(double metaval)
{
  double val = 0.0;

  switch (metatype_)
  {
    case INPAR::INVANA::stat_inv_meta_none:
    {
      val = 1.0;
      break;
    }
    case INPAR::INVANA::stat_inv_meta_quad:
    {
      val = metaval;
      break;
    }
    case INPAR::INVANA::stat_inv_meta_exp:
    {
      val = exp(metaval);
      break;
    }
    case INPAR::INVANA::stat_inv_meta_arctan:
    {
      val = metaval;
      val = 1. / PI / (val * val + 1.);
      break;
    }
    default:
      dserror("metaparams only implemented for none/quad/arctan");
      break;
  }

  return val;
}

/*----------------------------------------------------------------------*/
void INVANA::MetaParametrization::Meta2Material(
    Teuchos::RCP<Epetra_MultiVector> meta, Teuchos::RCP<Epetra_MultiVector> material)
{
  switch (metatype_)
  {
    case INPAR::INVANA::stat_inv_meta_none:
      material->Scale(1.0, *meta);
      break;
    case INPAR::INVANA::stat_inv_meta_quad:
    {
      material->PutScalar(0.1);
      material->Multiply(0.5, *meta, *meta, 1.0);
      break;
    }
    case INPAR::INVANA::stat_inv_meta_exp:
    {
      int numvecs = meta->NumVectors();

      double* metaval;
      double* matval;

      for (int i = 0; i < numvecs; i++)
      {
        (*meta)(i)->ExtractView(&metaval);
        (*material)(i)->ExtractView(&matval);

        for (int j = 0; j < material->MyLength(); j++) matval[j] = exp(metaval[j]);
      }
      break;
    }
    case INPAR::INVANA::stat_inv_meta_arctan:
    {
      material->PutScalar(0.5);
      for (int n = 0; n < meta->NumVectors(); ++n)
        for (int i = 0; i < meta->MyLength(); ++i)
        {
          (*(*material)(n))[i] += 1. / PI * atan((*(*material)(n))[i]);
        }
      break;
    }
    default:
      dserror("metaparams only implemented for none/quad/arctan");
      break;
  }

  return;
}
