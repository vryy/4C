/*----------------------------------------------------------------------------*/
/*!
 \file mc_mat_par_manager.cpp
\brief manages material parameters during UQ at some point this should
       be merged with the mat par manager in stat inv ana

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            089 - 28915276
</pre>

!*/
#ifdef HAVE_FFTW
/*----------------------------------------------------------------------------*/
/* headers */
#include "mc_mat_par_manager.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"
#include "../linalg/linalg_utils.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "randomfield_spectral.H"
#include "randomvariable.H"
#include "../drt_mat/maxwell_0d_acinus_Ogden.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
UQ::MCMatParManager::MCMatParManager(Teuchos::RCP<DRT::Discretization> discret)
    : discret_(discret),
      numstochparams_r_field_(0),
      numstochparams_r_var_(0),
      params_r_field_(Teuchos::null),
      params_r_field_median_(Teuchos::null),
      params_r_var_(Teuchos::null)
{
  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");

  // set up maps to link against materials, parameters and
  // materials/parameters for UQ
  InitStochParaMaps();

  if (numstochparams_r_field_)
    params_r_field_ = Teuchos::rcp(
        new Epetra_MultiVector(*(discret_->ElementColMap()), numstochparams_r_field_, true));

  if (numstochparams_r_field_)
    params_r_field_median_ = Teuchos::rcp(
        new Epetra_MultiVector(*(discret_->ElementColMap()), numstochparams_r_field_, true));

  // create map with one redundant entry on each proc
  std::vector<int> MyIds(1, 1);
  Teuchos::RCP<const Epetra_Map> my_map =
      Teuchos::rcp(new Epetra_Map(1, 1, &(MyIds[0]), 0, discret_->Comm()));

  if (numstochparams_r_var_)
    params_r_var_ = Teuchos::rcp(new Epetra_MultiVector(*(my_map), numstochparams_r_var_, true));

  // temp map to keep correspondence of parameter block position and eleids
  // used to build the mapextractor and the various maps to keep track of
  // parameters and elements to be filled with stochparaid to all GIDS of
  // elements having this parameter
  std::map<int, std::vector<int>> elemap;

  // number of elements with stochastic materials times
  // numstoch parameters of this material
  int nummyparams = 0;

  // loop all my elements
  for (int i = 0; i < discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    // get material ID
    int elematid = actele->Material()->Parameter()->Id();

    // if this material is not material stochastic skip rest of for loop ()
    if (stochparamap_r_field_.find(elematid) == stochparamap_r_field_.end()) continue;

    // if this material is stochastic get vector of ids of these parameters
    std::vector<int> actparapos = stochparaid_r_field_.at(elematid);
    std::vector<int>::const_iterator it;
    for (it = actparapos.begin(); it != actparapos.end(); it++)
    {
      // add GID of element to map
      elemap[*it].push_back(actele->Id());
      nummyparams++;
    }
  }
  // commpute median value once and for all
  ComputeMedianMatParamsMultivectorFromRandomFields();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::InitStochParaMaps()
{
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // the materials of the problem
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();

  if (discret_->Comm().MyPID() == 0)
  {
    IO::cout << "UQ::MatParManager ... SETUP" << IO::endl;
    IO::cout << "UQ for material with ids: ";
  }

  // uncertain parameters modeled as random field
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(mlmcp, "PARAMLIST_R_FIELD"));
  int matid;
  int actmatid = 0;
  char* pEnd;
  while (pstream >> word2)
  {
    matid = std::strtol(&word2[0], &pEnd, 10);
    if (*pEnd == '\0')  // if (matid != 0)
    {
      if (discret_->Comm().MyPID() == 0) std::cout << matid << " ";
      actmatid = matid;
      continue;
    }

    if (word2 != "none" && actmatid != 0)
    {
      // check whether this material exists in the problem
      if (mats.find(actmatid) == mats.end()) dserror("material %d not found in matset", actmatid);

      // check if this material has parameters that are stochastic
      // (This currently works only for AAAneohook material
      //  and can be checked using Optparams() )
      std::map<std::string, int> optparams;
      mats.at(actmatid)->Parameter()->OptParams(&optparams);
      if (optparams.find(word2) == optparams.end())
        dserror("parameter %s is not prepared to be optimized for mat %s", word2.c_str(),
            mats.at(actmatid)->Name().c_str());

      stochparamap_r_field_[actmatid].push_back(optparams.at(word2));
      stochparaid_r_field_[actmatid].push_back(numstochparams_r_field_);
      numstochparams_r_field_ += 1;
    }
    // else, it's ok maybe there are some stochastic parameters modeled as
    // random variables
  }

  // repeat the procedure for random variables
  // uncertain parameters modeled as random variables
  std::string word3;
  std::istringstream pstream2(Teuchos::getNumericStringParameter(mlmcp, "PARAMLIST_R_VAR"));
  actmatid = 0;
  while (pstream2 >> word3)
  {
    matid = std::strtol(&word3[0], &pEnd, 10);
    if (*pEnd == '\0')  // if (matid != 0)
    {
      if (discret_->Comm().MyPID() == 0) std::cout << matid << " ";
      actmatid = matid;
      continue;
    }

    if (word3 != "none" && actmatid != 0)
    {
      // check whether this material exists in the problem
      if (mats.find(actmatid) == mats.end()) dserror("material %d not found in matset", actmatid);

      // check if this material has parameters that are stochastic
      // (This currently works only for AAAneohook material
      //  and can be checked using Optparams() )
      std::map<std::string, int> optparams;
      mats.at(actmatid)->Parameter()->OptParams(&optparams);
      if (optparams.find(word3) == optparams.end())
        dserror("parameter %s is not prepared to be optimized for mat %s", word3.c_str(),
            mats.at(actmatid)->Name().c_str());

      stochparamap_r_var_[actmatid].push_back(optparams.at(word3));
      stochparaid_r_var_[actmatid].push_back(numstochparams_r_var_);
      numstochparams_r_var_ += 1;
    }
    // else, it's ok maybe there are some stochastic parameters modeled as
    // random fields
  }


  if (discret_->Comm().MyPID() == 0)
  {
    IO::cout << "" << IO::endl;
    IO::cout << "the number of stochastic material parameters modelled as random field is: "
             << numstochparams_r_field_ << IO::endl;
    IO::cout << "the number of stochastic material parameters modelled as random variable is: "
             << numstochparams_r_var_ << IO::endl;
  }
  if (numstochparams_r_field_ + numstochparams_r_var_ == 0)
    dserror("No uncertain quantities defined, fix your input file!");

  // safety check to prevent that the same material appears in both list
  // in the future we might check if the same parameter is modeled however
  // currently this is sufficient
  std::map<int, std::vector<int>>::const_iterator curr_r_var;
  std::map<int, std::vector<int>>::const_iterator curr_r_field;

  for (curr_r_field = stochparamap_r_field_.begin(); curr_r_field != stochparamap_r_field_.end();
       curr_r_field++)
  {
    if (stochparamap_r_var_.find(curr_r_field->first) != stochparamap_r_var_.end())
      dserror(
          "MAT %d must not contain parameters modeled as random field and random variables at the "
          "same time",
          curr_r_field->first);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::ComputeMatParamsMultivectorFromRandomFields(
    const double para_cont_parameter)
{
  Teuchos::RCP<Epetra_MultiVector> params = Teuchos::rcp(
      new Epetra_MultiVector(*(discret_->ElementRowMap()), numstochparams_r_field_, false));

  // must have EleRowMap layout
  params->PutScalar(0.0);
  for (int i = 0; i < (discret_->NumMyRowElements()); i++)
  {
    // check whether Element has this material
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();
    // if this material is not material stochastic skip rest of for loop ()
    if (stochparamap_r_field_.find(elematid) == stochparamap_r_field_.end()) continue;

    std::vector<double> ele_center;
    ele_center = DRT::UTILS::ElementCenterRefeCoords(actele);
    std::vector<int> actparapos = stochparaid_r_field_.at(elematid);

    std::vector<int>::const_iterator it;
    for (it = actparapos.begin(); it != actparapos.end(); it++)
    {
      std::vector<double> ele_center_temp;
      // get dim of field
      if (randomfields_[*it]->Dimension() == 2)
      {
        // special solution here assuming circular geometry with r=25 mm
        double phi = acos(ele_center[0] / 25);
        // compute x coord
        ele_center_temp.push_back(phi * 25);
        ele_center_temp.push_back(ele_center[2]);
        ele_center_temp.push_back(ele_center[2]);
        // compute random field values an store in params vector
        if (!randomfields_[*it]->HasSpatialMedian())
        {
          (*params)[*it][i] = randomfields_[*it]->EvalFieldAtLocation(
              ele_center_temp, para_cont_parameter, false, false);
        }
        else  // we should actuall never go in here because only spectral field has 2d
              // implementation
        {
          (*params)[*it][i] = randomfields_[*it]->EvalFieldAtLocation(ele_center_temp,
              (*params_r_field_median_)[*it][i], para_cont_parameter, false, false);
        }
      }
      else  // i.e. dim == 3
      {
        // compute random field values an store in params vector
        if (!randomfields_[*it]->HasSpatialMedian())
        {
          (*params)[*it][i] = randomfields_[*it]->EvalFieldAtLocation(
              ele_center, para_cont_parameter, false, false);
        }
        else
        {
          (*params)[*it][i] = randomfields_[*it]->EvalFieldAtLocation(
              ele_center, (*params_r_field_median_)[*it][i], para_cont_parameter, false, false);
        }
      }
    }  // eof loop over stoch material parameters
  }    // eof loop over RowElements

  discret_->Comm().Barrier();
  discret_->Comm().Barrier();

  // export to column layout to be able to run column elements
  LINALG::Export(*params, *params_r_field_);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::ComputeMatParamsMultivectorFromRandomVariables(
    const double para_cont_parameter)
{
  params_r_var_->PutScalar(0.0);

  // loop all materials in problem
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance(0)->Materials()->Map();
  std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
  for (curr = mats.begin(); curr != mats.end(); ++curr)
  {
    int matid = curr->second->Parameter()->Id();
    // if this material is not material stochastic skip rest of for loop
    if (stochparamap_r_var_.find(matid) == stochparamap_r_var_.end()) continue;

    std::vector<int> actparapos = stochparaid_r_var_.at(matid);

    std::vector<int>::const_iterator it;
    for (it = actparapos.begin(); it != actparapos.end(); it++)
    {
      // compute random field values an store in params vector
      (*params_r_var_)[*it][0] =
          randomvariables_[*it]->EvalVariable(para_cont_parameter, false, false);
    }  // eof loop over stoch material parameters
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::ComputeMedianMatParamsMultivectorFromRandomFields()
{
  // must have EleRowMap layout
  for (int i = 0; i < (discret_->NumMyRowElements()); i++)
  {
    // check whether Element has this material
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();
    // if this material is not material stochastic skip rest of for loop ()
    if (stochparamap_r_field_.find(elematid) == stochparamap_r_field_.end()) continue;

    std::vector<int> actparapos = stochparamap_r_field_.at(elematid);
    std::vector<int> actparapos2 = stochparaid_r_field_.at(elematid);

    // quick error check
    if (actparapos.size() != actparapos2.size()) dserror("Size mismatch cannot continue");

    // depending on forward problem, call either red_airways or structure specific function
    const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
    INPAR::MLMC::FWDProblem fwdprb =
        DRT::INPUT::IntegralValue<INPAR::MLMC::FWDProblem>(mlmcp, "FWDPROBLEM");

    if (fwdprb == INPAR::MLMC::structure)
    {
      std::vector<int>::const_iterator it;
      std::vector<int>::const_iterator it2;
      it2 = actparapos2.begin();
      for (it = actparapos.begin(); it != actparapos.end(); it++)
      {
        // IO::cout << "beta?" << actele->Material()->Parameter()->GetParameter(*it,actele->Id()) <<
        // IO::endl;
        (*params_r_field_median_)[*it2][i] =
            actele->Material()->Parameter()->GetParameter(*it, actele->Id());
        it2++;
      }  // eof loop over stoch material parameters
    }
    else if (fwdprb == INPAR::MLMC::red_airways)
    {
      switch (actele->Material()->MaterialType())
      {
        case INPAR::MAT::m_0d_maxwell_acinus_ogden:
        {
          Teuchos::RCP<MAT::Maxwell_0d_acinus_Ogden> mymat =
              Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus_Ogden>(actele->Material());

          // it iterator over number of RF per element, meaning stochastic parameters for element,
          // in this case now, only kappa is evaluated using RF, so only one RF per element.
          // loop the parameters to be optimized
          std::vector<int> actparams = stochparamap_r_field_.at(
              stochparamap_r_field_.find(actele->Material()->Parameter()->Id())->first);
          if (actparams.size() != 1)
            dserror("actparams.size() has to be one for red_airway problems");
          (*params_r_field_median_)[*actparams.begin()][i] = mymat->GetParams("kappa");
          break;
        }
        default:
        {
          dserror("Material not implemented yet for stochastic parameter setting");
          break;
        }

      }  // end switch
    }
    else
    {
      dserror("UQ is currently only implemented for structure and red_airways problems");
    }
  }  // eof loop over RowElements

  discret_->Comm().Barrier();
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::SetParamsStructure(double para_cont_parameter)
{
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance()->Materials()->Map();

  std::map<int, std::vector<int>>::const_iterator curr;

  // ************************************************************
  // deal with material parameters modeled as random field first
  // ************************************************************
  // do we have quantities modeled as random fields ?
  if (numstochparams_r_field_)
  {
    // this call updates params_r_field_
    ComputeMatParamsMultivectorFromRandomFields(para_cont_parameter);

    for (curr = stochparamap_r_field_.begin(); curr != stochparamap_r_field_.end(); curr++)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(curr->first);

      // loop the parameters to be optimized
      std::vector<int> actparams = stochparamap_r_field_.at(curr->first);
      std::vector<int>::const_iterator it;
      for (it = actparams.begin(); it != actparams.end(); it++)
      {
        actmat->Parameter()->SetParameter(*it,
            Teuchos::rcp(
                (*params_r_field_)(stochparaid_r_field_.at(curr->first).at(it - actparams.begin())),
                false));
      }
    }  // loop materials
  }

  // **************************************************************
  // now deal with material parameters modeled as random variables
  // **************************************************************
  // do we have quantities modeled as random variables ?
  if (numstochparams_r_var_)
  {
    ComputeMatParamsMultivectorFromRandomVariables(para_cont_parameter);
    for (curr = stochparamap_r_var_.begin(); curr != stochparamap_r_var_.end(); curr++)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(curr->first);
      // loop the parameters to be optimized
      std::vector<int> actparams = stochparamap_r_var_.at(curr->first);
      std::vector<int>::const_iterator it;
      for (it = actparams.begin(); it != actparams.end(); it++)
      {
        actmat->Parameter()->SetParameter(*it,
            Teuchos::rcp(
                (*params_r_var_)(stochparaid_r_var_.at(curr->first).at(it - actparams.begin())),
                false));
      }
    }  // loop materials
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::SetParamsRedAirways(double para_cont_parameter)
{
  // ************************************************************
  // deal with material parameters modeled as random field first
  // ************************************************************
  // do we have quantities modeled as random fields ?
  if (numstochparams_r_field_)
  {
    // updates params_r_field_
    ComputeMatParamsMultivectorFromRandomFields(para_cont_parameter);
    // loop over number of elements (on processor) -- the material class has element specific
    // parameters
    for (int i = 0; i < (discret_->NumMyColElements()); i++)
    {
      // get global element ID
      int eleGID = discret_->ElementColMap()->GID(i);
      // is this material a stochastic material, global material id needed
      if (stochparamap_r_field_.find(discret_->gElement(eleGID)->Material()->Parameter()->Id()) !=
          stochparamap_r_field_.end())
      {
        switch (discret_->gElement(eleGID)->Material()->MaterialType())
        {
          case INPAR::MAT::m_0d_maxwell_acinus_ogden:
          {
            Teuchos::RCP<MAT::Maxwell_0d_acinus_Ogden> mymat =
                Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus_Ogden>(
                    discret_->gElement(eleGID)->Material());

            // it iterator over number of RF per element, meaning stochastic parameters for element,
            // in this case now, only kappa is evaluated using RF, so only one RF per element.
            // loop the parameters to be optimized
            std::vector<int> actparams = stochparamap_r_field_.at(
                stochparamap_r_field_
                    .find(discret_->gElement(eleGID)->Material()->Parameter()->Id())
                    ->first);
            if (actparams.size() != 1)
              dserror("actparams.size() has to be one for red_airway problems");
            mymat->SetParams("kappa", (*params_r_field_)[*actparams.begin()][i]);
            break;
          }
          default:
          {
            dserror("Material not implemented yet for stochastic parameter setting");
            break;
          }
        }  // end switch
      }    // end if stochastic
    }      // end loop elements
  }        // end if

  // **************************************************************
  // now deal with material parameters modeled as random variables
  // **************************************************************
  // do we have quantities modeled as random variables ?
  if (numstochparams_r_var_)
  {
    dserror("Untested, check for correct implementation before use");
    ComputeMatParamsMultivectorFromRandomVariables(para_cont_parameter);
    // loop over number of elements (on processor) -- the material class has element specific
    // parameters
    for (int i = 0; i < (discret_->NumMyColElements()); i++)
    {
      int eleDID = discret_->ElementColMap()->GID(i);  // get global element ID
      int num_mymat =
          discret_->gElement(eleDID)->NumMaterial();  // get the number of element materials

      // loop over number of element materials
      for (int j = 0; j < num_mymat; j++)
      {
        // is this material a stochastic material, global material id needed
        if (stochparamap_r_field_.find(
                discret_->gElement(eleDID)->Material(j)->Parameter()->Id()) !=
            stochparamap_r_field_.end())
        {
          switch (discret_->gElement(eleDID)->Material(j)->MaterialType())
          {
            case INPAR::MAT::m_0d_maxwell_acinus_ogden:
            {
              Teuchos::RCP<MAT::Maxwell_0d_acinus_Ogden> mymat =
                  Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus_Ogden>(
                      discret_->gElement(eleDID)->Material(j));

              // it iterator over number of RF per element, meaning stochastic parameters for
              // element, in this case now, only kappa is evaluted using RF, so only one RF per
              // element. loop the parameters to be optimized
              std::vector<int> actparams = stochparamap_r_var_.at(
                  stochparamap_r_var_
                      .find(discret_->gElement(eleDID)->Material(j)->Parameter()->Id())
                      ->first);
              std::vector<int>::const_iterator it;
              for (it = actparams.begin(); it != actparams.end(); it++)
              {
                mymat->SetParams(0, (*params_r_var_)[*it][i]);
              }
              break;
            }
            default:
            {
              dserror("Material not implemented yet for stochastic parameter setting");
              break;
            }

          }  // end switch

        }  // end if stochastic
      }    // end loop element materials
    }      // end loop elements
  }        // end if
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::SetupAllRandomQuantities(unsigned int myseed)
{
  SetupRandomFields(myseed);
  SetupRandomVariables(myseed);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::SetupRandomFields(unsigned int myseed)
{
  // loop over num stoch parameter
  for (int i = 0; i < numstochparams_r_field_; i++)
  {
    randomfields_[i] = CreateRandomField(i + 1, myseed);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::SetupRandomVariables(unsigned int myseed)
{
  // loop over num stoch variables
  for (int i = 0; i < numstochparams_r_var_; i++)
  {
    randomvariables_[i] = CreateRandomVariable(i + 1, myseed);
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::WriteRandomVariablesToFile(std::string filename, int numrun)
{
  // loop over num stoch parameter
  for (int i = 0; i < numstochparams_r_field_; i++)
  {
    std::stringstream ss;  // create a stringstream
    ss << filename;
    ss << "_";
    ss << "rf_no_";
    ss << i + 1;
    ss << "_numbrun_";
    ss << numrun;  // add number to the stream
    ss << ".txt";
    randomfields_[i]->WriteRandomVariablesToFile(ss.str());
  }
  for (int i = 0; i < numstochparams_r_var_; i++)
  {
    std::stringstream ss;
    ss << filename << "_random_variable_" << i + 1 << ".txt";
    randomvariables_[i]->WriteRandomVariablesToFile(ss.str(), numrun);
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::CreateNewRealizationOfRandomQuantities(unsigned int myseed)
{
  // loop over num stoch parameter
  for (int i = 0; i < numstochparams_r_field_; i++)
  {
    // assumes that we never run more than 51200 samples
    randomfields_[i]->CreateNewSample(myseed + (i * 51200));
  }

  // loop over num stoch parameter
  for (int i = 0; i < numstochparams_r_var_; i++)
  {
    // assumes that we never run more than 51200 samples
    randomvariables_[i]->CreateNewSample(myseed + (i * 51200) + 23213);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::MCMatParManager::SetUpStochMats(
    unsigned int myseed, double para_cont_parameter, bool reuse_rf)
{
  // setup random fields
  if (!reuse_rf) CreateNewRealizationOfRandomQuantities(myseed);

  // depending on forward problem, call either red_airways or structure specific function
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  INPAR::MLMC::FWDProblem fwdprb =
      DRT::INPUT::IntegralValue<INPAR::MLMC::FWDProblem>(mlmcp, "FWDPROBLEM");

  switch (fwdprb)
  {
    case INPAR::MLMC::structure:
      SetParamsStructure(para_cont_parameter);
      break;
    case INPAR::MLMC::red_airways:
      SetParamsRedAirways(para_cont_parameter);
      break;
    default:
      dserror("Unknown forward problem type fix your input file");
      break;
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<UQ::RandomField> UQ::MCMatParManager::CreateRandomField(
    int random_field_id, unsigned int myseed)
{
  const Teuchos::ParameterList& rfp = DRT::Problem::Instance()->RandomFieldParams(random_field_id);
  // before calling the constructor make a quick safety check whether this
  // random field was activated in the input file or if the section contains
  // only the default parameters
  bool active = DRT::INPUT::IntegralValue<int>(rfp, "ACTIVE");
  if (!active) dserror("Trying to setup random field that is not active");

  // call constructor based on type of random field
  INPAR::MLMC::CalcMethod calcm =
      DRT::INPUT::IntegralValue<INPAR::MLMC::CalcMethod>(rfp, "CALC_METHOD");
  Teuchos::RCP<RandomField> test = Teuchos::null;

  switch (calcm)
  {
    case INPAR::MLMC::calc_m_fft:
      test = Teuchos::rcp(new RandomFieldSpectral(myseed, discret_, rfp));
      break;
    case INPAR::MLMC::calc_m_cos:
      test = Teuchos::rcp(new RandomFieldSpectral(myseed, discret_, rfp));
      break;
    case INPAR::MLMC::calc_m_fourier:
      test = Teuchos::rcp(new RandomFieldFourier(myseed, discret_, rfp));
      break;
    default:
      dserror("Unknown simulation method for RF choose fft or cos or fourier");
      break;
  }
  return test;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<UQ::RandomVariable> UQ::MCMatParManager::CreateRandomVariable(
    int random_variable_id, unsigned int myseed)
{
  const Teuchos::ParameterList& rfp =
      DRT::Problem::Instance()->RandomVariableParams(random_variable_id);
  // before calling the constructor make a quick safety check whether this
  // random variable was activated in the input file or if the section contains
  // only the default parameters

  bool active = DRT::INPUT::IntegralValue<int>(rfp, "ACTIVE");
  if (!active) dserror("Trying to setup random variable that is not active");

  // call constructor
  Teuchos::RCP<RandomVariable> temp =
      Teuchos::rcp(new RandomVariable(rfp, random_variable_id, myseed));
  return temp;
}

#endif /* #ifdef HAVE_FFTW */
