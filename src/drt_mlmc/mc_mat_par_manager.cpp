/*----------------------------------------------------------------------*/
/*!
\file matpar_manager.H
\brief manage material parameters during UQ at some point this should be merged with the
mat par manager in stat inv ana

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            089 - 28915276
</pre>

!*/

/*----------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------*/
/* constructor                                               jb 05/14   */
/*----------------------------------------------------------------------*/
STR::UQ::MCMatParManager::MCMatParManager(Teuchos::RCP<DRT::Discretization> discret):
discret_(discret),
//paramlayoutmap_(Teuchos::null),
//paramlayoutmapunique_(Teuchos::null),
numstochparams_(0),
params_(Teuchos::null)
{
  if (not discret_->Filled() || not discret_->HaveDofs())
      dserror("Discretisation is not complete or has no dofs!");

  // set up maps to link against materials, parameters and materials/parameters for UQ
  InitStochParaMap();

  params_ = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementColMap()),numstochparams_,true));


  // temp map to keep correspondence of parameter block position and eleids
    // used to build the mapextractor and the various maps to keep track of parameters and elements
     // to be filled with stochparaid to all GIDS of elementes having this parameter
    std::map< int,std::vector<int> > elemap;
    int nummyparams=0; // number of elememts with stochastic materials times numstochparameters of this material

    // loop all my elements
    for (int i=0; i<discret_->NumMyRowElements(); i++)
    {
      DRT::Element* actele;
      actele = discret_->lRowElement(i);
      // get material ID
      int elematid = actele->Material()->Parameter()->Id();

      // if this material is not material stochastic skip rest of for loop ()
      if (stochparamap_.find(elematid) == stochparamap_.end() )
        continue;

      // if this material is stochastic get vector of ids of these parameters
      std::vector<int> actparapos = stochparaid_.at( elematid );
      std::vector<int>::const_iterator it;
      for ( it=actparapos.begin(); it!=actparapos.end(); it++)
      {
        // add GID of element to map
        elemap[*it].push_back(actele->Id());
        nummyparams++;
      }
    }

//    // generate global ids plus build map paramsLIDtoeleGID_
//    // store <stochparameter id><some increasing integer which is unique across procs GID of elementwise parameter >
//    std::map<int, std::vector<int> > gids;
//    int count = 0;
//    for (int i=0; i<discret_->Comm().NumProc(); i++)
//    {
//      if (discret_->Comm().MyPID() == i)
//      {
//        // loop over all stochastic material parameters
//        for (int j=0; j<numstochparams_; j++)
//        {
//          // loop over all elements which have a material with this parameter
//          for (int k=0; k<(int)elemap[j].size(); k++)
//          {
//            gids[j].push_back(count);
//            // store element GIDs in a vector
//            std::vector<int> elements_with_stochpara_j=elemap[j];
//            // store all eleGID that have parameter j sequentially in this vector
//            // final layout will be like this [1 2 3 4 ... 1 2 3 4 ...]
//            // the length of this will be equal to the total number of elementwise constitutive parameters
//            paramsLIDtoeleGID_.push_back(elements_with_stochpara_j.at(k));
//            count++;
//          }
//        }
//      }
//      // broadcast count from proc i to all others
//      discret_->Comm().Broadcast(&count,1,i);
//    }
//
//    //build map eleGIDtoparamsLID_
//    for (int i=0; i<(int)paramsLIDtoeleGID_.size(); i++)
//    {
//      // the blocks are ordered so this just
//      int myElementGID=paramsLIDtoeleGID_[i];
//      // will hold <eleGID, paramsLID>
//      eleGIDtoparamsLID_[myElementGID].push_back(i);
//    }
//
//    // the full map of the vector layout (nummyparams contains number of locally store elementwise quantities)
//    paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(-1,nummyparams,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
//    paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));
//
//    // the partial maps:
//    std::vector< Teuchos::RCP<const Epetra_Map> > partials;
//    for (int i=0; i<numstochparams_; i++)
//    {
//      partials.push_back(Teuchos::rcp(new Epetra_Map(-1,gids[i].size(),gids[i].data(),0,discret_->Comm())));
//    }
//
//    // finally build the MapExtractor
//    paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmap_,partials));


}


/*----------------------------------------------------------------------*/
/* Setup map of material parameters that are stochastic       jb 05/14  */
/*----------------------------------------------------------------------*/
void STR::UQ::MCMatParManager::InitStochParaMap()
{
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // the materials of the problem
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  if (discret_->Comm().MyPID()==0)
  {
    IO::cout << "STR::INVANA::MatParManager ... SETUP" << IO::endl;
    IO::cout <<  "Optimizing material with ids: ";
  }

  // parameters to be optimized
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(mlmcp,"PARAMLIST"));
  int matid;
  int actmatid=0;
  char* pEnd;
  while (pstream >> word2)
  {
    matid = std::strtol(&word2[0],&pEnd,10);
    if (*pEnd=='\0') //if (matid != 0)
    {
      if (discret_->Comm().MyPID()==0) std::cout << matid << " ";
      actmatid = matid;
      continue;
    }

    if (word2!="none" && actmatid!=0)
    {
      //check whether this material exists in the problem
      if ( mats.find(actmatid) == mats.end() )
        dserror("material %d not found in matset", actmatid);

      //check if this material has parameters that are stochastic (This currently works only for AAAneohook material
      // and can be checked using Optparams() )
      std::map<std::string, int> optparams;
      mats.at(actmatid)->Parameter()->OptParams(&optparams);
      if ( optparams.find(word2) == optparams.end() )
        dserror("parameter %s is not prepared to be optimized for mat %s", word2.c_str(), mats.at(actmatid)->Name().c_str());

      stochparamap_[actmatid].push_back(optparams.at(word2));
      stochparaid_[actmatid].push_back(numstochparams_);
      numstochparams_ += 1;
      IO::cout << "numstochparams " << numstochparams_ << IO::endl;
    }
    else
      dserror("Give the parameters for the respective materials");
  }

  if (discret_->Comm().MyPID()==0)
  {
    IO::cout << "" << IO::endl;
    IO::cout << "the number of stochastic material parameters is: " << numstochparams_ << IO::endl;
  }

}

void STR::UQ::MCMatParManager::ComputeMatParamsMultivectorFromRandomFields(Teuchos::RCP<Epetra_MultiVector> params,double para_cont_parameter)
{
  // must have EleRowMap layout
  params->PutScalar(0.0);
    for (int i=0; i< (discret_->NumMyRowElements()); i++)
    {
      // check whether Element has this material
      DRT::Element* actele;
      actele = discret_->lRowElement(i);
      int elematid = actele->Material()->Parameter()->Id();
      // if this material is not material stochastic skip rest of for loop ()
      if (stochparamap_.find(elematid) == stochparamap_.end() )
        continue;

      std::vector<double> ele_center;
      ele_center = actele->ElementCenterRefeCoords();
      std::vector<int> actparapos = stochparaid_.at( elematid );

      std::vector<int>::const_iterator it;
      for ( it=actparapos.begin(); it!=actparapos.end(); it++)
      {
        std::vector<double> ele_center_temp;
        // get dim of field
        if(randomfields_[*it]->Dimension()==2)
        {
          //special hack here assuming circular geometry with r=25 mm
          double phi= acos(ele_center[0]/25);
          //compute x coord
          ele_center_temp.push_back(phi*25);
          ele_center_temp.push_back(ele_center[2]);
          ele_center_temp.push_back(ele_center[2]);
          // compute random field values an store in params vector
          (*params)[*it][i]=randomfields_[*it]->EvalFieldAtLocation(ele_center_temp,para_cont_parameter,false,false);
        }
        else
        {
          // compute random field values an store in params vector
          (*params)[*it][i]=randomfields_[*it]->EvalFieldAtLocation(ele_center,para_cont_parameter,false,false);
        }
      }// eof loop over stoch material parameters

    } // eof loop over RowElements
}

/*----------------------------------------------------------------------*/
/* Compute params vector and set Vectors of respective materials jb 05/14  */
/*----------------------------------------------------------------------*/
void STR::UQ::MCMatParManager::SetParams(double para_cont_parameter)
{
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  // get the actual set of elementwise material parameters from the derived classes
  Teuchos::RCP<Epetra_MultiVector> getparams = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementRowMap()),numstochparams_,false));
  ComputeMatParamsMultivectorFromRandomFields(getparams, para_cont_parameter);
  discret_->Comm().Barrier();
  discret_->Comm().Barrier();

  // export to column layout to be able to run column elements
  LINALG::Export(*getparams,*params_);

  //loop materials to be optimized
  std::map<int,std::vector<int> >::const_iterator curr;
  for (curr=stochparamap_.begin(); curr != stochparamap_.end(); curr++ )
  {
    Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(curr->first);

    // loop the parameters to be optimized
    std::vector<int> actparams = stochparamap_.at(curr->first);
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      actmat->Parameter()->SetParameter(*it,Teuchos::rcp((*params_)( stochparaid_.at(curr->first).at(it-actparams.begin()) ),false));
    }
  }//loop optimized materials
}


/*----------------------------------------------------------------------*
 |  Setup Random fields and store them in map                  jb 05/14  |
 *----------------------------------------------------------------------*/
void STR::UQ::MCMatParManager::SetupRandomFields(unsigned int myseed)
{
  // loop over num stoch parameter
  for(int i=0;i<numstochparams_; i++)
  {
    randomfields_[i]=CreateRandomField(i+1, myseed);
  }
}

/*----------------------------------------------------------------------*
 |  Write random variables to file                            jb 05/14  |
 *----------------------------------------------------------------------*/
void STR::UQ::MCMatParManager::WriteRandomVariablesToFile(std::string filename, int numrun)
{
  // loop over num stoch parameter
  for(int i=0;i<numstochparams_; i++)
  {
    std::stringstream ss; //create a stringstream
    ss << filename;
    ss<< "_";
    ss << "rf_no_" ;
    ss << i+1;
    ss << "_numbrun_";
    ss <<  numrun;//add number to the stream
    ss << ".txt";
    randomfields_[i]->WriteRandomVariablesToFile(ss.str());
  }
}


/*----------------------------------------------------------------------*
 |  Compute new realizations of random fields                 jb 05/14  |
 *----------------------------------------------------------------------*/
void STR::UQ::MCMatParManager::CreateNewRealizationOfRandomFields(unsigned int myseed)
{
  // loop over num stoch parameter
  for(int i=0;i<numstochparams_; i++)
  {
    randomfields_[i]->CreateNewSample(myseed+(i*51200));
  }
}

/*----------------------------------------------------------------------*
 |  Compute new realizations of random fields                 jb 05/14  |
 *----------------------------------------------------------------------*/
void STR::UQ::MCMatParManager::SetUpStochMats(unsigned int myseed, double para_cont_parameter, bool reuse_rf)
{
  if(!reuse_rf)
    CreateNewRealizationOfRandomFields(myseed);
  SetParams(para_cont_parameter);
}


/*----------------------------------------------------------------------*
 |  Create Random field based on input data                   jb 05/14  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<STR::UQ::RandomField> STR::UQ::MCMatParManager::CreateRandomField(int random_field_id, unsigned int myseed)
{
  const Teuchos::ParameterList& rfp = DRT::Problem::Instance()->RandomFieldParams(random_field_id);
  // before calling the construtor make a quick safety check whether this random field was activated in
  // the input file or if the section contains only the default parameters
  bool active = DRT::INPUT::IntegralValue<int>(rfp ,"ACTIVE");
  if (!active)
    dserror("Trying to setup random field that is not active");

  // call constructor based on type of random field
  INPAR::MLMC::CalcMethod calcm = DRT::INPUT::IntegralValue<INPAR::MLMC::CalcMethod>(rfp,"CALC_METHOD");
  Teuchos::RCP<RandomField> test = Teuchos::null;

  switch(calcm)
  {
    case INPAR::MLMC::calc_m_fft:
      test = Teuchos::rcp(new RandomFieldSpectral(myseed,discret_,rfp ));
      break;
    case INPAR::MLMC::calc_m_cos:
      test = Teuchos::rcp(new RandomFieldSpectral(myseed,discret_,rfp ));
      break;
    case INPAR::MLMC::calc_m_fourier:
      test = Teuchos::rcp(new RandomFieldFourier(myseed,discret_,rfp));
      break;
    default:
      dserror("Unknown simulation method for RF choose fft or cos or fourier");
      break;
  }
  return test;
}
