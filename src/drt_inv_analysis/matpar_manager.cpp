/*----------------------------------------------------------------------*/
/*!
\file matpar_manager.H
\brief manage material parameters during optimization

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"

/*----------------------------------------------------------------------*/
/* standard constructor */
STR::INVANA::MatParManager::MatParManager(Teuchos::RCP<DRT::Discretization> discret):
numparams_(0),
discret_(discret),
dofrowmap_(NULL),
elecolmap_(NULL),
params_(Teuchos::null),
params_o_(Teuchos::null)
{
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  if (not discret_->Filled() || not discret_->HaveDofs())
      dserror("Discretisation is not complete or has no dofs!");
  else
  {
    dofrowmap_ = discret_->DofRowMap();
    elecolmap_ = discret_->ElementColMap();
  }

  // set up maps to link against materials, parameters and materials/parameters for the optimization
  SetupMatOptMap();

  params_ = Teuchos::rcp(new Epetra_MultiVector(*elecolmap_,numparams_,true));
  params_o_ = Teuchos::rcp(new Epetra_MultiVector(*elecolmap_,numparams_,true));
  params_init_ = Teuchos::rcp(new Epetra_MultiVector(*elecolmap_,numparams_,true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();

  // set these to the elements
  SetParams();

  reg_weight_ = statinvp.get<double>("REG_WEIGHT");

}

void STR::INVANA::MatParManager::InitParams()
{
  const std::map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  /* ------------------------------*/
  /*    get initial parameters     */
  /* ------------------------------*/
  std::map<int,std::vector<std::string> >::const_iterator it;
  for (it=paramap_.begin(); it!=paramap_.end(); it++)
  {
    Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(it->first);
    switch(actmat->Parameter()->Type())
    {
      case INPAR::MAT::m_aaaneohooke_stopro:
      {
        std::vector<std::string>::const_iterator jt;
        for (jt = it->second.begin(); jt != it->second.end(); jt++)
        {
          (*params_)( parapos_.at(it->first).at(jt-it->second.begin()) )->PutScalar( actmat->GetDouble(*jt) );
          std::cout << *params_ << std::endl;
        }
      }
      break;
      default:
        dserror("Material not provided by the Material Manager for Optimization");
      break;
    }
  }

  params_init_->Update(1.0,*params_,0.0);

}

//--------------------------------------------------------------------------------------
// Setup map of parameters to be optimized
void STR::INVANA::MatParManager::SetupMatOptMap()
{
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  /* ------------------------------*/
  /*   parameters to be optimized  */
  /* ------------------------------*/
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(statinvp,"PARAMLIST"));
  int matid;
  int actmatid=0;
  char* pEnd;
  while (pstream >> word2)
  {
    matid = std::strtol(&word2[0],&pEnd,10);
    //std::cout << *pEnd << std::endl;
    if (*pEnd=='\0') //if (matid != 0)
    {
      std::cout <<  "matid " << matid << std::endl;
      actmatid = matid;
      std::cout << "this" << std::endl;
      continue;
    }

    std::cout << "testoutput: " << word2 << std::endl;
    if (word2!="none" && actmatid!=0)
    {
      paramap_[actmatid].push_back(word2);
      parapos_[actmatid].push_back(numparams_);
      numparams_ += 1;
    }
    else
      dserror("Give the parameters for the respective materials");
  }

  std::cout << "the number of parameters is: " << numparams_ << std::endl;

  /* ------------------------------*/
  /*    check input consistency    */
  /* ------------------------------*/
  const std::map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  //loop materials to be optimized
  std::map<int,std::vector<std::string> >::const_iterator curr;
  for (curr=paramap_.begin(); curr != paramap_.end(); curr++ )
  {
    //check whether this mat exists in the problem
    if ( mats.find(curr->first) == mats.end() )
      dserror("material %d not found in matset", curr->first);
    else
    {
      //check whether input params for this material are valid parameters for optimization
      RCP<MAT::PAR::Material> actmat = mats.at(curr->first);
      std::vector<std::string> actmatparams;
      actmat->Parameter()->OptParams(&actmatparams);
      std::set<std::string> actmatparamsset(actmatparams.begin(),actmatparams.end());
      for(std::vector<std::string>::const_iterator it = curr->second.begin(); it != curr->second.end(); ++it)
      {
        if ( actmatparamsset.find(*it) == actmatparamsset.end() )
          dserror("invalid optimization parameters for material %d", curr->first);
      }
    }
  }

}

//--------------------------------------------------------------------------------------
// update new material parameters from outside
void STR::INVANA::MatParManager::UpdateParams(Teuchos::RCP<Epetra_MultiVector> toadd)
{
  params_o_->Scale(1.0, *params_);
  params_->Update(1.0,*toadd,1.0);

  // bring updated parameters to the elements
  SetParams();
}

//--------------------------------------------------------------------------------------
// set new material parameters from outside
void STR::INVANA::MatParManager::ResetParams()
{
  params_->Scale(1.0, *params_o_);

  // bring updated parameters to the elements
  SetParams();
}
//--------------------------------------------------------------------------------------
// evaluate gradient based on dual solution
void STR::INVANA::MatParManager::Evaluate(Teuchos::RCP<Epetra_MultiVector> dfint, Teuchos::RCP<Epetra_Vector> disdual)
{

  for (int i=0; i< discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lColElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (paramap_.find(elematid) == paramap_.end() )
      continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;

    std::vector<std::string> actparams = paramap_.at( elematid );
    std::vector<std::string>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      p.set("action", "calc_struct_internalforce");
      p.set("matparderiv", *it);

      //initialize element vectors
      int ndof = actele->NumNode()*3;
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(1);
      actele->LocationVector(*discret_,la,false);

      actele->Evaluate(p,*discret_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      //std::cout << "this elevector " << elevector1 << std::endl;
      //std::cout << "dual solution  " << *disdual << std::endl;

      // derivativ w.r.t the actual parameter
      Epetra_Vector dfinti(*dofrowmap_,true);
      for (int l=0; l<(int)la[0].lm_.size(); l++)
      {
        int err = dfinti.SumIntoGlobalValue(la[0].lm_.at(l),0,elevector1[l]);
        if (err!=0) dserror("gid %d is not on this processor", la[0].lm_.at(l));
      }

      //std::cout << "dfinti  " << dfinti << std::endl;

      //wait for all processors to be able to perform the scalar multiplication
      discret_->Comm().Barrier();

      double val = 0.0;
      dfinti.Dot(*disdual,&val);

      int success = dfint->SumIntoGlobalValue(i,parapos_.at(elematid).at(it-actparams.begin()),val);
      if (success!=0) dserror("gid %d is not on this processor", i);


    }//loop this elements material parameters (only the one to be optimized)

  }//loop elements

  //simple tikhonov regularization on the parameter vector
  Epetra_MultiVector tmpvec=Epetra_MultiVector(*params_);
  tmpvec.Update(1.0,*params_,-1.0);

  dfint->Update(reg_weight_,tmpvec,1.0);


}

void STR::INVANA::MatParManager::SetParams()
{
  for (int i=0; i< discret_->NumMyColElements(); i++)
  {
    Teuchos::RCP<MAT::Material> actmat = discret_->lColElement(i)->Material();
    int actmatid = actmat->Parameter()->Id();
    if ( paramap_.find( actmatid ) != paramap_.end() )
    {
      std::vector<std::string> actparams = paramap_.at( actmatid );
      std::vector<std::string>::const_iterator it;
      for ( it=actparams.begin(); it!=actparams.end(); it++)
      {
        double val = (*(*params_)( parapos_.at(actmatid).at(it-actparams.begin()) ))[i];
        actmat->Init(val,*it);
      }

    }//is optimizable material?

  }//loop elements
}
