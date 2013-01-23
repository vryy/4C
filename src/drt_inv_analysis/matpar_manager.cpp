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
sizeparams_(0),
numparams_(0),
discret_(discret),
dofrowmap_(NULL),
elecolmap_(NULL),
dfint_(Teuchos::null),
params_(Teuchos::null)
{
  if (not discret_->Filled() || not discret_->HaveDofs())
      dserror("Discretisation is not complete or has no dofs!");
  else
  {
    dofrowmap_ = discret_->DofRowMap();
    elecolmap_ = discret_->ElementColMap();
  }

  // for now each element has its own set of parameters
  sizeparams_ = discret_->NumGlobalElements();

  // set up maps to links against materials, parameters and materials/parameters for the optimization
  SetupMatOptMap();

  // initialize "state" vectors
  params_ = Teuchos::rcp(new Epetra_MultiVector(*elecolmap_, numparams_, true));
  dfint_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,sizeparams_,true));

  // read initial parameters and sort into params_
  InitParams();

  // set initial parameters, so every elements has values
  SetParams();


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
          (*params_)(parapos_.at(it->first).at(jt-it->second.begin()))->PutScalar(actmat->GetDouble(*jt));
        }
      }
      break;
      default:
        dserror("Material not provided by the Material Manager for Optimization");
      break;
    }
  }

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
  int matid, actmatid;
  while (pstream >> word2)
  {
    matid = std::strtol(&word2[0],NULL,10);
    if (matid != 0)
    {
      actmatid = matid;
      continue;
    }

    if (word2!="none")
    {
      paramap_[actmatid].push_back(word2);
      parapos_[actmatid].push_back(numparams_);
      numparams_ += 1;
    }
    else
      dserror("Give the parameters for the respective materials");
  }

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
// Set new material parameters
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

    }//is optimizable material??

  }//loop elements
}

//--------------------------------------------------------------------------------------
// get new material parameters from outside
void STR::INVANA::MatParManager::UpdateParams(Teuchos::RCP<Epetra_MultiVector> newparams)
{
  params_->Update(1.0,*newparams,0.0);
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

      //cout << "this elevector " << elevector1 << endl;
      //cout << "dual solution  " << *disdual << endl;

      // derivativ w.r.t the actual parameter
      Epetra_Vector dfinti(*dofrowmap_,true);
      for (int l=0; l<(int)la[0].lm_.size(); l++)
        dfinti.SumIntoGlobalValue(la[0].lm_.at(l),0,elevector1[l]);

      //cout << "dfinti  " << dfinti << endl;

      double val = 0.0;
      dfinti.Dot(*disdual,&val);
      //cout << "val " << val << endl;
      dfint->SumIntoGlobalValue(i,parapos_.at(elematid).at(it-actparams.begin()),val);

    }//loop this elements material parameters (only the one to be optimized)

  }//loop elements

}





