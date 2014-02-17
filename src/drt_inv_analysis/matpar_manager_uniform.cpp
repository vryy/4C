/*----------------------------------------------------------------------*/
/*!

<pre>
Maintainer: Jonas Biehler
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"
#include "matpar_manager_uniform.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"
#include "smc_particle.H"

STR::INVANA::MatParManagerUniform::MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret){ ; }

void STR::INVANA::MatParManagerUniform::ComputeParamsMultiVectorFromSMCParticlePosition(Teuchos::RCP<Epetra_MultiVector>& params , std::vector<double> my_global_params)
{

  //check if size is a match
  if(my_global_params.size()!=(unsigned) numparams_)
    dserror("Size mismatch try again");

  //int numparams = my_particle.GetSizeOfPosition();
  params = Teuchos::rcp(new Epetra_MultiVector(*elecolmap_,numparams_,true));


  // get all materials from drt::problem
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  std::map<int,std::vector<std::string> >::const_iterator it;
  // loop over all mats that we have in the optimization
  for (it=paramap_.begin(); it!=paramap_.end(); it++)
  {
    // loop all mats in the DRT::Problem
    Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(it->first);
    switch(actmat->Parameter()->Type())
    {
    // for now this is only implemented for aaa_neohook_stopro
    case INPAR::MAT::m_aaaneohooke_stopro:
    {
      // get vector with mat parameter names
      std::vector<std::string> mat_param_names = it->second;
      int curr_mat_id = it->first;
      for (unsigned int i=0; i< mat_param_names.size();i++)
      {
        //get vector in which opt param id is stored
        std::vector<int> opt_param_ids = parapos_.at(curr_mat_id);
        (*params)( opt_param_ids.at(i))->PutScalar(my_global_params.at( opt_param_ids.at(i)) );
      }
    }
    break;
    default:
      dserror("Material not provided by the Material Manager for Optimization");
      break;
    }
  }

  //return params;

}


