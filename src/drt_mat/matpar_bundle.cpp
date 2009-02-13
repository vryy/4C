/*======================================================================*/
/*!
\file matpar_bundle.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_bundle.H"

/*----------------------------------------------------------------------*/
MAT::PAR::Bundle::Bundle()
  : materialreadfromproblem_(0)
{
}


/*----------------------------------------------------------------------*/
MAT::PAR::Bundle::~Bundle()
{
}


/*----------------------------------------------------------------------*/
void MAT::PAR::Bundle::Insert(
  int matid,
  Teuchos::RCP<MAT::PAR::Material> mat
  )
{
  matmap_.insert(std::pair<int,Teuchos::RCP<MAT::PAR::Material> >(matid,mat));
}

/*----------------------------------------------------------------------*/
int MAT::PAR::Bundle::Find(
  const int id
  ) const
{
  if (matmap_.find(id) == matmap_.end())
    return -1;
  else
    return matmap_.find(id)->first;
}

/*----------------------------------------------------------------------*/
void MAT::PAR::Bundle::MakeParameters()
{
  for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::iterator m=matmap_.begin();
       m!=matmap_.end();
       ++m)
  {
    int matid = m->first;
    // indirectly add quick access parameter members
    Teuchos::RCP<MAT::Material> elemat = MAT::Material::Factory(matid);
    // check if allocation was successful
    Teuchos::RCP<MAT::PAR::Material> mat = m->second;
    if (mat->Parameter()==NULL) dserror("Allocation of quick access parameters failed for material MAT %d", matid);
  }
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::PAR::Material> MAT::PAR::Bundle::ById(
  const int num
  ) const
{ 
  std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator m = matmap_.find(num);
        
  if (matmap_.size() == 0)
    dserror("No materials avialable");
        
  if (m == matmap_.end())
    dserror("Material 'MAT %d' could not be found", num);
  else
    return m->second;
        
  // catch up
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
int MAT::PAR::Bundle::FirstIdByType(
  const INPAR::MAT::MaterialType type
  ) const
{
  std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator m;
        
  int id = -1;
  for (m=matmap_.begin(); m!=matmap_.end(); ++m)
  {
    if (m->second->Type() == type)
    {
      id = m->first;
      break;
    }
  }

  return id;
}

#endif
