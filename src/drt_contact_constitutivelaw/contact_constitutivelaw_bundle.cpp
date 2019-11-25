/*----------------------------------------------------------------------*/
/*! \file

\brief This bundle is used to hold all contact constitutive laws from the input file

\level 3

\maintainer Nora Hagmeyer

*/
/*----------------------------------------------------------------------*/



/*----------------------------------------------------------------------*/
/* headers */
#include "contact_constitutivelaw_bundle.H"

#include "contactconstitutivelaw.H"
#include "contactconstitutivelaw_parameter.H"
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Bundle::Bundle() : readfromproblem_(0) {}
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Bundle::Insert(
    int matid, Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> mat)
{
  map_.insert(std::pair<int, Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container>>(matid, mat));
}

/*----------------------------------------------------------------------*/
int CONTACT::CONSTITUTIVELAW::Bundle::Find(const int id) const
{
  if (map_.find(id) == map_.end())
    return -1;
  else
    return map_.find(id)->first;
}

/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Bundle::MakeParameters()
{
  for (std::map<int, Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container>>::iterator m = map_.begin();
       m != map_.end(); ++m)
  {
    int lawid = m->first;

    // 1st try
    {
      // indirectly add quick access parameter members
      Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> law =
          CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(lawid);
      // check if allocation was successful
      Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> lawpar = m->second;
      if (law != Teuchos::null) continue;
    }
    dserror("Allocation of quick access parameters failed for contact constitutivelaw %d", lawid);
  }
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> CONTACT::CONSTITUTIVELAW::Bundle::ById(
    const int num) const
{
  std::map<int, Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container>>::const_iterator m =
      map_.find(num);

  if (map_.size() == 0) dserror("No contact constitutivelaws available, num=%d", num);

  if (m == map_.end())
    dserror("Contact Constitutive Law 'Law %d' could not be found", num);
  else
    return m->second;

  // catch up
  return Teuchos::null;
}
