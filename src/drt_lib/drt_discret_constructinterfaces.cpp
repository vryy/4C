/*----------------------------------------------------------------------*/
/*!
\file drt_constructinterfaces.cpp

\brief Construct interface within a discretization

<pre>
Maintainer: Fedderik van der Bos
            bos@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_ParameterListExceptions.hpp>

#include "drt_discret.H"

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#endif

int DRT::Discretization::ConstructInterfaces()
{
	if (not Filled()) dserror("FillComplete() must be called before call to ConstructInterface()\n");
	if (not RequiresInterfaces()) dserror("This discretization does not require interface, so no need to construct them");
		
	map<int,RefCountPtr<DRT::Element> >::iterator elecurr;
	for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
	{
		int nfacesl = elecurr->second->NumFace();
		for (int ifacel=0; ifacel<nfacesl; ifacel++)
		{
			DRT::Node* primnode = elecurr->second->FacePrimNode(ifacel);
			vector<int> nodeidsl = elecurr->second->FaceNodeIds(ifacel,DRT::Element::approachleft);
			
			int nelementsincloud = primnode->NumElement();
			DRT::Element** elementsincloud = primnode->Elements();
			for (int ielem=0; ielem<nelementsincloud; ielem++)
			{
				int nfacesr = elementsincloud[ielem]->NumFace();
				for (int ifacer=0; ifacer<nfacesr; ifacer++)
				{
					vector<int> nodeidsr = elementsincloud[ielem]->FaceNodeIds(ifacer,DRT::Element::approachright);
					
					if (nodeidsr==nodeidsl)
					{
						DRT::Element* newface = elecurr->second->CreateFace(ifacel,elementsincloud[ielem]->Id(),ifacer);
						if (newface==NULL) dserror("No interface has been created between element %d and %d\n",
								                   elecurr->second->Id(),elementsincloud[ielem]->Id());
						this->AddElement(rcp(newface));
					}
					
				}
			}
			
			
		}
		
		
	}
	
	return 0;
}

#endif
