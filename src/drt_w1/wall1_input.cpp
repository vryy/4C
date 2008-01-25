/*!----------------------------------------------------------------------
\file wall1_input.cpp
\brief

<pre>
Maintainer: Markus Gitterle 
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

extern "C" 
{
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         mgit 03/07
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
}
#include "wall1.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              mwgit 03/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Wall1::ReadElement()
{
  // read element's nodes
  int ierr=0;
  int nnode=0;
  int nodes[9]; 
  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    nnode = 4;
    frint_n("QUAD4",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("QUAD8",&ierr);
  if (ierr==1)
  {
    nnode = 8;
    frint_n("QUAD8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("QUAD9",&ierr);
  if (ierr==1)
  {
    nnode = 9;
    frint_n("QUAD9",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    nnode = 3;
    frint_n("TRI3",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI6",&ierr);
  if (ierr==1)
  {
    nnode = 6;
    frint_n("TRI6",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;
  
  SetNodeIds(nnode,nodes);
  
  // read number of material model
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of WALL1 element failed");
  
  // read wall thickness
  thickness_ = 1.0;
  frdouble("THICK",&thickness_,&ierr);
  if (ierr!=1) dserror("Reading of WALL1 element failed");

  // read gaussian points
  int ngp[2];
  frint_n("GP",ngp,2,&ierr);
  if (ierr!=1) dserror("Reading of WALL1 element failed");
  gaussrule_ = getGaussrule(ngp); 

  //read 2D problem type
  frchk("PLANE_STRESS",&ierr);
  if (ierr==1) wtype_ = plane_stress; 
  frchk("PLANE_STRAIN",&ierr);
  if (ierr==1) wtype_ = plane_strain;
  
  //read model (EAS or not)
  frchk("EAS_Model",&ierr);
  if (ierr==1)
  {
	
	iseas_=true;
    // EAS enhanced deformation gradient parameters
    Epetra_SerialDenseMatrix alpha(4,1);

//    alpha (0,0)=1;
//    alpha (1,0)=2;
//    alpha (2,0)=3;
//    alpha (3,0)=4;
                
    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    Epetra_SerialDenseVector feas(4);
    // EAS matrix K_{alpha alpha}, also called Dtilde
    Epetra_SerialDenseMatrix invKaa(2,2);
    // EAS matrix K_{d alpha}
    Epetra_SerialDenseMatrix Kda(2,2);
    
    // save EAS data into element container easdata_
    data_.Add("alpha",alpha);
    data_.Add("feas",feas);
    data_.Add("invKaa",invKaa);
    data_.Add("Kda",Kda);

  }

  //read lokal or global stresses
  char buffer [50];
  frchar("STRESSES",buffer,&ierr);
  if (ierr)
  {
     if      (strncmp(buffer,"XY",2)==0)       stresstype_ = w1_xy;
     else if (strncmp(buffer,"RS",2)==0)       stresstype_ = w1_rs;
     else dserror ("Reading of WALL1 element failed");
  }

//  //read TSI
// #if 0
//  tsi_couptyp = tsi_coup_none;  /* default */
//  char buffer[50];
//  frchar("TSI_COUPTYP",buffer,&ierr);
//  if (ierr)
//  {
//    if (strncmp(buffer,"None",4)==0) tsi_couptyp_ = tsi_coup_none;
//    if (strncmp(buffer,"Thermconf",9)==0) tsi_couptyp_ = tsi_coup_thermconf;
//    if (strncmp(buffer,"Thermcreate",11)==0) tsi_couptyp_ = tsi_coup_thermcreate;
//  }
//#endif
        
  return true;
} // Wall1::ReadElement()


//Get gaussrule on dependance of gausspoints 

DRT::UTILS::GaussRule2D DRT::ELEMENTS::Wall1::getGaussrule(int* ngp)
{
  DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;

  switch (Shape())
  {
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    {
       if ( (ngp[0]==1) && (ngp[1]==1) )
       {
         rule = DRT::UTILS::intrule_quad_1point;
       }
       else if ( (ngp[0]==2) && (ngp[1]==2) )
       {
         rule = DRT::UTILS::intrule_quad_4point;
       }
       else if ( (ngp[0]==3) && (ngp[1]==3) )
       {
         rule = DRT::UTILS::intrule_quad_9point;
       }
       else
         dserror("Unknown number of Gauss points for quad element");  
       break;
    }
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
        switch (ngp[0])
        {
            case 1:                   /* constant */
            {
                rule = DRT::UTILS::intrule_tri_1point;
                break;
            }
            case 3:  /* quadratic - type 1 and 2*/
            {
                //  GAUSS INTEGRATION 3 SAMPLING POINTS, DEG.OF PRECISION 2 
                if (ngp[1]-1 == 0)  // integration 1
                    {
                    rule = DRT::UTILS::intrule_tri_3point; 
                    }
                else if (ngp[1]-1 == 1)  // integration 2
                    {
                    rule = DRT::UTILS::intrule_tri_3point_on_corners; 
                    }
                else
                    {
                    dserror("Integration case %g is not available\n", ngp[1]);
                    }
     
                break;
            }
            case 6: 
            {
                rule = DRT::UTILS::intrule_tri_6point;
                break;
            }
            default:
            {
                cout << ngp[0] << "  " << ngp[1] << endl;
                dserror("Unknown number of Gauss points for tri element");
            }
        }
        break;        
    default:
       dserror("Unknown distype");
    }
  } 
  return rule;
}



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
