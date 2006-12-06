/*!----------------------------------------------------------------------
\file shell8.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "shell8.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_dserror.H"
extern "C" 
{
#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
}
#include "dstrc.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::Shell8(int id, int owner) :
DRT::Element(id,element_shell8,owner),
forcetype_(s8_none),
thickness_(0.0),
ngptri_(0),
nhyb_(0),
ans_(0),
sdc_(1.0),
material_(0)
{
  DSTraceHelper dst("Shell8::Shell8");  
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::Shell8(const DRT::Elements::Shell8& old) :
DRT::Element(old),
forcetype_(old.forcetype_),
thickness_(old.thickness_),
ngptri_(old.ngptri_),
nhyb_(old.nhyb_),
ans_(old.ans_),
sdc_(old.sdc_),
material_(old.material_)
{
  DSTraceHelper dst("Shell8::Shell8");  
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  for (int i=0; i<5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Shell8::Clone() const
{
  DSTraceHelper dst("Shell8::Clone");  
  DRT::Elements::Shell8* newelement = new DRT::Elements::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Elements::Shell8::Pack(int& size) const
{
  DSTraceHelper dst("Shell8::Pack");  
  const int sizeint    = sizeof(int);
  const int sizeforcetype = sizeof(enum ForceType);
  const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // pack base class
  int         basesize = 0;
  const char* basedata = Element::Pack(basesize);
   
  // calculate size of vector data
  size = 
  sizeint +        // size itself
  sizeint +        // type of this instance of ParObject, see top of ParObject.H
  basesize +       // base class data
  sizeforcetype +  // forcetype_
  sizedouble +     // thickness_
  3*sizeint +      // ngp_
  sizeint +        // ngptri_  
  sizeint +        // nhyb_
  5*sizeint +      // eas_
  sizeint +        // ans_
  sizedouble +     // sdc_
  sizeint +        // material_
  0;            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Shell8;
  AddtoPack(position,data,type);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // forcetype_
  AddtoPack(position,data,&forcetype_,sizeof(enum ForceType));
  // thickness_
  AddtoPack(position,data,thickness_);
  // ngp_
  AddtoPack(position,data,ngp_,3*sizeint);
  // ngptri_
  AddtoPack(position,data,ngptri_);
  // nhyb_
  AddtoPack(position,data,nhyb_);
  // eas_
  AddtoPack(position,data,eas_,5*sizeint);
  // ans_
  AddtoPack(position,data,ans_);
  // sdc_
  AddtoPack(position,data,sdc_);
  // material_
  AddtoPack(position,data,material_);
  // continue to add stuff here  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Elements::Shell8::Unpack(const char* data)
{
  DSTraceHelper dst("Shell8::Unpack");  
  const int sizeint    = sizeof(int);
  //const int sizeforcetype = sizeof(enum ForceType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  int position = 0;

  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_Shell8) dserror("Wrong instance type in data");
  
  // extract base class
  int basesize = SizePack(&data[position]);
  Element::Unpack(&data[position]);
  position += basesize;
  
  // forcetype_
  ExtractfromPack(position,data,forcetype_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeint);
  // ngptri_
  ExtractfromPack(position,data,ngptri_);
  // nhyb_
  ExtractfromPack(position,data,nhyb_);
  // eas_
  ExtractfromPack(position,data,eas_,5*sizeint);
  // ans_
  ExtractfromPack(position,data,ans_);
  // sdc_
  ExtractfromPack(position,data,sdc_);
  // material_
  ExtractfromPack(position,data,material_);
  // extract more stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::~Shell8()
{
  DSTraceHelper dst("Shell8::~Shell8");  
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::Print(ostream& os) const
{
  DSTraceHelper dst("Shell8::Print");  
  os << "Shell8 ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Shell8::Evaluate(ParameterList& params, 
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("Shell8::Evaluate");  
  DRT::Elements::Shell8::ActionType act = Shell8::none;
  
  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Shell8::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Shell8::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Shell8::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Shell8::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Shell8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = Shell8::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Shell8::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Shell8::calc_struct_fsiload;
  else dserror("Unknown type of action for Shell8");
  
  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);
  
  switch(act)
  {
    case calc_struct_linstiff:
      dserror("Case not yet implemented");
    break;
    case calc_struct_nlnstiff:
      dserror("Case not yet implemented");
    break;
    case calc_struct_internalforce:
      dserror("Case not yet implemented");
    break;
    case calc_struct_linstiffmass:
      dserror("Case not yet implemented");
    break;
    case calc_struct_nlnstiffmass: // do mass, stiffness and internal forces
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual");
      if (disp==null || res==null) dserror("Cannot get state vectors");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      
    }
    break;
    case calc_struct_stress:
      dserror("Case not yet implemented");
    break;
    case calc_struct_eleload:
      dserror("Case not yet implemented");
    break;
    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;
    default:
      dserror("Unknown type of action for Shell8");
  }
  
  
  return 0;
}














#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
