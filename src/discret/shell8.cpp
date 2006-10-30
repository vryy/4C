/*!----------------------------------------------------------------------
\file discret.cpp
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

#include "shell8.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Shell8::Shell8(int id) :
CCADISCRETIZATION::Element(id,element_shell8)
{

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Shell8::Shell8(const CCADISCRETIZATION::Shell8& old) :
CCADISCRETIZATION::Element(old)
{

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Shell8::Clone() const
{
  CCADISCRETIZATION::Shell8* newelement = new CCADISCRETIZATION::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* CCADISCRETIZATION::Shell8::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizetype   = sizeof(enum ElementType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);
   
  // calculate size of vector
  size = 
  sizeint +     // holds size itself
  sizeint +     // holds Id()
  sizetype +    // holds type of element
  0;            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  memcpy(&data[position],&size,sizeint);
  position += sizeint;
  
  // add Id()
  const int id = Id();
  memcpy(&data[position],&id,sizeint);
  position += sizeint;
  
  // add type
  ElementType etype = Type();
  memcpy(&data[position],&etype,sizetype);
  position += sizeint;
    
  // continue to add stuff here  
    
  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Shell8::Unpack(const char* data)
{
  const int sizeint    = sizeof(int);
  const int sizetype   = sizeof(enum ElementType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  int position = 0;

  // extract size
  int size = 0;
  memcpy(&size,&data[position],sizeint);
  position += sizeint;
  
  // extract id
  int id = 0;
  memcpy(&id,&data[position],sizeint);
  position += sizeint;
  id_ = id;
  
  // extract type
  ElementType etype = element_none;
  memcpy(&etype,&data[position],sizetype);
  position += sizeint;
  etype_ = etype;

  // extract more stuff here



  if (position != size)
  {
    cout << "CCADISCRETIZATION::Shell8::Unpack:\n"
         << "Mismatch in size of data " << size << " <-> " << position << endl
         << __FILE__ << ":" << __LINE__ << endl;
    exit(EXIT_FAILURE);
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Shell8::~Shell8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Shell8::Print() const
{
  cout << "Element " << Id() << " Shell8 " << endl;
  return;
}
















#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
