/*!----------------------------------------------------------------------
\file So3_volcoupl.cpp

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "so3_volcoupl.H"

#include "../drt_lib/drt_linedefinition.H"

template<class so3_ele, class coupltype>
DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::So3_volcoupl(int id, int owner) :
DRT::Element(id,owner),
so3_ele(id,owner),
coupltype (id,owner)
{
  return;
}

/*!
\brief Copy Constructor

Makes a deep copy of a Element

*/

template<class so3_ele, class coupltype>
DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::So3_volcoupl(const DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>& old) :
    DRT::Element(old),
    so3_ele(old),
    coupltype (old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
DRT::Element* DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::Clone() const
{
  DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>* newelement =
      new DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
void DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data,type);
  // add base class so3_ele Element
  so3_ele::Pack(data);
  // add base class coupling Element
  coupltype::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
void DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class so3_ele Element
  vector<char> basedata(0);
  so3_ele::ExtractfromPack(position,data,basedata);
  so3_ele::Unpack(basedata);
  // extract base class coupling Element
  basedata.clear();
  so3_ele::ExtractfromPack(position,data,basedata);
  coupltype::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::~So3_volcoupl()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
void DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::Print(ostream& os) const
{
  os << "So3_volcoupl ";
  coupltype::Print(os);
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
bool DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::ReadElement(const std::string& eletype,
                                             const std::string& eledistype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  return so3_ele::ReadElement(eletype, eledistype, linedef);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
int DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::NumDofPerNode(const unsigned nds,
                                              const DRT::Node& node) const
{
  if (nds==1)
    return coupltype::NumDofPerNode(nds, node);

  return so3_ele::NumDofPerNode(node);
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
DRT::ElementType & DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::ElementType() const
{
   return coupltype::ElementType();
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
int DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::UniqueParObjectId() const
{
   return coupltype::UniqueParObjectId();
};

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, class coupltype>
int DRT::ELEMENTS::So3_volcoupl< so3_ele, coupltype>::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  typename So3_volcoupl::ActionType act = So3_volcoupl::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_multidofsetcoupling")   act = So3_volcoupl::calc_struct_multidofsetcoupling;
  // what should the element do
  switch(act)
  {
  //==================================================================================
  // coupling terms in force-vector and stiffness matrix
  case So3_volcoupl::calc_struct_multidofsetcoupling:
  {
    coupltype::Evaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;
  //==================================================================================
  default:
  {
    //in some cases we need to write/change some data before evaluating
    coupltype::PreEvaluate(params,
                      discretization,
                      la);

    so3_ele::Evaluate(params,
                      discretization,
                      la[0].lm_,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);

    coupltype::Evaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  } // action

  return 0;
}

// TET4
template class DRT::ELEMENTS::So3_volcoupl<DRT::ELEMENTS::So_tet4,DRT::ELEMENTS::So3_Poro<DRT::Element::tet4> >;
template class DRT::ELEMENTS::So3_volcoupl<DRT::ELEMENTS::So_tet4,DRT::ELEMENTS::So3_Scatra<DRT::Element::tet4> >;
// HEX8
template class DRT::ELEMENTS::So3_volcoupl<DRT::ELEMENTS::So_hex8,DRT::ELEMENTS::So3_Poro<DRT::Element::hex8> >;
template class DRT::ELEMENTS::So3_volcoupl<DRT::ELEMENTS::So_hex8,DRT::ELEMENTS::So3_Scatra<DRT::Element::hex8> >;

