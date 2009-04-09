/*!----------------------------------------------------------------------
\file so_hex8.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/contchainnetw.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/anisotropic_balzani.H"

// inverse design object
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#include "inversedesign.H"
#endif

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(int id, int owner) :
DRT::Element(id,element_so_hex8,owner),
data_()
{
  kintype_ = soh8_totlag;
  eastype_ = soh8_easnone;
  neas_ = 0;
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
//  for (int i=0; i<NUMGPT_SOH8; ++i)
//    invJ_[i].Shape(3,3);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  prestress_ = rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8,NUMGPT_SOH8));
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  invdesign_ = rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOH8,NUMGPT_SOH8));
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(const DRT::ELEMENTS::So_hex8& old) :
DRT::Element(old),
kintype_(old.kintype_),
eastype_(old.eastype_),
neas_(old.neas_),
data_(old.data_),
detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (int i=0; i<(int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH8 x NUMDIM_SOH8?
    //invJ_[i].Shape(old.invJ_[i].M(),old.invJ_[i].N());
    invJ_[i] = old.invJ_[i];
  }

#if defined(PRESTRESS) || defined(POSTSTRESS)
  prestress_ = rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  invdesign_ = rcp(new DRT::ELEMENTS::InvDesign(*(old.invdesign_)));
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_hex8::Clone() const
{
  DRT::ELEMENTS::So_hex8* newelement = new DRT::ELEMENTS::So_hex8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_hex8::Shape() const
{
  return hex8;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // kintype_
  AddtoPack(data,kintype_);
  // eastype_
  AddtoPack(data,eastype_);
  // neas_
  AddtoPack(data,neas_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  // prestress_
  vector<char> tmpprestress(0);
  prestress_->Pack(tmpprestress);
  AddtoPack(data,tmpprestress);
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  // invdesign_
  vector<char> tmpinvdesign(0);
  invdesign_->Pack(tmpinvdesign);
  AddtoPack(data,tmpinvdesign);
#endif

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // kintype_
  ExtractfromPack(position,data,kintype_);
  // eastype_
  ExtractfromPack(position,data,eastype_);
  // neas_
  ExtractfromPack(position,data,neas_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  // prestress_
  vector<char> tmpprestress(0);
  ExtractfromPack(position,data,tmpprestress);
  prestress_->Unpack(tmpprestress);
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  // invdesign_
  vector<char> tmpinvdesign(0);
  ExtractfromPack(position,data,tmpinvdesign);
  invdesign_->Unpack(tmpinvdesign);
#endif

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size;
  ExtractfromPack(position,data,size);
  invJ_.resize(size);
  for (int i=0; i<size; ++i)
  {
    //invJ_[i].Shape(0,0);
    ExtractfromPack(position,data,invJ_[i]);
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::~So_hex8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Print(ostream& os) const
{
  os << "So_hex8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 02/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_expol
(
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>& stresses,
    LINALG::Matrix<NUMDOF_SOH8,1>& elevec1,
    LINALG::Matrix<NUMDOF_SOH8,1>& elevec2
)
{
  // static variables, that are the same for every element
  static LINALG::Matrix<NUMNOD_SOH8,NUMGPT_SOH8> expol;
  static bool isfilled;

  
  if (isfilled==false)
  {
    double sq3=sqrt(3.0);
    expol(0,0)=1.25+0.75*sq3;
    expol(0,1)=-0.25-0.25*sq3;
    expol(0,2)=-0.25+0.25*sq3;
    expol(0,3)=-0.25-0.25*sq3;
    expol(0,4)=-0.25-0.25*sq3;
    expol(0,5)=-0.25+0.25*sq3;
    expol(0,6)=1.25-0.75*sq3;
    expol(0,7)=-0.25+0.25*sq3;
    expol(1,1)=1.25+0.75*sq3;
    expol(1,2)=-0.25-0.25*sq3;
    expol(1,3)=-0.25+0.25*sq3;
    expol(1,4)=-0.25+0.25*sq3;
    expol(1,5)=-0.25-0.25*sq3;
    expol(1,6)=-0.25+0.25*sq3;
    expol(1,7)=1.25-0.75*sq3;
    expol(2,2)=1.25+0.75*sq3;
    expol(2,3)=-0.25-0.25*sq3;
    expol(2,4)=1.25-0.75*sq3;
    expol(2,5)=-0.25+0.25*sq3;
    expol(2,6)=-0.25-0.25*sq3;
    expol(2,7)=-0.25+0.25*sq3;
    expol(3,3)=1.25+0.75*sq3;
    expol(3,4)=-0.25+0.25*sq3;
    expol(3,5)=1.25-0.75*sq3;
    expol(3,6)=-0.25+0.25*sq3;
    expol(3,7)=-0.25-0.25*sq3;
    expol(4,4)=1.25+0.75*sq3;
    expol(4,5)=-0.25-0.25*sq3;
    expol(4,6)=-0.25+0.25*sq3;
    expol(4,7)=-0.25-0.25*sq3;
    expol(5,5)=1.25+0.75*sq3;
    expol(5,6)=-0.25-0.25*sq3;
    expol(5,7)=-0.25+0.25*sq3;
    expol(6,6)=1.25+0.75*sq3;
    expol(6,7)=-0.25-0.25*sq3;
    expol(7,7)=1.25+0.75*sq3;

    for (int i=0;i<NUMNOD_SOH8;++i)
    {
      for (int j=0;j<i;++j)
      {
        expol(i,j)=expol(j,i);
      }
    }
    isfilled = true;
  }

  LINALG::Matrix<NUMNOD_SOH8,NUMSTR_SOH8> nodalstresses;

  //nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);
  nodalstresses.Multiply(expol, stresses);

  for (int i=0;i<NUMNOD_SOH8;++i){
    elevec1(NODDOF_SOH8*i) = nodalstresses(i,0);
    elevec1(NODDOF_SOH8*i+1) = nodalstresses(i,1);
    elevec1(NODDOF_SOH8*i+2) = nodalstresses(i,2);
  }
  for (int i=0;i<NUMNOD_SOH8;++i){
    elevec2(NODDOF_SOH8*i) = nodalstresses(i,3);
    elevec2(NODDOF_SOH8*i+1) = nodalstresses(i,4);
    elevec2(NODDOF_SOH8*i+2) = nodalstresses(i,5);
  }
}

/*----------------------------------------------------------------------*
 |  allocate and return So_hex8Register (public)                maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_hex8::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Soh8Register(Type()));
}

  /*====================================================================*/
  /* 8-node hexhedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (r,s,t) of nodes
   * of biunit cube [-1,1]x[-1,1]x[-1,1]
   *  8-node hexahedron: node 0,1,...,7
   *                      t
   *                      |
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------s
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      r                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::VisNames(map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);

  if (Material()->MaterialType() == INPAR::MAT::m_contchainnetw){
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber3";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber4";
    names[fiber] = 3; // 3-dim vector
    fiber = "FiberCell1";
    names[fiber] = 3; // 3-dim vector
    fiber = "FiberCell2";
    names[fiber] = 3; // 3-dim vector
    fiber = "FiberCell3";
    names[fiber] = 3; // 3-dim vector
    fiber = "l1";
    names[fiber] = 1;
    fiber = "l2";
    names[fiber] = 1;
    fiber = "l3";
    names[fiber] = 1;
//    fiber = "l1_0";
//    names[fiber] = 1;
//    fiber = "l2_0";
//    names[fiber] = 1;
//    fiber = "l3_0";
//    names[fiber] = 1;
  }
  if ((Material()->MaterialType() == INPAR::MAT::m_artwallremod) ||
      (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic))
  {
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
  }
  if (Material()->MaterialType() == INPAR::MAT::m_anisotropic_balzani){
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisData(name,data);

  if (Material()->MaterialType() == INPAR::MAT::m_contchainnetw){
    RefCountPtr<MAT::Material> mat = Material();
    MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
    if (!chain->Initialized()){
      data[0] = 0.0; data[1] = 0.0; data[2] = 0.0;
    } else {
      RCP<vector<vector<double> > > gplis = chain->Getli();
      RCP<vector<vector<double> > > gpli0s = chain->Getli0();
      RCP<vector<LINALG::Matrix<3,3> > > gpnis = chain->Getni();

      vector<double> centerli (3,0.0);
      vector<double> centerli_0 (3,0.0);
      for (int i = 0; i < (int)gplis->size(); ++i) {
        LINALG::Matrix<3,1> loc(&(gplis->at(i)[0]));
        //Epetra_SerialDenseVector loc(CV,&(gplis->at(i)[0]),3);
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,gpnis->at(i),loc,0.0);
        glo.Multiply(gpnis->at(i),loc);
        // Unfortunately gpnis is a vector of Epetras, to change this
        // I must begin at a deeper level...
        centerli[0] += glo(0);
        centerli[1] += glo(1);
        centerli[2] += glo(2);

//        centerli[0] += gplis->at(i)[0];
//        centerli[1] += gplis->at(i)[1];
//        centerli[2] += gplis->at(i)[2];
//
        centerli_0[0] += gplis->at(i)[0];
        centerli_0[1] += gplis->at(i)[1];
        centerli_0[2] += gplis->at(i)[2];
      }
      centerli[0] /= gplis->size();
      centerli[1] /= gplis->size();
      centerli[2] /= gplis->size();

      centerli_0[0] /= gplis->size();
      centerli_0[1] /= gplis->size();
      centerli_0[2] /= gplis->size();

      // just the unit cell of the first gp
      int gp = 0;
//      Epetra_SerialDenseVector loc(CV,&(gplis->at(gp)[0]),3);
//      Epetra_SerialDenseVector glo(3);
//      glo.Multiply('N','N',1.0,gpnis->at(gp),loc,0.0);
      LINALG::Matrix<3,3> T(gpnis->at(gp).A(),true);
      vector<double> gpli =  chain->Getli()->at(gp);

      if (name == "Fiber1"){
        if ((int)data.size()!=3) dserror("size mismatch");
        data[0] = centerli[0]; data[1] = -centerli[1]; data[2] = -centerli[2];
      } else if (name == "Fiber2"){
        data[0] = centerli[0]; data[1] = centerli[1]; data[2] = -centerli[2];
      } else if (name == "Fiber3"){
        data[0] = centerli[0]; data[1] = centerli[1]; data[2] = centerli[2];
      } else if (name == "Fiber4"){
        data[0] = -centerli[0]; data[1] = -centerli[1]; data[2] = centerli[2];
      } else if (name == "FiberCell1"){
        LINALG::Matrix<3,1> e(true);
        e(0) = gpli[0];
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,T,e,0.0);
        glo.Multiply(T, e);
        data[0] = glo(0); data[1] = glo(1); data[2] = glo(2);
      } else if (name == "FiberCell2"){
        LINALG::Matrix<3,1> e(true);
        e(1) = gpli[1];
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,T,e,0.0);
        glo.Multiply(T, e);
        data[0] = glo(0); data[1] = glo(1); data[2] = glo(2);
      } else if (name == "FiberCell3"){
        LINALG::Matrix<3,1> e(true);
        e(2) = gpli[2];
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,T,e,0.0);
        glo.Multiply(T, e);
        data[0] = glo(0); data[1] = glo(1); data[2] = glo(2);
      } else if (name == "l1"){
        data[0] = centerli_0[0];
      } else if (name == "l2"){
        data[0] = centerli_0[1];
      } else if (name == "l3"){
        data[0] = centerli_0[2];
//      } else if (name == "l1_0"){
//        data[0] = centerli_0[0];
//      } else if (name == "l2_0"){
//        data[0] = centerli_0[1];
//      } else if (name == "l3_0"){
//        data[0] = centerli_0[2];
      } else if (name == "Owner"){
        if ((int)data.size()<1) dserror("Size mismatch");
        data[0] = Owner();
      } else {
        cout << name << endl;
        dserror("Unknown VisData!");
      }
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod){
    MAT::ArtWallRemod* art = static_cast <MAT::ArtWallRemod*>(Material().get());
    vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
    vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
    if (name == "Fiber1"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a1[0]; data[1] = a1[1]; data[2] = a1[2];
    } else if (name == "Fiber2"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a2[0]; data[1] = a2[1]; data[2] = a2[2];
    } else if (name == "Owner"){
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = Owner();
    } else {
      cout << name << endl;
      dserror("Unknown VisData!");
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic){
    MAT::ViscoAnisotropic* art = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
    vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
    if (name == "Fiber1"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a1[0]; data[1] = a1[1]; data[2] = a1[2];
    } else if (name == "Fiber2"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a2[0]; data[1] = a2[1]; data[2] = a2[2];
    } else if (name == "Owner"){
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = Owner();
    } else {
      cout << name << endl;
      dserror("Unknown VisData!");
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_anisotropic_balzani){
    MAT::AnisotropicBalzani* balz = static_cast <MAT::AnisotropicBalzani*>(Material().get());
    if (name == "Fiber1"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = balz->Geta1().at(0); data[1] = balz->Geta1().at(1); data[2] = balz->Geta1().at(2);
    } else if (name == "Fiber2"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = balz->Geta2().at(0); data[1] = balz->Geta2().at(1); data[2] = balz->Geta2().at(2);
    } else if (name == "Owner"){
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = Owner();
    } else {
      cout << name << endl;
      dserror("Unknown VisData!");
    }
  }


  return;
}



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register::Soh8Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register::Soh8Register(
                               const DRT::ELEMENTS::Soh8Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register* DRT::ELEMENTS::Soh8Register::Clone() const
{
  return new DRT::ELEMENTS::Soh8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Register::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register::~Soh8Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Register::Print(ostream& os) const
{
  os << "Soh8Register ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

