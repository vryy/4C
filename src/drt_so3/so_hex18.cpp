/*!----------------------------------------------------------------------
\file so_hex18.cpp
\brief 18-node hexahedral (bi-quadratic linear)
\level 1

<pre>
\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "so_hex18.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

DRT::ELEMENTS::So_hex18Type DRT::ELEMENTS::So_hex18Type::instance_;

DRT::ELEMENTS::So_hex18Type& DRT::ELEMENTS::So_hex18Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex18Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_hex18* object = new DRT::ELEMENTS::So_hex18(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex18Type::Create( const std::string eletype,
                                                                const std::string eledistype,
                                                                const int id,
                                                                const int owner )
{
  if ( eletype=="SOLIDH18" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex18(id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex18Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex18(id,owner));
  return ele;
}

void DRT::ELEMENTS::So_hex18Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_hex18Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_hex18Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH18"];

  defs["HEX18"]
    .AddIntVector("HEX18",18)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddOptionalNamedDoubleVector("RAD",3)
    .AddOptionalNamedDoubleVector("AXI",3)
    .AddOptionalNamedDoubleVector("CIR",3)
    .AddOptionalNamedDoubleVector("FIBER1",3)
    .AddOptionalNamedDoubleVector("FIBER2",3)
    .AddOptionalNamedDoubleVector("FIBER3",3)
    .AddOptionalNamedDouble("STRENGTH")
    .AddOptionalNamedDouble("HU")
    .AddOptionalNamedDouble("lambda")
    ;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex18::So_hex18(int id, int owner) :
So_base(id,owner)
{
  invJ_.resize(NUMGPT_SOH18, LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18>(true));
  detJ_.resize(NUMGPT_SOH18, 0.0);
  InitGp();

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex18::So_hex18(const DRT::ELEMENTS::So_hex18& old) :
So_base(old),
detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  // can this size be anything but NUMDIM_SOH27 x NUMDIM_SOH27?
  for (int i=0; i<(int)invJ_.size(); ++i)
    invJ_[i] = old.invJ_[i];
  InitGp();

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_hex18::Clone() const
{
  DRT::ELEMENTS::So_hex18* newelement = new DRT::ELEMENTS::So_hex18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  So_base::Pack(data);

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
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  So_base::Unpack(basedata);

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size,LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18>(true) );
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::Print(std::ostream& os) const
{
  os << "So_hex18 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)               seitz 11/14 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_hex18::Volumes()
{
  std::vector<Teuchos::RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                          seitz 11/14 |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_hex18::Surfaces()
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
|  get vector of lines (public)                            seitz 11/14 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_hex18::Lines()
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
|  Return names of visualization data (public)             seitz 11/14 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::VisNames(std::map<std::string,int>& names)
{
 SolidMaterial()->VisNames(names);

 return;
}

/*----------------------------------------------------------------------*
|  Return visualization data (public)                      seitz 11/14 |
*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex18::VisData(const std::string& name, std::vector<double>& data)
{
 // Put the owner of this element into the file (use base class method for this)
 if (DRT::Element::VisData(name,data))
   return true;

 return SolidMaterial()->VisData(name, data, NUMGPT_SOH18, this->Id());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex18::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);

  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)

  Teuchos::RCP<MAT::Material> mat = Material();

  SolidMaterial()->Setup(NUMGPT_SOH18, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
   //kintype_ = soh8_linear;
   dserror ("Only nonlinear kinematics for SO_SH8 implemented!");
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else dserror ("Reading SO_HEX18 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}

void DRT::ELEMENTS::So_hex18::InitGp()
{
  xsi_.resize(NUMGPT_SOH18,LINALG::Matrix<NUMDIM_SOH18,1>(true));
  wgt_.resize(NUMGPT_SOH18,0.);
  DRT::UTILS::IntPointsAndWeights<NUMDIM_SOH18> intpoints(DRT::UTILS::intrule_hex_18point);
  for (int gp=0; gp<NUMGPT_SOH18; ++gp)
  {
    wgt_.at(gp)=(intpoints.IP().qwgt)[gp];
    for (int idim=0; idim<NUMDIM_SOH18; idim++)
      xsi_.at(gp)(idim) = (intpoints.IP().qxg)[gp][idim];
  }
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           seitz 11/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex18::Evaluate(Teuchos::ParameterList&  params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,1> elevec2(elevec2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,1> elevec3(elevec3_epetra.A(),true);

  // start with "none"
  DRT::ELEMENTS::So_hex18::ActionType act = So_hex18::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act = So_hex18::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act = So_hex18::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act = So_hex18::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act = So_hex18::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act = So_hex18::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act = So_hex18::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act = So_hex18::calc_struct_stress;
  else if (action=="calc_struct_eleload")                         act = So_hex18::calc_struct_eleload;
  else if (action=="calc_struct_update_istep")                    act = So_hex18::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")                     act = So_hex18::calc_struct_reset_istep;
  else if (action=="calc_struct_reset_all")                       act = So_hex18::calc_struct_reset_all;
  else if (action=="postprocess_stress")                          act = So_hex18::postprocess_stress;
  else if (action=="calc_struct_recover")                         act = So_hex18::calc_recover;
  else dserror("Unknown type of action for So_hex8: %s", action.c_str());

  // what should the element do
  switch(act)
  {

    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

        nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                                    INPAR::STR::stress_none,INPAR::STR::strain_none);
      break;
    }


    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18> myemat(true);

      nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
        INPAR::STR::stress_none,INPAR::STR::strain_none);

      break;
    }

    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    case calc_struct_linstiffmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      // need current velocities and accelerations (for non constant mass matrix)
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");

      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

        nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none);

      if (act==calc_struct_nlnstifflmass) Lumpmass(&elemat2);

      break;
    }

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char> > stressdata = params.get<Teuchos::RCP<std::vector<char> > >("stress",Teuchos::null);
        Teuchos::RCP<std::vector<char> > straindata = params.get<Teuchos::RCP<std::vector<char> > >("strain",Teuchos::null);
        if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D> strain;
        INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

          nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);

        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
        }

        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
        }
      }
    }
    break;

    //==================================================================================
    // postprocess stresses/strains at gauss points
    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      std::string stresstype = params.get<std::string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D> gpstress(((*gpstressmap)[gid])->A(),true);
      Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
      if (poststress==Teuchos::null)
        dserror("No element stress/strain vector available");
      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        dserror("no node-based stress output");
      }
      else if (stresstype=="cxyz")
      {
        const Epetra_BlockMap& elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
          {
            double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH18; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH18;
          }
        }
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      SolidMaterial()->Update();
      Update();
    }
    break;

    //==================================================================================
    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    //==================================================================================
    case calc_struct_reset_all:
    {
      // Reset of history for materials
      SolidMaterial()->ResetAll(NUMGPT_SOH18);
    }
    break;

    case calc_recover:
    {
      Teuchos::RCP<const Epetra_Vector> res  =
          discretization.GetState("residual displacement");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      Recover(myres);
    }
      break;

  //==================================================================================
  default:
    dserror("Unknown type of action for So_hex18");
    break;
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 11/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex18::EvaluateNeumann(Teuchos::ParameterList&   params,
                                             DRT::Discretization&      discretization,
                                             DRT::Condition&           condition,
                                             std::vector<int>&         lm,
                                             Epetra_SerialDenseVector& elevec1,
                                             Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOH18)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOH18; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }


  // find out whether we will use time curves and get the factors
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  std::vector<double> curvefacs(NUMDIM_SOH18, 1.0);
  for (int i=0; i < NUMDIM_SOH18; ++i)
  {
    const int curvenum = (curve) ? (*curve)[i] : -1;
    if (curvenum>=0 && usetime)
      curvefacs[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }


  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int> >("funct");
  LINALG::Matrix<NUMDIM_SOH18,1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim=0; dim<NUMDIM_SOH18; dim++)
      if ((*funct)[dim] > 0)
        havefunct = havefunct or true;

  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH18; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH18; ++gp) {

    // shape function and derivatives
    LINALG::Matrix<NUMNOD_SOH18,1> shapefunct;
    DRT::UTILS::shape_function<DRT::Element::hex18>(xsi_[gp],shapefunct);
    LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> deriv;
    DRT::UTILS::shape_function_deriv1<DRT::Element::hex18>(xsi_[gp],deriv);

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> jac;
    jac.Multiply(deriv,xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct) {
      for (int dim=0; dim<NUMDIM_SOH18; dim++) {
        xrefegp(dim) = 0.0;
        for (int nodid=0; nodid<NUMNOD_SOH18; ++nodid)
          xrefegp(dim) += shapefunct(nodid) * xrefe(nodid,dim);
      }
    }

    // integration factor
    const double fac = wgt_[gp] * detJ;
    // distribute/add over element load vector
    for(int dim=0; dim<NUMDIM_SOH18; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac
          = (functnum>0)
          ? DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,xrefegp.A(),time,NULL)
          : 1.0;
        const double dim_fac = (*val)[dim] * fac * curvefacs[dim] * functfac;
        for (int nodid=0; nodid<NUMNOD_SOH18; ++nodid)
        {
          elevec1[nodid*NUMDIM_SOH18+dim] += shapefunct(nodid) * dim_fac;
        }
      }
    }


  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_hex18::EvaluateNeumann

int DRT::ELEMENTS::So_hex18::InitJacobianMapping()
{

  LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18> xrefe;
  for (int i=0; i<NUMNOD_SOH18; ++i)
  {
    Node** nodes=Nodes();
    if(!nodes) dserror("Nodes() returned null pointer");
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.reserve(NUMGPT_SOH18);
  detJ_.reserve(NUMGPT_SOH18);


  for (int gp=0; gp<NUMGPT_SOH18; ++gp)
  {
    // reset
    invJ_[gp].Clear();
    detJ_[gp]=0.;

    LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> deriv;
    DRT::UTILS::shape_function_deriv1<DRT::Element::hex18>(xsi_[gp],deriv);

    invJ_[gp].Multiply(deriv,xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp]<0.) return 1;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::nlnstiffmass(
    std::vector<int>& lm,  ///< location matrix
    std::vector<double>& disp,  ///< current displacements
    std::vector<double>& residual,  ///< current residual displ
    LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* stiffmatrix,  ///< element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* massmatrix,  ///< element mass matrix
    LINALG::Matrix<NUMDOF_SOH18,1>* force,  ///< element internal force vector
    LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,  ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  ///< stress output option
    const INPAR::STR::StrainType iostrain   ///< strain output option
    )
{
  LINALG::Matrix<NUMNOD_SOH18,3> xrefe(false);      // X, material coord. of element
  LINALG::Matrix<NUMNOD_SOH18,3> xcurr(false);      // x, current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH18; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH18+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH18+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH18+2];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> N_XYZ;
  // build deformation gradient wrt to material configuration
  LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> defgrd(false);
  // shape functions and their first derivatives
  LINALG::Matrix<NUMNOD_SOH18,1> shapefunct;
  LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> deriv;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH18; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    DRT::UTILS::shape_function<hex18>(xsi_[gp],shapefunct);
    DRT::UTILS::shape_function_deriv1<hex18>(xsi_[gp],deriv);

    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // (material) deformation gradient
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ);

    // calcualte total rcg
    LINALG::Matrix<3,3> cauchygreen(false);
    cauchygreen.MultiplyTN(defgrd,defgrd);
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> glstrain(false);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    // B-operator
    LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH18> bop(false);
    for (int i=0; i<NUMNOD_SOH18; ++i)
    {
      bop(0,NUMDIM_SOH18*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NUMDIM_SOH18*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NUMDIM_SOH18*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NUMDIM_SOH18*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NUMDIM_SOH18*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NUMDIM_SOH18*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NUMDIM_SOH18*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NUMDIM_SOH18*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NUMDIM_SOH18*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NUMDIM_SOH18*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NUMDIM_SOH18*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NUMDIM_SOH18*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NUMDIM_SOH18*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NUMDIM_SOH18*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NUMDIM_SOH18*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NUMDIM_SOH18*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NUMDIM_SOH18*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NUMDIM_SOH18*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> stress(true);
    params.set<int>("gp",gp);
    SolidMaterial()->Evaluate(&defgrd,&glstrain,params,&stress,&cmat,Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    double detJ_w = detJ*wgt_[gp];
    // update internal force vector
    if (force)
      force->MultiplyTN(detJ_w,bop,stress,1.);

    // update stiffness matrix
    if (stiffmatrix)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH18> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ_w); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH18; ++inod) {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
            + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
            + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
            + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<NUMNOD_SOH18; ++jnod) {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH18; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************
    }

    if (massmatrix) // evaluate mass matrix +++++++++++++++++++++++++
    {
      double density = Material()->Density();
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH18; ++inod)
      {
        ifactor = shapefunct(inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH18; ++jnod)
        {
          massfactor = shapefunct(jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH18*inod+0,NUMDIM_SOH18*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH18*inod+1,NUMDIM_SOH18*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH18*inod+2,NUMDIM_SOH18*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
  } // gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                              seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::Lumpmass(LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                               seitz 11/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex18Type::Initialize(DRT::Discretization& dis)
{
  // here we order the nodes such that we have a positive definite jacobian
  //       maybe the python script generating the hex18 elements would be a better place for this.
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex18* actele = dynamic_cast<DRT::ELEMENTS::So_hex18*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex18* failed");
    if (actele->InitJacobianMapping()==1)
      actele->FlipT();
  }
  dis.FillComplete(false,false,false);

  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex18* actele = dynamic_cast<DRT::ELEMENTS::So_hex18*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex18* failed");
    if (actele->InitJacobianMapping()==1)
      dserror("why");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  revert the 3rd parameter direction                      seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18::FlipT()
{
  if (NodeIds()==NULL) dserror("couldn't get node ids");
  // reorder nodes
  int new_nodeids[NUMNOD_SOH18];
  new_nodeids[0] = NodeIds()[9];
  new_nodeids[1] = NodeIds()[10];
  new_nodeids[2] = NodeIds()[11];
  new_nodeids[3] = NodeIds()[12];
  new_nodeids[4] = NodeIds()[13];
  new_nodeids[5] = NodeIds()[14];
  new_nodeids[6] = NodeIds()[15];
  new_nodeids[7] = NodeIds()[16];
  new_nodeids[8] = NodeIds()[17];

  new_nodeids[9]  = NodeIds()[0];
  new_nodeids[10] = NodeIds()[1];
  new_nodeids[11] = NodeIds()[2];
  new_nodeids[12] = NodeIds()[3];
  new_nodeids[13] = NodeIds()[4];
  new_nodeids[14] = NodeIds()[5];
  new_nodeids[15] = NodeIds()[6];
  new_nodeids[16] = NodeIds()[7];
  new_nodeids[17] = NodeIds()[8];

  SetNodeIds(NUMNOD_SOH18,new_nodeids);
  return;
}

LINALG::Matrix<18,3> DRT::ELEMENTS::So_hex18::NodeParamCoord()
{
  LINALG::Matrix<18,3> coord;
  for (int node=0; node<NUMNOD_SOH18; ++node)
  {
    LINALG::Matrix<3,1> nodeCoord = NodeParamCoord(node);
    for (int i=0; i<3; ++i)
      coord(node,i) = nodeCoord(i);
  }
  return coord;
}

LINALG::Matrix<3,1> DRT::ELEMENTS::So_hex18::NodeParamCoord(const int node)
{
  LINALG::Matrix<3,1> coord;

  switch (node%9)
  {
  case 0:
  case 3:
  case 7:
    coord(0)=-1.;break;
  case 4:
  case 6:
  case 8:
    coord(0)=+0.;break;
  case 1:
  case 2:
  case 5:
    coord(0)=+1.;break;
  }
  switch (node%9)
  {
  case 0:
  case 1:
  case 4:
    coord(1)=-1.;break;
  case 5:
  case 7:
  case 8:
    coord(1)=+0.;break;
  case 2:
  case 3:
  case 6:
    coord(1)=+1.;break;
  }

  if (node<9) coord(2)=-1.;
  else        coord(2)=+1.;

  return coord;
}

