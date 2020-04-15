/*-----------------------------------------------------------------------*/
/*! \file
\brief Contact interface capable of TSI

\maintainer Matthias Mayr

\level 3
*/
/*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    seitz 08/15|
 *----------------------------------------------------------------------*/
#include "contact_defines.H"
#include "contact_tsi_interface.H"
#include "contact_interface.H"
#include "contact_node.H"
#include "friction_node.H"
#include "contact_element.H"

#include "../drt_mortar/mortar_dofset.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include "../drt_inpar/inpar_mortar.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            seitz 08/15|
 *----------------------------------------------------------------------*/
CONTACT::CoTSIInterface::CoTSIInterface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    bool selfcontact)
    : CONTACT::CoInterface(interfaceData_ptr, id, comm, dim, icontact, selfcontact)
{
  return;
}
void CONTACT::CoTSIInterface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
    LINALG::SparseMatrix& linstickDISglobal, LINALG::SparseMatrix& linstickTEMPglobal,
    Epetra_Vector& linstickRHSglobal)
{
  CONTACT::CoInterface::AssembleLinStick(linstickLMglobal, linstickDISglobal, linstickRHSglobal);

#ifdef CONSISTENTSTICK

  // do the additional thermal linearizations
  // be aware, that there is another contribution to the linDIS part, since the temperature used
  // in the calculation of the friction coefficient may be the master side temperature, which has
  // a displacement derivative due to the projection


  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_, *slipnodes_);
  Teuchos::RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_, *slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements() == 0) return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");
  if (ftype != INPAR::CONTACT::friction_coulomb) dserror("only coulomb friction for CTSI");

  double frcoeff_in =
      InterfaceParams().get<double>("FRCOEFF");  // the friction coefficient from the input
  double cn = InterfaceParams().get<double>("SEMI_SMOOTH_CN");

  // some things that are not implemented
  bool gp_slip = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "FRLESS_FIRST");
  if (gp_slip || frilessfirst)
    dserror("this fancy option for the contact algorithm is not implemented for TSI");

  // consistent equation is:
  // mu(d,T)*(z_n-c_n*g)ut = 0

  typedef std::map<int, double>::const_iterator _CI;
  // loop over all stick nodes of the interface
  for (int i = 0; i < sticknodes->NumMyElements(); ++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

    const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
    const LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = n.Dot(lm);
    const double wgap = cnode->CoData().Getg();
    LINALG::Matrix<3, 1> txi(cnode->CoData().txi(), true);
    LINALG::Matrix<3, 1> teta(cnode->CoData().teta(), true);
    LINALG::Matrix<3, 1> jump(cnode->FriData().jump(), true);
    double jump_txi = jump.Dot(txi);
    double jump_teta = jump.Dot(teta);

    // row number of entries
    std::vector<int> row(Dim() - 1);
    if (Dim() == 2)
      dserror("2-Dimensional thermo-contact not implemented 'cause there's no 2-D thermo element");
    else if (Dim() == 3)
    {
      row[0] = stickt->GID(2 * i);
      row[1] = stickt->GID(2 * i) + 1;
    }
    else
      dserror("ERROR: AssemblelinStick: Dimension not correct");

    std::map<int, double> dfrdT, dfrdD;
    cnode->derivFrCoeffTemp(frcoeff_in, dfrdT, dfrdD);

    double fac = lm_n - cn * wgap;
    for (_CI p = dfrdT.begin(); p != dfrdT.end(); ++p)
    {
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); j++)
        {
          linstickTEMPglobal.Assemble(
              fac * jump_txi * p->second * txi(j), cnode->Dofs()[j], p->first);
          linstickTEMPglobal.Assemble(
              fac * jump_teta * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linstickTEMPglobal.Assemble(fac * jump_txi * p->second, row[0], p->first);
        linstickTEMPglobal.Assemble(fac * jump_teta * p->second, row[1], p->first);
      }
    }
    for (_CI p = dfrdD.begin(); p != dfrdD.end(); ++p)
    {
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); j++)
        {
          linstickDISglobal.Assemble(
              fac * jump_txi * p->second * txi(j), cnode->Dofs()[j], p->first);
          linstickDISglobal.Assemble(
              fac * jump_teta * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linstickDISglobal.Assemble(fac * jump_txi * p->second, row[0], p->first);
        linstickDISglobal.Assemble(fac * jump_teta * p->second, row[1], p->first);
      }
    }
  }
#else
  // if we don't do the inconsistent stick linearization, it reduces to a mesh-tying
  // problem and hence it is independent of the friction coefficient and therefore
  // independent of the temperature. So do nothing here.
#endif

  return;
}
void CONTACT::CoTSIInterface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
    LINALG::SparseMatrix& linslipDISglobal, LINALG::SparseMatrix& linslipTEMPglobal,
    Epetra_Vector& linslipRHSglobal)
{
  CONTACT::CoInterface::AssembleLinSlip(linslipLMglobal, linslipDISglobal, linslipRHSglobal);

  // do the additional thermal linearizations
  // be aware, that there is another contribution to the linDIS part, since the temperature used
  // in the calculation of the friction coefficient may be the master side temperature, which has
  // a displacement derivative due to the projection


#ifdef CONSISTENTSLIP
  dserror("CONSISTENTSLIP not implemented for CTSI");
#endif

  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> slipnodes = slipnodes_;
  Teuchos::RCP<Epetra_Map> slipt = slipt_;

  // nothing to do if no stick nodes
  if (slipnodes->NumMyElements() == 0) return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");
  if (ftype != INPAR::CONTACT::friction_coulomb) dserror("only coulomb friction for CTSI");

  if (Dim() != 3) dserror("CTSI only for 3D");

  double frcoeff_in =
      InterfaceParams().get<double>("FRCOEFF");  // the friction coefficient from the input
  double ct_input = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  double cn_input = InterfaceParams().get<double>("SEMI_SMOOTH_CN");

  // some things that are not implemented
  bool gp_slip = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "FRLESS_FIRST");
  if (gp_slip || frilessfirst)
    dserror("this fancy option for the contact algorithm is not implemented for TSI");

  // consistent equation is:
  // || z_t^tr || z_t - \mu (z_n -c_n*wgap) z_t^tr = 0

  typedef std::map<int, double>::const_iterator _CI;
  // loop over all stick nodes of the interface
  for (int i = 0; i < slipnodes->NumMyElements(); ++i)
  {
    int gid = slipnodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

    const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
    const LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = n.Dot(lm);
    const double wgap = cnode->CoData().Getg();
    LINALG::Matrix<3, 1> txi(cnode->CoData().txi(), true);
    LINALG::Matrix<3, 1> teta(cnode->CoData().teta(), true);
    LINALG::Matrix<3, 1> jump(cnode->FriData().jump(), true);
    const double jump_txi = jump.Dot(txi);
    const double jump_teta = jump.Dot(teta);
    const double lm_txi = lm.Dot(txi);
    const double lm_teta = lm.Dot(teta);

    // row number of entries
    std::vector<int> row(Dim() - 1);
    if (Dim() == 2)
      dserror("2-Dimensional thermo-contact not implemented 'cause there's no 2-D thermo element");
    else if (Dim() == 3)
    {
      row[0] = slipt->GID(2 * i);
      if (slipt->GID(2 * i) != cnode->Dofs()[1]) dserror("why");
      row[1] = slipt->GID(2 * i) + 1;
      if (slipt->GID(2 * i) + 1 != cnode->Dofs()[2]) dserror("why2");
    }
    else
      dserror("ERROR: AssemblelinStick: Dimension not correct");

    std::map<int, double> dfrdT, dfrdD;
    cnode->derivFrCoeffTemp(frcoeff_in, dfrdT, dfrdD);

    double cn = cn_input;
    double ct = ct_input;

    double fac = lm_n - cn * wgap;
    for (_CI p = dfrdT.begin(); p != dfrdT.end(); ++p)
    {
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); ++j)
        {
          linslipTEMPglobal.Assemble(
              -fac * (lm_txi + ct * jump_txi) * p->second * txi(j), cnode->Dofs()[j], p->first);
          linslipTEMPglobal.Assemble(
              -fac * (lm_teta + ct * jump_teta) * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linslipTEMPglobal.Assemble(-fac * (lm_txi + ct * jump_txi) * p->second, row[0], p->first);
        linslipTEMPglobal.Assemble(-fac * (lm_teta + ct * jump_teta) * p->second, row[1], p->first);
      }
    }
    for (_CI p = dfrdD.begin(); p != dfrdD.end(); ++p)
    {
      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); ++j)
        {
          linslipDISglobal.Assemble(
              -fac * (lm_txi + ct * jump_txi) * p->second * txi(j), cnode->Dofs()[j], p->first);
          linslipDISglobal.Assemble(
              -fac * (lm_teta + ct * jump_teta) * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linslipDISglobal.Assemble(-fac * (lm_txi + ct * jump_txi) * p->second, row[0], p->first);
        linslipDISglobal.Assemble(-fac * (lm_teta + ct * jump_teta) * p->second, row[1], p->first);
      }
    }
  }

  return;
}


void CONTACT::CoTSIInterface::AssembleLinConduct(LINALG::SparseMatrix& linConductDISglobal,
    LINALG::SparseMatrix& linConductTEMPglobal, LINALG::SparseMatrix& linConductThermoLMglobal,
    LINALG::SparseMatrix& linConductContactLMglobal)
{
  // nothing to do if no active contact nodes
  if (activenodes_->NumMyElements() == 0) return;

  // heat transfer factors
  const double alpha_s = InterfaceParams().get<double>("HEATTRANSSLAVE");
  const double alpha_m = InterfaceParams().get<double>("HEATTRANSMASTER");
  const double beta_bar = alpha_s * alpha_m / (alpha_s + alpha_m);
  const double delta_bar = alpha_s / (alpha_s + alpha_m);

  AssembleDualMassLumped(linConductThermoLMglobal, linConductDISglobal);

  AssembleLinDM_X(&linConductDISglobal, NULL, -delta_bar, LinDM_Diss, activenodes_);
  AssembleDM_linDiss(&linConductDISglobal, NULL, &linConductContactLMglobal, NULL, -delta_bar);

  AssembleDM_LMn(-beta_bar, &linConductTEMPglobal);
  AssembleLinLMnDM_Temp(-beta_bar, &linConductDISglobal, &linConductContactLMglobal);

  return;
}


void CONTACT::CoTSIInterface::AssembleDualMassLumped(
    LINALG::SparseMatrix& dualMassGlobal, LINALG::SparseMatrix& linDualMassGlobal)
{
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::CoNode* conode = dynamic_cast<CONTACT::CoNode*>(node);

    if (conode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDualMass: Node ownership inconsistency!");

    double thermo_lm = conode->CoTSIData().ThermoLM();
    std::map<int, std::map<int, double>>& derivDualMass = conode->CoData().GetDerivD();

    if (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(InterfaceParams(), "LM_QUAD") !=
        INPAR::MORTAR::lagmult_const)
    {
      /**********************************************dual mass matrix ******/
      if (conode->MoData().GetD().size() > 0)
      {
        const GEN::pairedvector<int, double>& dualMassmap = conode->MoData().GetD();
        GEN::pairedvector<int, double>::const_iterator colcurr;

        for (colcurr = dualMassmap.begin(); colcurr != dualMassmap.end(); ++colcurr)
        {
          DRT::Node* knode = Discret().gNode(colcurr->first);
          if (!knode) dserror("node not found");
          CoNode* kcnode = dynamic_cast<CoNode*>(knode);
          if (!kcnode) dserror("node not found");

          // create the mass matrix
          dualMassGlobal.FEAssemble(colcurr->second, conode->Dofs()[0], kcnode->Dofs()[0]);
        }
      }

      for (std::map<int, std::map<int, double>>::const_iterator a = derivDualMass.begin();
           a != derivDualMass.end(); ++a)
      {
        int sgid = a->first;
        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("ERROR: Cannot find node with gid %", sgid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        for (std::map<int, double>::const_iterator b = a->second.begin(); b != a->second.end(); ++b)
          // val               row                col
          linDualMassGlobal.FEAssemble(b->second * thermo_lm, csnode->Dofs()[0], b->first);
      }
    }

    else  // INPAR::MORTAR::lagmult_const
    {
      if (conode->NumElement() != 1)
        dserror(
            "some inconsistency: for lagmult_const every slave node may only have one element "
            "attached (it's the center-node!)");

      CONTACT::CoElement* coele = dynamic_cast<CONTACT::CoElement*>(conode->Elements()[0]);
      if (!coele) dserror("this should be a contact element");

      GEN::pairedvector<int, double> derivArea(coele->NumNode() * Dim());
      double area = coele->ComputeAreaDeriv(derivArea);

      dualMassGlobal.FEAssemble(area, conode->Dofs()[0], conode->Dofs()[0]);

      for (GEN::pairedvector<int, double>::const_iterator p = derivArea.begin();
           p != derivArea.end(); ++p)
        linDualMassGlobal.FEAssemble(p->second * thermo_lm, conode->Dofs()[0], p->first);
    }
  }
  return;
}

void CONTACT::CoTSIInterface::AssembleLinDM_X(LINALG::SparseMatrix* linD_X,
    LINALG::SparseMatrix* linM_X, const double fac, const LinDM_X_mode mode,
    const Teuchos::RCP<Epetra_Map> node_rowmap)
{
  // get out if there's nothing to do
  if (linD_X == NULL && linM_X == NULL) return;

  // no dissipation without friction
  if (!friction_ && mode == LinDM_Diss) return;

  const double dt = InterfaceParams().get<double>("TIMESTEP");

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < node_rowmap->NumMyElements(); ++j)
  {
    int gid = node_rowmap->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // Mortar matrix D and M derivatives
    std::map<int, std::map<int, double>>& dderiv = cnode->CoData().GetDerivD();
    std::map<int, std::map<int, double>>& mderiv = cnode->CoData().GetDerivM();

    // current Lagrange multipliers
    double lm = 0.;
    switch (mode)
    {
      case LinDM_ThermoLM:
        lm = cnode->CoTSIData().ThermoLM();
        break;
      case LinDM_Diss:
      {
        CONTACT::FriNode* frnode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (frnode == NULL)
          lm = 0.;
        else
        {
          double dval = 1.;
          if (cnode->MoData().GetD().size() == 0)
            continue;
          else
            cnode->MoData().GetD()[cnode->Id()];
          const LINALG::Matrix<3, 1> lmc(cnode->MoData().lm(), true);
          const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
          const LINALG::Matrix<3, 1> jump(frnode->FriData().jump(), true);
          double diss = (-lmc.Dot(jump) + lmc.Dot(n) * jump.Dot(n)) / (dt * dval);
          lm = diss;
        }
        break;
      }
      case linDM_ContactLMnormal:
      {
        const LINALG::Matrix<3, 1> contact_LM(cnode->MoData().lm(), true);
        const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
        lm = contact_LM.Dot(n);
        break;
      }
      default:
        dserror("mode not implemented");
        break;
    }

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    int mastersize = (int)mderiv.size();
    std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();
    std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    if (linD_X)
    {
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("ERROR: Cannot find node with gid %", sgid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        int row = csnode->Dofs()[0];
        std::map<int, double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = scolcurr->first;
          double val = lm * (scolcurr->second);
          ++scolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          // std::cout << "Assemble LinD: " << row << " " << col << " " << val << std::endl;
          if (abs(val) > 1.0e-12) linD_X->FEAssemble(val * fac, row, col);
        }

        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr != thisdderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
      }

      // check for completeness of DerivD-Slave-iteration
      if (scurr != dderiv.end())
        dserror("ERROR: AssembleLinDM: Not all DISP slave entries of DerivD considered!");
    } /******************************** Finished with LinDMatrix **********/

    /********************************************** LinMMatrix **********/
    if (linM_X)
    {
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %", mgid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->CoData().GetDerivM()[mgid];
        int mapsize = (int)(thismderiv.size());

        int row = cmnode->Dofs()[0];
        std::map<int, double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = mcolcurr->first;
          double val = lm * (mcolcurr->second);
          ++mcolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          // std::cout << "Assemble LinM: " << row << " " << col << " " << val << std::endl;
          if (abs(val) > 1.0e-12) linM_X->FEAssemble(-val * fac, row, col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
      }

      // check for completeness of DerivM-Master-iteration
      if (mcurr != mderiv.end())
        dserror("ERROR: AssembleLinDM: Not all master entries of DerivM considered!");
    } /******************************** Finished with LinMMatrix **********/
  }

  return;
}

void CONTACT::CoTSIInterface::AssembleDM_linDiss(LINALG::SparseMatrix* d_LinDissDISP,
    LINALG::SparseMatrix* m_LinDissDISP, LINALG::SparseMatrix* d_LinDissContactLM,
    LINALG::SparseMatrix* m_LinDissContactLM, const double fac)
{
  // get out if there's nothing to do
  if (d_LinDissDISP == NULL && m_LinDissDISP == NULL && d_LinDissContactLM == NULL &&
      m_LinDissContactLM == NULL)
    return;

  // there's no dissipation without friction
  if (!friction_) return;

  typedef std::map<int, double>::const_iterator _cim;
  typedef GEN::pairedvector<int, double>::const_iterator _cip;

  const double dt = InterfaceParams().get<double>("TIMESTEP");

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < activenodes_->NumMyElements(); ++j)
  {
    int gid = activenodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    FriNode* fnode = dynamic_cast<FriNode*>(cnode);
    // if this is not a friction node, there is no dissipation so go to the next one
    if (fnode == NULL) continue;

    // get nodal normal
    const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);

    // projection into tangential plane = 1 - n \otimes n
    LINALG::Matrix<3, 3> tang_proj(true);
    for (int i = 0; i < 3; ++i) tang_proj(i, i) = 1.;
    tang_proj.MultiplyNT(-1., n, n, 1.);

    // get D entry of this node
    // remember: D is diagonal
    int id = cnode->Id();
    const std::map<int, double>& derivD = cnode->CoData().GetDerivD(id);
    const double dval = cnode->MoData().GetD()[cnode->Id()];

    // get nodal values
    LINALG::Matrix<3, 1> jump(fnode->FriData().jump());
    const LINALG::Matrix<3, 1> lm(cnode->MoData().lm());
    LINALG::Matrix<3, 1> lm_t;
    lm_t.Multiply(tang_proj, lm);
    const double lm_n = lm.Dot(n);
    const double jump_n = jump.Dot(n);
    LINALG::Matrix<3, 1> jump_tan;
    jump_tan.Multiply(tang_proj, jump);


    // linearization w.r.t. displacements
    if (d_LinDissDISP != NULL || m_LinDissDISP != NULL)
    {
      // get nodal jump and deriv
      const std::vector<std::map<int, double>>& derivJump = fnode->FriData().GetDerivJump();
      const std::vector<GEN::pairedvector<int, double>>& derivN = cnode->CoData().GetDerivN();

      // calculate derivative of Dissipation *******************************
      std::map<int, double> derivDiss;
      for (int i = 0; i < 3; ++i)
      {
        for (_cim p = derivJump[i].begin(); p != derivJump[i].end(); ++p)
          derivDiss[p->first] -= ((lm(i) - lm_n * n(i)) * p->second) / (dt * dval);
        for (_cip p = derivN[i].begin(); p != derivN[i].end(); ++p)
          derivDiss[p->first] += ((lm_n * jump(i) + jump_n * lm(i)) * p->second) / (dt * dval);
        for (_cim p = derivD.begin(); p != derivD.end(); ++p)
          derivDiss[p->first] +=
              (-lm.Dot(jump) + lm.Dot(n) * jump.Dot(n)) / (dt * dval * dval) * (-p->second);
      }

      // put everything together*******************************************
      /**************************************************** D-matrix ******/
      if (d_LinDissDISP != NULL)
        if ((cnode->MoData().GetD()).size() > 0)
          for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
          {
            DRT::Node* knode = Discret().gNode(k->first);
            if (!knode) dserror("node not found");
            CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);
            for (_cim currDeriv = derivDiss.begin(); currDeriv != derivDiss.end(); ++currDeriv)
              d_LinDissDISP->FEAssemble(
                  fac * k->second * currDeriv->second, kcnode->Dofs()[0], currDeriv->first);
          }
      /**************************************************** M-matrix ******/
      if (m_LinDissDISP != NULL)
        if ((cnode->MoData().GetM()).size() > 0)
          for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
          {
            DRT::Node* knode = Discret().gNode(k->first);
            if (!knode) dserror("node not found");
            CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);
            for (_cim currDeriv = derivDiss.begin(); currDeriv != derivDiss.end(); ++currDeriv)
              m_LinDissDISP->FEAssemble(
                  fac * k->second * currDeriv->second, kcnode->Dofs()[0], currDeriv->first);
          }
    }  // linearization w.r.t. displacements

    // linearization wrt contact Lagrange multiplier
    if (d_LinDissContactLM != NULL)
      // put everything together*******************************************
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetD()).size() > 0)
        for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
        {
          DRT::Node* knode = Discret().gNode(k->first);
          if (!knode) dserror("node not found");
          CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);
          for (int d = 0; d < 3; ++d)
            d_LinDissContactLM->FEAssemble(
                -fac * k->second * jump_tan(d) / (dt * dval), kcnode->Dofs()[0], cnode->Dofs()[d]);
        }

    /**************************************************** M-matrix ******/
    if (m_LinDissContactLM)
      if ((cnode->MoData().GetM()).size() > 0)
        for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
        {
          DRT::Node* knode = Discret().gNode(k->first);
          if (!knode) dserror("node not found");
          CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);
          for (int d = 0; d < 3; ++d)
            m_LinDissContactLM->FEAssemble(
                -fac * k->second * jump_tan(d) / (dt * dval), kcnode->Dofs()[0], cnode->Dofs()[d]);
        }
  }  // loop over all active LM slave nodes (row map)
}

void CONTACT::CoTSIInterface::AssembleLinLMnDM_Temp(
    const double fac, LINALG::SparseMatrix* lin_disp, LINALG::SparseMatrix* lin_lm)
{
  // get out if there's nothing to do
  if (lin_disp == NULL) dserror("called to assemble something but didn't provide a matrix");

  typedef std::map<int, double>::const_iterator _cim;
  typedef GEN::pairedvector<int, double>::const_iterator _cip;
  typedef std::map<int, std::map<int, double>>::const_iterator _cimm;

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < activenodes_->NumMyElements(); ++j)
  {
    int gid = activenodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
    const LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = lm.Dot(n);

    for (_cimm k = cnode->CoData().GetDerivD().begin(); k != cnode->CoData().GetDerivD().end(); ++k)
    {
      DRT::Node* knode = Discret().gNode(k->first);
      if (!knode) dserror("ERROR: Cannot find node with gid %", gid);
      double temp_k = dynamic_cast<CoNode*>(knode)->CoTSIData().Temp();
      for (_cim l = k->second.begin(); l != k->second.end(); ++l)
        if (abs(l->second) > 1.e-12)
          lin_disp->FEAssemble(fac * lm_n * temp_k * l->second, cnode->Dofs()[0], l->first);
    }
    for (_cimm k = cnode->CoData().GetDerivM().begin(); k != cnode->CoData().GetDerivM().end(); ++k)
    {
      DRT::Node* knode = Discret().gNode(k->first);
      if (!knode) dserror("ERROR: Cannot find node with gid %", gid);
      double temp_k = dynamic_cast<CoNode*>(knode)->CoTSIData().Temp();
      for (_cim l = k->second.begin(); l != k->second.end(); ++l)
        if (abs(l->second) > 1.e-12)
          lin_disp->FEAssemble(-fac * lm_n * temp_k * l->second, cnode->Dofs()[0], l->first);
    }

    for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
    {
      DRT::Node* knode = Discret().gNode(k->first);
      if (!knode) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnodek = dynamic_cast<CoNode*>(knode);
      double temp_k = cnodek->CoTSIData().Temp();
      for (int d = 0; d < 3; ++d)
        for (_cip l = cnode->CoData().GetDerivN()[d].begin();
             l != cnode->CoData().GetDerivN()[d].end(); ++l)
          if (abs(l->second) > 1.e-12)
            lin_disp->FEAssemble(
                lm(d) * k->second * l->second * fac * temp_k, cnode->Dofs()[0], l->first);

      for (int d = 0; d < 3; ++d)
        lin_lm->FEAssemble(n(d) * temp_k * fac * k->second, cnode->Dofs()[0], cnode->Dofs()[d]);
    }

    for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
    {
      DRT::Node* knode = Discret().gNode(k->first);
      if (!knode) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnodek = dynamic_cast<CoNode*>(knode);
      double temp_k = cnodek->CoTSIData().Temp();
      for (int d = 0; d < 3; ++d)
        for (_cip l = cnode->CoData().GetDerivN()[d].begin();
             l != cnode->CoData().GetDerivN()[d].end(); ++l)
          if (abs(l->second) > 1.e-12)
            lin_disp->FEAssemble(
                -lm(d) * k->second * l->second * fac * temp_k, cnode->Dofs()[0], l->first);

      for (int d = 0; d < 3; ++d)
        lin_lm->FEAssemble(-n(d) * temp_k * fac * k->second, cnode->Dofs()[0], cnode->Dofs()[d]);
    }
  }

  return;
}

void CONTACT::CoTSIInterface::AssembleDM_LMn(const double fac, LINALG::SparseMatrix* DM_LMn)
{
  // get out if there's nothing to do
  if (DM_LMn == NULL) dserror("called to assemble something but didn't provide a matrix");

  typedef std::map<int, double>::const_iterator _cim;
  typedef GEN::pairedvector<int, double>::const_iterator _cip;

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < activenodes_->NumMyElements(); ++j)
  {
    int gid = activenodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    const LINALG::Matrix<3, 1> n(cnode->MoData().n(), true);
    const LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = lm.Dot(n);

    for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
      if (abs(k->second) > 1.e-12)
      {
        DRT::Node* knode = Discret().gNode(k->first);
        if (!knode) dserror("node not found");
        CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);
        DM_LMn->FEAssemble(fac * lm_n * k->second, cnode->Dofs()[0], kcnode->Dofs()[0]);
      }

    for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
      if (abs(k->second) > 1.e-12)
      {
        DRT::Node* knode = Discret().gNode(k->first);
        if (!knode) dserror("node not found");
        CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);
        DM_LMn->FEAssemble(-fac * lm_n * k->second, cnode->Dofs()[0], kcnode->Dofs()[0]);
      }
  }
  return;
}


void CONTACT::CoTSIInterface::AssembleInactive(LINALG::SparseMatrix* linConductThermoLM)
{
  // get out if there's nothing to do
  if (linConductThermoLM == NULL)
    dserror("called to assemble something but didn't provide a matrix");

  // inactive nodes
  Teuchos::RCP<Epetra_Map> inactivenodes = LINALG::SplitMap(*snoderowmap_, *activenodes_);

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < inactivenodes->NumMyElements(); ++j)
  {
    int gid = inactivenodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Active()) dserror("something with active and inactive maps is not in order");

    linConductThermoLM->FEAssemble(1., cnode->Dofs()[0], cnode->Dofs()[0]);
  }
}
/*----------------------------------------------------------------------*
 |  initialize / reset interface for tsi                     seitz 08/15|
 *----------------------------------------------------------------------*/
void CONTACT::CoTSIInterface::Initialize()
{
  CONTACT::CoInterface::Initialize();

  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
  {
    CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscret_->lColNode(i));
    node->InitializeTSIDataContainer(
        imortar_.get<double>("TEMP_REF"), imortar_.get<double>("TEMP_DAMAGE"));
    node->CoTSIData().Clear();
  }
  return;
}
