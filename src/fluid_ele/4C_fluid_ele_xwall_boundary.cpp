// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_boundary_parent_calc.hpp"
#include "4C_fluid_ele_xwall.hpp"


FOUR_C_NAMESPACE_OPEN

Discret::Elements::FluidXWallBoundaryType Discret::Elements::FluidXWallBoundaryType::instance_;

Discret::Elements::FluidXWallBoundaryType& Discret::Elements::FluidXWallBoundaryType::instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::Element> Discret::Elements::FluidXWallBoundaryType::create(
    const int id, const int owner)
{
  return nullptr;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::FluidXWallBoundary::FluidXWallBoundary(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::Elements::Fluid* parent,
    const int lsurface)
    : FluidBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::Elements::FluidXWallBoundary::FluidXWallBoundary(
    const Discret::Elements::FluidXWallBoundary& old)
    : FluidBoundary(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::FluidXWallBoundary::clone() const
{
  Discret::Elements::FluidXWallBoundary* newelement =
      new Discret::Elements::FluidXWallBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void Discret::Elements::FluidXWallBoundary::print(std::ostream& os) const
{
  os << "FluidXWallBoundary ";
  Element::print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void Discret::Elements::FluidXWallBoundary::location_vector(const Core::FE::Discretization& dis,
    Core::Elements::LocationArray& la, bool doDirichlet, const std::string& condstring,
    Teuchos::ParameterList& params) const
{
  // get the action required
  const auto act = Teuchos::getIntegralValue<FLD::BoundaryAction>(params, "action");

  switch (act)
  {
    case FLD::enforce_weak_dbc:
    case FLD::mixed_hybrid_dbc:
    case FLD::flow_dep_pressure_bc:
    case FLD::slip_supp_bc:
    case FLD::navier_slip_bc:
    case FLD::fpsi_coupling:
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      //       as input in the respective evaluate routines
      parent_element()->location_vector(dis, la, doDirichlet);
      break;
    case FLD::ba_none:
      FOUR_C_THROW("No action supplied");
      break;
    default:
      // standard case: element assembles into its own dofs only
      const int numnode = num_node();
      const Core::Nodes::Node* const* nodes = FluidXWallBoundary::nodes();

      la.clear();

      // we need to look at all DofSets of our discretization
      for (int dofset = 0; dofset < la.size(); ++dofset)
      {
        std::vector<int>& lm = la[dofset].lm_;
        std::vector<int>& lmdirich = la[dofset].lmdirich_;
        std::vector<int>& lmowner = la[dofset].lmowner_;
        std::vector<int>& lmstride = la[dofset].stride_;

        // fill the vector with nodal dofs
        if (nodes)
        {
          for (int i = 0; i < numnode; ++i)
          {
            const Core::Nodes::Node* node = nodes[i];

            const int owner = node->owner();
            std::vector<int> dofx;
            dis.dof(dofx, node, dofset, 0, this);
            std::vector<int> dof;
            // only take the first four dofs (the real dofs, but not the enriched ones)
            dof.push_back(dofx.at(0));
            dof.push_back(dofx.at(1));
            dof.push_back(dofx.at(2));
            dof.push_back(dofx.at(3));
            const int size = dof.size();
            if (size) lmstride.push_back(size);
            for (int j = 0; j < size; ++j)
            {
              lmowner.push_back(owner);
              lm.push_back(dof[j]);
            }

            if (doDirichlet)
            {
              const std::vector<int>* flag = nullptr;
              Core::Conditions::Condition* dirich = node->get_condition("Dirichlet");
              if (dirich)
              {
                if (dirich->type() != Core::Conditions::PointDirichlet &&
                    dirich->type() != Core::Conditions::LineDirichlet &&
                    dirich->type() != Core::Conditions::SurfaceDirichlet &&
                    dirich->type() != Core::Conditions::VolumeDirichlet)
                  FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
                flag = &dirich->parameters().get<std::vector<int>>("ONOFF");
              }
              for (unsigned j = 0; j < dof.size(); ++j)
              {
                if (flag && (*flag)[j])
                  lmdirich.push_back(1);
                else
                  lmdirich.push_back(0);
              }
            }
          }
        }

        // fill the vector with element dofs
        //      const int owner = Owner();
        std::vector<int> dofx = dis.dof(dofset, this);
        if (dofx.size()) FOUR_C_THROW("no element dofs expected");
        //      std::vector<int> dof;
        //      if(dofx.size())
        //      {
        //        //only take the first four dofs (the real dofs, but not the enriched ones)
        //        dof.push_back(dofx.at(0));
        //        dof.push_back(dofx.at(1));
        //        dof.push_back(dofx.at(2));
        //        dof.push_back(dofx.at(3));
        //      }

        //      if (dof.size()) lmstride.push_back(dof.size());
        //      for (unsigned j=0; j<dof.size(); ++j)
        //      {
        //        lmowner.push_back(owner);
        //        lm.push_back(dof[j]);
        //      }

        // fill the vector with face dofs
        if (this->num_dof_per_face(0) > 0)
          FOUR_C_THROW("set face_ from private to protected and uncomment");
        //      {
        //        for (int i=0; i<NumFace(); ++i)
        //        {
        //          const int owner = face_[i]->Owner();
        //          std::vector<int> dof = dis.Dof(dofset,face_[i]);
        //          if (dof.size())
        //            lmstride.push_back(dof.size());
        //          for (unsigned j=0; j<dof.size(); ++j)
        //          {
        //            lmowner.push_back(owner);
        //            lm.push_back(dof[j]);
        //          }
        //        }
        //      }

        //      if (doDirichlet)
        //      {
        //        const std::vector<int>* flag = nullptr;
        //        Core::Conditions::Condition* dirich = GetCondition("Dirichlet");
        //        if (dirich)
        //        {
        //          if (dirich->Type()!=Core::Conditions::geometry_type_pointDirichlet &&
        //              dirich->Type()!=Core::Conditions::geometry_type_lineDirichlet &&
        //              dirich->Type()!=Core::Conditions::geometry_type_surfaceDirichlet &&
        //              dirich->Type()!=Core::Conditions::geometry_type_volumeDirichlet)
        //            FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
        //          flag = dirich->get<std::vector<int> >("ONOFF");
        //        }
        //        for (unsigned j=0; j<dof.size(); ++j)
        //        {
        //          if (flag && (*flag)[j])
        //            lmdirich.push_back(1);
        //          else
        //            lmdirich.push_back(0);
        //        }
        //      }

      }  // for (int dofset=0; dofset<la.Size(); ++dofset)
      break;
  }


  return;
}

FOUR_C_NAMESPACE_CLOSE
