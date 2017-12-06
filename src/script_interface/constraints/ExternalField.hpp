/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SCRIPT_INTERFACE_CONSTRAINTS_EXTERNALPOTENTIAL_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_EXTERNALPOTENTIAL_HPP

#include "core/constraints/Constraint.hpp"
#include "core/constraints/ExternalField.hpp"

namespace ScriptInterface {
namespace Constraints {

class ExternalField : public Constraint {
public:
  ExternalField()
    : m_constraint(new ::Constraints::ExternalField()) {
    add_parameters({
                        {"field", [this](Variant const &v) { m_constraint->set_field( ??? get_value<Vector3d>(v)); },
                                  [this]()                 { return ??? m_constraint->potential(); }}
                   });
  }

  const std::string name() const override { return "Constraints::ExternalField"; }

  std::shared_ptr<::Constraints::Constraint> constraint() override {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<const ::Constraints::Constraint> constraint() const override {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<::Constraints::ExternalField> external_field() const {
    return m_constraint;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::ExternalField> m_constraint;

};

} /* namespace Constraints */
} /* namespace ScriptInterface */

#endif
