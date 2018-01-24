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
#include <boost/multi_array.hpp>
#include "integrate.hpp" //for skin
#include "grid.hpp"      //for box_l

namespace ScriptInterface {
namespace Constraints {

class ExternalField : public Constraint {
public:
  ExternalField() {
    add_parameters({
                        {"field", [this](Variant const &v) {
                            //std::cout << print_variant_types(v) << std::endl;
                            //field = (bins, data, order)
                            auto field = get_value<std::vector<Variant>>(v);
                            auto shape = get_value<std::vector<int>>(field[0]);
                            auto data = get_value<std::vector<double>>(field[1]);
                            auto order = get_value<int>(field[2]); 
                            //Increase halo to cover order + skin
                            //order = 0: ?
                            auto halo_x = order / 2 + int(skin / box_l[0]*shape[0]); 
                            auto halo_y = order / 2 + int(skin / box_l[1]*shape[1]); 
                            auto halo_z = order / 2 + int(skin / box_l[2]*shape[2]); 
                            auto array = boost::multi_array_ref<Vector3d, 3>(reinterpret_cast<Vector3d *>(data.data()), shape);

                            m_constraint = std::make_shared<::Constraints::ExternalField>(array, std::array<int, 3>{halo_x, halo_y, halo_z}); 
                            //m_constraint = std::make_shared<::Constraints::ExternalField<Vector3d,3>>(array, {halo, halo, halo}); 
                        }, [this]()                 { return std::vector<Variant>{}; }},
                        {"particle_weights", [this](Variant const &v) {
                                  std::unordered_map<int, double> weights;
                                  auto pairs = get_value<std::vector<Variant>>(v);
                                    for(auto const& e: pairs) {
                                        auto pair = get_value<std::vector<Variant>>(e);
                                        weights[get_value<int>(pair.front())] = get_value<double>(pair.back());
                                    }
                                    m_constraint->set_weights(weights);
                        }, []() {        return std::vector<Variant>{};   }           }
                                  
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
