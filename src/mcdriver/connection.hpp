//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include "unit_box.hpp"
#include "random_generator.hpp"

/* Implements combination event */

class Connection
{
public:

    ExponentialMap _align_bound(ExponentialMap &previous_rot, ExponentialMap &align_to, float delta) const;

    bool mono_branch(UnitBox &box, Polymer &activated, int32_t mon_index, Polymer &target, int32_t targ_moni);

    bool _is_intersecting_from_second_mon(Polymer &activated, Polymer &target);

    bool perform_scenario(UnitBox &box, std::vector<std::vector<int>> &con_stats, Polymer &activated, int32_t act_moni,
                          Branched_polymer &act_branpol, Polymer &target, int32_t targ_moni,
                          Branched_polymer &targ_branpol,
                          std::vector<Polymer> &pol_output, std::vector<Branched_polymer> &branpol_output);

    bool double_branch(UnitBox &box, Polymer &activated, int32_t act_mon, Polymer &target, int32_t targ_mon,
                       std::vector<Polymer> &result);

    bool _is_intersecting_from_first_mon(Polymer &activated, Polymer &target);

    bool _is_centers_too_close(Polymer &activated, int32_t amoni, Polymer &target);

    bool merge(UnitBox &box, Polymer &activated, Polymer &target, bool is_branched);

    bool branch_merge(UnitBox &box, Polymer &activated, int32_t act_mon, Branched_polymer &act_branpol, Polymer &target,
                      int32_t targ_mon, Branched_polymer &targ_branpol, Branched_polymer &result);

    bool branch_swap(UnitBox &box, Polymer &activated, int32_t act_mon, Branched_polymer &act_branpol, Polymer &target,
                     int32_t targ_mon, Branched_polymer &targ_branpol, std::vector<Branched_polymer> &bresults,
                     std::vector<Polymer> &presults);

    int _find_position(const Polymer &polymer, std::vector<Polymer> &polymers);

    Polymer _find_polymer(const std::string &label, std::vector<Polymer> &polymers);

    void _move_polymers(Polymer &target, int32_t act_mon, float fx, float fy, float fz, std::vector<Polymer> &polymers);
};
