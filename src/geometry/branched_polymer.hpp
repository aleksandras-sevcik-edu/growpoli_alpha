//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include <vector>
#include <optional>
#include <unordered_map>
#include <map>

#include "polymer.hpp"

/* class implements branched polymer structure */

unsigned long long static bid_counter = 1;

class Branched_polymer
{
public:
    unsigned long long _identifier;
    Polymer _parent;
    std::vector<Polymer> _branches;

    Branched_polymer();
    Branched_polymer(Polymer &parent, std::vector<Polymer> &branches);
    Branched_polymer(Polymer &parent, Polymer &branch);
    ~Branched_polymer();
    Branched_polymer(const Branched_polymer& source);
    Branched_polymer& operator=(const Branched_polymer& source);
    Branched_polymer(Branched_polymer&& source) noexcept;
    bool operator==(const Branched_polymer& other) const;
    bool operator!=(const Branched_polymer& other) const;
    friend std::ostream& operator<<(std::ostream& output, const Branched_polymer& source);

    unsigned long long _get_new_id();
    size_t get_id() const;
    Polymer get_parent() const;
    std::vector<Polymer> get_branches() const;
    void set_parent(Polymer &parent);
    void _set_branches(std::vector<Polymer> &branches);

    bool is_intersecting_with (std::vector<Polymer> &polymers);

    bool calculate_branches(std::vector<Polymer> &input, std::vector<Polymer> &branches,
                            std::vector<int32_t> &moni);

    bool calculate_all_branches(std::vector<Polymer> &input, std::vector<Polymer> &output);

    static void _move_element_to_back(Polymer &polymer, std::vector<Polymer> &polymers);

    size_t _find_position(const Polymer &polymer);

    bool make_division(Polymer &activated, int32_t moni, std::vector<Branched_polymer> &new_branpols,
                       std::vector<Polymer> &new_polymers);

    bool merge(Branched_polymer &other);

    void from_polymers(std::vector<Polymer> &polymers);

    std::vector<Polymer> to_polymers() const;

    void calculate_coords();

    bool division_validation_check(std::vector<Branched_polymer> &new_branpols, std::vector<Polymer> &new_polymers);

    bool add_branches(Polymer &activated, int32_t moni, std::vector<Polymer> &polymers);

    bool add_branch(Polymer &activated, int32_t moni, Polymer &polymer);

    bool is_intersecting_between() const;

    bool change(Polymer &old_polymer, Polymer &new_polymer);

    int32_t number_of_polymers() const;

    void make_longest_parent();

    bool _add_unique(std::vector<Polymer> &polymers, Polymer &polymer);

    static void _remove_same_element_from_other_vec(std::vector<Polymer> &a, std::vector<Polymer> &b);

    bool is_monomer_occupied(Polymer &activated, int32_t moni) const;

    void _add_uniques(std::vector<Polymer> &polymers, std::vector<Polymer> &others);

    void _move(float fx, float fy, float fz);

    bool make_division_for_swap(Polymer &activated, int32_t moni, std::vector<Branched_polymer> &new_branpols,
                                std::vector<Polymer> &new_polymers, std::vector<Polymer> &specials);

    bool calculate_all_branches_for_stats(std::vector<int> &number, std::vector<float> &length);

    bool _check_if_contains_already(std::vector<Polymer> &polymers, std::vector<Polymer> &others) const;
};
