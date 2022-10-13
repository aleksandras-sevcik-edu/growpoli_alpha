//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "division.hpp"

bool Division::perform_polymer_division_rotransl(Polymer &polymer, int32_t moni, float transl_dist,
                                                 std::vector<Polymer> &output)
{
    if (moni == 0) return false;

    Polymer first, second;
    polymer.divide(moni, first, second);

    Rotransl r;
    r.rotate_translate(first, transl_dist);
    r.rotate_translate(second, transl_dist);

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();
    if (!first.polymer_intersecting(second, &bd1, &bd2))
    {
        output.push_back(first);
        output.push_back(second);
        free_bd(&bd1);
        free_bd(&bd2);
        return true;
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}


bool Division::perform_branpol_division_rotransl(Branched_polymer &branpol, Polymer &polymer, int32_t moni,
                                                 float transl_dist,
                                                 std::vector<Branched_polymer> &output, std::vector<Polymer> &output_p)
{
    if (moni == 0) return false;

    if (branpol.make_division(polymer, moni, output, output_p))
    {
        for (auto &b: output) b.make_longest_parent();

        Rotransl r;
        for (auto &b: output) r.rotate_translate(b, transl_dist);
        for (auto &p: output_p) r.rotate_translate(p, transl_dist);

        for (auto &p: output_p) p.set_bid(0);

        for (auto &out: output)
        {
            auto pol = out.get_branches();
            pol.push_back(out.get_parent());
            std::cout << std::endl;
        }


        if (!is_intersecting(output, output_p)) return true;
    }
    return false;
}


bool Division::is_intersecting(std::vector<Branched_polymer> &output, std::vector<Polymer> &output_p)
{

    UnitBox check_box;
    for (auto &b: output) check_box.add_branpol(b);
    for (auto &p: output_p) check_box.add(p);

    if (check_box.is_intersecting_between_skip_same_branches())
    {
        return true;
    }

    for (auto &out: output)
    {
        if (out.is_intersecting_between())
        {
            return true;
        }
    }
    return false;
}
