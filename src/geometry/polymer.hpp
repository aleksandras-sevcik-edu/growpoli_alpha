//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <ostream>
#include <algorithm>

#include "monomer.hpp"
#include "exponential_map.hpp"
#include "gjk_intersection.hpp"
#include "poly_rotation.h"

/* class implements polymer object and associated functions */

unsigned long long static polymer_counter = 1;

enum class PolymerType
{
    valid, invalid
};

class Polymer
{
    std::vector<Monomer> _monomers;

    std::string _label_p;

    PolymerType _type;

    float _origin_x;

    float _origin_y;

    float _origin_z;

    std::vector<ExponentialMap> _rotations;

    unsigned long long _bid;

public:

    Polymer();

    Polymer(const std::vector<Monomer> &monomers,
            const std::string &label_p,
            const PolymerType &type,
            const float origin_x,
            const float origin_y,
            const float origin_z,
            const std::vector<ExponentialMap> &rotations,
            const unsigned long long bid);

    Polymer(const Polymer &source);

    Polymer(Polymer &&source) noexcept;

    ~Polymer();

    Polymer &operator=(const Polymer &source);

    Polymer &operator=(Polymer &&source) noexcept;

    bool operator==(const Polymer &other) const;

    bool operator!=(const Polymer &other) const;

    std::string _add_label();

    void add(Monomer &monomer,
             const ExponentialMap &rotation);

    void position(float &rx,
                  float &ry,
                  float &rz,
                  int32_t index) const;

    void position_center(float &rx,
                         float &ry,
                         float &rz,
                         int32_t index) const;

    bool intersecting_gjk(const Monomer &monom,
                          const ExponentialMap &rotation,
                          int32_t index,
                          bd *bd1,
                          bd *bd2) const;

    bool polymer_intersecting(Polymer &other) const;

    bool polymer_intersecting_connect(const Polymer &other) const;

    bool polymer_intersecting_from(const Polymer &other,
                                   const int32_t index_from) const;

    bool polymer_intersecting_return(Polymer &other, bd *bd1, bd *bd2, int32_t &im) const;

    void set_label(std::string label);

    std::string get_label() const;

    void set_type(PolymerType type);

    PolymerType get_type() const;

    int32_t number_of_monomers() const;

    void set_origins(const float originx,
                     const float originy,
                     const float originz);

    void set_rotation(const float rx,
                      const float ry,
                      const float rz,
                      const float teta,
                      const int32_t index);

    void get_rotation(float &rx,
                      float &ry,
                      float &rz,
                      float &teta,
                      int32_t index) const;

    ExponentialMap rotation(int32_t index) const;

    friend std::ostream &operator<<(std::ostream &output, const Polymer &source);

    void polymer_merge(Polymer &other);

    Monomer monomer(int n) const;

    void divide(int32_t monomer_index, Polymer &polyone, Polymer &polytwo);

    void add_for_saving(Monomer &monomer, ExponentialMap &exponentialMap);

    bool intersecting_itself();

    void fill_coords(int32_t index);

    void set_mcoord(std::vector<std::vector<double>> mcoord, int32_t index);

    std::vector<std::vector<double>> get_mcoord(int32_t index) const;

    void get_mcoord(int32_t index, bd *bd) const;

    bool polymer_intersecting(Polymer &other, bd *bd1, bd *bd2) const;

    void set_rotation(const ExponentialMap &new_rot, const int32_t index);

    void rotate_polymer(const Pt &centre, const ExponentialMap &newr);

    void get_mass_center(float &mx, float &my, float &mz) const;

    void reverse();

    void calculate_coords();

    bool polymer_intersecting_as_branch(Polymer &other, bd *bd1, bd *bd2) const;

    bool is_valid() const;

    bool is_branched() const;

    void shorten(int32_t nr);

    bool polymer_intersecting_from_second(Polymer &other, bd *bd1, bd *bd2) const;

    bool polymer_intersecting_return_branch(Polymer &other, struct bd *bd1, struct bd *bd2,
                                            int32_t &mon) const;
    void set_bid(unsigned long long bid);

    unsigned long long get_bid() const;

    bool has_coords() const;

    bool polymer_intersecting_return_connection(Polymer &other, bd *bd1, bd *bd2, int32_t &moni, int32_t &omoni) const;

    bool polymer_check_if_branch(Polymer &other, bd *bd1, bd *bd2) const;

    bool polymer_monomer_intersecting(int32_t moni, Polymer &other, bd *bd1, bd *bd2) const;

    void position_center_alternative(float &rx, float &ry, float &rz, const int32_t index) const;

    void get_coords_to_vec(float *ptr_coord, int index) const;

    void calculate_monomer_center(int32_t index, float x, float y, float z);

    void position_center_alternative_2(float &rx, float &ry, float &rz, const int32_t index) const;
};

