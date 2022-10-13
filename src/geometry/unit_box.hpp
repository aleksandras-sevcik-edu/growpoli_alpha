//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#ifndef unit_box_hpp
#define unit_box_hpp

#include <cstdint>
#include <string>
#include <ostream>
#include <algorithm>

#include "branched_polymer.hpp"
#include "coord_transform.hpp"

/* Class UnitBox implements unit box. */
class UnitBox
{
    std::vector<Polymer> _polymers;

    std::vector<Branched_polymer> _bran_polymers;

    std::string _box_label;

    float _box_origin_x;

    float _box_origin_y;

    float _box_origin_z;

    float _box_x;

    float _box_y;

    float _box_z;

    std::vector<std::vector<float>> _irradiation;

public:

    UnitBox();

    UnitBox(const std::vector<Polymer> &polymers,
            const std::vector<Branched_polymer> &bran_polymers,
            const std::string &label_box,
            const float box_origin_x,
            const float box_origin_y,
            const float box_origin_z,
            const float box_x,
            const float box_y,
            const float box_z,
            const std::vector<std::vector<float>> &irradiation);

    UnitBox(const UnitBox &source);

    UnitBox(UnitBox &&source) noexcept;

    ~UnitBox();

    UnitBox &operator=(const UnitBox &source);

    UnitBox &operator=(UnitBox &&source) noexcept;

    bool operator==(const UnitBox &other) const;

    bool operator!=(const UnitBox &other) const;

    void set_label(const std::string &box_label);

    void
    set_box_size(const float box_x,
                 const float box_y,
                 const float box_z);

    void add(const Polymer &polymer);

    int32_t number_of_polymers(const std::string &label) const;


    int32_t number_of_polymers() const;

    float origin_x() const;

    float origin_y() const;

    float origin_z() const;

    float size_x() const;

    float size_y() const;

    float size_z() const;


    std::string label() const;

    Polymer polymer(int n) const;

    friend std::ostream &operator<<(std::ostream &output, const UnitBox &source);

    bool check_polymer_intersection(Polymer &polymer);

    bool check_polymer_intersection_return(Polymer &polymer,
                                           int32_t &ip,
                                           int32_t &im);

    void remove_polymer(const Polymer &polymer);

    int find_polymer(const Polymer &polymer);

    bool find_nearest(int32_t polymer_index,
                      int32_t monomer_index,
                      float distance,
                      int32_t &nearest_polymer_index,
                      int32_t &nearest_monomer_index) const;

    bool check_boundaries(const Polymer &polymer) const;

    void save(std::string box_name);

    void load(std::string box_name);

    void irradiate(float dose_value, int32_t split_size);

    std::vector<std::vector<float>> &get_irradiation();

    void change_irradiation(int32_t index, float dose_value);

    float get_irradiation_max_value();

    void print_branches(int32_t polymer_index);

    void print_irradiation();

    void calculate_branching_alternative(std::vector<float> &average_number, std::vector<float> &average_length,
                                         std::vector<float> &correlation);

    void remove_branpol(const Branched_polymer &branpolymer);

    void add(std::vector<Polymer> &polymers);

    bool add_branpol(Branched_polymer &branpol);

    bool find_branpol(Branched_polymer &branpol) const;

    Branched_polymer branpol(int32_t index) const;

    bool check_branpol_intersection(Branched_polymer &branpol);

    int32_t number_of_polymers(const unsigned long long int bid) const;

    bool find_nearest_for_radical(float x, float y, float z, float distance, int32_t &nearest_polymer_index,
                                  int32_t &nearest_monomer_index) const;

    bool is_intersecting_between_skip_same_branches();

    bool check_group_intersection(std::vector<Polymer> &polymers);

    std::vector<Branched_polymer> get_branpols() const;

    std::vector<Polymer> get_polymers() const;

    bool find_nearest_skip_same_branch(int32_t polymer_index, int32_t monomer_index, float distance,
                                       int32_t &nearest_polymer_index, int32_t &nearest_monomer_index) const;

    void remove_polymer(int32_t index);


    Polymer polymer(const std::string &label) const;

    int32_t number_of_monomers() const;

    Polymer *polymer_ptr(int n);

    Branched_polymer find_branpol(unsigned long long int bid) const;

    int32_t number_of_branched_polymers() const;

    Branched_polymer *branpol_ptr(int n);

    float get_irradiation_min_value();

    void get_irradiation(int32_t index, float &x, float &y, float &z, float &value) const;
};

#endif /* UnitBox_hpp */
