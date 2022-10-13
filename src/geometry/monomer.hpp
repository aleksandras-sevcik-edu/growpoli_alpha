//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma  once

#include <cstdint>
#include <string>
#include <ostream>
#include <vector>

#include "rect_shape_3d.hpp"


/* Class Monomer implements monomer. */

class Monomer
{
    std::string _label_m;

    Rectangle3D _geometry;

    std::vector<std::vector<double>> _coords;

    float _center_x;

    float _center_y;

    float _center_z;

public:

    Monomer();

    Monomer(const std::string &label_m,
            const Rectangle3D &geometry,
            const std::vector<std::vector<double>> &coords,
            const float center_x,
            const float center_y,
            const float center_z);

    Monomer(const Monomer &source);

    Monomer(Monomer &&source) noexcept;

    ~Monomer();

    Monomer &operator=(const Monomer &source);

    Monomer &operator=(Monomer &&source) noexcept;

    bool operator==(const Monomer &other) const;

    bool operator!=(const Monomer &other) const;

    void set_label(const std::string &label_m);

    void set_length(const float length);

    void set_width(const float width_m);

    void set_height(const float height_m);

    std::string label() const;

    float length() const;

    float width() const;

    float height() const;

    float getSphereRadius() const;

    void set_coords (std::vector<std::vector<double>>& coords);

    std::vector<std::vector<double>>
    get_coords() const;

    friend std::ostream &operator<<(std::ostream &output, const Monomer &source);

    const double*
    get_coord(int32_t index) const;

    void set_center(float center_x, const float center_y, const float center_z);

    float get_center_x() const;

    float get_center_y() const;

    float get_center_z() const;
};

