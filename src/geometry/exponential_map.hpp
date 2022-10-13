//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include <cstdint>
#include <string>
#include <ostream>
#include <tuple>
#include <array>

/* Class exponential map implements axis-angle representation of rotation. */

class ExponentialMap
{
    float _coord_x;

    float _coord_y;

    float _coord_z;

    float _theta;

    void _normalize();

public:

    ExponentialMap();

    ExponentialMap(float coord_x,
                   float coord_y,
                   float coord_z,
                   float theta);

    ExponentialMap(const ExponentialMap& source);

    ExponentialMap(ExponentialMap&& source) noexcept;

    ~ExponentialMap();

    ExponentialMap& operator=(const ExponentialMap& source);

    ExponentialMap& operator=(ExponentialMap&& source) noexcept;

    bool operator==(const ExponentialMap& other) const;

    bool operator!=(const ExponentialMap& other) const;

    ExponentialMap&
    operator+=(const ExponentialMap &other);

    void
    rotate(float& x,
           float& y,
           float& z,
           float rx,
           float ry,
           float rz) const;

    void
    get_rotation(float &x,
                 float &y,
                 float &z,
                 float &teta) const;

    friend std::ostream& operator<<(      std::ostream&   output,
                                    const ExponentialMap& source);

    void to_quaternion(float &x, float &y, float &z, float &w) const;

    void from_quaternion(float &xcoord, float &ycoord, float &zcoord,
                         float &angle, float x, float y, float z, float w) const;

    void change_angle(float teta);

    float get_angle() const;
};
