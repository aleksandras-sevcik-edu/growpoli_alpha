//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#pragma once

#include <cstdint>
#include <string>
#include <ostream>
#include <array>
#include <tuple>
#include <vector>

/**
 Class Rectangle3D implements rectangle shape in 3D.
 */
class Rectangle3D
{
    float _cuboid_length;

    float _cuboid_width;

    float _cuboid_height;

public:

    Rectangle3D();

    Rectangle3D(const float cuboid_length,
                const float cuboid_width,
                const float cuboid_height);

    Rectangle3D(const Rectangle3D &source);

    Rectangle3D(Rectangle3D &&source) noexcept;

    ~Rectangle3D();

    Rectangle3D &operator=(const Rectangle3D &source);

    Rectangle3D &operator=(Rectangle3D &&source) noexcept;

    bool operator==(const Rectangle3D &other) const;

    bool operator!=(const Rectangle3D &other) const;

    void length(const float length);

    void width(const float width);

    void height(const float height);

    float length() const;

    float width() const;

    float height() const;

    std::tuple<float, float, float> getGeometricalCenter() const;

    float getSphereRadius() const;

    std::tuple<float, float, float> getReferencePoint() const;

    friend std::ostream &operator<<(std::ostream &output,
                                    const Rectangle3D &source);
};
