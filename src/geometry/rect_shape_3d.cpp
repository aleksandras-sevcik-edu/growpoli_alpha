//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "rect_shape_3d.hpp"
#include <cmath>
#include <fstream>
#include <algorithm>


Rectangle3D::Rectangle3D()

    : _cuboid_length(0.0f)

    , _cuboid_width(0.0f)

    , _cuboid_height(0.0f)
{

}

Rectangle3D::Rectangle3D(const float cuboid_length,
                         const float cuboid_width,
                         const float cuboid_height)

    : _cuboid_length(cuboid_length)

    , _cuboid_width(cuboid_width)

    , _cuboid_height(cuboid_height)
{
}

Rectangle3D::Rectangle3D(const Rectangle3D& source)

    : _cuboid_length(source._cuboid_length)

    , _cuboid_width(source._cuboid_width)

    , _cuboid_height(source._cuboid_height)
{

}

Rectangle3D::Rectangle3D(Rectangle3D&& source) noexcept

    : _cuboid_length(std::move(source._cuboid_length))

    , _cuboid_width(std::move(source._cuboid_width))

    , _cuboid_height(std::move(source._cuboid_height))
{

}

Rectangle3D::~Rectangle3D()
{

}

Rectangle3D&
Rectangle3D::operator=(const Rectangle3D& source)
{
    if (this == &source) return *this;

    _cuboid_length = source._cuboid_length;

    _cuboid_width = source._cuboid_width;

    _cuboid_height = source._cuboid_height;

    return *this;
}

Rectangle3D&
Rectangle3D::operator=(Rectangle3D&& source) noexcept
{
    if (this == &source) return *this;

    _cuboid_length = std::move(source._cuboid_length);

    _cuboid_width = std::move(source._cuboid_width);

    _cuboid_height = std::move(source._cuboid_height);

    return *this;
}

bool
Rectangle3D::operator==(const Rectangle3D& other) const
{

    if (std::fabs(_cuboid_length - other._cuboid_length) > 1.0e-6) return false;

    if (std::fabs(_cuboid_width - other._cuboid_width) > 1.0e-6) return false;

    if (std::fabs(_cuboid_height - other._cuboid_height) > 1.0e-6) return false;

    return true;
}

bool
Rectangle3D::operator!=(const Rectangle3D& other) const
{
    return !(*this == other);
}

void
Rectangle3D::length(const float length)
{
    _cuboid_length = length;
}

void
Rectangle3D::width(const float width)
{
    _cuboid_width = width;
}

void
Rectangle3D::height(const float height)
{
    _cuboid_height = height;
}

float
Rectangle3D::length() const
{
    return _cuboid_length;
}

float
Rectangle3D::width() const
{
    return _cuboid_width;
}

float
Rectangle3D::height() const
{
    return _cuboid_height;
}

std::tuple<float, float, float>
Rectangle3D::getGeometricalCenter() const
{
   return std::make_tuple(0.5 * _cuboid_length, 0.5 * _cuboid_width, 0.5 * _cuboid_height);
}

float
Rectangle3D::getSphereRadius() const
{
    const auto ccords = Rectangle3D::getGeometricalCenter();

    const auto rx = std::get<0>(ccords);

    const auto ry = std::get<1>(ccords);

    const auto rz = std::get<2>(ccords);

    return std::sqrt(rx * rx + ry * ry + rz * rz);
}

std::tuple<float, float, float>
Rectangle3D::getReferencePoint() const
{
    return std::make_tuple(0.5f * _cuboid_length, 0.5f * _cuboid_width, 0.0f);
}

std::ostream&
operator<<(std::ostream& output, const Rectangle3D& source)
{
    output << std::endl;

    output << "[Rectangle3D (Instance): " << &source << "]" << std::endl;

    output << "this->_cuboid_length: " << source._cuboid_length << std::endl;

    output << "this->_cuboid_width: " << source._cuboid_width << std::endl;

    output << "this->_cuboid_height: " << source._cuboid_height << std::endl;

    return output;
}



