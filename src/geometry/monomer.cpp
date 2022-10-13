//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include <algorithm>
#include "monomer.hpp"

Monomer::Monomer()

        :_label_m(std::string()),

        _geometry(Rectangle3D(0.0f, 0.0f, 0.0f)),

        _coords(std::vector<std::vector<double>>()),

        _center_x(0.0f),

        _center_y(0.0f),

        _center_z(0.0f)
{

}

Monomer::Monomer(const std::string &label_m, const Rectangle3D &geometry, const std::vector<std::vector<double>> &coords,
                 const float center_x, const float center_y, const float center_z)

        : _label_m(label_m), _geometry(geometry), _coords(coords), _center_x(center_x), _center_y(center_y), _center_z(center_z)
{
}

Monomer::Monomer(const Monomer &source)
        : _label_m(source._label_m),

        _geometry(source._geometry),

        _coords(source._coords),

        _center_x(source._center_x),

        _center_y(source._center_y),

        _center_z(source._center_z)
{
}

Monomer::Monomer(Monomer &&source) noexcept

        : _label_m(std::move(source._label_m)),

        _geometry(std::move(source._geometry)),

        _coords(std::move(source._coords)),

        _center_x(std::move(source._center_x)),

        _center_y(std::move(source._center_y)),

        _center_z(std::move(source._center_z))
{
}

Monomer::~Monomer()
{
}

Monomer &
Monomer::operator=(const Monomer &source)
{
    if (this == &source) return *this;

    _label_m = source._label_m;

    _geometry = source._geometry;

    _coords = source._coords;

    _center_x = source._center_x;

    _center_y = source._center_y;

    _center_z = source._center_z;

    return *this;
}

Monomer &
Monomer::operator=(Monomer &&source) noexcept
{
    if (this == &source) return *this;

    _label_m = std::move(source._label_m);

    _geometry = std::move(source._geometry);

    _coords = std::move(source._coords);

    _center_x = std::move(source._center_x);

    _center_y = std::move(source._center_y);

    _center_z = std::move(source._center_z);

    return *this;
}

bool
Monomer::operator==(const Monomer &other) const
{
    if (_label_m != other._label_m) return false;

    if (_geometry != other._geometry) return false;

    if (!std::equal(_coords.begin(), _coords.end(),
                    other._coords.begin(), other._coords.end()))
    {
        return false;
    }

    if (_center_x != other._center_x) return false;
    if (_center_y != other._center_y) return false;
    if (_center_z != other._center_z) return false;

    return true;
}

bool
Monomer::operator!=(const Monomer &other) const
{
    return !(*this == other);
}

void
Monomer::set_label(const std::string &label_m)
{
    _label_m = label_m;
}

void
Monomer::set_length(const float length)
{
    _geometry.length(length);
}

void
Monomer::set_width(const float width)
{
    _geometry.width(width);
}

void
Monomer::set_height(const float height)
{
    _geometry.height(height);
}

void
Monomer::set_center(const float center_x, const float center_y, const float center_z)
{
    _center_x = center_x;
    _center_y = center_y;
    _center_z = center_z;
}

float
Monomer::get_center_x() const
{
    return _center_x;
}

float
Monomer::get_center_y() const
{
    return _center_y;
}

float
Monomer::get_center_z() const
{
    return _center_z;
}

std::string
Monomer::label() const
{
    return _label_m;
}

float
Monomer::length() const
{
    return _geometry.length();
}

float
Monomer::width() const
{
    return _geometry.width();
}

float
Monomer::height() const
{
    return _geometry.height();
}

float
Monomer::getSphereRadius() const
{
    return _geometry.getSphereRadius();
}

std::ostream &
operator<<(std::ostream &output, const Monomer &source)
{
    output << std::endl;

    output << "[Monomer (Instance): " << &source << "]" << std::endl;

    output << "this->_label_m: " << source._label_m << std::endl;

    output << "this->_geometry: " << source._geometry << std::endl;

    output << "this->_coords(vector): " << std::endl;

    std::for_each(source._coords.begin(), source._coords.end(),
                  [&output](auto &tval)
                  {output << tval[0] << tval[1] << tval[2] << std::endl;});
    return output;
}

void
Monomer::set_coords(std::vector<std::vector<double>> &coords)
{
    _coords = coords;
}

std::vector<std::vector<double>>
Monomer::get_coords() const
{
    return _coords;
}


const double*
Monomer::get_coord(int32_t index) const
{
    return _coords[index].data();
}
