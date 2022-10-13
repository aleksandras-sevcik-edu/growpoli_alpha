//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "exponential_map.hpp"

#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include "math_operations.hpp"


ExponentialMap::ExponentialMap()

    : _coord_x(0.0f)

    , _coord_y(0.0f)

    , _coord_z(0.0f)

    , _theta(0.0f)
{
    
}

ExponentialMap::ExponentialMap(const float coord_x,
                               const float coord_y,
                               const float coord_z,
                               const float theta)
    : _coord_x(coord_x)

    , _coord_y(coord_y)

    , _coord_z(coord_z)

    , _theta(theta)
{
    _normalize();
}

ExponentialMap::ExponentialMap(const ExponentialMap& source)

    : _coord_x(source._coord_x)

    , _coord_y(source._coord_y)

    , _coord_z(source._coord_z)

    , _theta(source._theta)
{
    
}

ExponentialMap::ExponentialMap(ExponentialMap&& source) noexcept

    : _coord_x(std::move(source._coord_x))

    , _coord_y(std::move(source._coord_y))

    , _coord_z(std::move(source._coord_z))

    , _theta(std::move(source._theta))
{
    
}

ExponentialMap::~ExponentialMap()
{
    
}

ExponentialMap&
ExponentialMap::operator=(const ExponentialMap& source)
{
    if (this == &source) return *this;
    
    _coord_x = source._coord_x;
    
    _coord_y = source._coord_y;
    
    _coord_z = source._coord_z;
    
    _theta = source._theta;
    
    return *this;
}

ExponentialMap&
ExponentialMap::operator=(ExponentialMap&& source) noexcept
{
    if (this == &source) return *this;
    
    _coord_x = std::move(source._coord_x);
    
    _coord_y = std::move(source._coord_y);
    
    _coord_z = std::move(source._coord_z);
    
    _theta = std::move(source._theta);
    
    return *this;
}

bool
ExponentialMap::operator==(const ExponentialMap& other) const
{
    if (std::fabs(_coord_x - other._coord_x) > 1.0e-6) return false;
    
    if (std::fabs(_coord_y - other._coord_y) > 1.0e-6) return false;
    
    if (std::fabs(_coord_z - other._coord_z) > 1.0e-6) return false;
    
    if (std::fabs(_theta - other._theta) > 1.0e-6) return false;
    
    return true;
}

bool
ExponentialMap::operator!=(const ExponentialMap& other) const
{
    return !(*this == other);
}

ExponentialMap&
ExponentialMap::operator+=(const ExponentialMap& other)
{
    if (_theta > getPiValue())
    {
        _theta = 2 * getPiValue() - _theta;

        _coord_x = _coord_x * -1;

        _coord_y = _coord_y * -1;

        _coord_z = _coord_z * -1;
    }
    float x0, y0, z0, w0;

    to_quaternion(x0, y0, z0, w0);

    float x1, y1, z1, w1;

    other.to_quaternion(x1, y1, z1, w1);

    auto w2 = w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1;

    auto x2 = w0 * x1 + x0 * w1 + y0 * z1 - z0 * y1;

    auto y2 = w0 * y1 + y0 * w1 + z0 * x1 - x0 * z1;

    auto z2 = w0 * z1 + z0 * w1 + x0 * y1 - y0 * x1;

    from_quaternion(_coord_x, _coord_y, _coord_z, _theta, x2, y2, z2, w2);

    return *this;
}


void
ExponentialMap::_normalize()
{
    auto fnorm = 1.0f / std::sqrt(_coord_x * _coord_x + _coord_y * _coord_y + _coord_z * _coord_z);
    
    _coord_x *= fnorm;
    
    _coord_y *= fnorm;
    
    _coord_z *= fnorm;
}

void
ExponentialMap::to_quaternion (float& x,
                               float& y,
                               float& z,
                               float& w) const
{
    float s = std::sin(_theta / 2);

    x = _coord_x * s;

    y = _coord_y * s;

    z = _coord_z * s;

    w = std::cos(_theta / 2);
}

void
ExponentialMap::from_quaternion (float& xcoord,
                                 float& ycoord,
                                 float& zcoord,
                                 float& angle,
                                 float x,
                                 float y,
                                 float z,
                                 float w) const
{
    const auto m_pi = getPiValue();

    angle = acosf(w);

    float sinz = sinf(angle);

    if (fabsf(sinz) > 1e-4f) {
        sinz = 1.0f / sinz;

        xcoord = x * sinz;
        ycoord = y * sinz;
        zcoord = z * sinz;

        angle *= 2.0f;
        if (angle > m_pi)
            angle = 2 * m_pi - angle;
    } else {
        angle = 0.0f;
        xcoord = 1.0f;
        ycoord = 0.0f;
        zcoord = 0.0f;
    }
}


void
ExponentialMap::rotate (float& x,
                       float& y,
                       float& z,
                       const float rx,
                       const float ry,
                       const float rz) const
{
    const float tsin = std::sin(_theta);

    const float tcos = std::cos(_theta);

    x = rx * tcos;

    y = ry * tcos;

    z = rz * tcos;

    const float vx = _coord_y * rz - _coord_z * ry;

    const float vy = _coord_z * rx - _coord_x * rz;

    const float vz = _coord_x * ry  - _coord_y * rx;

    x += tsin * vx;

    y += tsin * vy;

    z += tsin * vz;

    const float fp = _coord_x * rx + _coord_y * ry + _coord_z * rz;

    const float tfact = (1.0f - tcos) * fp;

    x += tfact * _coord_x;

    y += tfact * _coord_y;

    z += tfact * _coord_z;
}





void
ExponentialMap::get_rotation (float& x,
                              float& y,
                              float& z,
                              float& teta) const
{
    x = _coord_x;

    y = _coord_y;

    z = _coord_z;

    teta = _theta;
}

void
ExponentialMap::change_angle (float teta)
{
    _theta = teta;
}


float
ExponentialMap::get_angle () const
{
    return _theta;
}


std::ostream &
operator<<(std::ostream &output, const ExponentialMap &source)
{
    output << "[Exponential Map (Instance): " << &source << "] : ";
    output << source._coord_x << " ";
    output << source._coord_y << " ";
    output << source._coord_z << " ";
    output << source._theta << std::endl;
    return output;
}


