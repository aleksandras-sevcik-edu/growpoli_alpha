//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "polymer.hpp"
#include "math_operations.hpp"

#include <omp.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>


Polymer::Polymer()

    : _monomers(std::vector<Monomer>())
	
	, _label_p(std::string())

	, _type(PolymerType::valid)
	
	, _origin_x(0.0f)

    , _origin_y(0.0f)

    , _origin_z(0.0f)
	
	, _rotations(std::vector<ExponentialMap>())

    , _bid (0)
{
    _label_p = _add_label();
}

Polymer::Polymer(const std::vector<Monomer> &monomers,
                 const std::string &label_p,
                 const PolymerType &type,
                 const float origin_x,
                 const float origin_y,
                 const float origin_z,
                 const std::vector<ExponentialMap> &rotations,
                 const unsigned long long bid)

    : _monomers(monomers)
	
	, _label_p(label_p)

    , _type(type)
	
	, _origin_x(origin_x)

    , _origin_y(origin_y)

    , _origin_z(origin_z)
	
	, _rotations(rotations)

    , _bid(bid)
{
    
}

Polymer::Polymer(const Polymer& source)

    : _monomers(source._monomers)
	
	, _label_p(source._label_p)

    , _type(source._type)
	
	, _origin_x(source._origin_x)

    , _origin_y(source._origin_y)

    , _origin_z(source._origin_z)
	
	, _rotations(source._rotations)

    , _bid(source._bid)
{
    
}

Polymer::Polymer(Polymer&& source) noexcept

    : _monomers(std::move(source._monomers))
	
    , _label_p(std::move(source._label_p))

    , _type(std::move(source._type))
	
    , _origin_x(std::move(source._origin_x))

    , _origin_y(std::move(source._origin_y))

    , _origin_z(std::move(source._origin_z))
	
	, _rotations(std::move(source._rotations))

    , _bid(std::move(source._bid))
{
    
}

Polymer::~Polymer()
{
    
}

Polymer&
Polymer::operator=(const Polymer& source)
{
    if (this == &source) return *this; 
    
	_monomers = source._monomers;	
	
	_label_p = source._label_p;

    _type = source._type;
	
	_origin_x = source._origin_x;

    _origin_y = source._origin_y;
    
    _origin_z = source._origin_z;
    	
	_rotations = source._rotations;

    _bid = source._bid;

    return *this;
}

Polymer&
Polymer::operator=(Polymer&& source) noexcept
{
    if (this == &source) return *this;
    
	_monomers = std::move(source._monomers);
	
	_label_p = std::move(source._label_p);

    _type = std::move(source._type);
	
	_origin_x = std::move(source._origin_x);

    _origin_y = std::move(source._origin_y);
 	
	_origin_z = std::move(source._origin_z);

    _rotations = std::move(source._rotations);

    _bid = std::move(source._bid);

    return *this;
}

bool
Polymer::operator==(const Polymer& other) const
{
    if (!std::equal(_rotations.begin(), _rotations.end(),
                    other._rotations.begin(), other._rotations.end()))
    {
        return false;
    }
    if (!std::equal(_monomers.begin(), _monomers.end(),
                    other._monomers.begin(), other._monomers.end()))
    {
        return false;
    }
    if (_label_p != other._label_p) return false;

    if (_type != other._type) return false;
	
	if (std::fabs(_origin_x - other._origin_x) > 1.0e-6) return false;
	
	if (std::fabs(_origin_y - other._origin_y) > 1.0e-6) return false;
    
	if (std::fabs(_origin_z - other._origin_z) > 1.0e-6) return false;

    if (_bid != other._bid) return false;

    return true;
}

bool
Polymer::operator!=(const Polymer& other) const
{
    return !(*this == other);
}

std::string
Polymer::_add_label()
{
    ++polymer_counter;
    return (std::to_string(polymer_counter));
}


void
Polymer::position(float& rx,float& ry, float& rz, const int32_t index) const
{
    rx = _origin_x + 0.5 * _monomers[0].length();
    ry = _origin_y + 0.5 * _monomers[0].width();
    rz = _origin_z;
    if (index > 0)
    {
        for (int32_t i = 0; i < index; i++)
        {
            float px, py, pz;

            _rotations[i].rotate(px, py, pz, 0.0f, 0.0f, _monomers[i].height());

            rx += px;
            ry += py;
            rz += pz;
        };
    }
}

void
Polymer::position_center(float& rx,float& ry, float& rz, const int32_t index) const
{
    rx = _origin_x + 0.5 * _monomers[0].length();
    ry = _origin_y + 0.5 * _monomers[0].width();
    rz = _origin_z;

    if (index > 0)
    {
        for (int32_t i = 0; i < index; i++)
        {
            float px, py, pz;
            _rotations[i].rotate(px, py, pz, 0.0f, 0.0f, _monomers[i].height());

            rx += px;
            ry += py;
            rz += pz;
        };
    }
    float px, py, pz;
    _rotations[index].rotate(px, py, pz, 0.0f, 0.0f, 0.5f * _monomers[index].height());
    rx += px;
    ry += py;
    rz += pz;
}


void
Polymer::position_center_alternative(float& rx,float& ry, float& rz, const int32_t index) const
{
    rx = (_monomers[index].get_coord(0)[0] + _monomers[index].get_coord(6)[0]) * 0.5f;
    ry = (_monomers[index].get_coord(0)[1] + _monomers[index].get_coord(6)[1]) * 0.5f;
    rz = (_monomers[index].get_coord(0)[2] + _monomers[index].get_coord(6)[2]) * 0.5f;
}

void
Polymer::position_center_alternative_2(float& rx,float& ry, float& rz, const int32_t index) const
{
    rx = _monomers[index].get_center_x();
    ry = _monomers[index].get_center_y();
    rz = _monomers[index].get_center_z();
}

void
Polymer::set_mcoord(std::vector<std::vector<double>> mcoord, int32_t index)
{
    _monomers[index].set_coords(mcoord);
}


std::vector<std::vector<double>>
Polymer::get_mcoord(int32_t index) const
{
    std::vector<std::vector<double>> mcoord;
    mcoord = _monomers[index].get_coords();
    return mcoord;
}


void
Polymer::get_mcoord(int32_t index, struct bd *bd) const
{
    for (int32_t i = 0; i < 8; i++)
    {
        for (int32_t j = 0; j < 3; j++)
        {
            bd->coord[i][j] = _monomers[index].get_coord(i)[j];
        }
    }
}


void
Polymer::get_coords_to_vec(float *ptr_coord, int index) const
{
    for (auto &m: _monomers)
    {
        for (int32_t i = 0; i < 8; i++)
        {
            for (int32_t j = 0; j < 3; j++)
            {
                ptr_coord[index + i*3 + j] = m.get_coord(i)[j];
            }
        }
    }
}

void
Polymer::calculate_monomer_center(int32_t index, float x, float y, float z)
{
    _monomers[index].set_center(x, y, z);
}

void
Polymer::fill_coords(int32_t index)
{
    std::vector<std::vector<double>> coords;

    auto rx = _origin_x + 0.5 * _monomers[0].length();
    auto ry = _origin_y + 0.5 * _monomers[0].width();
    auto rz = _origin_z;

    if (index > 0)
    {
        for (int32_t i = 0; i < index; i++)
        {
            float px, py, pz;

            _rotations[i].rotate(px, py, pz, 0.0f, 0.0f, _monomers[i].height());

            rx += px;
            ry += py;
            rz += pz;
        };
    }
    float v0x, v0y, v0z;
    _rotations[index].rotate(v0x, v0y, v0z, -0.5f * _monomers[index].length(), -0.5f * _monomers[index].width(),
                             0.0f * _monomers[index].height());
    std::vector<double> vertice_0 = {rx + v0x, ry + v0y, rz + v0z};
    coords.push_back(vertice_0);

    float v1x, v1y, v1z;
    _rotations[index].rotate(v1x, v1y, v1z, -0.5f * _monomers[index].length(), 0.5f * _monomers[index].width(),
                             0.0f * _monomers[index].height());
    std::vector<double> vertice_1 = {rx + v1x, ry + v1y, rz + v1z};
    coords.push_back(vertice_1);

    float v2x, v2y, v2z;
    _rotations[index].rotate(v2x, v2y, v2z, 0.5f * _monomers[index].length(), 0.5f * _monomers[index].width(), 0.0f * _monomers[index].height());
    std::vector<double> vertice_2 = {rx + v2x, ry+v2y, rz+v2z};
    coords.push_back(vertice_2);

    float v3x, v3y, v3z;
    _rotations[index].rotate(v3x, v3y, v3z, 0.5f * _monomers[index].length(), -0.5f * _monomers[index].width(), 0.0f * _monomers[index].height());
    std::vector<double> vertice_3 = {rx+v3x, ry+v3y, rz+v3z};
    coords.push_back(vertice_3);

    float v4x, v4y, v4z;
    _rotations[index].rotate(v4x, v4y, v4z, -0.5f * _monomers[index].length(), -0.5f * _monomers[index].width(), 1.0f * _monomers[index].height());
    std::vector<double> vertice_4 = {rx+v4x, ry+v4y, rz+v4z};
    coords.push_back(vertice_4);

    float v5x, v5y, v5z;
    _rotations[index].rotate(v5x, v5y, v5z, -0.5f * _monomers[index].length(), 0.5f * _monomers[index].width(), 1.0f * _monomers[index].height());
    std::vector<double> vertice_5 = {rx+v5x, ry+v5y, rz+v5z};
    coords.push_back(vertice_5);

    float v6x, v6y, v6z;
    _rotations[index].rotate(v6x, v6y, v6z, 0.5f * _monomers[index].length(), 0.5f * _monomers[index].width(), 1.0f * _monomers[index].height());
     std::vector<double> vertice_6 = {rx+v6x, ry+v6y, rz+v6z};
     coords.push_back(vertice_6);

    float v7x, v7y, v7z;
    _rotations[index].rotate(v7x, v7y, v7z, 0.5f * _monomers[index].length(), -0.5f * _monomers[index].width(), 1.0f * _monomers[index].height());
    std::vector<double> vertice_7 = {rx+v7x, ry+v7y, rz+v7z};
    coords.push_back(vertice_7);

    set_mcoord(coords, index);

    float center_x = (vertice_0[0] + vertice_6[0]) * 0.5f;
    float center_y = (vertice_0[1] + vertice_6[1]) * 0.5f;
    float center_z = (vertice_0[2] + vertice_6[2]) * 0.5f;

    calculate_monomer_center(index, center_x, center_y, center_z);
}

void
Polymer::calculate_coords()
{
    for (int32_t i = 0; i < number_of_monomers(); i++) fill_coords(i);
}


void
Polymer::add(Monomer &monomer,
             const ExponentialMap &rotation)
{
    if (number_of_monomers() > 0) fill_coords(number_of_monomers() - 1);
    auto nom = number_of_monomers() - 1;
    bool checker = false;

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();
    auto ptr_mon = &monomer;

    for (int32_t i = 0; i < nom; i++)
    {

        if (intersecting_gjk(*ptr_mon, rotation, i, &bd1, &bd2))
        {
            checker = true;
        }
    }
    free_bd(&bd1);
    free_bd(&bd2);

    if (!checker)
    {
        _monomers.push_back(monomer);
        _rotations.push_back(rotation);
        fill_coords(number_of_monomers() - 1);

    } else
    {
        return;
    }
}


bool
Polymer::intersecting_gjk(const Monomer &monom,
                          const ExponentialMap &rotation,
                          const int32_t index,
                          struct bd *bd1,
                          struct bd *bd2) const
{
    if (number_of_monomers() == 0) return false;

    get_mcoord(index, bd1);

    float mposx, mposy, mposz;
    position(mposx, mposy, mposz, number_of_monomers());

    float v0x, v0y, v0z;
    rotation.rotate(v0x, v0y, v0z, -0.5f * monom.length(), -0.5f * monom.width(),
                    0.0f * monom.height());
    bd2->coord[0][0] = mposx + v0x;
    bd2->coord[0][1] = mposy + v0y;
    bd2->coord[0][2] = mposz + v0z;

    float v1x, v1y, v1z;
    rotation.rotate(v1x, v1y, v1z, -0.5f * monom.length(), 0.5f * monom.width(),
                    0.0f * monom.height());
    bd2->coord[1][0] = mposx + v1x;
    bd2->coord[1][1] = mposy + v1y;
    bd2->coord[1][2] = mposz + v1z;

    float v2x, v2y, v2z;
    rotation.rotate(v2x, v2y, v2z, 0.5f * monom.length(), 0.5f * monom.width(),
                    0.0f * monom.height());
    bd2->coord[2][0] = mposx + v2x;
    bd2->coord[2][1] = mposy + v2y;
    bd2->coord[2][2] = mposz + v2z;

    float v3x, v3y, v3z;
    rotation.rotate(v3x, v3y, v3z, 0.5f * monom.length(), -0.5f * monom.width(),
                    0.0f * monom.height());
    bd2->coord[3][0] = mposx + v3x;
    bd2->coord[3][1] = mposy + v3y;
    bd2->coord[3][2] = mposz + v3z;

    float v4x, v4y, v4z;
    rotation.rotate(v4x, v4y, v4z, -0.5f * monom.length(), -0.5f * monom.width(),
                    1.0f * monom.height());
    bd2->coord[4][0] = mposx + v4x;
    bd2->coord[4][1] = mposy + v4y;
    bd2->coord[4][2] = mposz + v4z;

    float v5x, v5y, v5z;
    rotation.rotate(v5x, v5y, v5z, -0.5f * monom.length(), 0.5f * monom.width(),
                    1.0f * monom.height());
    bd2->coord[5][0] = mposx + v5x;
    bd2->coord[5][1] = mposy + v5y;
    bd2->coord[5][2] = mposz + v5z;

    float v6x, v6y, v6z;
    rotation.rotate(v6x, v6y, v6z, 0.5f * monom.length(), 0.5f * monom.width(),
                    1.0f * monom.height());
    bd2->coord[6][0] = mposx + v6x;
    bd2->coord[6][1] = mposy + v6y;
    bd2->coord[6][2] = mposz + v6z;

    float v7x, v7y, v7z;
    rotation.rotate(v7x, v7y, v7z, 0.5f * monom.length(), -0.5f * monom.width(),
                    1.0f * monom.height());
    bd2->coord[7][0] = mposx + v7x;
    bd2->coord[7][1] = mposy + v7y;
    bd2->coord[7][2] = mposz + v7z;

    if (check_gjk_intersection(bd1, bd2)) return true;

    return  false;
}


bool
Polymer::polymer_intersecting_connect(const Polymer &other) const
{
    auto nom = number_of_monomers();
    auto onom = other.number_of_monomers();

    bool checker = false;

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    for (int32_t i = 0; i < nom; i++)
    {
        if (!checker)
        {
            get_mcoord(i, &bd1);
            for (int32_t j = 0; j < onom; j++)
            {
                if ((i == number_of_monomers() - 1) && (j == 0)) continue;
                other.get_mcoord(j, &bd2);
                if (check_gjk_intersection(&bd1, &bd2))
                    checker = true;
                break;
            };
        }
    }
    free_bd(&bd1);
    free_bd(&bd2);

    return checker;
}


bool
Polymer::polymer_intersecting_from (const Polymer& other,
                                    const int32_t index_from ) const
{
    auto nom = number_of_monomers();
    auto onom = other.number_of_monomers();

    bool checker = false;

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    for (int32_t i = 0; i < nom; i++)
    {
        if (!checker)
        {
            get_mcoord(i, &bd1);

            for (int32_t j = index_from; j < onom; j++)
            {
                other.get_mcoord(j, &bd2);

                if (check_gjk_intersection(&bd1, &bd2))
                {
                    checker = true;
                    break;
                };
            }
        }
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return checker;
}


bool
Polymer::polymer_intersecting(Polymer &other, struct bd *bd1, struct bd *bd2) const
{
    const auto nom = number_of_monomers();
    const auto onom = other.number_of_monomers();

    for (int32_t i = 0; i < nom; i++)
    {
        get_mcoord(i, bd1);

        for (int32_t j = 0; j < onom; j++)
        {
            other.get_mcoord(j, bd2);
            if (check_gjk_intersection(bd1, bd2))
            {
                return true;
            }
        }
    }
    return false;
}


bool
Polymer::polymer_intersecting_from_second(Polymer &other, struct bd *bd1, struct bd *bd2) const
{
    const auto nom = number_of_monomers();
    const auto onom = other.number_of_monomers();

    for (int32_t i = 0; i < nom; i++)
    {
        get_mcoord(i, bd1);

        for (int32_t j = 1; j < onom; j++)
        {
            other.get_mcoord(j, bd2);
            if (check_gjk_intersection(bd1, bd2))
            {
                return true;
            }
        }
    }
    return false;
}

bool
Polymer::polymer_intersecting_return (Polymer &other, struct bd *bd1, struct bd *bd2, int32_t& im) const
{
    const auto nom = number_of_monomers();
    const auto onom = other.number_of_monomers();

    for (int32_t i = 0; i < nom; i++)
    {
        get_mcoord(i, bd1);
        for (int32_t j = 0; j < onom; j++)
        {
            other.get_mcoord(j, bd2);
            if (check_gjk_intersection(bd1, bd2))
            {
                im = i;
                return true;
            }
        }
    }
    return false;
}


bool
Polymer::polymer_intersecting_return_branch(Polymer &other, struct bd *bd1, struct bd *bd2,
                                            int32_t &mon) const
{
    other.get_mcoord(0, bd1);
    for (int32_t i = 1; i < (number_of_monomers()); i++)
    {
        get_mcoord(i, bd2);
        if (check_gjk_intersection(bd1, bd2))
        {
            mon = i;
            return true;
        }
    }
    //special cases
    get_mcoord(0, bd1);
    for (int32_t i = 0; i < (other.number_of_monomers()); i++)
    {
        other.get_mcoord(i, bd2);
        if (check_gjk_intersection(bd1, bd2))
        {
            mon = 0;
            return true;
        }
    }
    get_mcoord(number_of_monomers()-1, bd1);
    for (int32_t i = 0; i < (other.number_of_monomers()); i++)
    {
        other.get_mcoord(i, bd2);
        if (check_gjk_intersection(bd1, bd2))
        {
            mon = number_of_monomers() - 1;
            return true;
        }
    }
    return false;
}


bool
Polymer::polymer_intersecting_as_branch(Polymer &other, struct bd *bd1, struct bd *bd2) const
{
    const auto onom = other.number_of_monomers();

    get_mcoord(0, bd1);
    for (int32_t i = 0; i < onom; i++)
    {
        other.get_mcoord(i, bd2);
        if (check_gjk_intersection(bd1, bd2))
        {
            return true;
        }
    }
    return false;
}

bool
Polymer::polymer_intersecting_return_connection(Polymer &other, struct bd *bd1, struct bd *bd2, int32_t &moni, int32_t &omoni) const
{
    const auto nom = number_of_monomers();
    const auto onom = other.number_of_monomers();

    for (int32_t i = 0; i < nom; i++)
    {
        get_mcoord(i, bd1);

        for (int32_t j = 0; j < onom; j++)
        {
            other.get_mcoord(j, bd2);
            if (check_gjk_intersection(bd1, bd2))
            {
                moni = i;
                omoni = j;
                return true;
            }
        }
    }
    return false;
}

bool
Polymer::polymer_check_if_branch(Polymer &other, struct bd *bd1, struct bd *bd2) const
{
    const auto nom = number_of_monomers();
    const auto onom = other.number_of_monomers();

    for (int32_t i = 1; i < nom - 1; i++)
    {
        get_mcoord(i, bd1);
        for (int32_t j = 1; j < onom - 1; j++)
        {
            other.get_mcoord(j, bd2);
            if (check_gjk_intersection(bd1, bd2))
            {
//                std::cout<<"Kertasi "<<get_label()<< " vs  "<<other.get_label()<<std::endl;
//                std::cout<<i<< " vs   "<<j<<std::endl;
                if (j == 0) continue;
                return true;
            }
        }
    }
    return false;
}


bool
Polymer::polymer_monomer_intersecting(int32_t moni, Polymer &other, struct bd *bd1, struct bd *bd2) const
{
    get_mcoord(moni, bd1);
    other.get_mcoord(0, bd2);
    if (check_gjk_intersection(bd1, bd2))
    {
        return true;
    }
    return false;
}


void
Polymer::set_label(std::string label)
{
    _label_p = label;
}

std::string
Polymer::get_label() const
{
    return _label_p;
}

void
Polymer::set_type(PolymerType type)
{
    _type = type;
}


PolymerType
Polymer::get_type() const
{
    return _type;
}

void
Polymer::set_bid(unsigned long long bid)
{
    _bid = bid;
}

unsigned long long
Polymer::get_bid() const
{
    return _bid;
}


int32_t
Polymer::number_of_monomers() const
{
    return static_cast<int32_t>(_monomers.size());
}

void
Polymer::set_origins(const float originx,
                     const float originy,
                     const float originz)
{
    _origin_x = originx;
    _origin_y = originy;
    _origin_z = originz;
    for (int32_t i= 0; i < _monomers.size(); i++) fill_coords(i);
}


void
Polymer::set_rotation (const float rx,
                       const float ry,
                       const float rz,
                       const float teta,
                       const int32_t index)
{
    ExponentialMap rotation (rx, ry, rz, teta);
    _rotations[index] = rotation;
}


void
Polymer::set_rotation (const ExponentialMap& new_rot,
                       const int32_t index)
{
    _rotations[index] = new_rot;
}


void
Polymer::get_rotation(float &rx,
                      float &ry,
                      float &rz,
                      float &teta,
                      int32_t index) const
{
    _rotations[index].get_rotation(rx,ry,rz,teta);
}


ExponentialMap
Polymer::rotation(int32_t index) const
{
    return _rotations[index];
}


std::ostream&
operator<<(std::ostream& output, const Polymer& source)
{
    output << "[Polymer (Instance): " << &source << "]" << std::endl;
    output << "label_p: " << source._label_p <<  " ";
    static std::array<const char*, 3> enumNames = {"valid", "invalid"};
    output << "type: " << enumNames[static_cast<int>(source._type)]<< " ";
    output << "bid " << source._bid <<  " ";
    output << "nom " << source.number_of_monomers() <<  " ";
    output << "origin : " << source._origin_x << " " << source._origin_y <<  " " <<  source._origin_z <<  std::endl;
    output << "this->_rotations(vector): " << std::endl;
    std::for_each(source._rotations.begin(), source._rotations.end(),
                  [&output] (const auto &tval) {output << tval;});
    return output;
}


void
Polymer::polymer_merge(Polymer &other)
{
    _monomers.insert(_monomers.end(), other._monomers.begin(), other._monomers.end());
    _rotations.insert(_rotations.end(), other._rotations.begin(), other._rotations.end());
}

Monomer
Polymer::monomer(int n) const
{
    return _monomers.at(n);
}


void
Polymer::divide( int32_t monomer_index,
                 Polymer& first_half,
                 Polymer& second_half)
{
    std::vector<Monomer> first_half_m = {_monomers.begin(), _monomers.end() - (_monomers.size() - monomer_index)};
    std::vector<ExponentialMap> first_half_r = {_rotations.begin(), _rotations.end() - (_monomers.size() - monomer_index)};
    std::vector<Monomer> second_half_m= {_monomers.begin() + monomer_index, _monomers.end()};
    std::vector<ExponentialMap> second_half_r= {_rotations.begin() + monomer_index, _rotations.end()};

    first_half = Polymer(first_half_m, _add_label(), PolymerType::valid, _origin_x, _origin_y, _origin_z, first_half_r,0);
    float x,y,z;
    position(x, y, z, monomer_index);
    second_half = Polymer(second_half_m, _add_label(), PolymerType::valid, x - monomer(monomer_index).length() / 2,
                                                 y - monomer(monomer_index).length() / 2,
                                                 z, second_half_r, 0);
}


void Polymer::rotate_polymer(const Pt &centre, const ExponentialMap &newr)
{
    float mx, my, mz, mtheta;
    newr.get_rotation(mx, my, mz, mtheta);

    if ((mz == 1.0f) && (mtheta == -1.0f)) mtheta = -0.99999f;

    Pt axis = {mx, my, mz};
    Mat M = rotationMatrix(axis, mtheta);

    Pt start = {_origin_x, _origin_y, _origin_z};
    start = centre + M * (start - centre);
    _origin_x = start.x;
    _origin_y = start.y;
    _origin_z = start.z;

    for (int32_t i = 0; i < number_of_monomers(); i++)
    {
        float rx, ry, rz, theta;
        _rotations[i].get_rotation(rx, ry, rz, theta);

        Pt monomer_axis = {rx, ry, rz};
        Mat P = rotationMatrix(monomer_axis, theta);

        Mat R = M * P;

        axisAngle(R, monomer_axis, theta);

        ExponentialMap mr(monomer_axis.x, monomer_axis.y, monomer_axis.z, theta);
        set_rotation(mr, i);
        fill_coords(i);
    }
}


void
Polymer::shorten (int32_t nr)
{
    _monomers.erase(_monomers.begin() + (_monomers.size()-nr), _monomers.begin() + _monomers.size());

    _rotations.erase(_rotations.begin() + (_rotations.size()-nr), _rotations.begin() + _rotations.size());
}


void
Polymer::add_for_saving (Monomer& monomer, ExponentialMap& exponentialMap)
{
    _monomers.push_back(monomer);
    _rotations.push_back(exponentialMap);
}


bool
Polymer::intersecting_itself()
{
    if (number_of_monomers() == 0) return false;

    auto bd1 = allocate_bd();
    auto bd2 = allocate_bd();

    for (int32_t i = 0; i < _monomers.size(); i++)
    {
        get_mcoord(i, &bd1);

        for (int32_t j = 0; j < _monomers.size(); j++)
        {
            if ((j == i - 1) || (j == i) || (j == i + 1)) continue;

            get_mcoord(j, &bd2);

            if (check_gjk_intersection(&bd1, &bd2))
            {
                return true;
            }
        }
    }
    return  false;
}

void
Polymer::get_mass_center(float& mx, float& my, float& mz) const
{
    float monomer_mass = 1.0f;
    float sum_x = 0, sum_y = 0, sum_z = 0.0f;
    for (int32_t i = 0; i < number_of_monomers(); i++)
    {
        float px, py, pz;
        position_center(px, py, pz, i);

        sum_x += px * monomer_mass;
        sum_y += py * monomer_mass;
        sum_z += pz * monomer_mass;
    }
    mx = sum_x / number_of_monomers();
    my = sum_y / number_of_monomers();
    mz = sum_z / number_of_monomers();
}

void
Polymer::reverse()
{
    float px, py, pz;
    position(px, py, pz, number_of_monomers());
    py -= monomer(number_of_monomers() - 1).width() / 2;
    px -= monomer(number_of_monomers() - 1).length() / 2;

    Mat M = rotationMatrix(Pt{0.0, 1.0, 0.0}, getPiValue());

    std::reverse(_monomers.begin(), _monomers.end());
    std::reverse(_rotations.begin(), _rotations.end());

    for (int32_t i = 0; i < number_of_monomers(); i++)
    {
        float rx, ry, rz, theta;
        _rotations[i].get_rotation(rx, ry, rz, theta);
        Pt monomer_axis = {rx, ry, rz};

        Mat P = rotationMatrix(monomer_axis, theta);

        Mat R = P * M;

        axisAngle(R, monomer_axis, theta);

        ExponentialMap mr(monomer_axis.x, monomer_axis.y, monomer_axis.z, theta);
        set_rotation(mr, i);
    }
    set_origins(px, py, pz);
}

bool
Polymer::is_valid() const
{
    if (get_type() == PolymerType::valid)
        return true;
    return false;
}

bool
Polymer::is_branched() const
{
    if (_bid != 0)
        return true;
    return false;
}

bool
Polymer::has_coords() const
{
    for (auto i = 0; i < number_of_monomers(); i++)
    {
        auto coords = get_mcoord(i);
        if (coords.empty()) return false;
    }
    return true;
}