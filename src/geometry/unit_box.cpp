//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include <cmath>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <iterator>
#include <omp.h>

#include "unit_box.hpp"


UnitBox::UnitBox()

    : _polymers(std::vector<Polymer>())

    , _bran_polymers(std::vector<Branched_polymer>())

	, _box_label(std::string())
	
	, _box_origin_x(0.0f)

    , _box_origin_y(0.0f)

    , _box_origin_z(0.0f)
	
	, _box_x(0.0f)
	
	, _box_y(0.0f)
	
	, _box_z(0.0f)

    , _irradiation(std::vector<std::vector<float>>())

{
    
}

UnitBox::UnitBox(const std::vector<Polymer>& polymers,
                 const std::vector<Branched_polymer>& bran_polymers,
                 const std::string& box_label,
                 const float box_origin_x,
                 const float box_origin_y,
                 const float box_origin_z,
                 const float box_x,
                 const float box_y,
                 const float box_z,
                 const std::vector<std::vector<float>>& irradiation)

    : _polymers (polymers)

    , _bran_polymers(bran_polymers)

	, _box_label(box_label)
	
	, _box_origin_x(box_origin_x)

    , _box_origin_y(box_origin_y)

    , _box_origin_z(box_origin_z)
	
	, _box_x(0.0f)
	
	, _box_y(0.0f)
	
	, _box_z(0.0f)

	, _irradiation(irradiation)
	
{
}

UnitBox::UnitBox(const UnitBox& source)

    : _polymers(source._polymers)

    , _bran_polymers(source._bran_polymers)

	, _box_label(source._box_label)
	
	, _box_origin_x(source._box_origin_x)

    , _box_origin_y(source._box_origin_y)

    , _box_origin_z(source._box_origin_z)
	
	, _box_x(source._box_x)

    , _box_y(source._box_y)

    , _box_z(source._box_z)

    , _irradiation(source._irradiation)
{
    
}

UnitBox::UnitBox(UnitBox&& source) noexcept

    : _polymers(std::move(source._polymers))

    , _bran_polymers(std::move(source._bran_polymers))

	, _box_label(std::move(source._box_label))
		
	, _box_origin_x(std::move(source._box_origin_x))

    , _box_origin_y(std::move(source._box_origin_y))

    , _box_origin_z(std::move(source._box_origin_z))
	
	, _box_x(std::move(source._box_x))

    , _box_y(std::move(source._box_y))

    , _box_z(std::move(source._box_z))

    , _irradiation(std::move(source._irradiation))
{
    
}

UnitBox::~UnitBox()
{

}

UnitBox&
UnitBox::operator=(const UnitBox& source)
{
    if (this == &source) return *this;
    
    _polymers = source._polymers;
    _bran_polymers = source._bran_polymers;
	_box_label = source._box_label;
	_box_origin_x = source._box_origin_x;
    _box_origin_y = source._box_origin_y;
    _box_origin_z = source._box_origin_z;
	_box_x = source._box_x;
    _box_y = source._box_y;
    _box_z = source._box_z;
    _irradiation = source._irradiation;
    
    return *this;
}

UnitBox&
UnitBox::operator=(UnitBox&& source) noexcept
{
    if (this == &source) return *this;
    
    _polymers = std::move(source._polymers);
    _bran_polymers = std::move(source._bran_polymers);
    _box_label = std::move(source._box_label);
	_box_origin_x = std::move(source._box_origin_x);
    _box_origin_y = std::move(source._box_origin_y);
    _box_origin_z = std::move(source._box_origin_z);
	_box_x = std::move(source._box_x);
    _box_y = std::move(source._box_y);
    _box_z = std::move(source._box_z);
    _irradiation = std::move(source._irradiation);
        
    return *this;
}

bool
UnitBox::operator==(const UnitBox& other) const
{
    if (!std::equal(_polymers.begin(), _polymers.end(), other._polymers.begin(), other._polymers.end())) return false;
    if (!std::equal(_bran_polymers.begin(), _bran_polymers.end(), other._bran_polymers.begin(), other._bran_polymers.end())) return false;
    if (_box_label != other._box_label) return false;
	if (std::fabs(_box_origin_x - other._box_origin_x) > 1.0e-6) return false;
	if (std::fabs(_box_origin_y - other._box_origin_y) > 1.0e-6) return false;
	if (std::fabs(_box_origin_z - other._box_origin_z) > 1.0e-6) return false;
	if (std::fabs(_box_x - other._box_x) > 1.0e-6) return false;
	if (std::fabs(_box_y - other._box_y) > 1.0e-6) return false;
	if (std::fabs(_box_z - other._box_z) > 1.0e-6) return false;
    if (!std::equal(_irradiation.begin(), _irradiation.end(), other._irradiation.begin(), other._irradiation.end())) return false;
	
	return true;
}

bool
UnitBox::operator!=(const UnitBox& other) const
{
    return !(*this == other);
}

void
UnitBox::set_box_size(const float box_x, const float box_y, const float box_z)
{
}

void
UnitBox::set_label(const std::string& label)
{
    _box_label = label;
}

float
UnitBox::origin_x() const
{
    return _box_origin_x;
}

float
UnitBox::origin_y() const
{
    return _box_origin_y;
}

float
UnitBox::origin_z() const
{
    return _box_origin_z;
}

float
UnitBox::size_x() const
{
    return _box_x;
}

float
UnitBox::size_y() const
{
    return _box_y;
}

float
UnitBox::size_z() const
{
    return _box_z;
}

std::string
UnitBox::label() const
{
    return _box_label;
}

std::vector<Polymer>
UnitBox::get_polymers() const
{
    return _polymers;
}

std::vector<Branched_polymer>
UnitBox::get_branpols() const
{
    return _bran_polymers;
}

Branched_polymer
UnitBox::branpol(int32_t index) const
{
    return _bran_polymers.at(index);
}

void
UnitBox::add(const Polymer& polymer)
{
    _polymers.push_back(polymer);
}

void
UnitBox::add(std::vector<Polymer>& polymers)
{
    for (auto &polymer : polymers)
    {
        _polymers.push_back(polymer);
    }
}

int32_t
UnitBox::number_of_polymers() const
{
    return _polymers.size();
}

int32_t
UnitBox::number_of_polymers(const unsigned long long bid) const
{
    return static_cast<int32_t>(std::count_if(_polymers.begin(), _polymers.end(),
                                              [&] (const auto &tval) {return tval.get_bid() == bid;}));
}

int32_t
UnitBox::number_of_branched_polymers() const
{
    return static_cast<int32_t>(std::count_if(_polymers.begin(), _polymers.end(),
                                              [&] (const auto &tval) {return tval.get_bid() != 0;}));
}

int32_t
UnitBox::number_of_monomers() const
{
    auto nom = 0;
    for (auto &p : _polymers)
    {
        nom += p.number_of_monomers();
    }
    return  nom;
}


Polymer
UnitBox::polymer(int n) const
{
    return _polymers.at(n);
}

Polymer*
UnitBox::polymer_ptr(int n)
{
    return &_polymers.at(n);
}

Branched_polymer*
UnitBox::branpol_ptr(int n)
{
    return &_bran_polymers.at(n);
}

Polymer
UnitBox::polymer(const std::string& label) const
{
    for (auto &p : _polymers)
    {
        if (p.get_label() == label)
        {
            return p;
        }
    }
    Polymer p;
    return p;
}


std::ostream&
operator<<(std::ostream& output, const UnitBox& source)
{
    output << std::endl;
    
    output << "[UnitBox (Instance): " << &source << "]" << std::endl;
	output << "this->_polymers(vector): " << std::endl;
	std::for_each(source._polymers.begin(), source._polymers.end(),
                  [&output] (const auto &tval) {output << tval.get_label()<< " : " << tval.get_bid() << std::endl;});
    output << "this->_bran_polymers(vector): " << std::endl;
    std::for_each(source._bran_polymers.begin(), source._bran_polymers.end(),
                  [&output] (const auto &tval) {output << tval << std::endl;});

    return output;
}


int
UnitBox::find_polymer (const Polymer& polymer)
{
    auto iter = std::find(_polymers.begin(), _polymers.end(), polymer);
    return(std::distance(_polymers.begin(), iter));
}


void
UnitBox::remove_polymer (const Polymer& polymer)
{
    auto erPol = std::remove(_polymers.begin(), _polymers.end(), polymer);
    _polymers.erase(erPol, _polymers.end());
}


void
UnitBox::remove_polymer (int32_t index)
{
    auto it = _polymers.begin() + index;
    if (&(*it) == &(_polymers.back()))
    {
        _polymers.pop_back();
    }
    *it = std::move(_polymers.back());
    _polymers.pop_back();
}

void
UnitBox::remove_branpol (const Branched_polymer& branpolymer)
{
    auto bid = branpolymer.get_id();

    auto erPol = std::remove(_bran_polymers.begin(), _bran_polymers.end(), branpolymer);
    _bran_polymers.erase(erPol, _bran_polymers.end());

    auto it = std::remove_if(_polymers.begin(), _polymers.end(), [&](const auto &tval){ return tval.get_bid() == bid; });
    _polymers.erase(it, _polymers.end());
}


bool
UnitBox::check_polymer_intersection(Polymer &polymer)
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    for (int32_t i = 0; i < _polymers.size(); i++)
    {
        if (polymer.polymer_intersecting(_polymers[i], &bd1, &bd2))
        {
            free_bd(&bd1);
            free_bd(&bd2);
            return true;
        }
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}


bool
UnitBox:: check_polymer_intersection_return (Polymer& polymer,
                                             int32_t& ip,
                                             int32_t& im)
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    for (int32_t i = 0; i < _polymers.size(); i++)
    {
        if (_polymers[i].polymer_intersecting_return(polymer, &bd1, &bd2, im))
        {
            ip = i;
            free_bd(&bd1);
            free_bd(&bd2);
            return true;
        }
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}


bool UnitBox::is_intersecting_between_skip_same_branches()
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    for (auto &polymer : _polymers)
    {
        for (auto &other: _polymers)
        {
            if (polymer.get_bid() == other.get_bid()) continue;
            if (polymer.polymer_intersecting(other, &bd1, &bd2))
            {
                std::cout<<polymer<<std::endl;
                std::cout<<other<<std::endl;
                free_bd(&bd1);
                free_bd(&bd2);
                return true;
            }
        }
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}


bool
UnitBox::find_nearest_skip_same_branch(int32_t polymer_index,
                                       int32_t monomer_index,
                                       float distance,
                                       int32_t &nearest_polymer_index,
                                       int32_t &nearest_monomer_index) const
{
    if (number_of_polymers() < 2) return false;

    float d = 1.0f;
    while (d < distance)
    {
        auto sq_distance = d * d;
        float rx, ry, rz;
        _polymers[polymer_index].position_center_alternative_2(rx, ry, rz, monomer_index);
        for (int32_t i = 0; i < number_of_polymers(); i++)
        {
            if (i == polymer_index) continue;
            if (_polymers[polymer_index].get_bid() != 0 && (_polymers[i].get_bid() == _polymers[polymer_index].get_bid()))  continue;

            for (int32_t j = 0; j < _polymers[i].number_of_monomers(); j++)
            {
                float jx, jy, jz;
                _polymers[i].position_center_alternative_2(jx, jy, jz, j);

                if ((rx-jx)*(rx-jx) + (ry-jy)*(ry-jy) + (rz-jz)*(rz-jz) <= sq_distance)
                {
                    nearest_polymer_index = i;
                    nearest_monomer_index = j;
                    return true;
                }
            }
        }
        d += 1.0f;
    }
    return false;
}

bool
UnitBox::find_nearest_for_radical(float x,
                                  float y,
                                  float z,
                                  float distance,
                                  int32_t &nearest_polymer_index,
                                  int32_t &nearest_monomer_index) const
{
    auto npol = static_cast<int32_t>(_polymers.size());
    auto ptr_polymers = _polymers.data();
    const auto nthreads = omp_get_max_threads();
    const auto distance_sq = distance * distance;

    std::vector<float> vecdata_dist (nthreads, distance_sq);
    std::vector<int> vecdata_poli (nthreads, 0);
    std::vector<int> vecdata_moni (nthreads, 0);

    auto ptr_vecdata_dist = vecdata_dist.data();
    auto ptr_vecdata_poli = vecdata_poli.data();
    auto ptr_vecdata_moni = vecdata_moni.data();

#pragma omp parallel shared (npol,ptr_polymers, ptr_vecdata_dist, ptr_vecdata_poli, ptr_vecdata_moni, x, y, z, distance)
    {
#pragma omp single nowait
        {
            const auto batch_size = npol / nthreads;
            for (int32_t i = 0; i < nthreads; i++)
            {
                const auto bstart = batch_size * i;
                const auto bend = ((bstart + batch_size) > npol) ? npol : bstart + batch_size;

#pragma omp task firstprivate(i, bstart, bend)
                {
                    for (auto j = bstart; j < bend; j++)
                    {
                        for (int32_t l = 0; l < ptr_polymers[j].number_of_monomers(); l++)
                        {
                            float jx, jy, jz;
                            ptr_polymers[j].position_center_alternative(jx, jy, jz, l);

                            auto actual_distance_sq = (x - jx) * (x - jx) + (y - jy) * (y - jy) + (z - jz) * (z - jz);

                            if ((ptr_vecdata_dist[i] > actual_distance_sq))
                            {
                                ptr_vecdata_dist[i] = actual_distance_sq;
                                ptr_vecdata_poli[i] = j;
                                ptr_vecdata_moni[i] = l;
                            }
                        }
                    }

                }
            }
        }
    }

    auto minElement = std::min_element(vecdata_dist.begin(), vecdata_dist.end());

    if (*minElement >= distance_sq) return false;

    auto index = std::distance(vecdata_dist.begin(), minElement);
    nearest_polymer_index = vecdata_poli[index];
    nearest_monomer_index = vecdata_moni[index];

    return true;
}


bool
UnitBox::check_boundaries(const Polymer &polymer) const
{
    if (!polymer.is_valid())
        return false;

    for (int32_t i = 0; i < polymer.number_of_monomers(); i++)
    {
        float iposx, iposy, iposz;

        polymer.position_center(iposx, iposy, iposz, i);

        if ((iposx >= _box_x) || (iposy >= _box_y) || (iposz >= _box_z)
            || (iposx <= _box_origin_x) || (iposy <= _box_origin_y) || (iposz <= _box_origin_z))
            return true;
    }
    return false;
}

void
UnitBox::save(std::string box_name)
{
    std::ofstream myfile (box_name);

    myfile<<_box_label<<std::endl;
    myfile<<_box_origin_x<<" "<<_box_x<<" "<<_box_origin_y<<" "<<_box_y<<" "<<_box_origin_z<<" "<<_box_z<<std::endl;
    myfile<<_polymers.size()<<std::endl;

    for (int32_t i = 0; i < _polymers.size(); i++)
    {
        myfile<<"--------------------------------------"<<std::endl;
        myfile << _polymers[i].get_label()<< std::endl;
        myfile << _polymers[i].get_bid() << "bid " <<std::endl;
        float x, y, z;
        _polymers[i].position(x, y, z,0);
        myfile << x - _polymers[i].monomer(0).length() / 2 << " " << y - _polymers[i].monomer(0).width() / 2 << " " << z << std::endl;
        myfile << _polymers[i].number_of_monomers() << std::endl;
        for (int32_t j = 0; j < _polymers[i].number_of_monomers(); j++)
        {
            myfile << _polymers[i].monomer(j).label() + std::to_string(j)<< std::endl;;
            myfile << _polymers[i].monomer(j).length() << " " << _polymers[i].monomer(j).width() << " "
                   << _polymers[i].monomer(j).height() << std::endl;
            float rx, ry, rz, teta;
            _polymers[i].get_rotation(rx, ry, rz, teta, j);
            myfile << rx << " " << ry << " " << rz << " " << teta << std::endl;
        }
    }
    myfile.close();
}


void
UnitBox::load(std::string box_name)
{
    std::ifstream file( box_name );
    std::string line;
    std::getline(file, line);
    _box_label = line;

    std::getline(file, line);
    std::istringstream iss(line);
    std::vector<float> origins(std::istream_iterator<float>{iss},
                               std::istream_iterator<float>());
    _box_origin_x = origins[0];
    _box_x = origins[1];
    _box_origin_y = origins[2];
    _box_y = origins[3];
    _box_origin_z = origins[4];
    _box_z = origins[5];

    std::getline(file, line); // polymer number
    int32_t  polymer_nr = std::stoi(line);
    for (int32_t p=0; p < polymer_nr; p++) _polymers.push_back({{},{},{},0.0f, 0.0f, 0.0f, {},{}});
    for (int32_t i = 0; i < polymer_nr; i++)
    {
        std::getline(file, line);// ------------ line
        std::getline(file, line);
        _polymers[i].set_label(line);
        std::getline(file, line);
        _polymers[i].set_bid(static_cast<long >(std::stoi(line)));
        std::getline(file, line);
        std::istringstream isss(line);
        std::vector<float> porigins(std::istream_iterator<float>{isss},
                                   std::istream_iterator<float>());
        _polymers[i].set_origins(porigins[0], porigins[1], porigins[2]);
        std::getline(file, line);//monomer number
        int32_t  monomer_nr = std::stoi(line);
        for (int32_t m=0; m < monomer_nr; m++)
        {
            std::getline(file, line);
            std::string mo_name = line;
            std::getline(file, line);
            std::istringstream issz(line);
            std::vector<float> monsize(std::istream_iterator<float>{issz},
                                        std::istream_iterator<float>());

            std::getline(file, line);
            std::istringstream isr(line);
            std::vector<float> rot (std::istream_iterator<float>{isr},
                                       std::istream_iterator<float>());

            Monomer mo {mo_name, Rectangle3D (monsize[0], monsize[1], monsize[2]), {}, 0.0f, 0.0f, 0.0f};
            ExponentialMap exp {rot[0], rot[1], rot[2], rot[3]};
            _polymers[i].add_for_saving(mo, exp);
        }
    }

    std::vector<std::vector<std::vector<int>>> data(5, std::vector<std::vector<int>>(5, std::vector<int>(5, 0)));
    for (auto i = 0; i < _polymers.size();i++)
    {
        for (auto j = 0; j < _polymers[i].number_of_monomers(); j++)
        {
            float rx, ry, rz;
            _polymers[i].position_center(rx, ry, rz, j);
            auto indX = rx / 100;
            auto indY = ry / 100;
            auto indZ = rz / 100;
            data[indX][indY][indZ]++;
        }
    }
    for (auto& d1 : data)
    {
        for (auto& d2 : d1)
        {
            for (auto& d3: d2)
                std::cout << d3 << std::endl;
        }
    }
}


void
UnitBox::calculate_branching_alternative(std::vector<float> &average_number,
                                         std::vector<float> &average_length,
                                         std::vector<float> &correlation)
{
    average_number.emplace_back(float(_bran_polymers.size()) / float(_polymers.size()));
    float plen = 0.0f;
    for (auto &b : _bran_polymers) plen+=float(b.get_parent().number_of_monomers());
    average_length.emplace_back(plen / float(_bran_polymers.size()));

    std::vector<float> total_number_each_level (1000),total_length_each_level (1000);
    for (auto &bran : _bran_polymers)
    {
        std::vector<int32_t> number;
        std::vector<float> length;
        bran.calculate_all_branches_for_stats(number, length);

        for (auto i = 0; i < number.size(); i++)
        {
            total_number_each_level[i] += number[i];
            total_length_each_level[i] += length[i];
        }
        number.clear();
        length.clear();
    }

    for (auto i = 0; i < total_number_each_level.size(); i++)
    {
        if (total_number_each_level[i] == 0) break;
        float avlen = total_length_each_level[i] / total_number_each_level[i];
        average_length.push_back(avlen);
        average_number.emplace_back(float(total_number_each_level[i]) / float (_polymers.size()));
    }

    std::vector<float> correl;
    correl.reserve(average_number.size());
    for (float i : average_number) correl.push_back(i);
    if (correl.size() % 2 != 0) correl.pop_back();
    for (int32_t i = 0; i < correl.size(); i+=2) correlation.push_back(correl[i] * correl[i + 1]);

}


void
UnitBox::irradiate(float dose_value, int32_t split_size)
{
    std::vector<float> value;
    value.reserve(_box_x / split_size * _box_y / split_size * _box_z / split_size );
    for (int32_t x = split_size; x <= int(_box_x); x += split_size)
    {
        for (int32_t y = split_size; y <= int(_box_y); y += split_size)
        {
            for (int32_t z = split_size; z <= int(_box_z); z += split_size)
            {
                value = {static_cast<float>(x), static_cast<float>(y), static_cast<float>(z), dose_value};
                _irradiation.push_back(value);
            }
        }
    }
}

std::vector<std::vector<float>>&
UnitBox::get_irradiation ()
{
    return _irradiation;
}

void
UnitBox::get_irradiation (int32_t index, float &x, float &y, float &z, float &value) const
{
    x =  _irradiation[index][0];
    y =  _irradiation[index][1];
    z =  _irradiation[index][2];
    value =  _irradiation[index][3];
}

float
UnitBox::get_irradiation_max_value()
{
    auto max = 0.0f;
    for (auto i = 0; i < _irradiation.size(); i++)
    {
        if (max < _irradiation[i][3])
            max = _irradiation[i][3];
    }
    return max;
}

float
UnitBox::get_irradiation_min_value()
{
    auto min = 100.0f;
    for (auto i = 0; i < _irradiation.size(); i++)
    {
        if (min > _irradiation[i][3])
            min = _irradiation[i][3];
    }
    return min;
}

void
UnitBox::change_irradiation (int32_t index, float dose_value)
{
    _irradiation[index][3] = dose_value;
}

bool
UnitBox::find_branpol(Branched_polymer &branpol) const
{
    auto element_found = std::find(_bran_polymers.begin(), _bran_polymers.end(), branpol);
    if (element_found != _bran_polymers.end())
    {
        branpol = *element_found;
        return true;
    }
    return false;
}

Branched_polymer
UnitBox::find_branpol(unsigned long long bid) const
{
    for (auto &b : _bran_polymers)
    {
        if (bid == b.get_id())
            return b;
    }
    abort();
}


bool
UnitBox::add_branpol(Branched_polymer &branpol)
{
     auto branches = branpol.get_branches();
    branches.push_back(branpol.get_parent());
    for (auto & branch : branches)
        add(branch);
    _bran_polymers.push_back(branpol);
    return true;
}


bool
UnitBox::check_branpol_intersection(Branched_polymer &branpol)
{
    auto polymers = branpol.to_polymers();

    for (auto &branch: polymers)
    {
        if (check_polymer_intersection(branch) || check_boundaries(branch))
            return true;
    }
    return false;
}

bool
UnitBox::check_group_intersection(std::vector<Polymer> &polymers)
{
    for (auto &poli: polymers)
    {
        if (check_polymer_intersection(poli) || check_boundaries(poli))
            return true;
    }
    return false;
}
