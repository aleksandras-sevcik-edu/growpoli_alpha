//  Created by Aleksandras Sevcik and Zilvinas Rinkevicius on 2019-08-01.
//  Copyright © 2019 Zilvinas Rinkevicius and Aleksandras Sevcik. All rights reserved.

#include "branched_polymer.hpp"

Branched_polymer::Branched_polymer()

        : _identifier(0)
        , _parent(Polymer())
        , _branches(std::vector<Polymer>())
{
    _identifier = _get_new_id();
}

Branched_polymer::Branched_polymer(Polymer &parent,
                                   std::vector<Polymer> &branches)
     : _identifier(0)
     , _parent(parent)
     , _branches(branches)
{
    _identifier = _get_new_id();
    _parent.set_bid(_identifier);
    for (auto &b : _branches)
        b.set_bid(_identifier);
}

Branched_polymer::Branched_polymer(Polymer &parent,
                                   Polymer &branch)
        : _identifier(0)
        , _parent(parent)
        , _branches(std::vector<Polymer>())
{
    _identifier = _get_new_id();
    _parent.set_bid(_identifier);
    _branches.push_back(branch);
    _branches.back().set_bid(_identifier);
}

Branched_polymer::~Branched_polymer()
{
}

Branched_polymer &Branched_polymer::operator=(const Branched_polymer &source)
{
    if (this == &source) return *this;
    _identifier = source._identifier;
    _parent = source._parent;
    _branches = source._branches;
    return *this;
}


Branched_polymer::Branched_polymer(const Branched_polymer &source)
        :_identifier(source._identifier)
        ,_parent(source._parent)
        ,_branches(source._branches)
{
}

Branched_polymer::Branched_polymer(Branched_polymer &&source) noexcept
        :_identifier(std::move(source._identifier))
        ,_parent(std::move(source._parent))
        ,_branches(std::move(source._branches))
{
}

bool Branched_polymer::operator==(const Branched_polymer &other) const
{
    if (_identifier != other._identifier) return false;
    if (_parent != other._parent) return false;
    if (!std::equal(_branches.begin(), _branches.end(), other._branches.begin(), other._branches.end())) return false;
    return true;
}

bool Branched_polymer::operator!=(const Branched_polymer &other) const
{
    return !(*this == other);
}

unsigned long long Branched_polymer::_get_new_id()
{
    return ++bid_counter;
}


size_t Branched_polymer::get_id() const
{
    return _identifier;
}

Polymer Branched_polymer::get_parent() const
{
    return _parent;
}

std::vector<Polymer> Branched_polymer::get_branches() const
{
    return _branches;
}

void Branched_polymer::_set_branches(std::vector<Polymer> &branches)
{
    _branches = branches;
    for (auto &b : _branches)
        b.set_bid(_identifier);
}

void Branched_polymer::set_parent(Polymer &parent)
{
    _parent = parent;
    _parent.set_bid(_identifier);
}

bool Branched_polymer::add_branches(Polymer &activated, int32_t moni, std::vector<Polymer> &polymers)
{
    if (!is_intersecting_with(polymers))
    {
        _branches.insert(_branches.end(), polymers.begin(), polymers.end());
        for (auto &b: _branches) b.set_bid(_identifier);
        return  true;
    }
    return false;
}

bool Branched_polymer::add_branch(Polymer &activated, int32_t moni, Polymer &polymer)
{
    std::vector<Polymer> polymers = {polymer};
    if (!is_intersecting_with(polymers))
    {
        _branches.insert(_branches.end(), polymers.begin(), polymers.end());
        for (auto &b : _branches) b.set_bid(_identifier);
        return true;
    }
    return false;
}


bool Branched_polymer::change (Polymer &old_polymer, Polymer &new_polymer)
{
    if (_parent == old_polymer)
    {
        _parent = new_polymer;
        if (!is_intersecting_between()) return true;
    } else
    {
        auto index = _find_position(old_polymer);
        _branches[index] = new_polymer;
        if (!is_intersecting_between()) return true;
    }
    return false;
}


void Branched_polymer::from_polymers(std::vector<Polymer> &polymers)
{
    for (auto &p : polymers)
        p.set_bid(_identifier);

    _parent = polymers.back();
    polymers.pop_back();
    _branches = polymers;
}

void Branched_polymer::calculate_coords()
{
    _parent.calculate_coords();
    for (auto &b : _branches) b.calculate_coords();
}

std::vector<Polymer> Branched_polymer::to_polymers() const
{
    std::vector<Polymer> reserve;
    reserve.insert(reserve.end(), _branches.begin(), _branches.end());
    reserve.push_back(_parent);

    return reserve;
}


bool
Branched_polymer::calculate_branches(std::vector<Polymer> &input, std::vector<Polymer> &branches,
                                     std::vector<int32_t> &moni)
{
    if (input.empty()) return false;

    Polymer parent = input.back();
    input.pop_back();

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    for (auto & i : input)
    {
        int32_t  im;
        if (parent.polymer_intersecting_return_branch(i, &bd1, &bd2, im))
        {
            moni.push_back(im);
            branches.push_back(i);
        }
    }
    if (moni.empty())
    {
        free_bd(&bd1);
        free_bd(&bd2);
        return false;
    }
    free_bd(&bd1);
    free_bd(&bd2);
    return true;
}

/* input should contain the parent in the end */
bool
Branched_polymer::calculate_all_branches(std::vector<Polymer> &input, std::vector<Polymer> &output)
{
    if (input.empty()) return false;
    std::vector<Polymer> parents {input.back()};
    std::vector<Polymer> already_checked_as_parents;
    std::vector<std::vector<int32_t>> mon_ind;

    while (true)
    {
        std::vector<Polymer> temp;
        std::vector<Polymer> results;
        std::vector<int32_t> moni;

        for (auto &pol: parents)
        {
            _move_element_to_back(pol, input);
            if (calculate_branches(input, results, moni))
            {
                for (auto &p : results)
                {
                    for (auto &pr : parents)
                    {
                        if (p.get_label() == pr.get_label()) return false;
                    }
                }
                _add_uniques(temp, results);
                mon_ind.push_back(moni);
                results.clear();
                moni.clear();
            }
        }
        if (temp.empty()) break;

        already_checked_as_parents.insert(already_checked_as_parents.end(), parents.begin(), parents.end());
        parents.clear();
        parents = temp;
        if (_check_if_contains_already(already_checked_as_parents, parents)) return false;

        _add_uniques(output, temp);
    }
    if (output.empty())
        return false;
    return true;
}


bool
Branched_polymer::make_division(Polymer &activated, int32_t moni,
                                std::vector<Branched_polymer> &new_branpols,
                                std::vector<Polymer> &new_polymers)
{
    if (activated.number_of_monomers() < 2) return false;
    if (moni < 1) return false;
    std::vector<Polymer> initial_structure = _branches;
    Polymer first, second;
    activated.divide(moni, first, second);


    if (_parent == activated)
    {
        auto branpol_a_structure(initial_structure);
        branpol_a_structure.push_back(first);
        std::vector<Polymer> branpol_a_branches;
        bool c1 = calculate_all_branches(branpol_a_structure, branpol_a_branches);
        if (!c1) return false;
        if (branpol_a_branches.empty())
        {
            new_polymers.push_back(first);
        } else
        {
            Branched_polymer branpol_a(first, branpol_a_branches);
            new_branpols.push_back(branpol_a);
        }
        auto branpol_b_structure(initial_structure);
        branpol_b_structure.push_back(second);
        std::vector<Polymer> branpol_b_branches;
        bool c2 = calculate_all_branches(branpol_b_structure, branpol_b_branches);
        if (!c2) return false;
        if (branpol_b_branches.empty())
        {
            new_polymers.push_back(second);
        } else
        {
            Branched_polymer branpol_b(second, branpol_b_branches);
            new_branpols.push_back(branpol_b);
        }
        if (division_validation_check(new_branpols, new_polymers)) return true;
    }
    else if (std::find(_branches.begin(), _branches.end(), activated) != _branches.end())
    {
        auto erPol = std::remove(initial_structure.begin(), initial_structure.end(), activated);
        initial_structure.erase(erPol, initial_structure.end());
        initial_structure.push_back(_parent);
        auto branpol_a_structure(initial_structure);
        branpol_a_structure.push_back(first);
        std::vector<Polymer> branpol_a_branches;
        bool c3 = calculate_all_branches(branpol_a_structure, branpol_a_branches);
        if (!c3) return false;
        if (branpol_a_branches.empty())
        {
            new_polymers.push_back(first);
        } else
        {
            Branched_polymer branpol_a(first, branpol_a_branches);
            branpol_a.make_longest_parent();
            new_branpols.push_back(branpol_a);
        }
        auto branpol_b_structure(initial_structure);
        branpol_b_structure.push_back(second);
        std::vector<Polymer> branpol_b_branches;
        bool c4 = calculate_all_branches(branpol_b_structure, branpol_b_branches);
        if (!c4) return false;
        if (branpol_b_branches.empty())
        {
            new_polymers.push_back(second);
        } else
        {
            Branched_polymer branpol_b(second, branpol_b_branches);
            branpol_b.make_longest_parent();
            new_branpols.push_back(branpol_b);
        }
        if (division_validation_check(new_branpols, new_polymers)) return true;
    }
    else
        return false;
    return false;
}


bool
Branched_polymer::division_validation_check(std::vector<Branched_polymer> &new_branpols,
                                            std::vector<Polymer> &new_polymers)
{
    auto original_nom = _parent.number_of_monomers();
    for (auto &b : _branches) original_nom += b.number_of_monomers();

    auto npol_nom = 0;
    for (auto &nb : new_branpols)
    {
        npol_nom += nb.get_parent().number_of_monomers();
        auto nbr = nb.get_branches();
        for (auto &nbrp : nbr) npol_nom += nbrp.number_of_monomers();
    }
    for (auto &np :new_polymers) npol_nom += np.number_of_monomers();
    if (original_nom == npol_nom)
    {
        return true;
    }

    return false;
}

bool
Branched_polymer::merge (Branched_polymer &other)
{
    std::vector<Polymer> checking;
    auto other_branches = other.get_branches();
    checking.insert(checking.end(), other_branches.begin(), other_branches.end());
    checking.push_back(other.get_parent());
    for (auto &p : checking) p.calculate_coords();

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    _branches.push_back(_parent);
    for (auto & branch : _branches)
    {
        for (auto &polymer: checking)
        {
            if (branch.polymer_intersecting_from_second(polymer, &bd1, &bd2))
            {
                _branches.pop_back();
                free_bd(&bd1);
                free_bd(&bd2);
                return false;
            }
        }
    }
    _branches.pop_back();

    for (auto &c : checking)
        c.set_bid(_identifier);

    _branches.insert(_branches.end(), checking.begin(), checking.end());

    free_bd(&bd1);
    free_bd(&bd2);
    return true;
}


void  Branched_polymer::make_longest_parent()
{
    std::vector<Polymer> polymers = _branches;
    polymers.push_back(_parent);
    std::vector<int> noms;
    noms.reserve(polymers.size());

    for (auto &p : polymers)
    {
        noms.push_back(p.number_of_monomers());
    }
    auto maxElement = std::max_element(noms.begin(), noms.end());
    auto index = std::distance(noms.begin(), maxElement);
    Polymer parent = polymers[index];

    _parent = parent;
    polymers.erase(polymers.begin() + index);
    _branches = polymers;
}



bool Branched_polymer::is_monomer_occupied(Polymer &activated, int32_t moni) const
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();
    auto polymers = to_polymers();

    for (auto &branch : polymers)
    {
        if (activated.polymer_monomer_intersecting(moni,branch, &bd1, &bd2))
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


bool Branched_polymer::is_intersecting_with(std::vector<Polymer> &polymers)
{
    if (polymers.empty()) return false;

    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    _branches.push_back(_parent);
    for (auto & branch : _branches)
    {
        for (auto &polymer: polymers)
        {
            if (branch == polymer) continue;
            if(branch.polymer_intersecting_from_second(polymer, &bd1, &bd2))
            {
                _branches.pop_back();
                free_bd(&bd1);
                free_bd(&bd2);
                return true;
            }
        }
    }
    _branches.pop_back();
    free_bd(&bd1);
    free_bd(&bd2);
    return false;
}


bool Branched_polymer::is_intersecting_between() const
{
    struct bd bd1{}, bd2{};
    bd1 = allocate_bd();
    bd2 = allocate_bd();

    std::vector<Polymer> polymers = _branches;
    polymers.push_back(_parent);

    for (auto &p : polymers)
    {
        for (auto &o : polymers)
        {
            if (p == o) continue;
            if (p.polymer_check_if_branch(o, &bd1, &bd2))
            {
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


int32_t
Branched_polymer::number_of_polymers() const
{
    return _branches.size() + 1;
}

void
Branched_polymer::_move (float fx, float fy, float fz)
{
    float cxp, cyp, czp;
    _parent.position(cxp, cyp, czp, 0);
    _parent.set_origins(cxp - _parent.monomer(0).length() / 2 + fx,
                        cyp - _parent.monomer(0).width() / 2 + fy,
                        czp + fz);

    for (auto &polymer: _branches)
    {
        float cx, cy, cz;
        polymer.position(cx, cy, cz, 0);
        polymer.set_origins(cx - polymer.monomer(0).length() / 2 + fx,
                            cy - polymer.monomer(0).width() / 2 + fy,
                            cz + fz);
    }
}


size_t
Branched_polymer::_find_position (const Polymer& polymer)
{
    auto iter = std::find(_branches.begin(), _branches.end(), polymer);
    return(std::distance(_branches.begin(), iter));
}


void
Branched_polymer::_move_element_to_back (Polymer &polymer, std::vector<Polymer> &polymers)
{
    if (polymers.size() < 2) return;

    if (polymer == polymers.back()) return;

    auto it = std::find(polymers.begin(), polymers.end(), polymer);
    auto oldIndex = it - polymers.begin();
    auto newIndex = polymers.size() - 1;

    if (oldIndex > newIndex)
        std::rotate(polymers.rend() - oldIndex - 1, polymers.rend() - oldIndex, polymers.rend() - newIndex);
    else
        std::rotate(polymers.begin() + oldIndex, polymers.begin() + oldIndex + 1, polymers.begin() + newIndex + 1);

}


bool
Branched_polymer::_add_unique(std::vector<Polymer> &polymers, Polymer &polymer)
{
    if (std::find(polymers.begin(), polymers.end(), polymer) == polymers.end()) {
        polymers.push_back(polymer);
        return true;
    }
    return false;
}

void
Branched_polymer::_add_uniques(std::vector<Polymer> &polymers, std::vector<Polymer> &others)
{
    for (auto &o: others)
    {
        if (std::find(polymers.begin(), polymers.end(), o) == polymers.end())
        {
            polymers.push_back(o);
        }
    }
}


bool
Branched_polymer::_check_if_contains_already(std::vector<Polymer> &polymers, std::vector<Polymer> &others) const
{
    for (auto &o: others)
    {
        if (std::find(polymers.begin(), polymers.end(), o) != polymers.end())
        {
           return true;
        }
    }
    return false;
}


void
Branched_polymer::_remove_same_element_from_other_vec(std::vector<Polymer> &a, std::vector<Polymer> &b)
{
    a.erase( remove_if( begin(a),end(a),
                        [&](auto x){return find(begin(b),end(b),x)!=end(b);}), end(a) );
}


std::ostream&
operator<<(std::ostream& output, const Branched_polymer& source)
{
    output << std::endl;
    output << "[Branched_Polymer(Instance): " << &source << "]" << std::endl;
    output << "this->_nr: " << source._identifier <<  std::endl;
    output << "this->_branpolymer(parent): " <<std::endl;
    output << "label: "<<source.get_parent().get_label()<<" bid: "<<source.get_parent().get_bid()<<" len: "<<source.get_parent().number_of_monomers()<<std::endl;
    output << "this->_branpolymer(branches): " <<std::endl;
    for (auto &b : source.get_branches()) output<<"label: "<<b.get_label()<< " bid: " <<b.get_bid()<<" len: "<<b.number_of_monomers()<<std::endl;
    output <<std::endl;

    return output;
}


bool
Branched_polymer::make_division_for_swap(Polymer &activated, int32_t moni,
                                std::vector<Branched_polymer> &new_branpols,
                                std::vector<Polymer> &new_polymers, std::vector<Polymer> &specials)
{
    if (activated.number_of_monomers() < 2) return false;
    if (moni < 1) return false;

    std::vector<Polymer> initial_structure = _branches;
    Polymer first, second;
    activated.divide(moni, first, second);

    if (_parent == activated)
    {
        auto branpol_a_structure(initial_structure);
        branpol_a_structure.push_back(first);
        std::vector<Polymer> branpol_a_branches;
        bool c1 = calculate_all_branches(branpol_a_structure, branpol_a_branches);
        if (!c1) return false;
        if (branpol_a_branches.empty())
        {
            new_polymers.push_back(first);
        } else
        {
            Branched_polymer branpol_a(first, branpol_a_branches);
            new_branpols.push_back(branpol_a);
            specials.push_back(first);
        }

        auto branpol_b_structure(initial_structure);
        branpol_b_structure.push_back(second);
        std::vector<Polymer> branpol_b_branches;
        bool c2 = calculate_all_branches(branpol_b_structure, branpol_b_branches);
        if (!c2) return false;
        if (branpol_b_branches.empty())
        {
            new_polymers.push_back(second);
        } else
        {
            Branched_polymer branpol_b(second, branpol_b_branches);
            new_branpols.push_back(branpol_b);
            specials.push_back(second);
        }
        if (division_validation_check(new_branpols, new_polymers)) return true;
    }
    else if (std::find(_branches.begin(), _branches.end(), activated) != _branches.end())
    {
        auto erPol = std::remove(initial_structure.begin(), initial_structure.end(), activated);
        initial_structure.erase(erPol, initial_structure.end());
        initial_structure.push_back(_parent);

        auto branpol_a_structure(initial_structure);
        branpol_a_structure.push_back(first);
        std::vector<Polymer> branpol_a_branches;
        bool c3 = calculate_all_branches(branpol_a_structure, branpol_a_branches);
        if (!c3) return false;
        if (branpol_a_branches.empty())
        {
            new_polymers.push_back(first);
        } else
        {
            Branched_polymer branpol_a(first, branpol_a_branches);
            branpol_a.make_longest_parent();
            new_branpols.push_back(branpol_a);
            specials.push_back(first);
        }

        auto branpol_b_structure(initial_structure);
        branpol_b_structure.push_back(second);
        std::vector<Polymer> branpol_b_branches;
        bool c4 = calculate_all_branches(branpol_b_structure, branpol_b_branches);
        if (!c4) return false;
        if (branpol_b_branches.empty())
        {
            new_polymers.push_back(second);
        } else
        {
            Branched_polymer branpol_b(second, branpol_b_branches);
            branpol_b.make_longest_parent();
            new_branpols.push_back(branpol_b);
            specials.push_back(second);
        }
        if (division_validation_check(new_branpols, new_polymers)) return true;
    }
    else
        return false;
    return false;
}


bool
Branched_polymer::calculate_all_branches_for_stats(std::vector<int> &number, std::vector<float> &length)
{
    auto input = _branches;
    std::vector<Polymer> parents {_parent};
    input.push_back(_parent);

    while (true)
    {
        std::vector<Polymer> temp;
        std::vector<Polymer> results;
        std::vector<int32_t> moni;

        for (auto &pol: parents)
        {
            _move_element_to_back(pol, input);
            if (calculate_branches(input, results, moni))
            {
                for (auto &i : results) std::cout<<i.get_label()<<" ";
                std::cout<<std::endl;

                _add_uniques(temp, results);
                results.clear();
            }
        }
        _remove_same_element_from_other_vec(temp, parents);

        if (temp.empty()) break;
        parents.clear();
        parents = temp;
        number.emplace_back(temp.size());
        float total_length = 0;
        for (auto &t : temp) total_length+= t.number_of_monomers();
        length.push_back(total_length);
    }
    return true;
}