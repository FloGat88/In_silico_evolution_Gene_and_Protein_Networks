//
//  Preliminary.hpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 28/01/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#ifndef Preliminary_hpp
#define Preliminary_hpp


// #define NDEBUG            // must be written before the inclusion of assert.h??
#define _LIBCPP_DEBUG_LEVEL 1       // works only with 1, no larger values. Macro for STL that does e.g. explicit range checking for vector; for array, however, this macro was not build in in STL a priori. Therefore use DebugArray instead if _LIBCPP_DEBUG_LEVEL >= 1 (see include's below) which is an almost exact copy of <array> from STL except that in operator[] range checking (via _LIBCPP_ASSERT) is build just as in vector. In a similar way as for arrays this could also be done for deque and other STL container when needed. So, _LIBCPP_DEBUG_LEVEL >= 1 makes the compiler complain whenever the range of a vector or array is exceeded and is equvalent to replacing all a[x] by a.at(x) for all STL vectors and arrays a :-)


// #define PARALLELIZE
// #define RANDOM_SEED
#define IDENTICAL_INITIAL_PERTURBATION


#include "RandomNumberGenerator.h"
#include "Utility.hpp"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <vector>
#if _LIBCPP_DEBUG_LEVEL >= 1
#include "DebugArray.h"
#include "DebugDeque.h"
#else
#include <array>
#include <deque>        // needed for RecordAncestry in EvoAlg
#endif
#include <cfloat>
#include <list>
#include <algorithm>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/container/static_vector.hpp>
#include <assert.h>


using namespace std;

enum ElementTypes {Prot=1, Phos=2, Comp=3};
enum ArrowTypes {Trans=1, CompReact=2, Dec=3};
enum StateType {cytosolic = 0, membrane_bound = 1, arbitrary = -1};

struct TypesState  {
    boost::container::static_vector<short,3> types;
    short state;
    TypesState(boost::container::static_vector<short,3> types, short state) : types(types),state(state) {};
};

////////////////////////////////////// Parameters //////////////////////////////////////////////////////////////////////////
// Parameters that determine Mutations and general settings (relevant for Genome_Version2)
constexpr bool BackRatesOn = false;
constexpr bool HigherOrderComplexes = true;
constexpr float ConfigInterconnectionProb = 0.25;
constexpr bool OnlyOneEssFuncPerElem = true;        // only one essential function per element
constexpr short NumEssFunc = 1;                         // Number of Essential Functions
const array<TypesState, NumEssFunc> ess_props{TypesState({Phos, Comp}, membrane_bound)};    // remember to adapt the constructor of Genome_Version2 accordingly by creating the initial ess_elements manually!
constexpr bool TrackMutations = true;
constexpr float InitialRate_Linear_Max = 1.0;   // maximum rate that can be chosen in add_arrow for a rate constant that describes a linear process (arrows with number of sources = 1) (minimum rate is 0)
constexpr float InitialRate_Quadratic_Max = 1.0;  // maximum rate that can be chosen in 'add_arrow' for a rate constant that describes a second order  process (arrows with number of sources = 2) (minimum rate is 0)
constexpr float InitialConcentration_Max = 600.;  // [#/mu] maximum initial density that can be chosen in 'Add_Protein' for the initial homogeneous protein concentration (minimum concentration is 0)

// Parameters relevant for the simulation (relevant for Dynamics_and_Fitness and EvoAlg)
constexpr float L_small = 10.0;          // [mu] length of the small system   // alternatively, choose L_small = 1
constexpr float L_large = 30.0;          // length of the large system
//constexpr float D = 0.25;           // diffusion constant
constexpr float D_cytosolic = 10;     // [mu^2/sec] cytosolic diffusion constant  // alternatively, choose D_cytosolic = 1
constexpr float D_membrane = 0.1;      // [mu^2/sec] diffusion constnat on membrane
constexpr double T_max = 20.0;     // maximal time for the simulations of the dynamics to run: choose it such that T_max >> time scale for diffusion over the complete system (L^2/D)
constexpr short N_discretisation = 20;      // number of pieces with which the system sizes are discretized

#ifdef IDENTICAL_INITIAL_PERTURBATION
extern array<float, N_discretisation> Initial_Perturbations;
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct node {
    short type;
    short index;
    
    node() {};
    node(short type, short index) : type(type), index(index) {};
    node(short type, short species, short index) : type(type), index(index) {};
    
    friend bool operator==(const node& node1, const node& node2)  {
        return (node1.type==node2.type && node1.index==node2.index);
    };
    friend bool operator!=(const node& node1, const node& node2)  {
        return !(node1==node2);
    };
    friend ofstream& operator<<(ofstream& file, const node& no)  {
        file << left << setw(4) << no.type << setw(6) << no.index;
        return file;
    };
};

struct arrow {
    short type;
    short index;
    
    arrow() {};
    arrow(short type, short index) : type(type), index(index) {};
    
    friend bool operator==(const arrow& ar1, const arrow& ar2)  {
        return (ar1.type==ar2.type && ar1.index==ar2.index);
    };
    friend bool operator!=(const arrow& ar1, const arrow& ar2)  {
        return !(ar1==ar2);
    };
    friend ofstream& operator<<(ofstream& file, const arrow& ar)  {
        file << left << setw(4) << ar.type << setw(6) << ar.index;
        return file;
    };
};


// define element structs
struct Element {
    bool existent;
    const short type;           // these const members efford an explicit copy constructor and assignment operator. Making them non-constant would allow to remove the explict copy constructor and assignment operator and uns the implicit ones.
    const short index;
    short state;
    short ess_func;       // possible essential function of the element (-1 if none), otherwise ess_func describes the row-index within "essential_elems" in "Essential"
    float init_dens;           // homogeneous initial desity of the element
    node base_element;
    short dup_index;                            // duplicate index:
    array<vector<short>,2> phosphs;
    vector<arrow> ins;
    vector<arrow> outs;
    
    // Constructor
    Element(short type, short index) : type(type),index(index),dup_index(-1) {
        if(type!=Phos)            // if a base element is initialized (i.e. not a Phosphorylate), 'base_element' just references the element itself. This will remain so for the whole runtime of the program and does not need to be changed when reactivating the Element. Only when a Phosphorylate is initialized, the corresponding base element must be assigned when reactivating (this is done in 'add_Element').
            base_element = node(type,index);
        ins.reserve(5);
        outs.reserve(5);
    };
    
    // Copy Constructor (explicit definition needed because of the const members 'type' and 'index')
    Element(const Element& el) : existent(el.existent), type(el.type), index(el.index), state(el.state), ess_func(el.ess_func), init_dens(el.init_dens), base_element(el.base_element), dup_index(-1)
    {
        ins.reserve(5);
        outs.reserve(5);
        if(existent==true)  {
            phosphs = el.phosphs;
            ins = el.ins;
            outs = el.outs;
        }
    };
    
    // Assignment Operator (explicit definition needed because of the const members 'type' and 'index')
    Element& operator =(const Element& el)  {
        existent = el.existent;     // necessary to copy these vectors even if existent==false because otherwise these vectors do not become cleared in the assigned non-existent element (add_element does not clear the vectors and new entries will be push_backed subsequently).
        phosphs = el.phosphs;
        ins = el.ins;
        outs = el.outs;
        assert(type==el.type);
        assert(index==el.index);
        if(existent == true)  {     // no need to copy these elements if the element is not existent anyway
            state = el.state;
            ess_func = el.ess_func;
            init_dens = el.init_dens;
            base_element = el.base_element;
        }
        return *this;
    };
    
    friend bool operator==(const Element& el1, const Element& el2)  {       // operator== needed for Tests: Check_PrintRebuild to compare for equality a genome_version2 and a printed and rebuilt genome_version2.
        if(el1.existent==false && el2.existent == false)        // if both elements are not existent return true
            return true;
        else
            return (el1.existent == el2.existent &&
                    el1.type == el2.type &&
                    el1.index == el2.index &&
                    el1.state == el2.state &&
                    el1.ess_func == el2.ess_func &&
                    abs(el1.init_dens-el2.init_dens) <= 1.0e-5*el2.init_dens &&        // relative difference of 1e-5 or absolut difference of 1e-7 for the floating point type due to restricted precision in the data file
                    el1.base_element == el2.base_element &&
                    is_permutation(el1.phosphs[0].begin(), el1.phosphs[0].end(), el2.phosphs[0].begin()) &&   // phosphs and ins and outs of el1 only need to be permutations of the respective vectors of el2
                    is_permutation(el1.phosphs[1].begin(), el1.phosphs[1].end(), el2.phosphs[1].begin()) &&
                    is_permutation( el1.ins.begin(), el1.ins.end(), el2.ins.begin() ) &&
                    is_permutation( el1.outs.begin(), el1.outs.end(), el2.outs.begin() ));
    };
    
    void PrintElement(ofstream& file) {
        if(existent==false)
            file << left << setw(4) << existent << endl;
        else {
            file << left << setw(4) << existent << setw(6) << state;
            if(type == Prot)
                file << init_dens;
            else if(type == Phos)
                file << base_element;
            file << endl;
        }
    };
    
    friend ifstream& operator>>(ifstream& file, Element& El)  {
        file >> El.existent;
        El.ess_func = -1;       // set essential function to -1 first; maybe this is changed later when the essential elements are read in
        if(El.existent == true)  {
            file >> El.state;
            if(El.type == Prot)     // in case of Protein
                file >> El.init_dens;
            else  {                 // in case of Complex or Phosphorylate
                El.init_dens = 0;   // set initial density to 0
                if(El.type == Phos)
                    file >> El.base_element.type >> El.base_element.index;
            }
        }
        return file;
    };
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Arrow  {
    bool existent;
    const short type;            // these const members efford an explicit copy constructor and assignment operator. Making them non-constant would allow to remove the explict copy constructor and assignment operator and uns the implicit ones.
    const short index;
    float rate;
    bool back;
    float back_rate;
    bool duplicated;
    const short N_source, N_sink;
    array<node,2> source;
    array<node,2> sink;
    
    // Constructor
    Arrow(short type, short index, short N_source, short N_sink) : type(type), index(index), N_source(N_source), N_sink(N_sink), duplicated(false) {};
    
    // Copy Constructor (explicit definition needed because of the const members 'type', 'index', 'N_source' and 'N_sink')
    Arrow(const Arrow& Ar) : existent(Ar.existent), type(Ar.type), index(Ar.index), rate(Ar.rate), back(Ar.back), back_rate(Ar.back_rate), duplicated(false), N_source(Ar.N_source), N_sink(Ar.N_sink), source(Ar.source), sink(Ar.sink) {};
    
    // Assignment operator (explicit definition needed because of the const members 'type', 'index', 'N_source' and 'N_sink')
    Arrow& operator = (const Arrow& Ar)  {
        existent = Ar.existent;
        assert(type == Ar.type);
        assert(index == Ar.index);
        assert(N_source == Ar.N_source);
        assert(N_sink == Ar.N_sink);
        if(existent == true)  {
            rate = Ar.rate;
            back = Ar.back;
            back_rate = Ar.back_rate;
            duplicated = false;
            source = Ar.source;
            sink = Ar.sink;
        }
        return *this;
    };
    
    friend bool operator==(const Arrow& ar1, const Arrow& ar2)  {       // operator== needed for Tests: Check_PrintRebuild to compare for equality a genome_version2 and a printed and rebuilt genome_version2.
        if(ar1.existent==false && ar2.existent == false)        // if both arrows are not existent return true
            return true;
        else  {
            for(short i=0; i<ar1.N_source; ++i)
                if(ar1.source[i] != ar2.source[i])
                    return false;
            for(short i=0; i<ar1.N_sink; ++i)
                if(ar1.sink[i] != ar2.sink[i])
                    return false;
            return (ar1.existent == ar2.existent &&
                    ar1.type == ar2.type &&
                    ar1.index == ar2.index &&
                    abs(ar1.rate-ar2.rate) <= 1.0e-5*ar1.rate &&        // relative difference of 1e-5 or absolute difference of 1e-7 for the floating point type due to restricted precision in the data file
                    ar1.back == ar2.back &&
                    ( BackRatesOn==false || abs(ar1.back_rate-ar2.back_rate) <= 1.0e-5*ar1.back_rate ) &&
                    ar1.N_source == ar2.N_source &&
                    ar1.N_sink == ar2.N_sink);
        }
    };
    
    
    void PrintArrow(ofstream& file) {
        file << left << setw(4) << existent;
        if(existent==true) {
            file << setprecision(7) << setw(16) << rate;
            if(BackRatesOn)
                file << setprecision(7) << setw(16) << back_rate;
            //            else
            //                file << setw(8) << "0.0";
            file << source[0];
            if(N_source==2)
                file << source[1];
            file << sink[0];
            if(N_sink==2)
                file << sink[1];
        }
        file << endl;
    };
    
    
    friend ifstream& operator>>(ifstream& file, Arrow& Ar)  {
        file >> Ar.existent;
        if(Ar.existent == true)  {
            file >> Ar.rate;
            Ar.back = BackRatesOn;
            if(BackRatesOn)
                file >> Ar.back_rate;
            for(short i=0; i<Ar.N_source; ++i)
                file >> Ar.source[i].type >> Ar.source[i].index;
            for(short i=0; i<Ar.N_sink; ++i)
                file >> Ar.sink[i].type >> Ar.sink[i].index;
        }
        return file;
    };
};


// structs that are used to describe the basic interactions in the Ode System (OdeSyst1 and OdeSyst2, respectively)

struct ode1 {
    float rate;
    short ind;
    ode1(float r, short i): rate(r),ind(i) {};
    
    friend bool operator==(const ode1& od1, const ode1& od2)  {
        return (od1.rate==od2.rate && od1.ind==od2.ind);
    };
    friend bool operator!=(const ode1& od1, const ode1& od2)  {
        return !(od1==od2);
    };
};


struct ode2 {
    float rate;
    short ind1;
    short ind2;
    ode2(float r, short i1, short i2): rate(r),ind1(i1),ind2(i2) {};
    
    friend bool operator==(const ode2& od1, const ode2& od2)  {
        return (od1.rate==od2.rate && od1.ind1==od2.ind1 && od1.ind2==od2.ind2);
    };
    friend bool operator!=(const ode2& od1, const ode2& od2)  {
        return !(od1==od2);
    };
    
};


#endif /* Preliminary_hpp */
