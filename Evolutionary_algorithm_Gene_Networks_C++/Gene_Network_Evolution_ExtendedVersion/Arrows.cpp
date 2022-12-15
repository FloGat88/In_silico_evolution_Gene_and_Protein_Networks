//
//  Arrows.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 28/01/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"


Genome_Version2::Arrows::Arrows(short type, Genome_Version2& genome_v2, short reserve_size) : genome_v2(genome_v2), type(type) {
        empty.reserve(reserve_size);
        occupied.reserve(reserve_size);
        variable.reserve(reserve_size);
        rate_variable.reserve(reserve_size);
        arrows.reserve(reserve_size);
    };

const Genome_Version2::Arrows& Genome_Version2::Arrows::operator = (const Arrows& Arrs) {
    empty = Arrs.empty;
    occupied = Arrs.occupied;
    variable = Arrs.variable;
    rate_variable = Arrs.rate_variable;
    arrows = Arrs.arrows;
    return *this;
}

template<short N_source, short N_sink, class SourceArrayType, class SinkArrayType>
void Genome_Version2::Arrows::add_arrow(const SourceArrayType& source, const SinkArrayType& sink, float rate, float back_rate)  {
    short index;
    if(empty.size()>0)  {
        index = empty.back();
        empty.pop_back();
    }
    else  {
        index = arrows.size();
        arrows.push_back(Arrow(type,index,N_source,N_sink));
    }
    occupied.push_back(index);
    variable.push_back(index);
    rate_variable.push_back(index);
    Arrow& new_arrow = arrows[index];
    new_arrow.existent = true;
    //   new_arrow.type = type;
    if(rate==-1)  {
        if(N_source==1)
            rate = UnifRand() * InitialRate_Linear_Max;
        else if(N_source==2)
            rate = UnifRand() * InitialRate_Quadratic_Max;
    }
    new_arrow.rate = rate;
    new_arrow.back = BackRatesOn;
    if(BackRatesOn)  {
        if(back_rate==-1)  {
            if(N_sink==1)
                back_rate = UnifRand() * InitialRate_Linear_Max;
            else if(N_sink==2)
                back_rate = UnifRand() * InitialRate_Quadratic_Max;
        }
        new_arrow.back_rate = back_rate;
    }
    //        new_arrow.duplicated = false;      // one could also define an explicit constructor that sets the members 'type', 'duplicated', 'N_source' and 'N_sink' once and for all so that they do not need to be reassigned each time an arrow gets reactivated.
    //         new_arrow.N_source = N_source;
    //         new_arrow.N_sink = N_sink;
    arrow self_arrow = arrow(type,index);
    for(short i=0; i<N_source; ++i)  {
        new_arrow.source[i] = source[i];
        genome_v2.Build_OutConnection(self_arrow, source[i]);
    }
    for(short i=0; i<N_sink; ++i)  {
        new_arrow.sink[i] = sink[i];
        genome_v2.Build_InConnection(self_arrow, sink[i]);
    }
};
// exlicit template instantiation in order to have the definition of the template function exported to this cpp file
template void Genome_Version2::Arrows::add_arrow<1,1>(const array<node,1>&, const array<node,1>&, float, float);
template void Genome_Version2::Arrows::add_arrow<1,1>(const array<node,2>&, const array<node,2>&, float, float);
template void Genome_Version2::Arrows::add_arrow<1,2>(const array<node,1>&, const array<node,2>&, float, float);
template void Genome_Version2::Arrows::add_arrow<1,2>(const array<node,2>&, const array<node,2>&, float, float);
template void Genome_Version2::Arrows::add_arrow<2,1>(const array<node,2>&, const array<node,1>&, float, float);
template void Genome_Version2::Arrows::add_arrow<2,1>(const array<node,2>&, const array<node,2>&, float, float);


bool Genome_Version2::Arrows::delete_arrow(short index) {                 // delete the arrow with index 'index' from 'arrows', delete its reference from all sources and sinks via 'Cap_InConnection' and 'Cap_OutConnection' (these functions then decide whether the drop of the arrow-connection is lethal for the respective node, i.e. if the node must be destructed as a consequence). If no index is provided (default value -1), choose an index randomly of one of the variable arrows. Return 'false' if no arrow could have been deleted, otherwise return 'true'.
    if(index==-1)  {
        if(variable.size()==0)    // if there is no arrow that can be deleted (if variable.size()==0), then break up and return 'false'.
            return false;
        else
            index = variable[ChooseInd(variable.size())];
    }
    Arrow& ar = arrows[index];                  // create reference on the arrow that shall be deleted to make the following more readable
    if(ar.existent==false)                      // first check if the arrow has already been deleted earlier in a deletion cascade and if so return immediately. Nothing else to do then.
        return false;
    ar.existent = false;
    for(short i=0; i<ar.N_source; ++i)          // delete the references (in outs) from the source nodes to the arrow that is destroyed
        genome_v2.Cap_OutConnection(arrow(ar.type, index), ar.source[i]);
    for(short i=0; i<ar.N_sink; ++i)          // delete the references (in ins) from the sink nodes to the arrow that is destroyed
        genome_v2.Cap_InConnection(arrow(ar.type, index), ar.sink[i]);
    
    VecDel_element(occupied, index);            // cancel the deleted element from the vectors that register the elements and mark it as empty in 'empty'.
    VecDel_element(variable, index);
    VecDel_element(rate_variable, index);
    empty.push_back(index);
    return true;
};


bool Genome_Version2::Arrows::mutate_rate(short index) {
    if(index==-1)  {
        if(rate_variable.size()==0)
            return false;
        else
            index = ChooseInd(rate_variable.size());
    }
    if(BackRatesOn==false)
        arrows[index].rate = arrows[index].rate * UnifRand() * 2;
    else {
        if(ChooseInd(2))
            arrows[index].rate = arrows[index].rate * UnifRand() * 2;
        else
            arrows[index].back_rate = arrows[index].back_rate * UnifRand() * 2;
    }
    return true;
};


void Genome_Version2::Arrows::print_Arrows(ofstream& file)  {
    short num_print = (occupied.size()==0)  ?  0 : *max_element(occupied.begin(),occupied.end()) + 1;
    file << left << num_print << endl;
    for(short i=0; i<num_print; ++i)
        arrows[i].PrintArrow(file);
    file << endl;
};


ifstream& operator>>(ifstream& file, Genome_Version2::Arrows& Arrs)  {
    short index = Arrs.arrows.size();
    short N_source = (Arrs.type == CompReact)  ?  2 : 1;
    short N_sink = (Arrs.type == Dec)  ?  2 : 1;
    
    Arrs.arrows.push_back(Arrow(Arrs.type, index, N_source, N_sink));
    file >> Arrs.arrows.back();
    if(Arrs.arrows.back().existent==true)  {        // register read-in Arrow in the helper vectors (occupied, variable, rate_variable or empty)
        Arrs.occupied.push_back(index);
        Arrs.rate_variable.push_back(index);
        Arrs.variable.push_back(index);
    }
    else
        Arrs.empty.push_back(index);
    
    if(Arrs.arrows.back().existent==true)  {        // if the read-in Arrow is existent, reference it in the 'outs' and 'ins' of its sources and sinks. This can be done here only if all Elements have been read in earlier than the Arrows.
        arrow self = arrow(Arrs.type, index);
        for(short i=0; i<N_source; ++i)  {                  // reference the read-in Arrow in the 'outs' vector(s) of its sources
            Element& Elem = Arrs.genome_v2.Get_Element(Arrs.arrows[index].source[i]);
            Elem.outs.push_back(self);
        }
        for(short i=0; i<N_sink; ++i)   {                   // reference the read-in Arrow in the 'ins' vector(s) of its sinks
            Element& Elem = Arrs.genome_v2.Get_Element(Arrs.arrows[index].sink[i]);
            Elem.ins.push_back(self);
        }
    }
    return file;
};


bool Genome_Version2::Arrows::operator==(const Arrows& Arrs2)  {
    return (type==Arrs2.type  &&
            is_permutation(occupied.begin(), occupied.end(), Arrs2.occupied.begin())  &&
            is_permutation(variable.begin(), variable.end(), Arrs2.variable.begin())  &&
            is_permutation(rate_variable.begin(), rate_variable.end(), Arrs2.rate_variable.begin())  &&
            is_permutation(empty.begin(), empty.end(), Arrs2.empty.begin())  &&
            arrows == Arrs2.arrows);
};

