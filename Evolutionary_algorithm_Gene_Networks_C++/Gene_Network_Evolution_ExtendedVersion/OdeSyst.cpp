//
//  OdeSyst.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 24/03/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"

void Genome_Version2::RangeThrough(function<void(Element&)> func)   {
    auto RangeThrough_BaseElems = [&](vector<Element>& elements, function<void(Element&)> func) {
        for(short i=0; i<elements.size(); ++i)  {
            if(elements[i].existent==true)  {
                func(elements[i]);
                for(short j=0; j<elements[i].phosphs[0].size(); ++j)  {
                    short phosph_index = elements[i].phosphs[0][j];
                    func(Phosphorylates.elements[phosph_index]);
                }
                for(short j=0; j<elements[i].phosphs[1].size(); ++j)  {
                    short phosph_index = elements[i].phosphs[1][j];
                    func(Phosphorylates.elements[phosph_index]);
                }
            }
        }
    };
    
    RangeThrough_BaseElems(Proteins.elements, func);
    RangeThrough_BaseElems(Complexes.elements, func);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Genome_Version2::Make_OdeSyst(vector<short>& states) {
    
    auto Generate_OdeIndex_and_states_BaseElems = [&](vector<Element>& elements, vector<short>& OdeIndex, short& s)   {
        for(short i=0; i<elements.size(); ++i)  {
            if(elements[i].existent==true)  {
                states[s] = elements[i].state;
                OdeIndex[i] = s++;
                for(short j=0; j<elements[i].phosphs[0].size(); ++j)  {
                    states[s] = cytosolic;
                    short phosph_index = elements[i].phosphs[0][j];
                    this->phosph_OdeIndex[phosph_index] = s++;
                }
                for(short j=0; j<elements[i].phosphs[1].size(); ++j)   {
                    states[s] = membrane_bound;
                    short phosph_index = elements[i].phosphs[1][j];
                    this->phosph_OdeIndex[phosph_index] = s++;
                }
            }
        }
    };
    
    auto Generate_OdeIndex_and_states = [&]()   {
        short s=0;
        Generate_OdeIndex_and_states_BaseElems(Proteins.elements, prot_OdeIndex, s);
        Generate_OdeIndex_and_states_BaseElems(Complexes.elements, comp_OdeIndex, s);
    };
    
    
    auto Transfer_to_OdeSyst1 = [&](short N_sink, float rate, const array<node,2>& source, const array<node,2>& sink)  {
        short index_source = Get_OdeIndex(source[0]);
        OdeSyst1[index_source].push_back({ode1(-rate,index_source)});
        OdeSyst1[Get_OdeIndex(sink[0])].push_back({ode1(rate,index_source)});
        if(N_sink==2)
            OdeSyst1[Get_OdeIndex(sink[1])].push_back({ode1(rate,index_source)});
    };
    
    auto Transfer_to_OdeSyst2 = [&](float rate, const array<node,2>& source, const array<node,2>& sink)  {
        short index_source1 = Get_OdeIndex(source[0]);
        short index_source2 = Get_OdeIndex(source[1]);
        OdeSyst2[index_source1].push_back({ode2(-rate,index_source1,index_source2)});
        OdeSyst2[index_source2].push_back({ode2(-rate,index_source1,index_source2)});
        OdeSyst2[Get_OdeIndex(sink[0])].push_back({ode2(rate,index_source1,index_source2)});
    };
    
    
    auto Get_TransformArrow = [&](const Arrow& Ar) {
        Transfer_to_OdeSyst1(1, Ar.rate, Ar.source, Ar.sink);
        if(BackRatesOn)
            Transfer_to_OdeSyst1(1, Ar.back_rate, Ar.sink, Ar.source);
    };
    
    auto Get_CompReactionArrow = [&](const Arrow& Ar)  {
        Transfer_to_OdeSyst2(Ar.rate, Ar.source, Ar.sink);
        if(BackRatesOn)
            Transfer_to_OdeSyst1(2, Ar.back_rate, Ar.sink, Ar.source);
    };
    
    auto Get_DecayArrow = [&](const Arrow& Ar)  {
        Transfer_to_OdeSyst1(2, Ar.rate, Ar.source, Ar.sink);
        if(BackRatesOn)
            Transfer_to_OdeSyst2(Ar.back_rate, Ar.sink, Ar.source);
    };
    /*
     auto Get_Arrow = [&](const Arrow& Ar)  {
     if(Ar.N_source==1)
     Transfer_to_OdeSyst1(Ar.N_sink, Ar.rate, Ar.source, Ar.sink);
     else
     Transfer_to_OdeSyst2(Ar.rate, Ar.source, Ar.sink);
     if(BackRatesOn) {
     if(Ar.N_sink==1)
     Transfer_to_OdeSyst1(Ar.N_source, Ar.back_rate, Ar.sink, Ar.source);
     else
     Transfer_to_OdeSyst2(Ar.back_rate, Ar.sink, Ar.source);
     }
     };
     */
    auto Create_OdeSyst = [&]()  {
        for(short i=0; i<Transforms.arrows.size(); ++i)
            if(Transforms.arrows[i].existent==true)
                Get_TransformArrow(Transforms.arrows[i]);
        
        for(short i=0; i<CompReactions.arrows.size(); ++i)
            if(CompReactions.arrows[i].existent==true)
                Get_CompReactionArrow(CompReactions.arrows[i]);
        
        for(short i=0; i<Decays.arrows.size(); ++i)
            if(Decays.arrows[i].existent==true)
                Get_DecayArrow(Decays.arrows[i]);
    };
    // reserve space for states, OdeSyst1 and OdeSyst2
    short num_elements = Get_NumberElements();
    prot_OdeIndex.resize(Proteins.elements.size());
    comp_OdeIndex.resize(Complexes.elements.size());
    phosph_OdeIndex.resize(Phosphorylates.elements.size());
    states.resize(num_elements);
    OdeSyst1.clear();
    OdeSyst2.clear();
    OdeSyst1.resize(num_elements);
    OdeSyst2.resize(num_elements);
    for(short i=0; i<num_elements; ++i)  {
        OdeSyst1[i].reserve(12);
        OdeSyst2[i].reserve(12);
    }
    
    Generate_OdeIndex_and_states();
    Create_OdeSyst();
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Genome_Version2::Transfer_directly_to_OdeSystCompressed(vector<short>& states, vector<ode1>& OdeSyst1_compressed, vector<short>& OdeSyst1_numbers, vector<ode2>& OdeSyst2_compressed, vector<short>& OdeSyst2_numbers)   {
    
    auto Generate_OdeIndex_and_states_BaseElems = [&](vector<Element>& elements, vector<short>& OdeIndex, short& s)   {
        for(short i=0; i<elements.size(); ++i)  {
            if(elements[i].existent==true)  {
                states[s] = elements[i].state;
                OdeIndex[i] = s++;
                for(short j=0; j<elements[i].phosphs[0].size(); ++j)  {
                    states[s] = cytosolic;
                    short phosph_index = elements[i].phosphs[0][j];
                    this->phosph_OdeIndex[phosph_index] = s++;
                }
                for(short j=0; j<elements[i].phosphs[1].size(); ++j)   {
                    states[s] = membrane_bound;
                    short phosph_index = elements[i].phosphs[1][j];
                    this->phosph_OdeIndex[phosph_index] = s++;
                }
            }
        }
    };
    
    auto Generate_OdeIndex_and_states = [&]()   {
        short s=0;
        Generate_OdeIndex_and_states_BaseElems(Proteins.elements, prot_OdeIndex, s);
        Generate_OdeIndex_and_states_BaseElems(Complexes.elements, comp_OdeIndex, s);
    };
    
    auto Get_ode1 = [&](Arrow& Ar, short sign)  {
        float rate = (sign < 0)  ?  -Ar.rate : Ar.rate;
        return ode1(rate, Get_OdeIndex(Ar.source[0]));
    };
    
    auto Get_ode1_back = [&](Arrow& Ar, short sign)  {
        float back_rate = (sign < 0)  ?  -Ar.back_rate : Ar.back_rate;
        return ode1(back_rate, Get_OdeIndex(Ar.sink[0]));
    };
    
    auto Get_ode2 = [&](Arrow& Ar, short sign)  {
        float rate = (sign < 0)  ?  -Ar.rate : Ar.rate;
        return ode2(rate, Get_OdeIndex(Ar.source[0]), Get_OdeIndex(Ar.source[1]));
    };
    
    auto Get_ode2_back = [&](Arrow& Ar, short sign)  {
        float back_rate = (sign < 0)  ?  -Ar.back_rate : Ar.back_rate;
        return ode2(back_rate, Get_OdeIndex(Ar.sink[0]), Get_OdeIndex(Ar.sink[1]));
    };
    
    
    auto Add_ode_to_OdeSystCompressed = [&](Arrow& Ar, short sign)  {
        if(Ar.type == Trans)  {
            OdeSyst1_compressed.push_back( Get_ode1(Ar,sign) );
            OdeSyst1_numbers.back()++;
            if(BackRatesOn)  {
                OdeSyst1_compressed.push_back( Get_ode1_back(Ar,-sign) );
                OdeSyst1_numbers.back()++;
            }
        }
        else if(Ar.type == CompReact)  {
            OdeSyst2_compressed.push_back( Get_ode2(Ar,sign) );
            OdeSyst2_numbers.back()++;
            if(BackRatesOn)  {
                OdeSyst1_compressed.push_back( Get_ode1_back(Ar,-sign) );
                OdeSyst1_numbers.back()++;
            }
        }
        else  {
            OdeSyst1_compressed.push_back( Get_ode1(Ar,sign) );
            OdeSyst1_numbers.back()++;
            if(BackRatesOn)  {
                OdeSyst2_compressed.push_back( Get_ode2_back(Ar,-sign) );
                OdeSyst2_numbers.back()++;
            }
        }
    };
    
    
    function<void(Element& el)> Transfer_Element_to_OdeSystCompressed = [&](Element& el)  {
        OdeSyst1_numbers.push_back(0);
        OdeSyst2_numbers.push_back(0);
        for(short i=0; i<el.ins.size(); ++i)
            Add_ode_to_OdeSystCompressed(Get_Arrow(el.ins[i]),+1);
        for(short i=0; i<el.outs.size(); ++i)
            Add_ode_to_OdeSystCompressed(Get_Arrow(el.outs[i]),-1);
    };
    
    
    // reserve enough space for the vectors; this is not necessary if static_vectors are used, only reasonable for std::vector
    short NumberElements = Get_NumberElements();
    prot_OdeIndex.resize(Proteins.elements.size());
    comp_OdeIndex.resize(Complexes.elements.size());
    phosph_OdeIndex.resize(Phosphorylates.elements.size());
    states.resize(NumberElements);
    OdeSyst1_numbers.reserve(NumberElements);           // OdeSyst_numbers and OdeSyst_compressed must have size 0. This is assumed to be the case; they are not cleared here.
    OdeSyst2_numbers.reserve(NumberElements);
    OdeSyst1_compressed.reserve(Get_NumberOde1());
    OdeSyst2_compressed.reserve(Get_NumberOde2());
    
    Generate_OdeIndex_and_states();
    RangeThrough(Transfer_Element_to_OdeSystCompressed);
};
