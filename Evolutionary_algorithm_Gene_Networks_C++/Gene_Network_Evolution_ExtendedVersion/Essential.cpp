//
//  Essential.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 28/01/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"


Genome_Version2::Essential_class::Essential_class(Genome_Version2& genome_v2) : genome_v2(genome_v2) {};

const Genome_Version2::Essential_class& Genome_Version2::Essential_class::operator = (const Essential_class& Ess) {    // explicit assignment operator; is needed because of the reference genome_v2 in the struct
    ess_elements = Ess.ess_elements;
    return *this;
};

void Genome_Version2::Essential_class::declare_Essential(node nd, short ess_func)  {
    ess_elements[ess_func].push_back(nd);
    genome_v2.Get_Element(nd).ess_func = ess_func;
};


short Genome_Version2::Essential_class::add_Essential(short ess_func)  {
    auto choose_for_Essential = [&](short ess_Func)  {      // uniformly choose an element with the specific porperties of the function 'ess_Func' that will then be attributed this function. If none can be found return node(-1,-1).
        boost::container::static_vector<short,3> types = ess_props[ess_Func].types;
        short state = ess_props[ess_Func].state;
        random_shuffle(types.begin(), types.end());
        for(short i=0; i<types.size(); ++i)  {
            Elements& elems = genome_v2.ElementsDistr(types[i]);
            for(short j=0; j<elems.elements.size(); ++j)
                if(elems.elements[j].existent==true)
                    if(OnlyOneEssFuncPerElem==false || elems.elements[j].ess_func==-1)     
                        if(state==-1 || elems.elements[j].state==state)
                            return node(types[i],j);
        }
        return node(-1,-1);
    };
    /*
     if(types[0]==Prot) return node(type,genome_v2.Add_Protein(state));
     else if(types[0] == Comp) return node(type,genome_v2.Add_Binding_Domain(state, false));
     else if(types[0]==Phosph) return node(type,genome_v2.Add_Phosphorylate(node(-1,-1),state));
     */
    
    if(ess_func==-1)
        ess_func = ChooseInd(ess_elements.size());
    node nd = choose_for_Essential(ess_func);
    if(nd.type != -1)  {
        declare_Essential(nd, ess_func);
        return 1;
    }
    else
        return -1;
};

void Genome_Version2::Essential_class::undeclare_Essential(node nd, short ess_func)  {
    VecDel_element(ess_elements[ess_func], nd);
    genome_v2.Get_Element(nd).ess_func = -1;
    if(ess_elements[ess_func].size()==0)    // if no essential element is left for this function, add a new one for this function
        add_Essential(ess_func);
};

short Genome_Version2::Essential_class::delete_Essential(short ess_func)  {
    auto ChooseEssInd = [&]()  {        // choose uniformly one of the essential elements for deletion; return -1 if none can be found because no one is existing
        short NumEssElems = 0;
        for(short i=0; i<ess_elements.size(); ++i)
            NumEssElems += ess_elements[i].size();
        short d = ChooseInd(NumEssElems);
        for(short i=0; i<ess_elements.size(); ++i)  {
            d -= ess_elements[i].size();
            if(d<0)
                return i;
        }
    //    assert(0==1);
        cout << "No existing essential element" << endl;
        return (short)-1;
    };
    if(ess_func == -1)
        ess_func = ChooseEssInd();
    if(ess_func != -1)  {
        node nd = ess_elements[ess_func][ChooseInd(ess_elements[ess_func].size())];
        undeclare_Essential(nd, ess_func);
        return 1;
    }
    else
        return -1;
};

void Genome_Version2::Essential_class::print_Essential(ofstream& file)  {
    for(short i=0; i<NumEssFunc; ++i)  {
        if(ess_elements[i].size()==0)
            file << -1;
        else
            for(short j=0; j<ess_elements[i].size(); ++j)
                file << ess_elements[i][j];
        file << endl;
    }
};

ifstream& operator>>(ifstream& file, Genome_Version2::Essential_class& Ess)  {
    short type, index;
    string line;
    
    for(short i=0; i<NumEssFunc; ++i)  {
        getline(file, line);
        istringstream iss(line);
        while (iss >> type)  {
            if(type!=-1) {                  // if no essential element exists, the line starts with -1 to indicate empty line
                iss >> index;
                Ess.ess_elements[i].push_back(node(type, index));
            }
        }
    }
    return file;
};

void Genome_Version2::reference_essFunc()  {
    for(short i=0; i<NumEssFunc; ++i)
        for(short j=0; j<Essential.ess_elements[i].size(); ++j)
            Get_Element(Essential.ess_elements[i][j]).ess_func = i;
};





