//
//  Elements.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 28/01/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"

// constructor
Genome_Version2::Elements::Elements(short type, Genome_Version2& genome_v2, short reserve_size) : genome_v2(genome_v2), type(type) {
        empty.reserve(reserve_size);
        occupied.reserve(reserve_size);
        variable.reserve(reserve_size);
        elements.reserve(reserve_size);
};
// assignment operator. Explicit definition needed??
const Genome_Version2::Elements& Genome_Version2::Elements::operator = (const Elements& Elems) {
    empty = Elems.empty;
    occupied = Elems.occupied;
    variable = Elems.variable;
    elements = Elems.elements;
    return *this;
}

short Genome_Version2::Elements::add_element(short state, float init_dens, const node& base_element, bool SuppressReferencing_in_BaseElem)  {     // 'base_element' and 'SuppressReferencing_in_BaseElem' are only relevant for Phosphorylates. In that case, if 'SuppressReferencing_in_BaseElem' is true, the index of the Phosphorylate is not added to 'phosphs' of the base element. This is useful e.g. in 'Add_Phosphorylate' for practical reasons, however should not be used e.g. in 'Duplicate_Protein'.
    short index;
    if(empty.size()>0)  {
        index = empty.back();
        empty.pop_back();
    }
    else  {
        index = elements.size();
        elements.push_back(Element(type,index));
    }
    occupied.push_back(index);
    variable.push_back(index);
    Element& new_element = elements[index];
    new_element.existent = true;
    //    new_element.type = type;
    new_element.state = state;
    new_element.ess_func = -1;          // set ess_func to -1. 'ess_func' is assigned if needed by 'declareEssential' of the class 'Essential'
    new_element.init_dens = init_dens;
    if(type==Phos)  {      // if the element is a Phosphorylate, also assign 'base_element' (for Phosphorylates, the argument 'base_element' is compulsory!). The reference to the index of the Phosphorylate in 'phosphs' of the base_element is esablished only if 'SuppressReferencing_in_BaseElem' is set to false.
        new_element.base_element = base_element;
        if(SuppressReferencing_in_BaseElem == false)
            genome_v2.Build_PhosphRef(state, index, base_element);
    }
    assert(new_element.ins.size()==0);          // check if all vectors are cleared!
    assert(new_element.outs.size()==0);
    assert(new_element.phosphs[0].size()==0);
    assert(new_element.phosphs[1].size()==0);
    
    /*        if(type!=Phos)
     new_element.base_element = node(type,index);
     else
     new_element.base_element = base_element;  */
    //    new_element.dup_index = -1;
    //    new_element.ins.reserve(5);
    //    new_element.outs.reserve(5);
    return index;
};
    
short Genome_Version2::Elements::delete_element(short index) {        // delete the element of index 'index' from 'elements'. If no index is provided (default value -1), choose an index randomly of one of the variable elements. Return 'false' if no element could have been deleted, otherwise return 'true'.
    if(index==-1)   {      // if no index has been provided (default argument -1) choose one index at random.
        if(variable.size()==0)    // if there is no element that can be deleted (if variable.size()==0), then break up and return 'false'.
            return false;
        else
            index = variable[ChooseInd(variable.size())];
    }
    
    if(elements[index].existent==false)                  // first check if the arrow has already been deleted earlier in a deletion cascade and if so return immediately. Nothing else to do then.
        return -1;
    
    Element& el = elements[index];          // create reference on the respective element for better handability
    el.existent = false;               // first thing to do is set 'existent' to 'false'
    
    vector<arrow> ins_copy = el.ins;            // make a copy of the ins and outs vectors first because without making a copy it's hard to make sure that you really got all in and out arrows or you might end up in an infinite loop since the original vectors may be modified while the loop is running by the deletion cascade.
    for(short i=ins_copy.size()-1; i>=0; --i)     // delete all influx-arrows (stored in 'ins_copy')
        genome_v2.Delete_Arrow(ins_copy[i]);
    vector<arrow> outs_copy = el.outs;
    for(short i=outs_copy.size()-1; i>=0; --i)     // delete all outflux-arrows (stored in 'outs_copy')
        genome_v2.Delete_Arrow(outs_copy[i]);
    
    if(type == Prot || type == Comp)           // if the element to be deleted is either a Protein or a Complex, delete all Phosphorylates thereof (whoes indices are stored in 'Phosphorylates') as these Phosphorylates can't have any further influxes and thus cannot participate in the dynamics any more with a constant concentration of 0 (zombies).
        for(short i=0; i<2; ++i)
            while(el.phosphs[i].size()>0)
                genome_v2.Delete_Node(node(Phos,el.phosphs[i].back()));
    if(type == Phos)                      // on the other hand, if the element to be deleted is a Phosphorylate, then delete also its reference from the base element (stored in 'phosphs'). Its own Phosphorylates, however, must not be destroyed here as they might still be interlinked with the rest of the network.
        genome_v2.Delete_PhosphRef(el.state, index, el.base_element);
    el.ins.clear();         // delete the vectors of the element
    el.outs.clear();
    el.phosphs[0].clear();
    el.phosphs[1].clear();
    
    VecDel_element(occupied, index);      // remove reference on the element from 'occupied' and 'variable' and register the new empty entry in 'empty'.
    VecDel_element(variable, index);
    empty.push_back(index);
    if(el.ess_func >= 0)            // if the deleted element had an essential function (as indicated by 'ess_func' being >= 0), undeclare this function.
        genome_v2.Essential.undeclare_Essential(node(type,index), el.ess_func);
    /*            if(index != elements.size()-1)
     empty.push_back(index);   */
    return index;           // return index of the deletet element
};

void Genome_Version2::Elements::cap_InConnection(const arrow ar, short index)  {  // delete arrow 'ar' from 'ins' of the element with index 'index'. If a CompReaction is destroyed or the last element from 'ins' is being destroyed, trigger deletion of the element as it would otherwise be a zombie node.
    VecDel_element(elements[index].ins, ar);        // delete the corresponding reference on the arrow from 'ins' of the respective element.
    if( ar.type == CompReact || (elements[index].ins.size() == 0 && type!=Prot) )  // in case a CompReaction arrow is destroyed as InConnection to the node, delete the whole node as then there will be no flow into that node (which must hence be a complex) anyways. Also destroy the node if no in-fluxes are left and the element is not an Protein (that has an initial concentration and thus could still influence the dynamics). By these deletions avoid the emergence of zombies (unconnected network parts that do not contribute to the dynamics (b/c all concentrations remain zero))
        delete_element(index);
};

void Genome_Version2::Elements::cap_OutConnection(const arrow ar, short index) {     // delete arrow 'ar' from 'outs' of the element with index 'index'
    VecDel_element(elements[index].outs, ar);
};

void Genome_Version2::Elements::build_InConnection(const arrow ar, short index)  {  // add reference to arrow 'ar' to 'ins' of the element with index 'index'
    elements[index].ins.push_back(ar);
};

void Genome_Version2::Elements::build_OutConnection(const arrow ar, short index)  {   // add reference to arrow 'ar' to 'outs' of the element with index 'index'
    elements[index].outs.push_back(ar);
};

void Genome_Version2::Elements::delete_PhosphRef(const short ref_state, const short ref_index, const short elem_index)  {   // delete reference to the Phosphorylate with index 'ref_index' from the element with index 'elem_index' (stored in 'phosphs').
    VecDel_element(elements[elem_index].phosphs[ref_state], ref_index);
};

void Genome_Version2::Elements::print_Elements(ofstream& file)  {
    short num_print = (occupied.size()==0)  ?  0 : *max_element(occupied.begin(),occupied.end()) + 1;
    file << left << num_print << endl;
    for(short i=0; i<num_print; ++i)
        elements[i].PrintElement(file);
    file << endl;
};


bool Genome_Version2::Elements::operator==(const Elements& Elems2)  {
    return (type==Elems2.type  &&
            is_permutation(occupied.begin(), occupied.end(), Elems2.occupied.begin())  &&
            is_permutation(variable.begin(), variable.end(), Elems2.variable.begin())  &&
            is_permutation(empty.begin(), empty.end(), Elems2.empty.begin())  &&
            elements == Elems2.elements);
};


ifstream& operator>>(ifstream& file, Genome_Version2::Elements& Elems)  {
    short index = Elems.elements.size();
    Elems.elements.push_back(Element(Elems.type, index));
    file >> Elems.elements.back();
    if(Elems.elements.back().existent==true)  {
        Elems.occupied.push_back(index);
        Elems.variable.push_back(index);
    }
    else
        Elems.empty.push_back(index);
    
    return file;
};


void Genome_Version2::reference_Phosphorylates()  {
    for(Element& El : Phosphorylates.elements)
        if(El.existent==true)  {
            Element& BaseElem = Get_Element(El.base_element);
            BaseElem.phosphs[El.state].push_back(El.index);
        }
};
