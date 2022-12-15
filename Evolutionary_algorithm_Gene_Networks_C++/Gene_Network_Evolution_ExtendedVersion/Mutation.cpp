//
//  Mutation.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 29/03/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"



bool Genome_Version2::Mutate_Rate() {
    short weights[3] = {short(Transforms.rate_variable.size()), short(CompReactions.rate_variable.size()), short(Decays.rate_variable.size())};
    short d = ChooseInd_respective_weights(weights, 3);
    if(d==-1)
        return -1;
    
    switch(d)  {
        case 0:
            Transforms.mutate_rate();
            Save_Mutation(1,1);
            break;
        case 1:
            CompReactions.mutate_rate();
            Save_Mutation(1,2);
            break;
        case 2:
            Decays.mutate_rate();
            Save_Mutation(1,3);
            break;
    }
    return +1;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Add_Phosphorylate(node node_parent, short state, bool link_nodes_randomly, vector<node> nodes_linked)   {     // default: node node_parent = node(-1,-1), short state = -1, bool link_nodes_randomly = true, vector<node> nodes_linked={}
    if(Proteins.occupied.size()==0)              // check if the network is non-trivial (at least one Protein exists). Otherwise break up here and indicate that the mutation has not been performed by returning 'false'.
        return -1;
    if(node_parent.type==-1) {        // if no parent element is provided (default value -1) determine the parent uniformly among all elements
        short weights[3] = {short(Proteins.occupied.size()), short(Phosphorylates.occupied.size()), short(Complexes.occupied.size())};
        short type_parent = node_parent.type = 1 + ChooseInd_respective_weights(weights, 3);
        node_parent.index = ElementsDistr(type_parent).occupied[ ChooseInd(ElementsDistr(type_parent).occupied.size()) ];
    }
    
    node node_BaseElement = Get_Element(node_parent).base_element;    // get the corresponding base element from the parent
    
    if(state==-1)          // if no state for the phosphorylate has been provided (default value -1) choose the state randomly as either  cytosolic (0) or membrane_bound (1).
        state = ChooseInd(2);
    
    short index_phosph = Phosphorylates.add_element(state, 0, node_BaseElement, true);   // make a node object for the new Phosphorylate to call the add_arrow functions.
    node node_phosph(Phos,index_phosph);
    
    Transforms.add_arrow<1,1>( array<node,1>({node_parent}), array<node,1>({node_phosph}) );  // create Transform-Arrow for Phosphorylation
    if(BackRatesOn == false)   // if back reaction is not automatically included because BackRatesOn is set to false, also add a Transform-Arrow for the back reaction (Dephosphorylation)
        Transforms.add_arrow<1,1>( array<node,1>({node_phosph}), array<node,1>({node_parent}) );
    
    array<vector<short>,2>& phosphs = Get_Element(node_BaseElement).phosphs;  // create a reference on the phosphs vector of the base element to make the following more readable.
    
    if(link_nodes_randomly)   {     // if the nodes that are linked to the new Phosphorylate shall be chosen randomly, first clear the vector 'nodes_linked' (just to be sure) and then decide for each configuration of the element (all phosphorylates (saved in 'phosphs') and also the base element) EXCEPT for the parent (which is already linked) whether the shall be linked to the new Phosphorylate (establish a link to each of these elements with probability 'ConfigInterconnectionProb'). Check first whether the Phosphorylate or BaseElement is the parent!
        nodes_linked.clear();
        nodes_linked.reserve(6);         // 6 is just a guess for the close-to-maximum number of linked nodes
        for(short i=0; i<2; ++i)
            for(short j=0; j<phosphs[i].size(); ++j)
                if(node_parent.type != Phos || node_parent.index != phosphs[i][j])    // check first whether the respective Phosphorylate is not the parent which is already linked.
                    if(UnifRand() < ConfigInterconnectionProb)    // if not, link the Phosphorylate to the new Phosphorylate with probablity ConfigInterconnectionProb.
                        nodes_linked.push_back( node(Phos,phosphs[i][j]) );
        
        if(node_BaseElement != node_parent && UnifRand() < ConfigInterconnectionProb)    // also for the BasisElement check if it is not the parent and if not link it with prob ConfigInterconnectionProb
            nodes_linked.push_back(node_BaseElement);
    }
    
    for(short i=0; i<nodes_linked.size(); ++i)   {      // create the Transform and back-transform arrows to link the new Phosphorylate to the elements in nodes_linked
        Transforms.add_arrow<1,1>( array<node,1>({node_phosph}), array<node,1>({nodes_linked[i]}) );
        if(BackRatesOn == false)      // if BackRatesOn = false add the back Arrow manually, otherwise it is already included in the Tranform arrow
            Transforms.add_arrow<1,1>( array<node,1>({nodes_linked[i]}), array<node,1>({node_phosph}) );
    }
    
    phosphs[state].push_back(node_phosph.index);    // finally, add the index of the new Phosphorylate to the referenced vector 'phosphs' of the Base Element. Note that state is either cytosolic = 0  or membrane_bound = 1 and hence can be used directly as index.
    Save_Mutation(index_phosph,7);
    return index_phosph;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Add_Phosphorylate_randomized()   {    // the same as 'Add_Phosphorylate' when called without arguments (all arguments are on their default value and are thus chosen randomly)
    if(Proteins.occupied.size()==0)              // check if the network is non-trivial (at least one Protein exists). Otherwise break up here and indicate that the mutation has not been performed by returning 'false'.
        return -1;
    short weights[3] = {short(Proteins.occupied.size()), short(Phosphorylates.occupied.size()), short(Complexes.occupied.size())};
    short type_parent = 1 + ChooseInd_respective_weights(weights, 3);
    short index_parent = ElementsDistr(type_parent).occupied[ ChooseInd(ElementsDistr(type_parent).occupied.size()) ];
    node node_parent(type_parent, index_parent);
    
    node node_BaseElement = Get_Element(node_parent).base_element;
    
    short state = ChooseInd(2);
    
    short index_phosph = Phosphorylates.add_element(state, 0, node_BaseElement, true);
    node node_phosph(Phos,index_phosph);
    
    Transforms.add_arrow<1,1>( array<node,1>({node_parent}), array<node,1>({node_phosph}) );
    if(BackRatesOn == false)
        Transforms.add_arrow<1,1>( array<node,1>({node_phosph}), array<node,1>({node_parent}) );
    
    array<vector<short>,2>& phosphs = Get_Element(node_BaseElement).phosphs;
    for(short i=0; i<2; ++i)
        for(short j=0; j<phosphs[i].size(); ++j)
            if(node_parent.type != Phos || node_parent.index != phosphs[i][j])
                if(UnifRand() < ConfigInterconnectionProb)  {
                    node node_linked = node(Phos,phosphs[i][j]);
                    Transforms.add_arrow<1,1>( array<node,1>({node_phosph}), array<node,1>({node_linked}) );
                    if(BackRatesOn == false)
                        Transforms.add_arrow<1,1>( array<node,1>({node_linked}), array<node,1>({node_phosph}) );
                }
    if(node_BaseElement != node_parent && UnifRand()<ConfigInterconnectionProb)  {
        Transforms.add_arrow<1,1>( array<node,1>({node_phosph}), array<node,1>({node_BaseElement}) );
        if(BackRatesOn == false)
            Transforms.add_arrow<1,1>( array<node,1>({node_BaseElement}), array<node,1>({node_phosph}) );
    }
    phosphs[state].push_back(node_phosph.index);
    Save_Mutation(index_phosph,7);
    return index_phosph;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Add_Protein()  {
    float init_dens = UnifRand()*InitialConcentration_Max;
    short index_prot = Proteins.add_element(cytosolic, init_dens);
    Save_Mutation(index_prot, 11);
    return index_prot;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Genome_Version2::Duplicate_SubNetwork(short dup_Prot, vector<short>& dup_Comps, vector<short>& dup_Phosphs)  {   // function very similar to 'Duplicate_Protein' but can duplicate an arbitrary subnetwork consisting of the elements given in 'dup_BaseElems' and 'dup_Phophs' (the corresponding connecting arrows are duplicated automatically). Relys on that each element appears at most once in the two vectors. Function duplicates the listed elements and all arrows that connect to those. For the details of the implementation of the nested functions look in the description of 'Duplicate_Protein'.
    vector<arrow> dup_arrows;
    dup_arrows.reserve(10);
    
    auto collect_arrows_for_duplication = [&](vector<arrow>& arrow_vec)   {
        for(short i=0; i<arrow_vec.size(); ++i)   {
            Arrow& Ar = Get_Arrow(arrow_vec[i]);
            if(Ar.duplicated==false)  {
                dup_arrows.push_back(arrow_vec[i]);
                Ar.duplicated = true;
            }
        }
    };
    
    auto Duplicate_Element = [&](node nd) {
        Element& el = Get_Element(nd);
        if(el.dup_index!=-1)
            return;
        if(nd.type == Prot)           // nd.type == el.type, but with nd.type there is the possibility for compiler optimization as nd.type is already known at compiletime
            el.dup_index = Proteins.add_element(el.state,el.init_dens);
        else if(nd.type == Comp)
            el.dup_index = Complexes.add_element(el.state, 0);
        else if(nd.type == Phos)   {
            Element& base_el = Get_Element(el.base_element);
            el.dup_index = Phosphorylates.add_element(el.state, 0, node(base_el.type ,base_el.dup_index) );
        }
        if(el.ess_func >= 0)
            Essential.declare_Essential(node(el.type,el.dup_index), el.ess_func);
        collect_arrows_for_duplication(el.ins);
        collect_arrows_for_duplication(el.outs);
    };
    
    auto Duplicate_Arrow = [&](arrow ar)  {
        Arrow& Ar = Get_Arrow(ar);      // can use reference here because add_arrow is called at the very end and not intermediately, so no problem that Ar is used after reallocation has shifted the vector to another position.
        array<node,2> source = Ar.source;
        array<node,2> sink = Ar.sink;
        for(short i=0; i<Ar.N_source; ++i)   {
            short dup_index = Get_Element(source[i]).dup_index;
            if(dup_index != -1)
                source[i].index = dup_index;
        }
        for(short i=0; i<Ar.N_sink; ++i)   {
            short dup_index = Get_Element(sink[i]).dup_index;
            if(dup_index != -1)
                sink[i].index = dup_index;
        }
        if (Ar.type == Trans)
            Transforms.add_arrow<1,1>(source, sink, Ar.rate, Ar.back_rate);
        else if(Ar.type == CompReact)
            CompReactions.add_arrow<2,1>(source, sink, Ar.rate, Ar.back_rate);
        else if(Ar.type == Dec)
            Decays.add_arrow<1,2>(source, sink, Ar.rate, Ar.back_rate);
    };
    
    auto Reset_All_dup_Indices = [&]() {     // reset all 'dup_indices' of Elements stored in dup_BaseElems and dup_Phosphs to -1 and 'duplicated' of Arrows stored in dup_arrows to 'false' to mark them as unduplicated for the next mutation.
        
        Proteins.elements[dup_Prot].dup_index = -1;
        for(short i=0; i<dup_Comps.size(); ++i)
            Complexes.elements[dup_Comps[i]].dup_index = -1;
        for(short i=0; i<dup_Phosphs.size(); ++i)
            Phosphorylates.elements[dup_Phosphs[i]].dup_index = -1;
        
        for(short i=0; i<dup_arrows.size(); ++i)
            Get_Arrow(dup_arrows[i]).duplicated = false;
    };
    
    // reserve enough space in the 'elements' vectors in order to duplicate all needed elements
    Proteins.elements.reserve(Proteins.occupied.size()+1);
    Complexes.elements.reserve(Complexes.occupied.size()+dup_Comps.size());
    Phosphorylates.elements.reserve(Phosphorylates.occupied.size()+dup_Phosphs.size());
    
    // duplicate elements and register the to-be-duplicated arrows in dup_arrows
    Duplicate_Element( node(Prot,dup_Prot) );
    
    for(short i=0; i<dup_Comps.size(); ++i)
        Duplicate_Element( node(Comp,dup_Comps[i]) );
    
    for(short i=0; i<dup_Phosphs.size(); ++i)
        Duplicate_Element( node(Phos,dup_Phosphs[i]) );
    
    // duplicate the arrows
    for(short i=0; i<dup_arrows.size(); ++i)
        Duplicate_Arrow(dup_arrows[i]);
    
    // reset all 'Element.dup_index' to -1 and 'Arrow.duplicated' to 'false'
    Reset_All_dup_Indices();
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Duplicate_Protein_V1(short index)  {    // default: index = -1;  // duplicates the Protein with index 'index' or chooses a random index if no index is provided (default value -1) and duplicates all elements that form or transform from the Protein (the builders).
    if(Proteins.occupied.size() == 0)  // if no protein exists in the network, break up and return 'false'.
        return -1;
    
    vector<node> dup_BaseElems;   // saves references on the base elements (Proteins and Complexes) that must be duplicated in the process (the Protein itself and all its or its Phosphorylates' Complexes and so on)
    vector<node> dup_Phosphs; // saves references on the Phosphorylates that must be duplicated in the process (all Phosphorylates of the respective Protein and of its or the Phosphorylates' Complexes and so on)
    vector<arrow> dup_arrows;  // saves references on the Arrows that must be duplicated in the process (all ins- and outs- Arrows of all the duplicated elements)
    
    dup_BaseElems.reserve(6);  // rough estimate of the needed space for the above vectors in order to avoid unnecessary reallocations.
    dup_Phosphs.reserve(15);
    dup_arrows.reserve(20);
    
    // make sure enough space is available in the 'elements' vectors so that reallocation is not necessary while the function performs. Reallocation would cause trouble with temporary references and pointers in the lambda functions that would point to deallocated memory when the vectors are shifted during the lifetime of the pointers. For Complexes and Phosphorylates, reserve the double amount of space so that with certainty all duplicated element can be accomodated.
    Proteins.elements.reserve(Proteins.occupied.size()+1);
    Complexes.elements.reserve(2*Complexes.occupied.size());
    Phosphorylates.elements.reserve(2*Phosphorylates.occupied.size());
    
    function<void(Element& el)> Duplicate_Element; // declare function objects for all the nested functions needed in the following in order to be able to handle them more easily (in particular necessary in order to be able to call 'Duplicate_BaseEl_and_Dependants' and 'Duplicate_Phosph_and_Dependants' iteratively and mutually one another)
    function<void(const arrow ar)> Duplicate_Arrow;
    function<void(vector<arrow>& arrow_vec)> Collect_Arrows_for_Duplication;
    function<void(const node no)> Duplicate_BaseEl_and_Dependants;
    function<void(const node no)> Duplicate_Phosph_and_Dependants;
    function<void(void)> Reset_All_dup_Indices;
    
    Duplicate_Element = [&](Element& el) {     // takes reference on an Element and creates a duplicate of that element (with the same state) and saves the index of the duplicate in dup_index of the original element. If a Phosphorylate is duplicated, also save the duplicate of the base element as 'base_element' in the new duplicated Phosphorylate (relies on that the base element has been duplicated before and its dup_index contains the index of the duplicated base element). Copies the essential function of the original element (in case there is one) also to the duplicate. Furthermore, registers all connecting Arrows ('ins' and 'outs') in 'dup_arrows'.
        if(el.dup_index!=-1)   // if the element is already marked as duplicated (dup_index > 0) return immediatly
            return;
        if(el.type == Prot)
            el.dup_index = Proteins.add_element(el.state, el.init_dens);
        else if(el.type == Comp)
            el.dup_index = Complexes.add_element(el.state, 0);
        else if(el.type == Phos)   {
            Element& base_el = Get_Element(el.base_element);
            el.dup_index = Phosphorylates.add_element(el.state, 0, node(base_el.type ,base_el.dup_index) );
        }
        if(el.ess_func >= 0)        // if the element has an essential function, as indicated by 'ess_func' bearing an index larger or equal 0, declare the same function also for the duplicate.
            Essential.declare_Essential(node(el.type,el.dup_index), el.ess_func);
        Collect_Arrows_for_Duplication(el.ins);   // register the connecting in-and out-arrows in 'dup_arrows'.
        Collect_Arrows_for_Duplication(el.outs);
    };
    
    Collect_Arrows_for_Duplication = [&](const vector<arrow>& arrow_vec)   {    // takes a vector of arrows (as reference), usually either 'ins' or 'outs' of an element and adds the arrow references to the vector 'dup_arrows' after checking that they have not already been added ('duplicated' == false) (otherwise, each Arrow would be added so many times in the process as it has sources and sinks). If the address has succesfully been transferred to 'dup_arrows', set 'duplicated' to true and thus mark the Arrow as duplicated.
        for(short i=0; i<arrow_vec.size(); ++i)   {
            Arrow& Ar = Get_Arrow(arrow_vec[i]);
            if(Ar.duplicated==false)  {
                dup_arrows.push_back(arrow_vec[i]);
                Ar.duplicated = true;
            }
        }
    };
    
    Duplicate_Arrow = [&](const arrow ar)  {    // takes an arrow and duplicates the corresponding Arrow such that its duplicate connects to the duplicates of the sources and sinks if these have been duplicated and to the original sources and sinks otherwise (efforts that the duplication process for the elements has already been finished and each (original) element contains the index of its duplicated version in 'dup_index').
        Arrow& Ar = Get_Arrow(ar);
        array<node,2> _source = Ar.source;   // copy (by value) the arrays 'source' and 'sink' of the respective Arrow. Then check for each entry if a duplicate exists (Element.dup_index != -1) in which case the respective 'node.index' in source or sink is modified to the duplicate, otherwise, leave it as it is.
        array<node,2> _sink = Ar.sink;
        for(short i=0; i<Ar.N_source; ++i)   {
            short dup_index = Get_Element(_source[i]).dup_index;
            if(dup_index != -1)
                _source[i].index = dup_index;
        }
        for(short i=0; i<Ar.N_sink; ++i)   {
            short dup_index = Get_Element(_sink[i]).dup_index;
            if(dup_index != -1)
                _sink[i].index = dup_index;
        }
        if (Ar.type == Trans)   // finally distinguish the type of the arrow and accordingly call the template function 'add_arrow' with the respective template parameters N_source and N_sink and with the modified 'source' and 'sink'.
            Transforms.add_arrow<1,1>(_source, _sink, Ar.rate, Ar.back_rate);        // back_rate will be ignored if BackRates is flase
        else if(Ar.type == CompReact)
            CompReactions.add_arrow<2,1>(_source, _sink, Ar.rate, Ar.back_rate);
        else if(Ar.type == Dec)
            Decays.add_arrow<1,2>(_source, _sink, Ar.rate, Ar.back_rate);
    };
    
    Duplicate_BaseEl_and_Dependants = [&](const node no)  {   // the function takes a node that MUST be a base element. If not already duplicated (dup_index == -1), this base element is duplicated and registered in 'dup_BaseElems' (by calling 'Duplicate_Element' also all its in- and out-arrows (connecting arrows) are registered in the vector dup_arrows). Afterwards, the function iteratively calls itself and 'Duplicate_Phosph_and_Dependants' (which is similar but takes a Phosphorylate as input), respectively, for the dependants, i.e. the Phosphorylates and Complexes of the base element, which leads to a duplication also of those elements. This therefore triggers a cascade in who's process all the builders (all elements that form or transform from the starting base element) are finally duplicated (i.e. also the Phosphorylates and Complexes of the Complexes of the base element and the Complexes of the Phosphorylates of the base element in the subsequent round of iteration and so on...)).
        Element& el = Get_Element(no);
        if(el.dup_index != -1)      // check if the elment has already been duplicated (dup_index != -1) and if so, break up immediately. This is important as otherwise Complexes that consist of two elements of the same family (e.g. dimers or Complexes of an Element and a Phosphorylate of that Element) would be duplicated twice.
            return;
        Duplicate_Element(el);
        dup_BaseElems.push_back(no);
        
        for(short i=0; i<el.outs.size(); ++i)   // iterative call of the function for the Complexes of the Element (which are base elements themselves)
            if(el.outs[i].type == CompReact)
                Duplicate_BaseEl_and_Dependants(Get_Arrow(el.outs[i]).sink[0]);
        for(short i=0; i<2; ++i)
            for(short j=0; j<el.phosphs[i].size(); ++j)   // iterative call for the Phosphorylates of the Element via call of 'Duplicate_Phosph_and_Dependants'.
                Duplicate_Phosph_and_Dependants( node(Phos, el.phosphs[i][j]) );
    };
    
    Duplicate_Phosph_and_Dependants = [&](const node node)  {    // similar to the function above except that it takes a Phosphorylate as input for which it only duplicates the Complexes via iterative function call but not the Phosphorylates (all Phosphorylates of the family are already duplicated in 'Duplicate_BaseEl_and_Dependants' when called with the respective base element).
        Element& el = Get_Element(node);
        if(el.dup_index != -1)
            return;
        Duplicate_Element(el);
        dup_Phosphs.push_back(node);
        
        for(short i=0; i<el.outs.size(); ++i)  // iterative call of 'Duplicate_BaseEl_and_Dependants' for the Complexes of the Element (which are base elements again)
            if(el.outs[i].type == CompReact)
                Duplicate_BaseEl_and_Dependants(Get_Arrow(el.outs[i]).sink[0]);
    };
    
    Reset_All_dup_Indices = [&]() {     // reset all 'dup_indices' of Elements stored in dup_BaseElems and dup_Phosphs to -1 and 'duplicated' of Arrows stored in dup_arrows to 'false' to mark them as unduplicated for the next mutation.
        for(short i=0; i<dup_BaseElems.size(); ++i)
            Get_Element(dup_BaseElems[i]).dup_index = -1;
        
        for(short i=0; i<dup_Phosphs.size(); ++i)
            Get_Element(dup_Phosphs[i]).dup_index = -1;
        
        for(short i=0; i<dup_arrows.size(); ++i)
            Get_Arrow(dup_arrows[i]).duplicated = false;
    };
    
    if(index==-1)     // if no index has been provided (default value = -1) choose a random index for a Protein to be duplicated
        index = Proteins.occupied[ ChooseInd(Proteins.occupied.size()) ];
    
    Duplicate_BaseEl_and_Dependants(node(Prot,index));   // trigger the cascade of iterative duplications of Complexes and Phosphorylates that form or transform from the respective Protein by calling 'Duplicate_BaseEl_and_Dependants' for the to-be-duplicated Protein.
    
    for(short i=0; i<dup_arrows.size(); ++i)   // duplicate also the Arrows that connect to at least one duplicated element that are now all stored (via arrow-references) in 'dup_arrows'.
        Duplicate_Arrow(dup_arrows[i]);
    
    Reset_All_dup_Indices();   // finally, reset all 'Element.dup_indices' to -1 and 'Arrow.duplicated' to 'false'.
    Save_Mutation(index, 10);
    return index;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Duplicate_Protein_V2(short index)  {    // duplicates the Protein with index 'index' or chooses a random index if no index is provided (default value -1) and duplicates all elements that form or transform from the Protein (the builders).
    if(Proteins.occupied.size() == 0)  // if no protein exists in the network, break up and return 'false'.
        return -1;
    
    vector<short> dup_Comps;   // saves references on the Complexes that must be duplicated in the process (all the Protein's or its Phosphorylates' Complexes and so on)
    vector<short> dup_Phosphs; // saves references on the Phosphorylates that must be duplicated in the process (all Phosphorylates of the respective Protein and of its or the Phosphorylates' Complexes and so on)
    
    dup_Comps.reserve(10);  // rough estimate of the needed space for the above vectors in order to avoid unnecessary reallocations.
    dup_Phosphs.reserve(15);
    
    function<void(const node no)> Collect_BaseEl_and_Dependants_for_Duplication;
    function<void(const short index)> Collect_Phosphs_and_Dependants_for_Duplication;
    
    Collect_BaseEl_and_Dependants_for_Duplication = [&](const node no)  {   // the function takes a node that MUST be a base element. If not already duplicated (dup_index == -1), this base element is registered in 'dup_BaseElems' (also all its in- and out-arrows (connecting arrows) are registered in the vector dup_arrows). Afterwards, the function iteratively calls itself and 'Collect_Phosph_and_Dependants_for_Duplication' (which is similar but takes a Phosphorylate as input), respectively, for the dependants, i.e. the Phosphorylates and Complexes of the base element, which leads to a duplication also of those elements. This therefore triggers a cascade in who's process all the builders (all elements that form or transform from the starting base element) are finally duplicated (i.e. also the Phosphorylates and Complexes of the Complexes of the base element and the Complexes of the Phosphorylates of the base element in the subsequent round of iteration and so on...)).
        Element& el = Get_Element(no);
        if(el.dup_index != -1)      // check if the elment has already been duplicated (dup_index != -1) and if so, break up immediately. This is important as otherwise Complexes that consist of two elements of the same family (e.g. dimers or Complexes of an Element and a Phosphorylate of that Element) would be duplicated twice.
            return;
        if(el.type==Comp)
            dup_Comps.push_back(no.index);
        
        for(short i=0; i<el.outs.size(); ++i)   // iterative call of the function for the Complexes of the Element (which are base elements themselves)
            if(el.outs[i].type == CompReact)
                Collect_BaseEl_and_Dependants_for_Duplication(Get_Arrow(el.outs[i]).sink[0]);
        for(short i=0; i<2; ++i)
            for(short j=0; j<el.phosphs[i].size(); ++j)   // iterative call for the Phosphorylates of the Element via call of 'Duplicate_Phosph_and_Dependants'.
                Collect_Phosphs_and_Dependants_for_Duplication( el.phosphs[i][j] );
    };
    
    Collect_Phosphs_and_Dependants_for_Duplication = [&](const short index)  {    // similar to the function above except that it takes the index of a Phosphorylate as input for which it only duplicates the Complexes via iterative function call but not the Phosphorylates (all Phosphorylates of the family are already duplicated in 'Duplicate_BaseEl_and_Dependants' when called with the respective base element).
        Element& el = Phosphorylates.elements[index];
        if(el.dup_index != -1)
            return;
        
        dup_Phosphs.push_back(index);
        
        for(short i=0; i<el.outs.size(); ++i)  // iterative call of 'Duplicate_BaseEl_and_Dependants' for the Complexes of the Element (which are base elements again)
            if(el.outs[i].type == CompReact)
                Collect_BaseEl_and_Dependants_for_Duplication(Get_Arrow(el.outs[i]).sink[0]);
    };
    
    if(index==-1)     // if no index has been provided (default value = -1) choose a random index for a Protein to be duplicated
        index = Proteins.occupied[ ChooseInd(Proteins.occupied.size()) ];
    
    Collect_BaseEl_and_Dependants_for_Duplication(node(Prot,index));   // trigger the cascade of iterative duplications of Complexes and Phosphorylates that form or transform from the respective Protein by calling 'Duplicate_BaseEl_and_Dependants' for the to-be-duplicated Protein.
    
    Duplicate_SubNetwork(index, dup_Comps, dup_Phosphs);
    Save_Mutation(index, 10);
    return index;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
node Genome_Version2::Choose_Config(Element& baseEl, short state, node exclude_node)  {      // chooses randomly a configuration of state 'state' of the base element 'baseEl' which is not the configuration 'exclude_node'. For that purpose, the function modifies temporarily the vector 'phosphs' of the base element (add base element and/or delete 'exclude_node') but resets all relevant changes before returning.
    array<vector<short>,2>& phosphs = baseEl.phosphs;     // create reference on the 'phosphs' vector of 'baseEl' for better handability
    short index;
    short Add_BaseEl_state = -1, Del_ExcludeNo_state = -1;   // shorts that indicate whether 'BaseEl' was added to or 'exclude_node' was deleted from 'phosphs', respectively, and in this case the state to which state the respective element has been added or deleted as the last one in 'phosphs' (either 0 (cytosolic) or 1 (membrane_bound)), otherwise set them to -1.
    bool Failed = false;     // bool that indicates whether an acceptable configuaration could or could not be found.
    
    if(exclude_node.type == Phos) {     // check whether exclude_node is a Phosphorylate and in that case delete it from the 'phosphs' vector and save its sate in 'Del_ExcludeNo_state'.
        Del_ExcludeNo_state = Get_Element(exclude_node).state;
        VecDel_element( phosphs[Del_ExcludeNo_state], exclude_node.index);
    }
    if(exclude_node.type != baseEl.type)  {    // if it's not the base element that shall be excluded push back '-1' to 'phosphs' which represents the base element and save its sate in 'Add_BaseEl_state'.
        Add_BaseEl_state = baseEl.state;
        phosphs[Add_BaseEl_state].push_back(-1);
    }
    
    if(state == arbitrary)  {  // first determine the state of the Configuration which is only necessary if no particular target state has been provided.
        short weights[2] = {short(phosphs[0].size()), short(phosphs[1].size())};
        state = ChooseInd_respective_weights(weights, 2);
    }
    
    if(state == cytosolic && phosphs[0].size() > 0)      // according to the required state check whether at least one possible configuration exists (all possible configurations are now contained in 'phosphs'; the base element as '-1') and if so, randomly choose one among them, otherwise do nothing and set 'Failed' to 'true'.
        index = phosphs[0][ ChooseInd(phosphs[0].size()) ];
    else if(state == membrane_bound && phosphs[1].size() > 0)
        index = phosphs[1][ ChooseInd(phosphs[1].size()) ];
    else
        Failed = true;
    
    if(Add_BaseEl_state >= 0)       // undo the changes done to 'phosphs' (if they have been done).
        phosphs[Add_BaseEl_state].pop_back();
    if(Del_ExcludeNo_state >= 0)
        phosphs[Del_ExcludeNo_state].push_back(exclude_node.index);
    
    if(Failed)       // if 'Failed' has been set to 'true', return 'node(-1,-1)', otherwise check if the choosen index is -1 in which case it represents the base element, otherwise return the respective Phosphorylate.
        return node(-1,-1);
    else
        return (index == -1) ? baseEl.base_element : node(Phos, index);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Add_Binding_Domain(bool EnzymaticReaction) {      // alternative name: add_Complex // chooses first two proteins (or if 'HigherOrderComplexes' == 'true' also a complex) as base element and then randomly pick a configuration of each one that define the two components of a new complex. Add that new complex as well as an CompReaction Arrow that connects the two components to the complex the a reverse Decay Arrow back into the components. If 'EnzymaticReaction' == 'true', add an additional Decay Arrow. One component (preferentially component1) is chosen as the substrate whereas the other one is chosen as the enzyme. The decay products of the additional Decay Arrow are then the enzyme (usually component2) and a different configuration of the substrate (usually component1). Only if this is not possible component1 instead of component2 may be chosen as the enzyme or the enzyme may also undergo a change in configuration.
    if(Proteins.occupied.size()==0)    // check whether the network is non-trivial and an element exists to which a binding domain can be added. If the network is trivial break up and return 'false'.
        return -1;
    
    auto Get_Number_Configs = [&](const Element& base_el, short state, const node& exclude_node = node(-1,-1))  {   // return the number of possible configurations of Element 'base_el' of state 'state' where 'exclude_node' is excluded as possible configuration.
        const array<vector<short>,2>& phosphs = base_el.phosphs;
        short n_configs = (state == arbitrary) ? phosphs[0].size()+phosphs[1].size() : phosphs[state].size();            // starting point is the number of Phosphorylates of the respective state or the sum of both states if state==arbitrary
        if(state == arbitrary || state == base_el.state)     // if the base element has also the desired state or if state==arbitrary, the number of possible configurations increases about 1
            n_configs += 1;
        if(exclude_node.type != -1)    // if on the other hand an element shall be excluded and it has the relevant state or if state==arbitrary, the number decreases about 1.
            if(state == arbitrary || state == Get_Element(exclude_node).state)
                n_configs -= 1;
        return n_configs;
    };
    
    auto Choose_BaseEl = [&](short type = -1)  {    // Randomly choose a base element of type 'type' (either Prot = 1, Comp = 2, or arbitrary = -1)
        if(type == -1)
            type = (ChooseInd(Proteins.occupied.size()+Complexes.occupied.size()) < Proteins.occupied.size())  ?  Prot : Comp;
        short index = (type == Prot)  ?  Proteins.occupied[ChooseInd(Proteins.occupied.size())] : Complexes.occupied[ChooseInd(Complexes.occupied.size())];
        return node(type,index);
    };
    
    auto Choose_Component = [&]()  {     // choose a component either as some configuration of a (randomly picked) Protein or allow also for a configuration of a Complex (if HigherOrderComplexes == true).
        short type_baseEl  =  HigherOrderComplexes  ?  -1 : Prot;
        node baseEl = Choose_BaseEl(type_baseEl);
        return Choose_Config(Get_Element(baseEl));
    };
    
    node no_component1 = Choose_Component();   // choose component1 and component2
    node no_component2 = Choose_Component();
    
    short state1 = Get_Element(no_component1).state;
    short state2 = Get_Element(no_component2).state;
    short state_complex  =  (state1==cytosolic && state2==cytosolic)  ?  cytosolic : ChooseInd(2);   // determine state of the new complex as cytosolic if the states of both component1 and 2 are cytosolic or randomly otherwise.
    /*
     short state_complex;            // same as the line above just more detailed
     if(state1==cytosolic && state2==cytosolic)
     state_complex = cytosolic;
     else if(state1==membrane_bound && state2==membrane_bound)
     state_complex = ChooseInd(2);
     else
     state_complex = ChooseInd(2);
     */
    short index_complex = Complexes.add_element(state_complex, 0);     // add the respective Arrows to 'CompReactions' and 'Decays'
    node no_complex = node(Comp, index_complex);
    CompReactions.add_arrow<2, 1>(array<node, 2>({no_component1, no_component2}), array<node, 1>({no_complex}));
    if(BackRatesOn==false)
        Decays.add_arrow<1, 2>(array<node, 1>({no_complex}), array<node, 2>({no_component1, no_component2}));
    
    if(EnzymaticReaction == false)  {
        Save_Mutation(index_complex, 8);
        return index_complex;
    }
    else  {     // if EnzymaticReaction == 'true' add a second decay reaction
        node no_enzyme, no_product;
        short decay_state = (state_complex == cytosolic) ? cytosolic : arbitrary;     // the state of both decay products must be either cytosolic if the complex is cytosolic or each of them is choosen randomly and independent otherwise.
        Element& baseEl1 = Get_Element( Get_Element(no_component1).base_element );
        Element& baseEl2 = Get_Element( Get_Element(no_component2).base_element );
        /*
         no_substrate = Choose_Config(no_prot1, decay_state, no_component1);    // no_substrate = (-1,-1) if not feasible because no additional configuration exists.
         no_enzyme =  (decay_state == arbitrary || Get_Element(no_component2).state == decay_state)   ?   no_component2 : Choose_Config(no_prot2, decay_state);     // change type if cytosolic -> membrane
         
         if(no_substrate.type == -1 || no_enzyme.type == -1)  {   // if it didn't work above interchange the role of component1 and 2 as substrate and enzyme
         no_substrate = Choose_Config(no_prot2, decay_state, no_component2);
         no_enzyme = (decay_state == arbitrary || Get_Element(no_component1).state == decay_state) ? no_component1 : Choose_Config(no_prot1, decay_state);
         if(no_substrate.type == -1 || no_enzyme.type == -1)   // if it still did not work return.
         return;
         }
         */
        if(Get_Number_Configs(baseEl1, decay_state, no_component1) > 0 && Get_Number_Configs(baseEl2, decay_state) > 0)   {   // check if component1 is viable as the substrate and component2 as enzyme (check whether enough configurations of them exist)
            no_product = Choose_Config(baseEl1, decay_state, no_component1);
            no_enzyme =  (decay_state == arbitrary || Get_Element(no_component2).state == decay_state)   ?   no_component2 : Choose_Config(baseEl2, decay_state);     // change type of enzyme if its state does not correspond with the prescribed decay_state (e.g. if cytosolic -> membrane + ?)
        }
        else if(Get_Number_Configs(baseEl2, decay_state, no_component2) > 0 && Get_Number_Configs(baseEl1, decay_state) > 0)  {   // if this was not viable, check for the possiblility to have component2 as substrate and component 1 as enzyme instead.
            no_product = Choose_Config(baseEl2, decay_state, no_component2);
            no_enzyme = (decay_state == arbitrary || Get_Element(no_component1).state == decay_state) ? no_component1 : Choose_Config(baseEl1, decay_state);
        }
        else  {          // if this is also impossible, return and do not add an additional Decay.
            Save_Mutation(index_complex, 8);
            return index_complex;
        }
        
        Decays.add_arrow<1, 2>( array<node, 1>({no_complex}), array<node, 2>({no_product, no_enzyme}) );        // if some constellation was possible add the additional Decay arrow.
        Save_Mutation(index_complex, 9);
        return index_complex;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Mutate_Concentration()  {
    short index = Proteins.occupied[ChooseInd(Proteins.occupied.size())];
    if(index >= 0)
        Proteins.elements[index].init_dens *= UnifRand()*1.5 + 0.5;      // multiply by a random number between 0.5 and 2.
    Save_Mutation(index, 4);
    return index;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Mutate_Function()  {
    short weights[2] = {3,2};
    short d = ChooseInd_respective_weights(weights,2);
    if(d==0)  {
        short index = Essential.delete_Essential();
        Save_Mutation(index, 6);
        return index;
    }
    else  {
        short index = Essential.add_Essential();
        Save_Mutation(index, 5);
        return index;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Delete_Phosphorylate()  {
    short index = Phosphorylates.delete_element();
    Save_Mutation(index, 12);
    return index;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Delete_Complex()  {
    short index = Complexes.delete_element();
    Save_Mutation(index, 13);
    return index;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Genome_Version2::Delete_Protein()  {
    short index = Proteins.delete_element();
    Save_Mutation(index, 14);
    return index;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Genome_Version2::Mutate_constant(array<short,9> weights)  {     // default: weights = {8,3,2,3,2,3,2,2,1}
    short Executed = -1;
    while(Executed < 0)  {
        short d=ChooseInd_respective_weights(&weights[0],9);
        switch(d)  {
            case 0:
                Executed = Mutate_Rate();
                break;
            case 1:
                Executed = Mutate_Concentration();
                break;
            case 2:
                Executed = Mutate_Function();
                break;
            case 3:
                Executed = Delete_Phosphorylate();
                break;
            case 4:
                Executed = Add_Phosphorylate_randomized();
                break;
            case 5:
                Executed = Delete_Complex();    // delete binding site
                break;
            case 6:
                Executed = Add_Binding_Domain(ChooseInd(2));
                break;
            case 7:
                Executed = Delete_Protein();
                break;
            case 8:
                Executed = Duplicate_Protein_V2();        // use Version 2 wich uses 'Duplicate_SubNetwork' and only allocates the necessary space in the 'elements' vectors.
                break;
        }
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Genome_Version2::Mutate_Rate_Concentration_Function()  {
    if(Phosphorylates.occupied.size()+Complexes.occupied.size() == 0)  {     // if the network is empty except for some Proteins, add a new element instead of mutating a rate (would not be possible anyway), a concentration or a function.
        Mutate_AddElement();
        return;
    }
    short weights[3] = {5,2,1};
    short d=ChooseInd_respective_weights(weights, 3);
    switch(d)  {
        case 0:
            Mutate_Rate();
            break;
        case 1:
            Mutate_Concentration();
            break;
        case 2:
            Mutate_Function();
            break;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
short Heaviside(size_t arg)  {           // Heaviside step function for use in 'Add_Element'
    return (arg <= 0)  ?  0 : 1;
};

void Genome_Version2::Mutate_AddElement()  {
    constexpr short N_Prot = 4;             // numbers that indicate mamimum numbers of elements. This is not strict for Phosphorylates and Complexes since by 'Duplicate_Element' these maximum numbers can be exceeded. If no element can be added without exceeding the maximum number an element is deleted instead.
    constexpr short N_Phosph = 20;
    constexpr short N_Comp = 15;
    short weights[3] = {short(2*Heaviside(N_Phosph-Phosphorylates.occupied.size())), short(2*Heaviside(N_Comp-Complexes.occupied.size())), Heaviside(N_Prot-Proteins.occupied.size())};
    short d=ChooseInd_respective_weights(weights, 3);
    if(d==-1)  {
        Mutate_DeleteElement();
        return;
    }
    switch(d)  {
        case 0:
            Add_Phosphorylate();
            break;
        case 1:
            Add_Binding_Domain(ChooseInd(2));    // delete binding site
            break;
        case 2:
            Duplicate_Protein_V2();
            break;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Genome_Version2::Mutate_DeleteElement()  {
    short weights[3] = {short(2*Phosphorylates.variable.size()), short(2*Complexes.variable.size()), short(Proteins.variable.size())};
    short d=ChooseInd_respective_weights(weights, 3);
    if(d==-1)  {
        Mutate_AddElement();
        return;
    }
    switch(d)  {
        case 0:
            Delete_Phosphorylate();
            break;
        case 1:
            Delete_Complex();    // delete binding site
            break;
        case 2:
            Delete_Protein();
            break;
    }
};


















