//
//  Tests.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 04/03/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"

// #ifndef NDEBUG

template<class VecClass, class T>
bool CheckPresence(const VecClass& vec, const T& val, short VecSize=-1)
{
    if(VecSize==-1)
        VecSize=vec.size();
    bool Present=false;
    for(short i=0; i<VecSize; ++i)
        if(vec[i]==val)
            Present=true;
    return (Present == false)  ?  false : true;
};


void Genome_Version2::Test_Elements(Elements& Elems)  {
    for(short i=0; i<Elems.empty.size(); ++i)  {      // check if all elements listed in 'empty' have 'existent' == false
        Element& el = Elems.elements[Elems.empty[i]];
        assert(el.existent == false);
        assert(el.ins.size()==0);           // check if all vectors 'ins', 'outs', and 'phosphs' have been cleared (have size 0)
        assert(el.outs.size()==0);
        assert(el.phosphs[0].size()==0);
        assert(el.phosphs[1].size()==0);
    }
    
    for(short i=0; i<Elems.occupied.size(); ++i)  {    // check all elements listed in 'occupied'
        short index = Elems.occupied[i];
        Element& el = Elems.elements[index];
        node self_node = node(Elems.type,index);
        assert(el.existent == true);         // check if 'existent' = ture
        assert(el.type == Elems.type);       // check if el.type coincides with Elems.type
        if(el.type == Prot)                   // for Proteins check if state is cytosolic
            assert(el.state == cytosolic);
        
        if(el.type == Prot || el.type == Comp)  {
            assert(el.base_element == self_node);    // for base elements (Proteins and Complexes) check if base_element == self_node
            for(short i=0; i<2; ++i)        // for base elements (Proteins and Complexes) check if all entries in 'phosph' reference an existent Phosphorylate
                for(short j=0; j<el.phosphs[i].size(); ++j)
                    assert(Get_Element(node(Phos,el.phosphs[i][j])).existent == true);
        }
        if(el.type == Phos)    // for Phosphorylates check whether its index is contained in 'phosphs' of the base element
            assert(CheckPresence(Get_Element(el.base_element).phosphs[el.state], index) == true);
        assert(el.dup_index == -1);      // dup_index should always be -1
        
        short num_Trans = 0, num_CompReact = 0, num_Decay = 0;     // count the number of Transform and CompReact in-connections
        for(short j=0; j<el.ins.size(); ++j)  {
            Arrow& Ar = Get_Arrow(el.ins[j]);
            if(Ar.type==Trans) ++num_Trans;
            if(Ar.type==CompReact) ++num_CompReact;
            if(Ar.type==Dec) ++num_Decay;
            assert(CheckPresence(Ar.sink, self_node, Ar.N_sink) == true);     // check if each in-connecting arrow has the same self_node as a sink.
        }
        if(el.type==Phos)    // Phosphorylates should have at least one in-connection (Transform in-connection or Decay in-connection)
            assert(num_Trans+num_Decay!=0);
        if(el.type==Comp)    // Complexes should have exactly one CompReact in-connection
            assert(num_CompReact==1);
        
        for(short j=0; j<el.outs.size(); ++j)  {
            Arrow& Ar = Get_Arrow(el.outs[j]);
            assert(CheckPresence(Ar.source, self_node, Ar.N_source) == true);  // check if each out-connecting arrow has the self_node as a source.
        }
    }
};


void Genome_Version2::Test_Arrows(Arrows& Arrs)  {
    for(short i=0; i<Arrs.empty.size(); ++i)      // check if all arrows listed in 'empty' have 'existent' == false
        assert(Arrs.arrows[Arrs.empty[i]].existent == false);
    
    for(short i=0; i<Arrs.occupied.size(); ++i)  {    // check all arrows listed in 'occupied'
        short index = Arrs.occupied[i];
        Arrow& ar = Arrs.arrows[index];
        arrow self_arrow = arrow(Arrs.type,index);
        assert(ar.existent == true);         // check if 'existent' = ture
        assert(ar.type == Arrs.type);       // check if ar.type coincides with Arrs.type
        
        assert(ar.duplicated == false);      // duplicated should always be false
        if(ar.type==Trans)          // check if N_source and N_sink correspond with the type of the arrow
            assert(ar.N_source==1 && ar.N_sink==1);
        if(ar.type==CompReact)
            assert(ar.N_source==2 && ar.N_sink==1);
        if(ar.type==Dec)
            assert(ar.N_source==1 && ar.N_sink==2);
        
        for(short j=0; j<ar.N_source; ++j)  {          // check whether the respective arrow is referenced in the 'outs' vector of the source element(s).
            Element& el = Get_Element(ar.source[j]);
            assert(CheckPresence(el.outs, self_arrow) == true);
        }
        for(short j=0; j<ar.N_sink; ++j)  {            // check whether the respective arrow is referenced in the 'ins' vector of the sink element(s).
            Element& el = Get_Element(ar.sink[j]);
            assert(CheckPresence(el.ins, self_arrow) == true);
        }
    }
};

void Genome_Version2::Test_Essential()  {        // test if all elements listed in Essential.ess_elements are actually existent
    for(short i=0; i<Essential.ess_elements.size(); ++i)
        for(short j=0; j<Essential.ess_elements[i].size(); ++j)  {
            assert(Get_Element(Essential.ess_elements[i][j]).existent==true);
            assert(Get_Element(Essential.ess_elements[i][j]).ess_func==i);
        }
};
/*
 template<typename T>
 void check_identical_content(T* pt_lower1, T* pt_upper1, T* pt_lower2, T* pt_upper2)  {
 for(T* pt1=pt_lower1; pt1!=pt_upper1; ++pt1)  {
 for(T* pt2=pt_lower2; pt2!=pt_upper2; ++pt2)  {
 if(*pt2 == *pt1)
 break;
 }
 assert(0==1);
 }
 }
 */

void Genome_Version2::Test_OdeSyst()  {
    vector<ode1> OdeSyst1_compressed;
    vector<short> OdeSyst1_numbers;
    vector<ode2> OdeSyst2_compressed;
    vector<short> OdeSyst2_numbers;
    vector<short> states1;                 // create two states vectors to call the two differeent methods each with one of them and then compare the two outcomes for consistency
    vector<short> states2;
    
    Make_OdeSyst(states1);
    Transfer_directly_to_OdeSystCompressed(states2,OdeSyst1_compressed,OdeSyst1_numbers,OdeSyst2_compressed,OdeSyst2_numbers);
    assert(states1==states2);                              // the produced states vectors should be exactly identical
    assert(OdeSyst1.size() == OdeSyst2.size());            // OdeSsyt1 and OdeSyst2 must have the same size
    assert(OdeSyst1.size()==OdeSyst1_numbers.size());      // compare sizes of OdeSyst1 and OdeSyst1_numbers
    assert(OdeSyst2.size()==OdeSyst2_numbers.size());
    
    for(short i=0; i<OdeSyst1.size(); ++i)                // check whether the size of each OdeSyst[i] is equal to the corresponding entry in OdeSyt_numbers, so the number of contributions to each element is identical for both methods.
        assert(OdeSyst1[i].size()==OdeSyst1_numbers[i]);
    for(short i=0; i<OdeSyst2.size(); ++i)
        assert(OdeSyst2[i].size()==OdeSyst2_numbers[i]);
    
    // check consistency of OdeSyst1_compressed with OdeSyst1
    {
        vector<ode1>::iterator it_lower;
        vector<ode1>::iterator it_upper = OdeSyst1_compressed.begin();
        for(short i=0; i<OdeSyst1.size(); ++i)  {
            it_lower = it_upper;
            it_upper = it_upper + OdeSyst1_numbers[i];
            assert( is_permutation(it_lower, it_upper, OdeSyst1[i].begin()) );
        }
    }
    // check consistency of OdeSyst2_compressed with OdeSyst2
    {
        vector<ode2>::iterator it_lower;
        vector<ode2>::iterator it_upper = OdeSyst2_compressed.begin();
        for(short i=0; i<OdeSyst2.size(); ++i)  {
            it_lower = it_upper;
            it_upper = it_upper + OdeSyst2_numbers[i];
            assert( is_permutation(it_lower, it_upper, OdeSyst2[i].begin()) );
        }
    }
    
    // check mass-conservation of OdeSyst1: check if each ode1 in OdeSyst1 with positive rate is opposed by a corresponding ode1 with negative rate for the ind of this ode1
    float rate;
    for(short i=0; i<OdeSyst1.size(); ++i)
        for(short j=0; j<OdeSyst1[i].size(); ++j)
            if((rate = OdeSyst1[i][j].rate) > 0)  {
                short ind = OdeSyst1[i][j].ind;
                auto p = find(OdeSyst1[ind].begin(),OdeSyst1[ind].end(),ode1(-rate,ind));
                assert(p!=OdeSyst1[ind].end());
            }
    // check mass-conservation of OdeSyst2: check if each ode2 in OdeSyst2 with positive rate is opposed by two corresponding ode2 with negative rate for the ind1 and ind2 of this ode2
    for(short i=0; i<OdeSyst2.size(); ++i)
        for(short j=0; j<OdeSyst2[i].size(); ++j)
            if((rate = OdeSyst2[i][j].rate) > 0)  {
                short ind1 = OdeSyst2[i][j].ind1;
                short ind2 = OdeSyst2[i][j].ind2;
                auto p = find(OdeSyst2[ind1].begin(),OdeSyst2[ind1].end(),ode2(-rate,ind1,ind2));
                auto q = find(OdeSyst2[ind2].begin(),OdeSyst2[ind2].end(),ode2(-rate,ind1,ind2));
                assert(p!=OdeSyst2[ind1].end() && q!=OdeSyst2[ind2].end());
            }
};


void Genome_Version2::Test_All()  {
#ifndef NDEBUG
    Test_Elements(Proteins);
    Test_Elements(Phosphorylates);
    Test_Elements(Complexes);
    Test_Essential();
    Test_Arrows(Transforms);
    Test_Arrows(CompReactions);
    Test_Arrows(Decays);
    Test_OdeSyst();
    Test_PrintRebuild();
#endif
};


void Genome_Version2::Make_SpeciesObject()  {
    function<void(Element& el, vector<short>& affiliations)> Get_SpeciesAffiliation;
    Get_SpeciesAffiliation = [this, &Get_SpeciesAffiliation](Element& el, vector<short>& affiliations)  {
        switch(el.type)  {
            case Prot:
                affiliations.push_back(el.index);
                return;
            case Phos:
                if(el.base_element.type==Prot)    // shortcut; could also call Get_SpeciesAffiliation() recursively also in this first case of the if condition.
                    affiliations.push_back(el.base_element.index);
                else              // in case base_element is a complex, which can only happen if HigherOrderComplexes is set to true
                    Get_SpeciesAffiliation(Get_Element(el.base_element), affiliations);
                return;
            case Comp:
                for(short i=0; i<el.ins.size(); ++i)
                    if(el.ins[i].type==CompReact)  {
                        auto source = Get_Arrow(el.ins[i]).source;
                        Get_SpeciesAffiliation(Get_Element(source[0]), affiliations);
                        Get_SpeciesAffiliation(Get_Element(source[1]), affiliations);
                    }
                return;
        }
    };
    auto Affiliate_Phosphs = [this](Element& el, short species_ind)  {
        for(short phosph_ind : el.phosphs[0])
            SpeciesObject[species_ind].push_back(Get_OdeIndex(node(Phos,phosph_ind)));
        for(short phosph_ind : el.phosphs[1])
            SpeciesObject[species_ind].push_back(Get_OdeIndex(node(Phos,phosph_ind)));
    };
    
    SpeciesObject.clear();
    SpeciesObject.resize(Proteins.elements.size());
    TotalMasses.clear();
    TotalMasses.reserve(Proteins.occupied.size());
    
    for(short i=0; i<Proteins.elements.size(); ++i)         // assign Proteins and their Phosphorylates to species
        if(Proteins.elements[i].existent == true)  {
            TotalMasses.push_back(Proteins.elements[i].init_dens);
            SpeciesObject[i].push_back(Get_OdeIndex(node(Prot,i)));
            Affiliate_Phosphs(Proteins.elements[i], i);
        }
    
    for(Element& el : Complexes.elements)
        if(el.existent==true)  {
            vector<short> affiliations;
            Get_SpeciesAffiliation(el, affiliations);
            for(short species_ind : affiliations)  {
                SpeciesObject[species_ind].push_back(Get_OdeIndex(node(Comp,el.index)));
                Affiliate_Phosphs(el, species_ind);
            }
        }
    // Compress SpeciesObject by removing all vectors that have size 0 (that belong to non-existent Proteins)
    auto pend = remove_if (SpeciesObject.begin(), SpeciesObject.end(), [](vector<short>& species) -> bool {return species.size()==0;});
    assert(Proteins.occupied.size() == pend-SpeciesObject.begin());         // check if SpeciesObject is now of size Proteins.occupied.size().
    SpeciesObject.resize(Proteins.occupied.size());
};


void Genome_Version2::Save_InitialConcentrations(vector<float>& concentrations)  {       // recalculate TotalMasses by saving the initial total protein concentrations if initial conditions with perturbatins are applied
    TotalMasses.clear();
    TotalMasses.reserve(Proteins.elements.size());
    for(short i=0; i<Proteins.elements.size(); ++i)
        if(Proteins.elements[i].existent==true)  {
            auto begin_it = concentrations.begin() + Get_OdeIndex(node(Prot,i))*N_discretisation;
            TotalMasses.push_back( accumulate(begin_it, begin_it+N_discretisation, 0.) );
        }
};


void Genome_Version2::Check_MassConservation(vector<float>& concentrations)  {
    for(short i=0; i<SpeciesObject.size(); ++i)  {
        double total_mass = 0;
        for(short ode_ind : SpeciesObject[i])  {
            auto begin_it = concentrations.begin() + N_discretisation * ode_ind;
            total_mass += accumulate(begin_it, begin_it + N_discretisation, 0.);
        }
        assert(abs(total_mass-TotalMasses[i])/TotalMasses[i] < 5.0e-3);    // check if relative difference between initial and final concentration is smaller than 1.0e-4
     //   cout << total_mass/N_discretisation << "  " << InitialConcentrations[i]/N_discretisation << "  " << Get_Element(node(Prot,i)).init_dens << "  " << N_discretisation * Get_Element(node(Prot,i)).init_dens / InitialConcentrations[i] << endl;
    }
};


void Genome_Version2::cut_empty()  {
    auto Elements_cut_empty = [this](Elements& Elems)  {
        while(Elems.elements.size()>0 && Elems.elements.back().existent == false)  {
            Elems.elements.pop_back();
            VecDel_element(Elems.empty, short(Elems.elements.size()) );
        }
    };
    auto Arrows_cut_empty = [this](Arrows& Arrs)  {
        while(Arrs.arrows.size()>0 && Arrs.arrows.back().existent == false)  {
            Arrs.arrows.pop_back();
            VecDel_element(Arrs.empty, short(Arrs.arrows.size()) );
        }
    };
    Elements_cut_empty(Proteins);
    Elements_cut_empty(Phosphorylates);
    Elements_cut_empty(Complexes);
    Arrows_cut_empty(Transforms);
    Arrows_cut_empty(CompReactions);
    Arrows_cut_empty(Decays);
};

void Genome_Version2::Test_PrintRebuild()  {
    ofstream file("PrintRebuild_TEST.txt");
    assert(file.is_open());
    PrintGenome(file);
    file.close();
    
    Genome_Version2 genome_rebuild;
    ifstream File("PrintRebuild_TEST.txt");
    assert(File.is_open());
    for(short i=1; i<=5; ++i)
        File.ignore(numeric_limits<streamsize>::max(), File.widen('\n'));
    genome_rebuild.RebuildGenome(File);
    cut_empty();
    assert(Proteins == genome_rebuild.Proteins);
    assert(Phosphorylates == genome_rebuild.Phosphorylates);
    assert(Complexes == genome_rebuild.Complexes);
    assert(Transforms == genome_rebuild.Transforms);
    assert(CompReactions == genome_rebuild.CompReactions);
    assert(Decays == genome_rebuild.Decays);
    assert(Essential.ess_elements == genome_rebuild.Essential.ess_elements);
    
    File.close();
}

// #endif

