//
//  Print_and_Rebuild.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 06/04/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#include "Genome_Version2.hpp"



void Genome_Version2::PrintGenome(ofstream& file)   // takes a file reference and and writes all relevnat data of Genome into it (all tab arrays as well as the rates and indices of OdeSyst1 and OdeSyst2)
{
    /*
     ofstream file (name_of_file);
     */
    file << endl;
    file << "Constants\n";
    file << left << setw(4) << BackRatesOn << setw(4) << HigherOrderComplexes << setw(8) << ConfigInterconnectionProb << setw(4) << OnlyOneEssFuncPerElem << setw(4) << NumEssFunc << setw(8) << InitialRate_Linear_Max << setw(8) << InitialRate_Quadratic_Max << setw(8) << InitialConcentration_Max << endl;      // constants that regulate mutations
    file << left << setw(4) << L_small << setw(4) << L_large << setw(4) << N_discretisation << setw(12) << D_cytosolic << setw(12) << D_membrane << setw(12) << T_max << "\n\n";        // constants of the simulation
    
    file << "Proteins\n";
    Proteins.print_Elements(file);
    
    file << "Phosphorylates\n";
    Phosphorylates.print_Elements(file);
    
    file << "Complexes\n";
    Complexes.print_Elements(file);
    
    file << "Transforms\n";
    Transforms.print_Arrows(file);
    
    file << "CompReactions\n";
    CompReactions.print_Arrows(file);
    
    file << "Decays\n";
    Decays.print_Arrows(file);
    
    file << "Essential Elements\n";
    Essential.print_Essential(file);
};


void Genome_Version2::Print_Mutations(ofstream& file)  {
    for(short i=8; i<Mutations.size(); ++i)       // start at index 8 because the first 8 Mutations are from the initial setup.
        file << left << setw(5) << Mutations[i];
};

void Genome_Version2::RebuildGenome(ifstream& file)  {
    // assumes that all vectors are empty (cleared) which they should be if the class Genome has been created with the usual constructor (but not applied to "Create_Random_Initial_Network").
    short Proteins_size, Phosphorylates_size, Complexes_size, Transforms_size, CompReactions_size, Decays_size;;
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));      // sarting point where the file pointer is assumed to point to is "Proteins";
    file >> Proteins_size;
    for(short i=0; i<Proteins_size; ++i)
        file >> Proteins;
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // rest of line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // empty line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // "Phosphorylates"
    
    file >> Phosphorylates_size;
    for(short i=0; i<Phosphorylates_size; ++i)
        file >> Phosphorylates;
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // rest of line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // empty line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // "Complexes"
    
    file >> Complexes_size;
    for(short i=0; i<Complexes_size; ++i)
        file >> Complexes;
    
    reference_Phosphorylates();         // after all three Elements types have been read in reference the phosphorylates in 'phosphs' in their respective (Protein or Complex) base elements.
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // rest of line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // empty line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // "Transforms"
    
    file >> Transforms_size;
    for(short i=0; i<Transforms_size; ++i)
        file >> Transforms;
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // rest of line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // empty line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // "CompReactions"
    
    file >> CompReactions_size;
    for(short i=0; i<CompReactions_size; ++i)
        file >> CompReactions;
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // rest of line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // empty line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // "Decays"
    
    file >> Decays_size;
    for(short i=0; i<Decays_size; ++i)
        file >> Decays;
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // rest of line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // empty line
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));    // "Essential Elements"
    
    file >> Essential;
    reference_essFunc();            // set the ess_func in the essential elements to the index of the respective essential function
};

void Genome_Version2::RebuildGenome(const string filename)  {
    ifstream file(filename);
    if(file.is_open())  {
        float init_concentration_max, l1, l2, d_cytosolic, d_membrane, t_max;
        short n_discr;
        for(short i=0; i<4; ++i)
            file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
        float temp;
        for(short i=0; i<7; ++i)
            file >> temp;
        file >> init_concentration_max;
        file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
        file >> l1;
        file >> l2;
        file >> n_discr;
        file >> d_cytosolic;
        file >> d_membrane;
        file >> t_max;
        if(init_concentration_max!=InitialConcentration_Max || l1!=L_small || l2!=L_large || N_discretisation!=n_discr || d_cytosolic!=D_cytosolic || d_membrane!=D_membrane || t_max!=T_max)
            cout << "Warning: Parameters do not coincide!\n";
        file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
        file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
        RebuildGenome(file);
    }
    file.close();
};

void Genome_Version2::Print_OdeSyst(ofstream& file)  {
    vector<short> states;
    Make_OdeSyst(states);
    
    short OdeSyst1_size = 0, OdeSyst2_size = 0;
    for(short i=0; i<OdeSyst1.size(); ++i)          // caluclate OdeSyst1_size and OdeSyst2_size as the number of all entries in the respective OdeSyst.
        OdeSyst1_size += OdeSyst1[i].size();
    for(short i=0; i<OdeSyst2.size(); ++i)
        OdeSyst2_size += OdeSyst2[i].size();
    
    file << "Proteins\n";
    file << Proteins.occupied.size() << endl;
    for(short i=0; i<Proteins.elements.size(); ++i)
        if(Proteins.elements[i].existent == true)
            file << left << Get_OdeIndex(node(Prot,i)) << '\t' << Proteins.elements[i].init_dens << endl;
    file << endl;
    
    file << "Essential Elements (Ode Index)\n";
    for(short i=0; i<Essential.ess_elements.size(); ++i)  {
        for(short j=0; j<Essential.ess_elements[i].size(); ++j)
            file << left << setw(4) << Get_OdeIndex(Essential.ess_elements[i][j]);
        file << endl;
    }
    file << endl;
    
    file << "initial concentrations\n";
    vector<float> concentrations(Get_NumberElements() * N_discretisation, 0);     // define the state type that the Ode Solver deals with as a vector that is initialized with 0
    FillIn_InitialHomogeneousProteinConcentrations(concentrations, N_discretisation);   // set initial protein concentrations
    file << concentrations.size() << endl;
    for(short i=0; i<concentrations.size(); ++i)
        file << left << setprecision(7) << concentrations[i] << "  ";
    file << "\n\n";
    
    file << "states\n";
    file << states.size() << endl;
    for(short i=0; i<states.size(); ++i)
        file << left << setw(4) << states[i];
    file << "\n\n";
    
    file << "OdeSyst1\n";
    file << OdeSyst1_size << endl;
    for(short i=0; i<OdeSyst1.size(); ++i)
        for(short j=0; j<OdeSyst1[i].size(); ++j)
            file << left << i  << '\t' << setw(16) << OdeSyst1[i][j].rate << OdeSyst1[i][j].ind << endl;
    file << endl;
    
    file << "OdeSyst2\n";
    file << OdeSyst2_size << endl;
    for(short i=0; i<OdeSyst2.size(); ++i)
        for(short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << i << '\t' << setw(16) << OdeSyst2[i][j].rate << OdeSyst2[i][j].ind1 << '\t' << OdeSyst2[i][j].ind2 << endl;
    file << endl;
    
    file << "prot_OdeIndex\n";
    file << prot_OdeIndex.size() << endl;
    for(short i=0; i<prot_OdeIndex.size(); ++i)
        file << left << setw(4) << prot_OdeIndex[i];
    file << "\n\n";
    
    file << "comp_OdeIndex\n";
    file << comp_OdeIndex.size() << endl;
    for(short i=0; i<comp_OdeIndex.size(); ++i)
        file << left << setw(4) << comp_OdeIndex[i];
    file << "\n\n";
    
    file << "phosph_OdeIndex\n";
    file << phosph_OdeIndex.size() << endl;
    for(short i=0; i<phosph_OdeIndex.size(); ++i)
        file << left << setw(4) << phosph_OdeIndex[i];
    file << "\n\n";
    
    Make_SpeciesObject();
    file << "Species Affiliations\n";
    file << SpeciesObject.size() << endl;
    for(short i=0; i<SpeciesObject.size(); ++i)  {
        for(short j=0; j<SpeciesObject[i].size(); ++j)
            file << left << setw(4) << SpeciesObject[i][j];
        file << endl;
    }
    file << endl;
};


short count_smaller(short value, vector<short>& vec)  {
    short count = 0;
    for(short i : vec)
        if(i<value)
            count++;
    return count;
};

void Genome_Version2::Print_Biograph(ofstream& file)  {
    constexpr bool SimpleNodeNames = true;
//    vector<string> nodes_name;
//    vector<Element*> waiting_list;      // saves pionters on Complexes whose name generation must be postponed due to the fact that a complex subcomponent appears with higher OdeIndex and therefore the name for the subcomponent has not yet been created.
    vector<short> nodes_type;
    vector<short> edges_foot;
    vector<short> edges_top;
    vector<short> edges_type;
    vector<float> edges_weight;
    short num_edges = Transforms.occupied.size() + 2*CompReactions.arrows.size() + 2*Decays.arrows.size();
//    nodes_name.reserve(Get_NumberElements());
    nodes_type.reserve(num_edges);
    edges_foot.reserve(num_edges);
    edges_top.reserve(num_edges);
    edges_type.reserve(num_edges);
    edges_weight.reserve(num_edges);
    
    short control_OdeIndex = 0;
    auto print_NodeNameSimple_BaseEl = [&](Element& el)  {       // input: reference on a base element
        if(el.existent==true)  {
            nodes_type.push_back(el.type);
            char graph_type_specifier = (el.type == Prot) ? 'P' : 'C';
            short graph_ind = el.index+1-count_smaller(el.index, Proteins.empty);
            file << graph_type_specifier << graph_ind << ((el.state==cytosolic) ? 'c' : 'm') << endl;
            assert(Get_OdeIndex({el.type, el.index})==control_OdeIndex++);
            for(short i=0; i<el.phosphs[0].size(); ++i)  {
                nodes_type.push_back(Phos);
                file << graph_type_specifier << graph_ind << '*' << 'c' << endl;
                assert(Get_OdeIndex({Phos,el.phosphs[0][i]})==control_OdeIndex++);
            }
            for(short i=0; i<el.phosphs[1].size(); ++i)  {
                nodes_type.push_back(Phos);
                file << graph_type_specifier << graph_ind << '*' << 'm' << endl;
                assert(Get_OdeIndex({Phos,el.phosphs[1][i]})==control_OdeIndex++);
            }
        }
    };
/*
    auto make_NodeName_BaseEl = [&](Element& el)  {       // input: reference on a base element
        string name_base;
        if(el.existent==true)  {
            if(el.type == Prot)
                name_base = 'P' + to_string(el.index+1-count_smaller(el.index, Proteins.empty));
            else if(el.type == Comp)  {
                Arrow* Ar;
                for(arrow ar : el.ins)
                    if(ar.type==CompReact)
                        Ar = &Get_Arrow(ar);
                short OdeIndex_component1 = Get_OdeIndex(Ar->source[0]);
                short OdeIndex_component2 = Get_OdeIndex(Ar->source[1]);
                if(OdeIndex_component1<nodes_name.size() && OdeIndex_component2<nodes_name.size())
                    name_base = '('+ nodes_name[OdeIndex_component1] + ':' + nodes_name[OdeIndex_component2] + ')';
                else
                    waiting_list.push_back(&el);        // if the complex name cant be generated yet because of a still lacking subcomponent's name put it on the waiting list.
            }
                
            nodes_name.push_back(name_base+((el.state==cytosolic) ? 'c' : 'm'));
            assert(Get_OdeIndex({el.type,el.index})==nodes_name.size()-1);
            for(short i=0; i<el.phosphs[0].size(); ++i)  {
                nodes_name.push_back(name_base + '*' + 'c');
                assert(Get_OdeIndex({Phos,el.phosphs[0][i]})==nodes_name.size()-1);
            }
            for(short i=0; i<el.phosphs[1].size(); ++i)  {
                nodes_name.push_back(name_base + '*' + 'm');
                assert(Get_OdeIndex({Phos,el.phosphs[1][i]})==nodes_name.size()-1);
            }
        }
    };
*/
    function<string(Element& el)> make_NodeName;       // preliminary declaration of function object in order for the function to be able to call itself iteratively
    make_NodeName = [&](Element& el)  {       // input: reference on a base element
        if(el.type == Prot)
            return ('P' + to_string(el.index+1-count_smaller(el.index, Proteins.empty)));
        else if(el.type == Phos)
            return make_NodeName(Get_Element(el.base_element)) + '*';
        else    {     // (el.type == Comp)
            Arrow* Ar;
            for(arrow ar : el.ins)
                if(ar.type==CompReact)
                    Ar = &Get_Arrow(ar);
            string component1 = make_NodeName(Get_Element(Ar->source[0]));
            string component2 = make_NodeName(Get_Element(Ar->source[1]));
            return '('+ component1 + ':' + component2 + ')';
        }
    };
    
    function<void(Element& el)> print_NodeName = [&](Element& el)  {
        nodes_type.push_back(el.type);
        file << make_NodeName(el) + (el.state==cytosolic ? 'c' : 'm') << endl;
    };
    
    auto add_edge_to_graph = [&](node source, node sink, short type, float weight)   {
        edges_foot.push_back(Get_OdeIndex(source));
        edges_top.push_back(Get_OdeIndex(sink));
        edges_type.push_back(type);
        edges_weight.push_back(weight);
    };

    
// create biograph data from network
// create edges_foot, edges_top, edges_weight, edges_type; node names are printed directly below
    for(Arrow& Ar : Transforms.arrows)
        if(Ar.existent==true)  {
            add_edge_to_graph(Ar.source[0], Ar.sink[0], Trans, Ar.rate);
            if(Ar.back==true)
                add_edge_to_graph(Ar.sink[0], Ar.source[0], Trans, Ar.back_rate);
        }
    for(Arrow& Ar : CompReactions.arrows)
        if(Ar.existent==true)  {
            add_edge_to_graph(Ar.source[0], Ar.sink[0], CompReact, Ar.rate);
            add_edge_to_graph(Ar.source[1], Ar.sink[0], CompReact, Ar.rate);
            if(Ar.back==true) {
                add_edge_to_graph(Ar.sink[0], Ar.source[0], Dec, Ar.back_rate);
                add_edge_to_graph(Ar.sink[0], Ar.source[1], Dec, Ar.back_rate);
            }
        }
    for(Arrow& Ar : Decays.arrows)
        if(Ar.existent==true)  {
            add_edge_to_graph(Ar.source[0], Ar.sink[0], Dec, Ar.rate);
            add_edge_to_graph(Ar.source[0], Ar.sink[1], Dec, Ar.rate);
            if(Ar.back==true)  {
                add_edge_to_graph(Ar.sink[0], Ar.source[0], CompReact, Ar.back_rate);
                add_edge_to_graph(Ar.sink[1], Ar.source[0], CompReact, Ar.back_rate);
            }
        }
    
// print biograph data
    file << "nodes_name" << endl;
    file << Get_NumberElements() << endl;
    if(SimpleNodeNames)  {
        for(Element& el : Proteins.elements)
            print_NodeNameSimple_BaseEl(el);
        for(Element& el : Complexes.elements)
            print_NodeNameSimple_BaseEl(el);
    }
    else
        RangeThrough(print_NodeName);
    file << endl;
    
    file << "nodes_type" << endl;
    file << nodes_type.size() << endl;
    for(short type : nodes_type)
        file << left << setw(4) << type;
    file << "\n\n";
    
    file << "edges_foot" << endl;
    file << edges_foot.size() << endl;
    for(short foot : edges_foot)
        file << left << setw(4) << foot;
    file << "\n\n";
    
    file << "edges_top" << endl;
    file << edges_top.size() << endl;
    for(short top : edges_top)
        file << left << setw(4) << top;
    file << "\n\n";
    
    file << "edges_type" << endl;
    file << edges_type.size() << endl;
    for(short type : edges_type)
        file << left << setw(4) << type;
    file << "\n\n";
    
    file << "edges_weight" << endl;
    file << edges_weight.size() << endl;
    for(float weight : edges_weight)
        file << left << setprecision(7) << setw(15) << weight;
    file << "\n\n";
};


/////////////////////////////////////////////////// Matlab Interface //////////////////////////////////////////////////////////////////
void Genome_Version2::Append_OdeSyst_Biograph(string path_from, string path_to)  {
    ifstream file(path_from);
    if(file.is_open()==false)  {
        cout << "Could not open file\n";
        return;
    }
    for(short i=0; i<7; ++i)
        file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    RebuildGenome(file);
    file.seekg(0);
    
    ofstream file_copy(path_to);
    if(file.is_open()==false)  {
        cout << "Could not open file_copy\n";
        return;
    }
    string line;
    for(short i=0; i<7; ++i)  {
        getline(file, line);
        file_copy << line << endl;
    }
    Print_OdeSyst(file_copy);
    Print_Biograph(file_copy);
    while(getline(file,line))  {
        file_copy << line << endl;
    }
    file.close();
    file_copy.close();
};

