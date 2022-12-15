//
//  Genome_Version2.hpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 20/12/16.
//  Copyright © 2016 Florian Gartner. All rights reserved.
//

#ifndef Genome_Version2_hpp
#define Genome_Version2_hpp

#include "Preliminary.hpp"



class Genome_Version2 {
    friend int main();     // make main a friend for debugging purposes
    
// define structs that save all existent elements of a type (or of a type and species in the case of Phosphorylates) or interactions of a type, respectively
    
///////////////////////////////////////////////////// Elements //////////////////////////////////////////////////////////////////////////////////
    struct Elements  {      // implementation in "Elements.cpp"
        Genome_Version2& genome_v2;
        const short type;
        vector<short> empty;           
        vector<short> occupied;
        vector<short> variable;
        vector<Element> elements;
            
        Elements(short type, Genome_Version2& genome_v2, short reserve_size);           // constructor
        const Elements& operator = (const Elements& Elems);            // explicit assignment operator; is needed because of the reference genome_v2 in the struct
        
        short add_element(short state, float init_dens, const node& base_element = node(-1,-1), bool SuppressReferencing_in_BaseElem = false);        // 'base_element' and 'SuppressReferencing_in_BaseElem' are only relevant for Phosphorylates. In that case, if 'SuppressReferencing_in_BaseElem' is true, the index of the Phosphorylate is not added to 'phosphs' of the base element. This is useful e.g. in 'Add_Phosphorylate' for practical reasons, however should not be used e.g. in 'Duplicate_Protein'.
        short delete_element(short index = -1);        // delete the element of index 'index' from 'elements'. If no index is provided (default value -1), choose an index randomly of one of the variable elements. Return 'false' if no element could have been deleted, otherwise return 'true'.
        void cap_InConnection(const arrow ar, short index);    // delete arrow 'ar' from 'ins' of the element with index 'index'. If a CompReaction is destroyed or the last element from 'ins' is being destroyed, trigger deletion of the element as it would otherwise be a zombie node.
        void cap_OutConnection(const arrow ar, short index);     // delete arrow 'ar' from 'outs' of the element with index 'index'
        void build_InConnection(const arrow ar, short index);  // add reference to arrow 'ar' to 'ins' of the element with index 'index'
        void build_OutConnection(const arrow ar, short index);   // add reference to arrow 'ar' to 'outs' of the element with index 'index'
        void delete_PhosphRef(const short ref_state, const short ref_index, const short elem_index);   // delete reference to the Phosphorylate with index 'ref_index' from the element with index 'elem_index' (stored in 'phosphs').
        
        void print_Elements(ofstream& file);
 //       friend ifstream& operator>>(ifstream& file, Elements& Elems);
        bool operator==(const Elements& Elems2);
    };
    friend ifstream& operator>>(ifstream& file, Elements& Elems);       // maybe reasonable to make this a member function of "Elements" instead and call it "Rebuild_Element" or something. Operator >> can only be defined as friend or normal function. Since "Elements" is a private member of Genome_Version2 it must be defined outside the scope of "Elements" as a friend of Genome_Version2.
    void reference_Phosphorylates();       // references each Phosphorylate in the corresponding base element in 'phosphs'. To be called in 'Rebuild_Genome' after all Elements have been read in.
    
    
///////////////////////////////////////////////////// Arrows //////////////////////////////////////////////////////////////////////////////////
    struct Arrows {         // implementation in "Arrows.cpp"
        Genome_Version2& genome_v2;         // refenrence from to the nested class to the mother class in order to call member functions from the mother class within the nested class, for example in DeleteArrow Cap_OutConnection and Cap_InConnection.
        const short type;                   // type of the arrows saved in the object (either Trans, CompReact or Dec)
        vector<short> empty;           // indices of elements (smaller than arrows.size()) that are at the moment unoccupied in 'arrows'
        vector<short> occupied;        // indices of elements (smaller than arrows.size()) that are at the moment occupied in 'arrows'
        vector<short> variable;       // indices of elements (smaller than arrows.size()) that are at the moment occupied and allowed to be modified, both structurally and by their rates, by the algorithm (via mutations).
        vector<short> rate_variable; // indices of elements (smaller than arrows.size()) that are occupied and allowed to be modified only by their rate.
        vector<Arrow> arrows;   // vector that saves all arrows.
        
        Arrows(short type, Genome_Version2& genome_v2, short reserve_size);         // constructor
        const Arrows& operator = (const Arrows& Arrs);            // explicit assignment operator; is needed because of the reference genome_v2 in the struct
        
        template<short N_source, short N_sink, class SourceArrayType, class SinkArrayType>
        void add_arrow(const SourceArrayType& source, const SinkArrayType& sink, float rate = -1, float back_rate = -1);

        bool delete_arrow(short index=-1);                 // delete the arrow with index 'index' from 'arrows', delete its reference from all sources and sinks via 'Cap_InConnection' and 'Cap_OutConnection' (these functions then decide whether the drop of the arrow-connection is lethal for the respective node, i.e. if the node must be destructed as a consequence). If no index is provided (default value -1), choose an index randomly of one of the variable arrows. Return 'false' if no arrow could have been deleted, otherwise return 'true'.
        bool mutate_rate(short index=-1);
        
        void print_Arrows(ofstream& file);
 //       friend ifstream& operator>>(ifstream& file, Arrows& Arrs);
        bool operator==(const Arrows& Arrs2);
    };
    friend ifstream& operator>>(ifstream& file, Arrows& Arrs);        // maybe reasonable to make this a member function of "Arrows" instead and call it "Rebuild_Arrow" or something. Operator >> can only be defined as friend or normal function. Since "Arrows" is a private member of Genome_Version2 it must be defined outside the scope of "Arrows" as a friend of Genome_Version2
    
///////////////////////////////////////////////////// Essential //////////////////////////////////////////////////////////////////////////////////
    struct Essential_class  {
        
  //    const short NumEssElems;
        Genome_Version2& genome_v2;
        array<vector<node>, NumEssFunc> ess_elements;
        
        Essential_class(Genome_Version2& genome_v2);                  // constructor
        const Essential_class& operator = (const Essential_class& Ess);     // explicit assignment operator; is needed because of the reference genome_v2 in the struct
        
        void declare_Essential(node nd, short ess_func);
        
        short add_Essential(short ess_func = -1);
        
        void undeclare_Essential(node nd, short ess_func);
        
        short delete_Essential(short ess_func = -1);
        
        void print_Essential(ofstream& file);
    };
    friend ifstream& operator>>(ifstream& file, Essential_class& Ess);
    void reference_essFunc();
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Elements Proteins;
    Elements Phosphorylates;
    Elements Complexes;
    
    Arrows Transforms;
    Arrows CompReactions;
    Arrows Decays;

    Essential_class Essential;
    vector<short> Mutations;
    
public:
    // Constructor
    Genome_Version2() :
        Proteins(Prot, *this, 10),
        Phosphorylates(Phos, *this, 18),
        Complexes(Comp, *this, 18),
        Transforms(Trans, *this, 50),
        CompReactions(CompReact, *this, 30),
        Decays(Dec, *this, 40),
        Essential(*this)
    {
    if(TrackMutations)
        Mutations.reserve(100);
    }
    
    void Create_Random_Initial_Network()        // function to create a random initial network; to be called directly after the Constructor. Formerly contained in the constructor but better to separate these two functions in order to use the constructor also in combination with "Rebuild_Genome"
    {
        //array<vector<node>, NumEssFunc>& ess_elements = Essential.ess_elements;
        Add_Protein();
        Add_Phosphorylate(node(-1,-1), membrane_bound, false);
        
        Essential.declare_Essential(node(Phos,Phosphorylates.occupied[0]), 0);
        
        for(short i=0; i<2; ++i)  {
            Add_Protein();
            Add_Binding_Domain();
            Add_Phosphorylate();
        }
#ifndef NDEBUG
        Test_All();
#endif
    };
    
    // explicit copy constructor: do not assign all the OdeSyst memebers as these need to be recalculated anyways when Fitness is evaluated the next time
    Genome_Version2(const Genome_Version2& genome_V2) : Proteins(genome_V2.Proteins), Phosphorylates(genome_V2.Phosphorylates), Complexes(genome_V2.Complexes), Transforms(genome_V2.Transforms), CompReactions(genome_V2.CompReactions), Decays(genome_V2.Decays), Essential(genome_V2.Essential), Mutations(genome_V2.Mutations)  {};
    
    // explicit assignment operator: do not assign all the OdeSyst memebers as these need to be recalculated anyways when Fitness is evaluated the next time
    const Genome_Version2& operator = (const Genome_Version2& genome_V2)  {
        Proteins = genome_V2.Proteins;
        Phosphorylates = genome_V2.Phosphorylates;
        Complexes = genome_V2.Complexes;
        Transforms = genome_V2.Transforms;
        CompReactions = genome_V2.CompReactions;
        Decays = genome_V2.Decays;
        Essential = genome_V2.Essential;
        Mutations = genome_V2.Mutations;
        return *this;
    };
    
private:
    inline Elements& ElementsDistr(short type) {
        switch(type)   {
            case Prot:
                return Proteins;
            case Phos:
                return Phosphorylates;
            case Comp:
                return Complexes;
            default:
                cerr << "ERROR in ElementsDistr()\n";
                abort();
        }
    };
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline Arrows& ArrowsDistr(short type) {
        switch(type)   {
            case Trans:
                return Transforms;
            case CompReact:
                return CompReactions;
            case Dec:
                return Decays;
            default:
                cerr << "ERROR in ArrowsDistributor()\n";
                abort();
        }
    };

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// these functions perform directly without calling to the memeber functions of 'Elements'. Use the below set of functions instead if 'Elements' shall be promoted to a data-capsulated class (instead of a public struct).
    inline Element& Get_Element(const node& node) {
        switch(node.type)   {
            case Prot:
                return Proteins.elements[node.index];
            case Phos:
                return Phosphorylates.elements[node.index];
            case Comp:
                return Complexes.elements[node.index];
            default:
                cerr << "ERROR in Get_Element()\n";
                abort();
        }
    };
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline void Delete_Node(const node node)  {
        ElementsDistr(node.type).delete_element(node.index);
    };
    
    inline void Build_InConnection(const arrow ar,const node node)  {
        Get_Element(node).ins.push_back(ar);
    };
    
    inline void Build_OutConnection(const arrow ar,const node node)  {
        Get_Element(node).outs.push_back(ar);
    };
    
    // only for Cap_InConnection use the respective member function of 'Elements' which performs additional tests and avoids the emergence of zombie-nodes.
    inline void Cap_InConnection(const arrow ar, const node node)  {
        ElementsDistr(node.type).cap_InConnection(ar, node.index);
    };
    
    inline void Cap_OutConnection(const arrow ar,const node node)  {
        VecDel_element(Get_Element(node).outs, ar);
    };
    
    inline void Build_PhosphRef(short ref_state, short ref_index, const node node)  {
        Get_Element(node).phosphs[ref_state].push_back(ref_index);
    };
    
    inline void Delete_PhosphRef(short ref_state, short ref_index, const node node)  {
        VecDel_element(Get_Element(node).phosphs[ref_state], ref_index);
    };
    
    inline Arrow& Get_Arrow(const arrow arrow) {
        switch(arrow.type)   {
            case Trans:
                return Transforms.arrows[arrow.index];
            case CompReact:
                return CompReactions.arrows[arrow.index];
            case Dec:
                return Decays.arrows[arrow.index];
            default:
                cerr << "ERROR in ArrowsDistributor/n";
                exit(1);
        }
    };
    
    inline void Delete_Arrow(const arrow arrow)  {
        ArrowsDistr(arrow.type).delete_arrow(arrow.index);
    };
 
// this set of functions calls the corresponding member functions of 'Elements' that perform the desired actions. They are thus more suitable if 'Elements' is to be promoted to a data-capsulated class, otherwise the above set is more directly and maybe more efficient.
    /*
     inline Element& Get_Element(const node& node) {
     return ElementsDistr(node.type).elements[node.index];
     };
     
     inline void Delete_Node(const node& node)  {
     ElementsDistr(node.type).delete_element(node.index);
     };
     
     inline void Build_InConnection(const arrow& ar, node& node)  {
     ElementsDistr(node.type).build_InConnection(ar, node.index);
     };
     
     inline void Build_OutConnection(const arrow& ar, node& node)  {
     ElementsDistr(node.type).build_OutConnection(ar, node.index);
     };
     
     inline void Cap_InConnection(const arrow& ar, node& node)  {
     ElementsDistr(node.type).cap_InConnection(ar, node.index);
     };
     
     inline void Cap_OutConnection(const arrow& ar, node& node)  {
     ElementsDistr(node.type).cap_OutConnection(ar, node.index);
     };
     
     inline void Delete_PhosphRef(short ref_state, short ref_index, node& node)  {
     ElementsDistr(node.type).delete_PhosphRef(ref_state, ref_index, node.index);
     };
     
     inline Arrow& Get_Arrow(const arrow& arrow) {
     return ArrowsDistr(arrow.type).arrows[arrow.index];
     };
     
     inline void Delete_Arrow(const arrow& arrow)  {
     ArrowsDistr(arrow.type).delete_arrow(arrow.index);
     };
     */
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public:
    size_t Get_NumberElements(short type=-1)  {
        if(type==-1)
            return Proteins.occupied.size() + Complexes.occupied.size() + Phosphorylates.occupied.size();
        else if(type==Prot)
            return Proteins.occupied.size();
        else if(type==Comp)
            return Complexes.occupied.size();
        else if(type==Phos)
            return Phosphorylates.occupied.size();
        else
            abort();
    };
    
    size_t Get_NumberOde1()  {
        if(BackRatesOn)
            return 4*Transforms.occupied.size() + 3*Decays.occupied.size() + 3*CompReactions.occupied.size();
        else
            return 2*Transforms.occupied.size() + 3*Decays.occupied.size();
    };
    size_t Get_NumberOde2()  {
        if(BackRatesOn)
            return 3*CompReactions.occupied.size() + 3*Decays.occupied.size();
        else
            return 3*CompReactions.occupied.size();
    };
    
    template<class vector_type>         // template this for the vector type that the Ode solver uses. For example ublas::vector for the implicit solvers of OdeInt
    void FillIn_InitialHomogeneousProteinConcentrations(vector_type& concentrations, short N_disc, bool WithPerturbation=true)  {
        for(short i=0; i<Proteins.elements.size(); ++i)
            if(Proteins.elements[i].existent == true)  {
       //         fill_n(concentrations.begin()+(prot_OdeIndex[i]*N_disc), N_disc, Proteins.elements[i].init_dens);   // use this to create unperturbed initial concentrations
                short j = 0;
                auto it_begin = concentrations.begin()+(prot_OdeIndex[i]*N_disc);   // fill in initial concentrations with a small perturbation of 10% of 'init_dens'
                for(auto it=it_begin; it!=it_begin+N_disc; ++it)  {
                    if(WithPerturbation)
#ifdef IDENTICAL_INITIAL_PERTURBATION
                        *it = Proteins.elements[i].init_dens * Initial_Perturbations[j++];
#else
                        *it = Proteins.elements[i].init_dens * (UnifRand()*0.2 + 0.9);      // multiply initial concentration by a random factor between 0.9 and 1.1, so add a small perturbation of 10% of 'init_dens' to the initial concentration
#endif
                    else
                        *it = Proteins.elements[i].init_dens;
                }
            }
    };
    
    template<typename T>
    void Make_ProbeProfiles(T& ProbeProfiles, float* concentrations_begin, short N_disc)  {
        for(short i=0; i<NumEssFunc; ++i)
            for(short j=0; j<Essential.ess_elements[i].size(); ++j)
                ProbeProfiles[i].push_back(concentrations_begin + Get_OdeIndex(Essential.ess_elements[i][j]) * N_disc);
    };
    
    void Save_Mutation(short executed, short identifier)  {
        if(TrackMutations && executed >= 0)
            Mutations.push_back(identifier);
    };

//////////////////////////////////////////// Mutation ///////////////////////////////////////////////////////////////////////////////////////
    bool Mutate_Rate();
    
    short Add_Phosphorylate(node node_parent = node(-1,-1), short state = -1, bool link_nodes_randomly = true, vector<node> nodes_linked={});
    short Add_Phosphorylate_randomized();    // the same as 'Add_Phosphorylate' when called without arguments (all arguments are on their default value and are thus chosen randomly)
    short Add_Protein();
    
    void Duplicate_SubNetwork(short dup_Prot, vector<short>& dup_Comps, vector<short>& dup_Phosphs);   // function very similar to 'Duplicate_Protein' but can duplicate an arbitrary subnetwork consisting of the elements given in 'dup_BaseElems' and 'dup_Phophs' (the corresponding connecting arrows are duplicated automatically). Relys on that each element appears at most once in the two vectors. Function duplicates the listed elements and all arrows that connect to those. For the details of the implementation of the nested functions look in the description of 'Duplicate_Protein'.
    short Duplicate_Protein_V1(short index = -1);       // duplicates the Protein with index 'index' or chooses a random index if no index is provided (default value -1) and duplicates all elements that form or transform from the Protein (the builders).
    short Duplicate_Protein_V2(short index = -1);       // duplicates the Protein with index 'index' or chooses a random index if no index is provided (default value -1) and duplicates all elements that form or transform from the Protein (the builders).
    
    node Choose_Config(Element& baseEl, short state = arbitrary, node exclude_node = node(-1,-1));        // chooses randomly a configuration of state 'state' of the base element 'baseEl' which is not the configuration 'exclude_node'. For that purpose, the function modifies temporarily the vector 'phosphs' of the base element (add base element and/or delete 'exclude_node') but resets all relevant changes before returning.
    short Add_Binding_Domain(bool EnzymaticReaction = ChooseInd(2));       // alternative name: add_Complex // chooses first two proteins (or if 'HigherOrderComplexes' == 'true' also a complex) as base element and then randomly pick a configuration of each one that define the two components of a new complex. Add that new complex as well as an CompReaction Arrow that connects the two components to the complex the a reverse Decay Arrow back into the components. If 'EnzymaticReaction' == 'true', add an additional Decay Arrow. One component (preferentially component1) is chosen as the substrate whereas the other one is chosen as the enzyme. The decay products of the additional Decay Arrow are then the enzyme (usually component2) and a different configuration of the substrate (usually component1). Only if this is not possible component1 instead of component2 may be chosen as the enzyme or the enzyme may also undergo a change in configuration.
    
    short Mutate_Concentration();
    
    short Mutate_Function();
    
    short Delete_Phosphorylate();
    short Delete_Complex();
    short Delete_Protein();
    
    void Mutate_constant(array<short,9> weights = {8,3,2,3,2,3,2,2,1});
    
    void Mutate_Rate_Concentration_Function();
    void Mutate_AddElement();
    void Mutate_DeleteElement();
    void Mutate_constant_inflationary()  {    // mutate with constant rates (except in 'Delete_Elements' where probability is proportional to number of elements) such that addition of an element is twice more frequent as deletion mutations however with upper limits (N_Prot, N_Phosph and N_Comp) for the number of elements as defined in 'Add_Element'. Therefore this mutation functions is infaltionary as it tends to increase the number of elements in the network.
        short weights[3] = {5,3,2};
        short d = ChooseInd_respective_weights(weights, 3);
        switch(d)  {
            case 0:
                Mutate_Rate_Concentration_Function();
                break;
            case 1:
                Mutate_AddElement();
                break;
            case 2:
                Mutate_DeleteElement();
                break;
        }
    };

    void Mutate(short times)  {            // the Mutate function that is called from all the Evolver objects and calls a specific mutate function
        for(short i=1; i<=times; ++i)  {
            Mutate_constant_inflationary();
        }
    };
    
//////////////////////////////////////////////////////// OdeSyst ///////////////////////////////////////////////////////////////////////////
    
    vector< vector<ode1> > OdeSyst1;
    vector< vector<ode2> > OdeSyst2;
    
    vector<short> prot_OdeIndex;
    vector<short> comp_OdeIndex;
    vector<short> phosph_OdeIndex;
    
    short Get_OdeIndex(const node no) {
        if(no.type == Phos)
            return phosph_OdeIndex[no.index];
        else if(no.type == Comp)
            return comp_OdeIndex[no.index];
        else
            return prot_OdeIndex[no.index];
    };
    
    void RangeThrough(function<void(Element&)> func);
    
    void Make_OdeSyst(vector<short>& states);
    
    void Transfer_directly_to_OdeSystCompressed(vector<short>& states, vector<ode1>& OdeSyst1_compressed, vector<short>& OdeSyst1_numbers, vector<ode2>& OdeSyst2_compressed, vector<short>& OdeSyst2_numbers);
    
////////////////////////////////////////////////////////  Print & Read  /////////////////////////////////////////////////////////////////////////
    void PrintGenome(ofstream& file);   // takes a file reference and and writes all relevnat data of Genome into it (all tab arrays as well as the rates and indices of OdeSyst1 and OdeSyst2)
    
    void Print_Mutations(ofstream& file);
    
    void RebuildGenome(ifstream& file);
        // assumes that all vectors are empty (cleared) which they should be if the class Genome has been created with the usual constructor (but not applied to "Create_Random_Initial_Network").
    
    void RebuildGenome(const string filename);
    
    void Print_OdeSyst(ofstream& file);              // Prints the OdeSystem into the delivered stream
    void Print_Biograph(ofstream& file);            // prints the data necessary for matlab to make a biograph object into the delivered stream.

    void Append_OdeSyst_Biograph(string path_from, string path_to);      // Matlab Interface. Copies a data file 'path_from' to another location (usually the matlab folder) and appends the copy by the OdeSystem (via call to 'Print_OdeSystem') and Biograph data (via call to 'Print_Biograph').
    
////////////////////////////////////////////////////// Tests //////////////////////////////////////////////////////////////////////////
//#ifndef NDEBUG      // SpeciesObject is needed also for StabilityAnalysis()
public:
    // for the implementation see 'Tests.cpp'
    void Test_Elements(Elements& Elems);
    void Test_Arrows(Arrows& Arrs);
    void Test_Essential();
    void Test_OdeSyst();
    
    void Test_All();

    vector<vector<short>> SpeciesObject;
    vector<float> TotalMasses;
    
    void Make_SpeciesObject();
    void Save_InitialConcentrations(vector<float>& concentrations);
    void Check_MassConservation(vector<float>& concentrations);
    
    void cut_empty();
    void Test_PrintRebuild();
//#endif
    
};



// TODO:
//
//
//

// maybe TODO:
// rename: arrow -> edge
// make source and sink a union of array<node,1> and array<node,2> each
// change array<short,2> to static_vector<short,2> for sink and source in 'Arrow'. Use typedef boost::container::static_vector<short,2> static_vector_2; in Preliminary!
// put definitions of Genome_Version2.hpp in cpp file

// Ask Fridtjof:
// reasonable values for initial concentrations in comparison to system length, Diffusion constatns, etc...
// when duplicating a protein, both proteins with half the initial concentration of the original protein??

#endif /* Genome_Version2_hpp */
