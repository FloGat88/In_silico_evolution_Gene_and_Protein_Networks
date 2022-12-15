#ifndef ____Genome__
#define ____Genome__


#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <cfloat>
#include <list>
#include <string>
#include <fstream>
#include <sstream>
// #include <boost/variant.hpp>
// #include <boost/numeric/odeint.hpp>

using namespace std;

constexpr short Number_Essential_Proteins = 3;
constexpr short N_P = 6;        // max number of genes or proteins
constexpr short N_GPC = 6;      // max number of gene-protein complexes
constexpr short N_PPC = 6;      // max number of protein-protein complexes
constexpr short N_Phosph = 6;   // max number of phosphorylates
constexpr short N_PPR = 6;      // max number of reactions
constexpr short N_tot = 2*N_P+N_GPC+N_PPC+N_Phosph;
constexpr float InitialRateMax = 1.0;

constexpr size_t OdeSyst1_MaxSize = 2*N_P + 4*N_GPC + 3*N_PPC + 7*N_Phosph;     //  constant equal to the maximal number of entries in Genome.OdeSyst1. The integer numbers in its definition arise from the number of entries that an element adds to Genome.OdeSyst1 when created (corresponding add_...-functions in Genome). constant relevant for Print_Genome and furthermore for reaction_diffusion_dynamics_fast as it contains the dimensionality of reaction_diffusion_dynamics_fast.OdeSyst1.
constexpr size_t OdeSyst2_MaxSize = 3*N_GPC + 3*N_PPC + 4*N_PPR;        // same as above just for OdeSyst2



union intfloat {
    short i;
    float f;
    intfloat(): i(0) {}
    intfloat(short d): i(d) {}
    intfloat(int d): i(d) {}
    intfloat(float fl): f(fl) {}
    intfloat(double fl): f(fl) {}
};


template<short N, short N_Col>
struct GeneticInstance
{
    array< array<intfloat, N_Col>, N > tab;
    vector<short> empty;
    vector<short> cons0;
    vector<short> cons01;
    vector<short> all;
    array< vector< array<short,2> >, N > builders;
    GeneticInstance();
    void Print_Instance(ofstream& file, vector<short>);
    void Rebuild_Instance(ifstream& file, short N_read, vector<short> float_entries);
};

struct TabAddress {
    short type;
    short index;
    short tab_pos;
    TabAddress(short type, short index, short tab_pos): type(type),index(index),tab_pos(tab_pos) {}
};

struct ode1 {
    float rate;
    short ind;
    short* add;
    ode1(float r, short i, short* a): rate(r),ind(i),add(a) {}
    ode1(float r, short i): rate(r),ind(i),add(NULL) {}
};

struct ode2 {
    float rate;
    short ind1;
    short ind2;
    short* add;
    ode2(float r, short i1, short i2, short* a): rate(r),ind1(i1),ind2(i2),add(a) {}
    ode2(float r, short i1, short i2): rate(r),ind1(i1),ind2(i2),add(NULL) {}
};

enum GeneticInfo {G=0, P=1, GPC=2, PPC=3, Phosph=4, PPR=5};



class Genome {
    
    friend int main();              // for debugging purposes declare main temporarily as friend; remove this later on again!

// Constructors and Addenda
public:
    Genome();                               // Constructor
    void Diversify();                       // Constructor Addendum1
    void Rebuild_Genome(ifstream& file);    // Constructor Addendum2
    void OdeSyst_ShiftAddresses();          // Copy Constructor Addendum

private:
// generally needed data and functions
    array<short,5> OdeIndexSummands;
    array<void (Genome::*)(short),4> DestroyInstances;
    
    array<intfloat, 15>&        Instances_tab(short, short);
    vector<array<short,2>>&     Instances_builders(short, short);
    vector<short>&              Instances_all(short);
    short&                      Instances_num_phos(short, short);
    inline bool                 Exists(short type, short index);
    
    template<class T>
    void VecDel_element(vector<T>&, const T&);
    
    short OdeIndex(const short&, const short&);
    
    void DestroyBuilders(vector<array<short,2>>&);
    
    template<size_t T>
    void SetComponent(short&, short&, array<short,T>, array<short,2>);
    
    short number_protein_components(short, short);
    
    short number_subunits(short, short);
    
    void ChooseDecayProduct(short, short, short, short&, short&);
    
    inline void get_Phosphorylates(vector<short>&, short, short);
    
    inline intfloat* TabAddress_to_Pt(TabAddress& tab_address);
            
    void PhosphorylationSwitch(short&, short&, TabAddress);
    
    
// Proteins (and Genes) related data and functions
    struct GeneticInstance<N_P , 3> Proteins;
public: vector<short> Essential_Proteins;
    
private:
    void Add_Protein(short = -1);
    void Proteins_MutateRate();
    void Delete_Protein(short = -1);
    
// GPComp related data and functions
    struct GeneticInstance<N_GPC , 13> GPComp;
    
    void Add_GPComp();
    void GPComp_MutateRate();
    void Delete_GPComp(short = -1);

// PPComp related data and functions
    struct GeneticInstance<N_PPC , 15> PPComp;
    
    void Add_PPComp();
    void PPComp_MutateRate();
    void Delete_PPComp(short = -1);

// Phosphorylates related data and functions
    struct GeneticInstance<N_Phosph , 15> Phosphorylates;
    array< vector<TabAddress>, N_Phosph> PhosphorylationSwitches;   // for each Phosphorylate create a vector that saves TabAddress tupels (vectors are filled by function PhosphorylationSwitch) that contain connection to ode_pos entries in PPComp, Phosphorylates or in PPReact whose decays have been Phnosphorylation-switched to the respective Phosphorylate. At deletion of the Phosphorylate, these ode_pos entries must be set to -1 because a Complex decaying to Phosphorylation-switched Phsphorylates may escape the deletion cascade of that Phosphorylate.
    void Print_PhosphorylationSwitches(ofstream& file);
    void Rebuild_PhosphorylationSwitches(ifstream& file, short N_Phosph_read);
    
    void Add_Phosphorylate();
    void Phosphorylates_MutateRate();
    void Delete_Phosphorylate(short = -1);
    
// PPReact related data and functions
    struct GeneticInstance<N_PPR, 15> PPReact;
    
    void Add_PPReact();
    void PPReact_MutateRate();
    void Delete_PPReact(short = -1);
 
    
// Functions that control mutation
    void MutateRate();
    void Add_Element();
    void Delete_Element();
    
public:
    void Mutate(int);

    
// ode_syst related data and functions
public:
    array< vector<ode1> , N_tot > OdeSyst1;
    array< vector<ode2> , N_tot > OdeSyst2;
    
private:
    inline void OdeSyst1_add(ode1, short);
    
    inline void OdeSyst2_add(ode2, short);
    
    inline void OdeSyst1_SetRate(const short&, const short&, const short&, float);
    
    inline void OdeSyst2_SetRate(const short&, const short&, const short&, float);
    
    inline void OdeSyst1_SetAddress(const short& type, const short& ind, short& ode_pos);
    
    inline void OdeSyst2_SetAddress(const short& type, const short& ind, short& ode_pos);
    
    template<size_t T>
    void OdeSyst_SetComplex(array<intfloat,T>&, short, short);
    
    inline void OdeSyst1_delete(const short&, const short&, const short&);
    
    inline void OdeSyst1_delete(const short&, const short&, const short&, const short&, const short&);
    
    inline void OdeSyst2_delete(const short&, const short&, const short&);
    
    inline void OdeSyst2_delete(const short&, const short&, const short&, const short&, const short&);
  
// Functions that regulate the writing of the Genome data to a file and saving of data
public:
    void PrintGenome(ofstream& file);
    
};


#endif /* defined(____Genome__) */
