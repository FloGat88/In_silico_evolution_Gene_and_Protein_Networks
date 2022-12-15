#include "RandomNumberGenerator.h"
#include "Genome.h"

using namespace std;

template<short N, short N_Col>
GeneticInstance<N,N_Col>::GeneticInstance()    // Constructor of the struct GeneticInstance. For cons0, cons01 and all reserve space for N entries and fill empty with the indices from N-1 to 0 in decreasing order. Corresponds to all the initial situatin where there are 0 elements of the instance and therefore all entries in tab are empty.
{
    for (size_t i=0; i<N; ++i)      // set the zeroth entry in the tab for all elements to -1.0 indicating that the element is empty (this is used e.g. in the function 'Exists')
        tab[i][0] = -1.0;
    
    cons0.reserve(N);
    cons01.reserve(N);
    all.reserve(N);
    
    empty.resize(N);
    for (short i = 0; i<N; ++i)
        empty[i] = (N-1)-i;
}


Genome::Genome()     // Konstructor
{
    OdeIndexSummands = {0, N_P, 2*N_P, 2*N_P+N_GPC, 2*N_P+N_GPC+N_PPC};
    DestroyInstances = {&Genome::Delete_GPComp, &Genome::Delete_PPComp, &Genome::Delete_Phosphorylate, &Genome::Delete_PPReact};
    
    for (short i=0; i<Number_Essential_Proteins; i++)   // for each Essential Protein in the System fill one index in the vector Essential_Proteins. The Proteins remain however unoccupied for now.
        Essential_Proteins.push_back(i);
}


void Genome::Diversify()
{
    for (short i=0; i<Number_Essential_Proteins; i++)    // for each Essential Protein in the System add one Protein.
        Add_Protein();
        
    // add some Elements to create some initial diversity
    if (Proteins.empty.size() != 0)
        Add_Protein();
    for (size_t i=1; i<=2; ++i) {
        if (GPComp.empty.size() != 0)
            Add_GPComp();
        if (PPComp.empty.size() != 0)
            Add_PPComp();
        if (Phosphorylates.empty.size() != 0)
            Add_Phosphorylate();
        if (PPReact.empty.size() != 0)
            Add_PPReact();
            
    }
}


short Genome::OdeIndex(const short& type, const short& ind)          // function that translates the table index 'ind' ind of type 'type' in the corresponding index in OdeSyst by adding the corresponding 'offset' that are stored in array OdeIndexSummands
{
    return OdeIndexSummands[type] + ind;
}


array<intfloat, 15>& Genome::Instances_tab(short type, short index)
{
    switch(type)
    {
        case PPC:
            return PPComp.tab[index];
        case Phosph:
            return Phosphorylates.tab[index];
        case PPR:
            return PPReact.tab[index];
        default:
            cerr << "ERROR: no match in Instances_tab for type=" << type << '\n';
            exit(2);
    }
}


vector<array<short,2>>& Genome::Instances_builders(short type, short index)
{
    switch(type)
    {
        case G:
        case P:
            return Proteins.builders[index];
        case GPC:
            return GPComp.builders[index];
        case PPC:
            return PPComp.builders[index];
        case Phosph:
            return Phosphorylates.builders[index];
        default:
            cerr << "ERROR: no match in Instances_builders" << '\n';
            exit(2);
    }
}


vector<short>& Genome::Instances_all(short type)
{
    switch(type)
    {
        case G:
        case P:
            return Proteins.all;
        case GPC:
            return GPComp.all;
        case PPC:
            return PPComp.all;
        case Phosph:
            return Phosphorylates.all;
        case PPR:
            return PPReact.all;
        default:
            cerr << "ERROR: no match in Instances_all" << '\n';
            exit(2);
    }
}


short& Genome::Instances_num_phos(short type, short index)
{
    switch(type)
    {
        case P:
            return Proteins.tab[index].back().i;
        case PPC:
            return PPComp.tab[index].back().i;
        case Phosph:
            return Phosphorylates.tab[index].back().i;
        default:
            cerr << "ERROR: no match in Instances_num_phos" << '\n';
            exit(2);
    }
}


inline bool Genome::Exists(short type, short index)            // function that checks whether the element (type, index) still exists by probing if the zeroth entry of its .tab array is >= 0 (at deletion the zeroth entry is set to -1.0)
{
    switch(type)
    {
        case G:
        case P:
            return (Proteins.tab[index][0].f>=0);
        case GPC:
            return (GPComp.tab[index][0].f>=0);
        case PPC:
            return (PPComp.tab[index][0].f>=0);
        case Phosph:
            return (Phosphorylates.tab[index][0].f>=0);
        case PPR:
            return (PPReact.tab[index][0].f>=0);
        default:
            cerr << "ERROR: no match in 'Exists' for " << type << ' ' << index << '\n';
            exit(2);
    }
}


template<class T>
void Genome::VecDel_element(vector<T>& vec, const T& val)          // VectorDelete: in vector vec find (the last) element that is identical to val and delete it, i.e. transfer the last vector element on its position and pop_back
{
    typename vector<T>::iterator it;
    it = vec.end()-1;
    while (*it != val)                              // search vector from behind until we find an element that is identical to val
    {
        --it;
        if (it<vec.begin()) {cout << "ERROR" << '\n'; break;}   // !!! Control mechanism: to be removed later !!!
    }
    *it = vec.back();                       // transfer last element to that position
    vec.pop_back();                         // delete last element via pop_back
}

void Genome::DestroyBuilders(vector<array<short,2>>& builders)      // Function that destroys all builders contained in a builders vector
{
    for (short i=0;i<builders.size();++i)
    {
        if ( Exists(builders[i][0],builders[i][1]) )                // check if the element still exists and has not already been destroyed in the deletion cascade
            (this->*DestroyInstances[builders[i][0]-2])(builders[i][1]);   // feeding the elements of builders to DestroyInstances that is an array of function pointers that point to the corresponding delete-commands; builders[i][0]-2 is necessary because DestroyInstances has only 4 entries because G and P don't appear in any builders vector (therefore DeleteGPComp has index 0 and so on)
    }
}


template<size_t T>
inline void Genome::SetComponent(short& type_selected, short& index_selected, array<short, T> Types, array<short,2> Builder)   // template function that takes an array Types that contains the instances among which the function randomly picks type and index of one instance element that can form one component of a complex or a reaction. The chosen type and index are saved to the given references type_selected and index_selected, respectively. The identity of the complex or reaction (handed over to the function in the input array 'Builder') is added to the .builder vector of the chosen component.
{
    short d;
    array<short,T> numbers;    // array of equal size as the array Types that in a moment will contain the numbers of the available elements of the instance types listed in Types (.all.size())
    short sum = 0;              // short that will save the total number of elements of all instance types
    
    for (short i=0; i<Types.size(); ++i)  {                  // fill the array numbers and calculate sum
        numbers[i] = Instances_all(Types[i]).size();
        sum += numbers[i];
    }

    d = ChooseInd(sum);                   // use ChooseInd to randomly pick an integer number between 0 and (sum - 1)
    
    for (short i=0; i<Types.size(); i++)   {               // we subtract numbers[.] from d so long until the condition d<numbers[i] is satisfied for some i, then we pick the i-th type and the index as d (attention that d can take values between 0 and (sum of all - 1)! therefore check for '<' in the if-statement instead '<=' and the calculation of index_selected is straight forward)
        if (d < numbers[i])   {
            type_selected = Types[i];
            index_selected = Instances_all(type_selected)[d];
            break;
        }
        else
             d -= numbers[i];
    }
    
    Instances_builders(type_selected,index_selected).push_back(Builder);    // add the Builder to the .builders vector of the chosen component
}


short Genome::number_protein_components(short type, short index)      // function that returns the number of Proteins that the element (type, index) consists of (= 1 if type == Proteins)
{
    while (type==Phosph) {                                              // if the complex is a phosphorylate, find the corresponding basis element (change type and index to the type and index of the basis element that is then either P or PPComp)
        type = Phosphorylates.tab[index][1].i;
        index = Phosphorylates.tab[index][2].i;
    }
    if (type == P)                                                      // for Proteins return 1, for PP Complexes return the corresponding entry in .tab containing the number of Protein components
        return 1;
    else
        return PPComp.tab[index][7].i;
}


void Genome::ChooseDecayProduct(short Complex_type, short Complex_index, short d, short& DecProd_type, short& DecProd_index)   // Function that chooses a decay product for a Complex (Complex_type, Complex_index) (either a PPComp or a Phosphorylate (phosphorylated PPComp, not Protein)) using a random number d between 2 and 2N-1 where N is the number of Proteins that build the complex. 2N-1 is the number of nodes if you draw the tree of the complex starting from the corresponding proteins. d=1 would correspond to the complex itself and is therefore excluded, then we count the branch of component 1 (with 2*N_1-1 subunits) and proceed with the branch of component 2 in this numeration. Thus d=2 corresponds to component 1, whereas d = (2*N_1-1) + 2 (the number 2 counts the complex and component 2) corresponds to component 2. If none of these is the case, the function recursively calls itself handing over either component 1 (if d<= 2*N_1-1 + 1) or component 2 otherwise and adapts d accordingly.
{
    while (Complex_type==Phosph) {                                              // if the complex is a phosphorylate, find the corresponding basis element (change Complex_type-and index to the type and index of the basis element that is then either P or PPComp) and replace R by the .tab vector of that basis element
        Complex_type = Phosphorylates.tab[Complex_index][1].i;
        Complex_index = Phosphorylates.tab[Complex_index][2].i;
    }
    array<intfloat,15>& R = PPComp.tab[Complex_index];        // create Reference on the concerning .tab vector which now points to a PPComp.tab vector.
    
    short CountLeft = 2*number_protein_components(R[1].i,R[2].i)-1;     // 2*N_1-1 where N_1 is the number of Protein components of the first component (1 for Proteins and R[7] for PP Complexes; Phosphorylates must first be reduced to the basis element)
    if(d==2)  {                                                                         // if d=2 decay product is component 1
        DecProd_type = R[1].i;
        DecProd_index = R[2].i;
    }
    else if(d==(CountLeft + 2))  {                                                  // if d=CountLeft + 2 decay product is component 2
        DecProd_type = R[4].i;
        DecProd_index = R[5].i;
    }
    else if (d  <= CountLeft+1)                                                     // otherwise recursive function call delivering either component1 or 2 and adapting d accordingly
        ChooseDecayProduct(R[1].i, R[2].i, d-1, DecProd_type, DecProd_index);
    else
        ChooseDecayProduct(R[4].i, R[5].i, d-(CountLeft+1), DecProd_type, DecProd_index);
}


inline void Genome::get_Phosphorylates(vector<short>& configurations, short type, short index)    // recursive function that checks for Phosphorylated versions of the delivered agent (type, index) (by searching its .builders vector) and saves the indices of all those phosphorylated versions that it finds in the reference vector configurations before the function calls itself with those phosphorylated versions such that at the end of the recursion 'configurations' contains the indices of all Phosphorylates of the initial element (ususally a basis element)
{
    short num_phos = Instances_num_phos(type,index);                    // short that saves the number of direct Phosphorylates in association to the respective agent (type, index)
    if (num_phos > 0) {                                                     // not necessary
        vector<array<short,2>>& B = Instances_builders(type,index);       // create reference to the .builders vector of the corresponding element for faster access and better readability
        for (short i=0; /*i<B.size(),*/num_phos>0; i++)                     // search the builders vector for all elements of type Phosph and save these in 'configurations' (stop if num_phos, which counts the remaining Ponsphorylates in builders, becomes 0)
            if (B[i][0]==Phosph) {
                num_phos--;
                configurations.push_back(B[i][1]);
                get_Phosphorylates(configurations, B[i][0], B[i][1]);      // for each builder of type Phosph call get_Phosphorylates to iterate and also capture the Phosphoryates of the Phosphorylate itslef
            }
    }
}


inline intfloat* Genome::TabAddress_to_Pt(TabAddress& tab_address)    // Function that converts a TabAddress (tuple of three shorts describing an element in a tab array) of either PPComp, PPreact or Phosphorylates to an intfloat pointer pointing at this tab-entry. Function is used for the vector PhosphorylationSwitches to convert its entries to pointers.
{
    return &Instances_tab(tab_address.type, tab_address.index)[tab_address.tab_pos];
}


void Genome::PhosphorylationSwitch(short& type, short& index, TabAddress tab_address)    // Function that randomly and with equal prob. picks one out of ALL configurations (also twice, ect. phosphorylated...) of an agent with type 'type' and index 'index' and exchanges type and index with those of the chosen configuration.
{
    vector<short> configurations;                       // vector that is filled with the indices of all the phosphorylates to a certain basis type.
    short basetype = type;
    short baseindex = index;
    
    while (basetype == Phosph) {                            // reduce the delivered element to its basis element (the unphosphorylated version of the element) if it is not already basic.
        basetype = Phosphorylates.tab[baseindex][1].i;
        baseindex = Phosphorylates.tab[baseindex][2].i;
    }
    
    get_Phosphorylates(configurations, basetype, baseindex);   // with type and index now pointing to an basis element (unphosphorylated element) call get_Phosphorylates that fills the indices of all phosphorylated versions of the basis element into the vector 'configurations'
    
    short d = ChooseInd(configurations.size()+1);       // finally, choose one of all configurations, where d=configurations.size() results to the choice of the unphosphorylated version (if-statement false, therfore no further change of type and index) and d<configurations.size() chooses the Phosphorylate at position d in configurations.
    if ( d<configurations.size() && (type!=Phosph || index!=configurations[d]) )   {     // if decision is for a phosphorylated version that differs from the submitted element
        type = Phosph;
        index = configurations[d];
        PhosphorylationSwitches[index].push_back(tab_address);
    }
    else if ( d==configurations.size() && type!=basetype )   {       // if decision is for base type that differs from the submitted element
        type = basetype;
        index = baseindex;
        // no need to add something to PhosphorylationSwitches because here we have a switch from phosphorylated version to base type which is save b/c if the base type gets destroyed also the Complex will be killed in its deletion cascade via the Phosphorylate (that is also killed)
    }
}



void Genome::Proteins_MutateRate()      // mutates either the expression or the decay rate of a randomly picked Protein by multiplying the rate with a random factor in [0,2]
{
    short p = ChooseInd(2);                                   // decide whether to change expression rate (p=0) or decay rate (p=1)
    short ind = Proteins.cons01[ChooseInd(Proteins.cons01.size())];     // choose index of Protein whose rate will be changed
    float newrate = Proteins.tab[ind][p].f *= UnifRand() * 2;                                      // change rate by multiplying it with a unif. random number in [0,2]
    if (p==0)                                        // modify that rate also in OdeSyst1
        OdeSyst1_SetRate( P, ind, 0, newrate);      // if decision for expression rate (p=0)
    else
        OdeSyst1_SetRate( P, ind, 1, -newrate);     // if decision for decay rate (p=1)
}


void Genome::Add_Protein(short ind)             // adds a new protein to the system at index ind (default value = -1). If no index is delivered the function picks a random index
{
    if (ind<0)                      //  if no index has been deliverd (default = -1) determine the index of the new protein where the Proteins matrix is empty right now
        ind = Proteins.empty.back();
    
    Proteins.tab[ind] = {UnifRand()*InitialRateMax, UnifRand()*InitialRateMax, 0};     // randomly choose the rates in [0,1] and enter the data in Proteins
    
    OdeSyst1_add({Proteins.tab[ind][0].f , OdeIndex(G,ind)}, OdeIndex(P,ind));  // modify OdeSyst1 for the expression (column 0) of the new Protein
    OdeSyst1_add({-Proteins.tab[ind][1].f , OdeIndex(P,ind)}, OdeIndex(P,ind));     // modify OdeSyst1 for the decay (column 1) of the new Protein
    
    Proteins.cons0.push_back(ind);        // register the new Protein in cons0, cons01, all and delete its index from empty
    Proteins.cons01.push_back(ind);
    Proteins.all.push_back(ind);
    Proteins.empty.pop_back();
}


void Genome::Delete_Protein(short index)    // deletes Protein with index 'index'; if the function is called without an index (default value -1) it will randomly pick a Protein
{
    if (index<0)                                        // if no index has been submitted randomly choose index of a Protein that shall be deleted
        index = Proteins.cons0[ChooseInd(Proteins.cons0.size())];
    
    OdeSyst1[OdeIndex(G,index)].clear();                // clear vector of the concerning Gene from OdeSyst1...
    OdeSyst2[OdeIndex(G,index)].clear();                // ... and from OdeSyst2    !!! noch prüfen ob wirklich nötig dh. ob GPComp wirklich an PPReact teilnehmen können
    OdeSyst1[OdeIndex(P,index)].clear();             // clear vector of concerning Protein from OdeSyst1...
    OdeSyst2[OdeIndex(P,index)].clear();             // ... and from OdeSyst2
    
    Proteins.tab[index][0] = -1.0;                                               // 4. deletion of .tab entry of the Phosphorylate is not necessary, instead only the zeroth entry is set to -1.0 to indicate that this element has been deleted (the function Destroy_Builders checks for negative zeroth element and destroys the element only if it is >= 0 indicating that the element has not already been deleted, similarly for OdeSyst_Delete)
    
    VecDel_element(Proteins.cons0,index);                 // modify the status vectors: remove the index of the destroyed Protein from cons0, cons01 and all and add it to empty
    VecDel_element(Proteins.cons01,index);
    VecDel_element(Proteins.all,index);
    Proteins.empty.push_back(index);
    
    DestroyBuilders(Proteins.builders[index]);      // call the DestroyBuilders function that calls the deletion commands for all builders on the Proteins stored in the vector Protein.builders
    
    Proteins.builders[index].clear();              // clear the builders vector of the deleted Protein
    
    for (short i=0; i<Essential_Proteins.size(); ++i)     // check if one of the essential Proteins has been deleted. If so, add a new Protein at that index. Note that with this choice (in contrast to the one below) the Essential Proteins will always be at fixed indices. Note also that if Essential Proteins are declared as cons1 in the Constructor, they will never be deleted and this passage will never be executed.
        if (Essential_Proteins[i] == index)  {
                Add_Protein(index);
//                break;          // we can break here if all Essential Proteins have different (fixed) indices
        }
    
    
/*    // this choice for Essential_Protein handling can set more essential Proteins on one and the same index which turns out to be indeed some attractor for the mutational dynamics, however, not the original intention. So rather choose the above alternative instead!
 
    for (short i=0; i<Essential_Proteins.size(); ++i)     // check if one of the essential Proteins has been deleted. If so, either add a new Protein at that index (especially do this if the Protein number is now 0) or choose some other already existing Protein as new essential Protein (note that in this case another essential Protein can be picked another time so that it fulfills two or more properties/functions of the essential Proteins). Alternatively one could declare the Essential Proteins as cons1 in the Constructor, then they will never be deleted
        if (Essential_Proteins[i] == index)  {
            if (Proteins.all.size()==0 || ChooseInd(2) == 0)   {
                Add_Protein(index);
                break;          // break so that maximally one new Protein is added at the index. note that other essential proteins sharing this index are served as well which is not the case if the else branch is realised, therefore we have no break there.
            }
            else
                Essential_Proteins[i] = Proteins.all[ChooseInd(Proteins.all.size())];
        }
 */
}


void Genome::GPComp_MutateRate()
{
    short p = ChooseInd(3);                                                 // decide whether to change reaction rate (p=0), expression rate (p=1) or decay rate (p=2)
    short index = GPComp.cons01[ChooseInd(GPComp.cons01.size())];           // choose index of GPComp whose rate will be changed
    array<intfloat,13>& R = GPComp.tab[index];                              // create reference on the corresponding row in GPComp.tab to make the following more readable
    float newrate;
    
    if (p==0)                                                               // decision for reaction rate
    {
        newrate = R[0].f = 2*UnifRand() * R[0].f;       // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst2_SetRate( R[1].i, R[2].i, R[3].i, -newrate);     // Modify rate for first component in OdeSyst2 using the function OdeSyst2_SetRate
        OdeSyst2_SetRate( R[4].i, R[5].i, R[6].i, -newrate);       // Modify rate for second component in OdeSyst2
        OdeSyst2_SetRate( GPC, index, 0, newrate);                  // Modify rate for formation of the GP Complex in OdeSyst2
    }
    else if (p==1)                                                          // decision for expression rate
    {
        newrate = R[7].f = 2*UnifRand() * R[7].f;       // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst1_SetRate( P, R[8].i, R[9].i, newrate);   // modify rate for the expression of the corresponding gene in OdeSyst1
        
    }
    else                                                                    // decision for decay rate
    {
        newrate = R[10].f = 2*UnifRand() * R[10].f;         // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst1_SetRate( R[1].i, R[2].i, R[11].i, newrate);        // modify rate for first decay product (= first component) in OdeSyst1
        OdeSyst1_SetRate( R[4].i, R[5].i, R[12].i, newrate);        // modify rate for second decay product (= second component) in OdeSyst1
        OdeSyst1_SetRate( GPC, index, 0, -newrate);                 // modify rate for decay process for the GP Complex in OdeSyst1
    }
}


void Genome::Add_GPComp()
{
    short index = GPComp.empty.back();              // determine the index of the new GPComp where GPComp.tab is empty right now
    array<intfloat,13>& R = GPComp.tab[index];      // create reference on the corresponding row in GPComp.tab to make the following more readable
    R[0] = UnifRand() * InitialRateMax;                              // randomly choose the reaction rate in [0,1)
    R[7] = UnifRand() * InitialRateMax;                              // randomly choose the expression rate in [0,1)
    R[10] = UnifRand() * InitialRateMax;                             // randomly choose the decay rate in [0,1)
    
    SetComponent(R[1].i,R[2].i,array<short, 2> {G,GPC},{GPC,index});                  // set type and index of the first and second component by using the function SetComponent (that also adds the identity of the GPComplex to the builders vector of the chosen component)
    SetComponent(R[4].i,R[5].i,array<short, 3> {P,PPC,Phosph},{GPC,index});
    
    if (R[1].i==G)                                    // to fill in the index of the expressed Gene (R[8]) check if first component is a Gene...
        R[8] = R[2].i;                                // ... then the corresponding index is the same as the index of the component
    else                                            // or if it is a GP Complex...
        R[8] = GPComp.tab[R[2].i][8].i;               // ... then the index can be taken from the corresponding entry of GPComp.tab of the corrsponding GP complex
            
    OdeSyst_SetComplex(R,GPC,index);             // call OdeSyst_SetComplex that sets up all the dynamics in OdeSyst and fills the referrals in R thereon
    
    GPComp.cons0.push_back(index);                  // register the new GP Complex in cons0, cons01, all and delete its index from empty
    GPComp.cons01.push_back(index);
    GPComp.all.push_back(index);
    GPComp.empty.pop_back();
}


void Genome::Delete_GPComp(short index)    // deletes GP Complex with index 'index'; if the function is called without an index (default value -1) it will randomly pick a GP Comp
{
    if (index<0)                                        // 0. if no index has been submitted randomly choose index of a GPComp that shall be deleted
        index = GPComp.cons0[ChooseInd(GPComp.cons0.size())];
    array<intfloat,13>& R = GPComp.tab[index];
    
    // 1./2. delete coresponding entries from OdeSyst1 and OdeSyst2 (one to each Ode-Column-Index in the GPComp.tab-row) and for Component1 and 2 also destroy the Complex from the .builders vectors. Use the function OdeSyst_delete (two overloaded methods)
    OdeSyst2_delete(R[1].i,R[2].i,R[3].i,GPC,index);          // 1./2. clear for the complex forming reaction the dynamic contributions to component 1 and 2 (and earase the GPcomp from their builders)
    OdeSyst2_delete(R[4].i,R[5].i,R[6].i,GPC,index);
    
    OdeSyst1_delete(P,R[8].i,R[9].i);                       // 1./2. earase dynamic contribution due to gene expression
    
    OdeSyst1_delete(R[1].i,R[2].i,R[11].i);                // 1./2. earase contribution to dynamics of decay product due to decay of the GPComp
    OdeSyst1_delete(R[4].i,R[5].i,R[12].i);
    
    OdeSyst1[OdeIndex(GPC,index)].clear();                // 3. clear vector of the concerning PPComp from OdeSyst1 and from OdeSyst2
    OdeSyst2[OdeIndex(GPC,index)].clear();
    
    R[0] = -1.0;                                               // 4. deletion of .tab entry of the Phosphorylate is not necessary, instead only the zeroth entry is set to -1.0 to indicate that this element has been deleted (the function Destroy_Builders checks for negative zeroth element and destroys the element only if it is >= 0 indicating that the element has not already been deleted, similarly for OdeSyst_Delete)
    
    VecDel_element(GPComp.cons0,index);                 // 5. modify the status vectors: remove the index of the destroyed GPComp from cons0, cons01 and all and add it to empty
    VecDel_element(GPComp.cons01,index);
    VecDel_element(GPComp.all,index);
    GPComp.empty.push_back(index);
    
    DestroyBuilders(GPComp.builders[index]);      // 6. call the DestroyBuilders function that calls the deletion commands for all builders on the GPComp stored in the vector GPComp.builders
    
    GPComp.builders[index].clear();              // 7. clear the builders vector of the deleted GPComp
}


void Genome::PPComp_MutateRate()
{
    short p = ChooseInd(2);                                                 // decide whether to change reaction rate (p=0) or decay rate (p=1)
    short index = PPComp.cons01[ChooseInd(PPComp.cons01.size())];           // choose index of PPComp whose rate will be changed
    array<intfloat,15>& R = PPComp.tab[index];                              // create reference on the corresponding row in PPComp.tab to make the following code more readable
    float newrate;
    
    if (p==0)                                                               // decision for reaction rate
    {
        newrate = R[0].f = 2*UnifRand() * R[0].f;               // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst2_SetRate( R[1].i, R[2].i, R[3].i, -newrate);     // Modify rate for first component in OdeSyst2 using the function OdeSyst2_SetRate
        OdeSyst2_SetRate( R[4].i, R[5].i, R[6].i, -newrate);       // Modify rate for second component in OdeSyst2
        OdeSyst2_SetRate( PPC, index, 0, newrate);                  // Modify rate for formation of the GP Complex in OdeSyst2
    }
        else                                                                    // decision for decay rate
    {
        newrate = R[8].f = 2*UnifRand() * R[8].f;                   // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst1_SetRate( PPC, index, 0, -newrate);                 // modify rate for decay process for the PP Complex in OdeSyst1
        if (R[9].i==1)                                                // if decay_mode == 1, i.e. decay into a single decay product
            OdeSyst1_SetRate( R[10].i, R[11].i, R[12].i, newrate);        // modify rate for decay product in OdeSyst1
        if (R[9].i==2) {                                                   // if decay_mode == 2, i.e. decay into the two components of the PPComp
            OdeSyst1_SetRate( R[1].i, R[2].i, R[12].i, newrate);        // modify rate for first decay product (= first component) in OdeSyst1
            OdeSyst1_SetRate( R[4].i, R[5].i, R[13].i, newrate);        // modify rate for second decay product (= second component) in OdeSyst1
        }
    }
}


void Genome::Add_PPComp()
{
    short index = PPComp.empty.back();              // determine the index of the new PPComp where PPComp.tab is empty right now
    array<intfloat,15>& R = PPComp.tab[index];      // create reference on the corresponding row in PPComp.tab to make the following more readable
    R[0] = UnifRand() * InitialRateMax;                              // randomly choose the reaction rate in [0,1)
    R[8] = UnifRand() * InitialRateMax;                             // randomly choose the decay rate in [0,1)
    
    SetComponent(R[1].i,R[2].i,array<short, 3>{P,PPC,Phosph},{PPC,index});           // set type and index of the first and second component by using the function SetComponent (that also adds the identity of the GPComplex to the builders vector of the chosen component)
    SetComponent(R[4].i,R[5].i,array<short, 3>{P,PPC,Phosph},{PPC,index});
    R[7] = number_protein_components(R[1].i,R[2].i) + number_protein_components(R[4].i,R[5].i);
    
// determin decay mode and product(s)
    R[9] = ChooseInd(3);                            // determine the decay mode as either 0,1 or 2
    if (R[9].i==1)   {                                   // if decay mode == 1 use ChooseDecayProduct to choose one protein or subcomplex of the PP Complex as its decay product
        ChooseDecayProduct(PPC,index,ChooseInd(2*R[7].i-2)+2, R[10].i,R[11].i);    // call ChooseDecayProduct for the complex delivering a random number d between 2 and 2N-1 (N being R[7], the nmber of involved proteins)
        PhosphorylationSwitch(R[10].i,R[11].i,{PPC,index,10});     // randomly switch to some other phosphorylate (configuration) of the chosen decay product if existent
    }
    else if (R[9].i==2) {                             // if decay mode == 2 just copy type and index of the first component to the corresponding fields of the first decay product
        R[10] = R[1].i;
        R[11] = R[2].i;
    }
    R[14] = 0;                                      // set number of Phosphorylates to 0
    
    OdeSyst_SetComplex(R,PPC,index);             // call OdeSyst_SetComplex that sets up all the dynamics in OdeSyst and fills the referrals in R thereon
    
    PPComp.cons0.push_back(index);                  // register the new PP Complex in cons0, cons01, all and delete its index from empty
    PPComp.cons01.push_back(index);
    PPComp.all.push_back(index);
    PPComp.empty.pop_back();
}

    
void Genome::Delete_PPComp(short index)    // deletes PP Complex with index 'index'; if the function is called without an index (default value -1) it will randomly pick a PP Comp
{
    if (index<0)                                        // 0. if no index has been submitted randomly choose index of a PPComp that shall be deleted
        index = PPComp.cons0[ChooseInd(PPComp.cons0.size())];
    array<intfloat,15>& R = PPComp.tab[index];
    
    // 1./2. delete coresponding entries from OdeSyst1 and OdeSyst2 (one to each Ode-Column-Index in the PPComp.tab-row) and for Component1 and 2 also destroy the Complex from the .builders vectors. Use the function OdeSyst_delete (two overloaded methods)
    OdeSyst2_delete(R[1].i,R[2].i,R[3].i,PPC,index);          // 1./2. clear for the complex forming reaction the dynamic contributions to component 1 and 2 (and earase the PPcomp from their builders)
    OdeSyst2_delete(R[4].i,R[5].i,R[6].i,PPC,index);
    
    if (R[9].i != 0)                                         // 1./2. earase contribution to dynamics of decay product due to decay of the PPComp (pay attention on the respective decay mode)
        OdeSyst1_delete(R[10].i,R[11].i,R[12].i);
    if (R[9].i == 2)
        OdeSyst1_delete(R[4].i,R[5].i,R[13].i);
    
    OdeSyst1[OdeIndex(PPC,index)].clear();                // 3. clear vector of the concerning PPComp from OdeSyst1 and from OdeSyst2
    OdeSyst2[OdeIndex(PPC,index)].clear();
    
    R[0] = -1.0;                                               // 4. deletion of .tab entry of the Phosphorylate is not necessary, instead only the zeroth entry is set to -1.0 to indicate that this element has been deleted (the function Destroy_Builders checks for negative zeroth element and destroys the element only if it is >= 0 indicating that the element has not already been deleted, similarly for OdeSyst_Delete )
    
    VecDel_element(PPComp.cons0,index);                 // 5. modify the status vectors: remove the index of the destroyed GPComp from cons0, cons01 and all and add it to empty
    VecDel_element(PPComp.cons01,index);
    VecDel_element(PPComp.all,index);
    PPComp.empty.push_back(index);
    
    DestroyBuilders(PPComp.builders[index]);      // 6. call the DestroyBuilders function that calls the deletion commands for all builders on the PPComp stored in the vector PPComp.builders
    
    PPComp.builders[index].clear();              // 7. clear the builders vector of the deleted PPComp
}


void Genome::PPReact_MutateRate()
{
    short index = PPReact.cons01[ChooseInd(PPReact.cons01.size())];           // choose index of Reaction whose rate will be changed
    array<intfloat,15>& R = PPReact.tab[index];                              // create reference on the corresponding row in PPReact.tab to make the following code more readable
    float newrate;
    newrate = R[0].f = 2*UnifRand() * R[0].f;                               // change rate by multiplying it with a unif. random number in [0,2]
    
    if (R[4].i != -1)                                                       // if reaction is a two component reaction (in that case type of second component is >=0, otherwise it would be -1)
    {
        OdeSyst2_SetRate( R[1].i, R[2].i, R[3].i, -newrate);     // Modify rate for first component in OdeSyst2 using the function OdeSyst2_SetRate
        OdeSyst2_SetRate( R[4].i, R[5].i, R[6].i, -newrate);       // Modify rate for second component in OdeSyst2
        if (R[8].i!=0)                                                // if mode is not 0, i.e. there is at least a single product
            OdeSyst2_SetRate( R[9].i, R[10].i, R[11].i, newrate);        // modify rate for product in OdeSyst2
        if (R[8].i==2)                                                    // if decay_mode == 2, i.e. there are two products
            OdeSyst2_SetRate( R[12].i, R[13].i, R[14].i, newrate);        // modify rate for second product in OdeSyst2
    }
    else                                                                    // if reaction is a one component reaction (a decay reaction) (in that case type of second component is -1)
    {
        OdeSyst1_SetRate( R[1].i, R[2].i, R[3].i, -newrate);     // Modify rate for component in OdeSyst1
        if (R[8].i!=0)                                                // if (decay) mode == 1, i.e. decay into a single decay product
            OdeSyst1_SetRate( R[9].i, R[10].i, R[11].i, newrate);        // modify rate for decay product in OdeSyst1
        if (R[8].i==2)                                                    // if decay_mode == 2, i.e. decay into two components of the PPComp
            OdeSyst1_SetRate( R[12].i, R[13].i, R[14].i, newrate);        // modify rate for second decay product in OdeSyst1
    }
}


void Genome::Add_PPReact()           // at the moment, function only creates true reactions but no decays, i.e. only reactions with two components. That should be extended also to decay reactions when time permits
{
    short index = PPReact.empty.back();              // determine the index of the new PPReact where PPReact.tab is empty right now
    array<intfloat,15>& R = PPReact.tab[index];      // create reference on the corresponding row in PPReact.tab to make the following more readable
    short d;                                           // used as random number later
    float rate;
    rate = R[0].f = UnifRand() * InitialRateMax;                              // randomly choose the reaction rate in [0,1)
    
        SetComponent(R[1].i,R[2].i,array<short, 3>{P,PPC,Phosph},{PPR,index});           // set type and index of the first and second component by using the function SetComponent (that also adds the identity of the PP Reaction to the builders vector of the chosen component)
        SetComponent(R[4].i,R[5].i,array<short, 3>{P,PPC,Phosph},{PPR,index});
        R[7] = number_protein_components(R[1].i,R[2].i) + number_protein_components(R[4].i,R[5].i);     // determin the number of involved proteins for the reaction as the sum of protein components for both components
  
        // determin mode of reaction and product(s)

        /*        if (R[7].i>=2)                                  // if the number of protein components is larger of equal 2, i.e. there is at least one complex involved, choose mode as 0,1, or 2. Otherwise mode = 2 is excluded because it would be without effect.
            R[8] = ChooseInd(3);
        else
            R[8] = ChooseInd(2);      */
        
        R[8] = ChooseInd(3);                                // choose mode as either 0,1 or 2. Thanks to Phosphorylation Switch, also the reaction of two Proteins to themselves can be nontrivial and therefore mode = 2 need not be excluded (as out-commented above)
        
        if (R[8].i!=0) {                                 // if mode != 0, i.e. there is at least one product needed for the reaction, create random number d between 1 and 2*N_1-1 (N_1 being the number of protein components in the first reactant). for d==1 choose the first reactant as product, otherwise randomly pick some subunit of the first reactant by calling ChooseDecayProduct. Finally, perform a PhosphorylationSwitch for the chosen reaction product.
            d = ChooseInd( 2 * number_protein_components(R[1].i,R[2].i) - 1) + 1;
            if (d==1) {
                R[9] = R[1].i;
                R[10] = R[2].i;
            }
            else
                ChooseDecayProduct(R[1].i,R[2].i,d,R[9].i,R[10].i);
            PhosphorylationSwitch(R[9].i,R[10].i,{PPR,index,9});
        }
        
        if (R[8].i == 2)   {                            // if mode ==2, i.e. a second product is needed, perform same as above now for the second reactant
            d = ChooseInd( 2 * number_protein_components(R[4].i,R[5].i) - 1) + 1;
            if (d==1) {
                R[12] = R[4].i;
                R[13] = R[5].i;
            }
            else
                ChooseDecayProduct(R[4].i,R[5].i,d, R[12].i,R[13].i);
            PhosphorylationSwitch(R[12].i,R[13].i,{PPR,index,12});
        }
        // register the PP Reaction in OdeSyst and fill the references in PPReact.tab thereon
        short OdeIndex_R1 = OdeIndexSummands[R[1].i] + R[2].i;      // determine the OdeIndex for the first Rectant with help of the array OdeIndexSummands
        short OdeIndex_R2 = OdeIndexSummands[R[4].i] + R[5].i;      // and also for the second Reactant
        
        OdeSyst2_add({-rate,OdeIndex_R1,OdeIndex_R2,&(R[3].i)}, OdeIndex_R1);    // update OdeSyst2 for Reactant 1
        OdeSyst2_add({-rate,OdeIndex_R1,OdeIndex_R2,&(R[6].i)}, OdeIndex_R2);    // update OdeSyst2 for Reactant 2
        
        if (R[8].i!=0)                                       // if mode is 1 or 2, i.e. there is at least one product
            OdeSyst2_add({rate,OdeIndex_R1,OdeIndex_R2,&(R[11].i)}, OdeIndex(R[9].i, R[10].i) );    // update OdeSyst2 for Product 1
        if (R[8].i==2)                                                       // if mode == 2, i.e. there is a second product
            OdeSyst2_add({rate,OdeIndex_R1,OdeIndex_R2,&(R[14].i)}, OdeIndex(R[12].i, R[13].i) );   // update OdeSyst2 for Product 2
    
    PPReact.cons0.push_back(index);                  // register the new PP Reaction in cons0, cons01, all and delete its index from empty
    PPReact.cons01.push_back(index);
    PPReact.all.push_back(index);
    PPReact.empty.pop_back();
}


void Genome::Delete_PPReact(short index)    // deletes PP Reaction with index 'index'; if the function is called without an index (default value -1) it will randomly pick a PP Reaction
{
    if (index<0)                                        // 0. if no index has been submitted randomly choose index of a PPReact that shall be deleted
        index = PPReact.cons0[ChooseInd(PPReact.cons0.size())];
    array<intfloat,15>& R = PPReact.tab[index];
    
    // 1./2. delete coresponding entries from OdeSyst1 and OdeSyst2 (one to each Ode-Column-Index in the PPReact.tab-row) and for Reactant 1 and 2 also destroy the Complex from the .builders vectors. Use the function OdeSyst_delete (two overloaded methods)
    OdeSyst2_delete(R[1].i,R[2].i,R[3].i,PPR,index);          // 1./2. clear for the reaction the dynamic contributions to reactant 1 and 2 (and earase the Reaction from their builders)
    OdeSyst2_delete(R[4].i,R[5].i,R[6].i,PPR,index);
    
    if (R[8].i != 0)                                         // 1./2. earase contribution to dynamics of product(s) (pay attention on the respective mode)
        OdeSyst2_delete(R[9].i,R[10].i,R[11].i);
    if (R[8].i == 2)
        OdeSyst2_delete(R[12].i,R[13].i,R[14].i);
    
    R[0] = -1.0;                                               // 4. deletion of .tab entry of the Phosphorylate is not necessary, instead only the zeroth entry is set to -1.0 to indicate that this element has been deleted (the function Destroy_Builders checks for negative zeroth element and destroys the element only if it is >= 0 indicating that the element has not already been deleted, similarly for OdeSyst_Delete )
    

    VecDel_element(PPReact.cons0,index);                 // 5. modify the status vectors: remove the index of the destroyed PPReact from cons0, cons01 and all and add it to empty
    VecDel_element(PPReact.cons01,index);
    VecDel_element(PPReact.all,index);
    PPReact.empty.push_back(index);
}


void Genome::Phosphorylates_MutateRate()
{
    short p = ChooseInd(3);                                                 // decide whether to change phosphorylation rate (p=0), dephosphorylation rate (p=1) or decay rate (p=2)
    short index = Phosphorylates.cons01[ChooseInd(Phosphorylates.cons01.size())];           // choose index of Phosphorylate whose rate will be changed
    array<intfloat,15>& R = Phosphorylates.tab[index];                              // create reference on the corresponding row in Phosphorylates.tab to make the following more readable
    float newrate;
    
    if (p==0)                                                               // decision for phosphorylation rate
    {
        newrate = R[0].f = 2*UnifRand() * R[0].f;       // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst1_SetRate( R[1].i, R[2].i, R[3].i, -newrate);     // Modify rate for basis element in OdeSyst1
        OdeSyst1_SetRate( Phosph, index, 0, newrate);       // Modify rate for Phosphorylate in OdeSyst1
    }
    else if (p==1)                                                          // decision for dephosphorylation rate
    {
        newrate = R[4].f = 2*UnifRand() * R[4].f;       // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst1_SetRate( R[1].i, R[2].i, R[5].i, newrate);     // Modify rate for basis element
        OdeSyst1_SetRate( Phosph, index, 2, -newrate);       // Modify rate for Phosphorylate
    }
    else                                                                    // decision for decay rate
    {
        newrate = R[6].f = 2*UnifRand() * R[6].f;         // change rate by multiplying it with a unif. random number in [0,2]
        OdeSyst1_SetRate( Phosph, index, 1, -newrate);                 // modify rate for decay process for the Phosphorylate in OdeSyst1
        if (R[7].i!=0)                                                // if decay_mode != 0, i.e. there is at least one decay product
            OdeSyst1_SetRate( R[8].i, R[9].i, R[10].i, newrate);        // modify rate for first decay product in OdeSyst1
        if (R[7].i==2)                                                    // if decay_mode == 2, i.e. there are two decay products
            OdeSyst1_SetRate( R[11].i, R[12].i, R[13].i, newrate);        // modify rate for second decay product in OdeSyst1
    }
}


void Genome::Add_Phosphorylate()
{
    short index = Phosphorylates.empty.back();              // determine the index of the new Phosphorylate where Phosphorylates.tab is empty right now
    short num_components;
    array<intfloat,15>& R = Phosphorylates.tab[index];      // create reference on the corresponding row in Phosphorylates.tab to make the following more readable
    
    R[0] = UnifRand() * InitialRateMax;                              // randomly choose the phosphorylation rate in [0,1)
    R[4] = UnifRand() * InitialRateMax;                             // randomly choose the dephosphorylation rate in [0,1)
    R[6] = UnifRand() * InitialRateMax;                               // randomly choose the decay rate in [0,1)
    
    SetComponent(R[1].i,R[2].i,array<short, 3>{P,PPC,Phosph},{Phosph,index});           // set type and index of the element that is phosphorylated by using the function SetComponent (that also adds the identity of the Phosphorylate to the builders vector of the chosen element)
    Instances_num_phos(R[1].i,R[2].i)++;        // increase the number of phosphorylates in the basis element
   
    // determin decay mode and product(s)
    num_components = number_protein_components(R[1].i,R[2].i);
    if (num_components==1)                          // if the Phosphorylate consists only of one component, the decay mode is always 0, for larger Phosphorylates it is chosen as either 0,1 or 2
        R[7] = 0;
    else
        R[7] = ChooseInd(3);
    
    if (R[7].i==1)   {                                   // if decay mode == 1 use ChooseDecayProduct to choose one protein or subcomplex of the PP Complex as its decay product
        ChooseDecayProduct(R[1].i,R[2].i,ChooseInd(2 * num_components - 2) + 2, R[8].i,R[9].i);    // call ChooseDecayProduct for the complex delivering a random number d between 2 and 2N-1 (N being the nmber of proteins of the basis element and hence also the Phosphorylate)
        PhosphorylationSwitch(R[8].i,R[9].i,{Phosph,index,8});     // randomly switch to some other phosphorylate (configuration) of the chosen decay product if existent
    }
    else if (R[7].i==2) {                             // if decay mode == 2 the two decay products are the two components of the related unphosphorylated entity of the Phosphorylate. To find these, first reduce the Phosphorylate to its unphosphorylated element (save this in the reference R2) and then extract its respective components.
        short type = Phosph;
        short index2 = index;
        while (type==Phosph) {                                              // if the complex is a phosphorylate, find the corresponding basis element (change type and index to the type and index of the basis element that in that case can only be PPComp)
            type = Phosphorylates.tab[index2][1].i;
            index2 = Phosphorylates.tab[index2][2].i;
        }
        if (type!=PPC) { cerr << "ERROR in Add_Phosphorylate" << '\n';  }
        array<intfloat,15>& R2 = PPComp.tab[index2];                             // create reference on the .tab vector of the PPComp that the Phosphorylate has been rduced to (note that it can only be PPComp and no Protein!)
        
        R[8] = R2[1].i;
        R[9] = R2[2].i;
        R[11] = R2[4].i;
        R[12] = R2[5].i;
    }
    
    R[14] = 0;                                      // set number of Phosphorylates to 0
    
    // register the new Phosphorylate in OdeSyst
    short OdeIndex_Basis = OdeIndexSummands[R[1].i] + R[2].i;      // determine the OdeIndex for the basis element with help of the array OdeIndexSummands
    short OdeIndex_Phosph = OdeIndexSummands[Phosph] + index;      // and same for the Phosphorylate
    // Phosphorylation (rate R[0])
    OdeSyst1_add({-R[0].f,OdeIndex_Basis,&(R[3].i)}, OdeIndex_Basis);    // update OdeSyst1 for basis element for Phosphorylation
    OdeSyst1_add({R[0].f,OdeIndex_Basis}, OdeIndex_Phosph);           // update OdeSyst1 for Phosphorylate for Phosphorylation
    // Decay (rate R[6])
    OdeSyst1_add({-R[6].f,OdeIndex_Phosph}, OdeIndex_Phosph);             // update OdeSyst1 for Phosphorylate and Decay
    if (R[7].i!=0)                                                      // if mode is 1 or 2, i.e. there is at least one decay product
        OdeSyst1_add({R[6].f,OdeIndex_Phosph,&(R[10].i)}, OdeIndex(R[8].i, R[9].i) );    // update OdeSyst1 for Product 1 and Decay
    if (R[7].i==2)                                                       // if mode == 2, i.e. there is a second product
        OdeSyst1_add({R[6].f,OdeIndex_Phosph,&(R[13].i)}, OdeIndex(R[11].i, R[12].i) );   // update OdeSyst1 for Product 2 and Decay
    // Dephosphorylation (rate R[4])
    OdeSyst1_add({R[4].f,OdeIndex_Phosph,&(R[5].i)}, OdeIndex_Basis);    // update OdeSyst1 for basis element and Dephosphorylation
    OdeSyst1_add({-R[4].f,OdeIndex_Phosph}, OdeIndex_Phosph);    // update OdeSyst1 for Phosphorylate and Dephosphorylation
    
    Phosphorylates.cons0.push_back(index);                  // register the new Phosphorylate in cons0, cons01, all and delete its index from empty
    Phosphorylates.cons01.push_back(index);
    Phosphorylates.all.push_back(index);
    Phosphorylates.empty.pop_back();
}


void Genome::Delete_Phosphorylate(short index)    // deletes Phosphorylate with index 'index'; if the function is called without an index (default value -1) it will randomly pick an index
{
    if (index<0)                                        // 0. if no index has been submitted randomly choose index of a Phosphorylate that shall be deleted
        index = Phosphorylates.cons0[ChooseInd(Phosphorylates.cons0.size())];
    array<intfloat,15>& R = Phosphorylates.tab[index];
    // 1./2. delete coresponding entries from OdeSyst1 (one to each Ode-Column-Index in the Phosphorylate.tab-row) and for the basis element also destroy the Phosphorylate from the .builders vectors. Use the function OdeSyst_delete (two overloaded methods)
    OdeSyst1_delete(R[1].i,R[2].i,R[3].i,Phosph,index);          // 1./2. clear for the phosphorylation process the dynamic contributions to the basis element (and earase the Phosphorylate from its builders)
    OdeSyst1_delete(R[1].i,R[2].i,R[5].i);                       // clear for dephosphorylation and the basis element
    
    if (R[7].i != 0)                                         // 1./2. earase contribution to dynamics of decay product(s) due to decay of the Phosphorylate (pay attention on the respective decay mode)
        OdeSyst1_delete(R[8].i,R[9].i,R[10].i);
    if (R[7].i == 2)
        OdeSyst1_delete(R[11].i,R[12].i,R[13].i);
    
    OdeSyst1[OdeIndex(Phosph,index)].clear();                // 3. clear vector of the concerning Phosphorylate from OdeSyst1 and from OdeSyst2
    OdeSyst2[OdeIndex(Phosph,index)].clear();
    
    for (size_t i=0; i<PhosphorylationSwitches[index].size(); ++i)  {  // 3.5. For all entries in PhosphorylationSwitches check if they still point to a reference to this phosphorylate (to exclude the possibility that the respective complex that has originally undergone the posphorylation switch in its decay product has been destroyed in the meantime) and if so set the corresponding ode_pos to -1 to indicate that this decay product is no longer active and must not be deleted from OdeSyst when the Complex is destroyed later on.
        intfloat* Pt = TabAddress_to_Pt(PhosphorylationSwitches[index][i]);
        if (Pt[0].i==Phosph && Pt[1].i==index)
            Pt[2].i = -1;
    }
    PhosphorylationSwitches[index].clear();         // 3.5. finally clear the corresponding PosphorylationSwitches vector.
    
    Instances_num_phos(R[1].i,R[2].i)--;        // 3.5. decrease the number of phosphorylates in the basis element

    R[0] = -1.0;                                               // 4. deletion of .tab entry of the Phosphorylate is not necessary, instead only the zeroth entry is set to -1.0 to indicate that this element has been deleted (the function Destroy_Builders checks for negative zeroth element and destroys the element only if it is >= 0 indicating that the element has not already been deleted, similarly for OdeSyst_Delete )
    
    VecDel_element(Phosphorylates.cons0,index);                 // 5. modify the status vectors: remove the index of the destroyed Phosphorylate from cons0, cons01 and all and add it to empty
    VecDel_element(Phosphorylates.cons01,index);
    VecDel_element(Phosphorylates.all,index);
    Phosphorylates.empty.push_back(index);
    
    DestroyBuilders(Phosphorylates.builders[index]);      // 6. call the DestroyBuilders function that calls the deletion commands for all builders on the Phosphorylate stored in the vector Phosphorylates.builders
    
    Phosphorylates.builders[index].clear();              // 7. clear the builders vector of the deleted Phosphorylate

/*
 // no longer valid. Was used before direct tracking of Phosphorylation Switch via PhosphorylationSwitches was introduced
 
    while (OdeSyst1[OdeIndex(Phosph,index)].size() > 3) {         // before clearing concerning Phosphorylate from OdeSyst1 and from OdeSyst2 check if in OdeSyst1 more than 3 (phosphorylation, dephosphorylation and decay) and in OdeSyst2 more than 0 elements have survived direct cancellation in the deletion cascade that can therefore only be associated with a decay or reaction whoes pruducts have undergone phosphorylation switch from a minor phosphorylate or basis element to this one to be deleted here, so that the respective complex or reaction is not cleared in the deletion cascade and still continues to decay in this posphorylate. We set the corresponding Ode_position in the complex or reaction to -1 to indicate that the dynamics of the porduct does no longer exist (when the complex or reaction is deleted later on, it does not delete the dynamic contribution of the phosphorylate)
        *OdeSyst1[OdeIndex(Phosph,index)].back().add = -1;
        OdeSyst1[OdeIndex(Phosph,index)].pop_back();
    }
    while (OdeSyst2[OdeIndex(Phosph,index)].size() > 0) {
        *OdeSyst2[OdeIndex(Phosph,index)].back().add = -1;
        OdeSyst2[OdeIndex(Phosph,index)].pop_back();
    }
 
    OdeSyst1[OdeIndex(Phosph,index)].clear();                // 3. finally clear vector of the concerning Phosphorylate from OdeSyst1 and from OdeSyst2
    OdeSyst2[OdeIndex(Phosph,index)].clear();  
 */
}



template<short N, short N_Col>
void GeneticInstance<N,N_Col>::Print_Instance(ofstream& file, vector<short> float_entries)    // Prints the content of the tab array to a file. Input: the file hat the function shall write to as an ofstream object and a vector that contains the row indices who's entries are float type (instead of short in the union intfloat)
{
    float_entries.push_back(N_Col);            // add N_Col to the float_entries vector as a 'stop-codon': j never reaches N_Col in the loop below so that once the last real float_entry has passed k becomes trapped at the last entry N_Col and is not incremented any more as the if-condition below is always negative henceforth.
    
    for (short i=0; i<N; ++i)  {
        short k=0;
        if (tab[i][0].f<0)                    // check if tab array is empty (the 0th element is set to -1 if it is empty)
            file << -1 << '\n';
        else  {
            for (short j=0; j<N_Col; ++j)
                if (j==float_entries[k])  {
                    file << left << setprecision(7) << setw(16) << tab[i][j].f;// << "  \t";
                    ++k;
                }
                else
                    file << left << setw(4) << tab[i][j].i;// << '\t';
            file << '\n';
        }
    }
    file << '\n';
    
    for (short i=0; i<N; ++i)   {
        if (builders[i].size()==0)                    // check if builders vector is empty (builders.size()==0)
            file << -1 << '\n';
        else  {
            for (short j=0; j<builders[i].size(); ++j)
                file << left << setw(4) << builders[i][j][0] << setw(6) << builders[i][j][1];
            file << '\n';
        }
    }
}


void Genome::Print_PhosphorylationSwitches(ofstream& file)  {
    for (short i=0; i<N_Phosph; ++i)   {
        if (PhosphorylationSwitches[i].size()==0)                    // check if builders vector is empty (builders.size()==0)
            file << -1 << '\n';
        else  {
            for (short j=0; j<PhosphorylationSwitches[i].size(); ++j)
                file << left << setw(4) << PhosphorylationSwitches[i][j].type << setw(4) << PhosphorylationSwitches[i][j].index << setw(6) << PhosphorylationSwitches[i][j].tab_pos;
            file << '\n';
        }
    }
}


void Genome::PrintGenome(ofstream& file)   // takes a file reference and and writes all relevnat data of Genome into it (all tab arrays as well as the rates and indices of OdeSyst1 and OdeSyst2)
{
    /*
     ofstream file (name_of_file);
     */
    
    file << "max number elements\n";
    file << setw(4) << N_P << setw(4) << N_GPC << setw(4) << N_PPC << setw(4) << N_Phosph << setw(4) << N_PPR << "\n\n";
    
    file << "Essential Proteins (Indices)\n";
    for (short i=0; i<Essential_Proteins.size(); ++i)
        file << setw(4) << N_P + Essential_Proteins[i];
    file << "\n\n";

    // print OdeSyst1 (as 1d vertical vectors)
    file << "OdeSyst1:\n";
    
    short sum=0;
    for (short i=0; i<N_tot; ++i)
        sum += OdeSyst1[i].size();
    
    file << sum << '\n';
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst1[i].size(); ++j)
            file << left << setw(5) << i << setprecision(7) << setw(18) << OdeSyst1[i][j].rate << setw(5) << OdeSyst1[i][j].ind << '\n';
    file << '\n';
    
    // print OdeSyst2 (as 1d vertical vectors)
    file << "OdeSyst2:\n";
    sum=0;
    for (short i=0; i<N_tot; ++i)
        sum += OdeSyst2[i].size();
    
    file << sum << '\n';
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setw(5) << i << setprecision(7) << setw(18) << OdeSyst2[i][j].rate << setw(5) << OdeSyst2[i][j].ind1 << OdeSyst2[i][j].ind2 << '\n';
    file << '\n';
    
    // print OdeSyst1 (as 1d horizontal vectors)
/*        vector<short> OdeSyst1_rec;
    OdeSyst1_rec.reserve((OdeSyst1_MaxSize));
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst1[i].size(); ++j)
            OdeSyst1_rec.push_back(i);
    
    file << "OdeSyst1:\n";
    
    file << OdeSyst1_rec.size() << '\n';
    
    for (short i=0; i<OdeSyst1_rec.size(); ++i)
        file << left << setw(5) << OdeSyst1_rec[i];
    file << '\n';
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst1[i].size(); ++j)
            file << left << setprecision(7) << setw(13) << OdeSyst1[i][j].rate;
    file << '\n';
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst1[i].size(); ++j)
            file << left << setw(5) << OdeSyst1[i][j].ind;
    file << "\n\n";
    
    // print OdeSyst2 (as 1d-vectors)
    vector<short> OdeSyst2_rec;
    OdeSyst2_rec.reserve((OdeSyst2_MaxSize));
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            OdeSyst2_rec.push_back(i);
    
    file << "OdeSyst2:\n";
    
    file << OdeSyst2_rec.size() << '\n';
    
    for (short i=0; i<OdeSyst2_rec.size(); ++i)
        file << left << setw(5) << OdeSyst2_rec[i];
    file << '\n';
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setprecision(7) << setw(13) << OdeSyst2[i][j].rate;
    file << '\n';
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setw(5) << OdeSyst2[i][j].ind1;
    file << '\n';
    
    for (short i=0; i<N_tot; ++i)
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setw(5) << OdeSyst2[i][j].ind2;
    file << "\n\n";
*/
 
//alternative representation of OdeSyst as matrices
/*
    file << "OdeSyst1: Rates\n";
    for (short i=0; i<N_tot; ++i)  {
        file << i << '\t';
        for (short j=0; j<OdeSyst1[i].size(); ++j)
            file << left << setprecision(7) << setw(13) << OdeSyst1[i][j].rate;
        file << '\n';
    }
    file << '\n';
    
    file << "OdeSyst1: Index\n";
    for (short i=0; i<N_tot; ++i)  {
        file << i << '\t';
        for (short j=0; j<OdeSyst1[i].size(); ++j)
            file << left << setw(5) << OdeSyst1[i][j].ind;
        file << '\n';
    }
    file << '\n';
    
    file << "OdeSyst2: Rates\n";
    for (short i=0; i<N_tot; ++i)  {
        file << i << '\t';
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setprecision(7) << setw(13) << OdeSyst2[i][j].rate;
        file << '\n';
    }
    file << '\n';
    
    file << "OdeSyst2: Index1\n";
    for (short i=0; i<N_tot; ++i)  {
        file << i << '\t';
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setw(5) << OdeSyst2[i][j].ind1;
        file << '\n';
    }
    file << '\n';
    
    file << "OdeSyst2: Index2\n";
    for (short i=0; i<N_tot; ++i)  {
        file << i << '\t';
        for (short j=0; j<OdeSyst2[i].size(); ++j)
            file << left << setw(5) << OdeSyst2[i][j].ind2;
        file << '\n';
    }
    file << '\n';
*/
    
    file << "Proteins\n";
    Proteins.Print_Instance(file, {0,1});
    file << '\n';
    
    file << "GPComp\n";
    GPComp.Print_Instance(file, {0,7,10});
    file << '\n';
    
    file << "PPComp\n";
    PPComp.Print_Instance(file, {0,8});
    file << '\n';
    
    file << "Phosphorylates\n";
    Phosphorylates.Print_Instance(file, {0,4,6});
    file << '\n';
    
    file << "PhosphorylationSwitches\n";
    Print_PhosphorylationSwitches(file);
    file << '\n';
    
    file << "PPReact\n";
    PPReact.Print_Instance(file, {0});
    file << '\n';
}


template<short N, short N_Col>
void GeneticInstance<N,N_Col>::Rebuild_Instance(ifstream& file, short N_read, vector<short> float_entries)    // Prints the content of the tab array to a file. Input: the file hat the function shall write to as an ofstream object and a vector that contains the row indices who's entries are float type (instead of short in the union intfloat)
{
    float_entries.push_back(N_Col);            // add N_Col to the float_entries vector as a 'stop-codon': j never reaches N_Col in the loop below so that once the last real float_entry has passed k becomes trapped at the last entry N_Col and is not incremented any more as the if-condition below is always negative henceforth.
    empty.clear();
    cons0.clear();
    cons01.clear();
    all.clear();
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    for (short i=0; i<N_read; ++i)  {
        short k=1;
        file >> tab[i][0].f;
        if (tab[i][0].f < 0)                    // check if tab array is empty (the 0th element is set to -1 if it is empty)
            empty.push_back(i);
        else  {
            cons0.push_back(i);
            cons01.push_back(i);
            all.push_back(i);
            for (short j=1; j<N_Col; ++j)
                if (j==float_entries[k])  {
                    file >> tab[i][j].f;
                    ++k;
                }
                else
                    file >> tab[i][j].i;
        }
    }
    for (short i=N_read; i<N; ++i)    // if N > N_read register the additional sites as empty sites in empty
        empty.push_back(i);
    
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    
    for (short i=0; i<N_read; ++i)   {
        short type;
        short index;
        string line;
        getline(file, line);
        istringstream iss(line);
        while ( iss >> type >> index)
            builders[i].push_back({type,index});
    }
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
}

 
void Genome::Rebuild_PhosphorylationSwitches(ifstream& file, short N_Phosph_read)  {
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    for (short i=0; i<N_Phosph_read; ++i)   {
        short type;
        short index;
        short ode_pos;
        string line;
        getline(file, line);
        istringstream iss(line);
        while ( iss >> type >> index >> ode_pos)
            PhosphorylationSwitches[i].push_back({type,index,ode_pos});
    }
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
}


void Genome::Rebuild_Genome(ifstream& file)   // takes a file reference from which it reads all relevant information to rebuild the individual stored in the file.
{
     short N_P_read, N_GPC_read, N_PPC_read, N_Phosph_read, N_PPR_read;
     
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     
     file >> N_P_read >> N_GPC_read >> N_PPC_read >> N_Phosph_read >> N_PPR_read;
     if (N_P_read > N_P || N_GPC_read > N_GPC || N_PPC_read > N_PPC || N_Phosph_read > N_Phosph || N_PPR_read > N_PPR)  {
         cerr << "ERROR: Max Number of elements provided not large enough to accomodate all elements\n";
         exit(11);
         //cout << "Warning: Max Number of elements provided not large enough to accomodate all elements\n";
     }
     
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     
     // read essential Proteins
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     
    Essential_Proteins.clear();       // clear existing vector in case the indices of the essential Proteins are different and reset it
    short index;
    string line;
    getline(file, line);
    istringstream iss(line);
    while (iss >> index)
        Essential_Proteins.push_back(index - N_P_read);
    if (Essential_Proteins.size() != Number_Essential_Proteins)       // check if the constant of the number of essential Proteins is identical to the number of read indices.
        cerr << "ERROR: Numbers of Essential Proteins do not coincide\n";
    
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    
     // read OdeSyst1 and OdeSyst2
    short OdeSyst1_size, OdeSyst2_size;
    short index_rec, index1, index2;
    float rate;
    
    array<short,5> OdeIndexSummands_read = {0, N_P_read, static_cast<short>(2*N_P_read), static_cast<short>(2*N_P_read+N_GPC_read), static_cast<short>(2*N_P_read+N_GPC_read+N_PPC_read)};
    auto OdeSystAdaptInd = [this, OdeIndexSummands_read](short index_read)  {
        for (short i=4; i>=0; --i)
            if (OdeIndexSummands_read[i] <= index_read)
                return static_cast<short>(index_read - OdeIndexSummands_read[i] + OdeIndexSummands[i]);
        cerr << "Error in OdeSystAdapt\n";
        exit(4);
    };
    
    
     file >> OdeSyst1_size;
     for (short i=0; i!=OdeSyst1_size; ++i) {
         file >> index_rec >> rate >> index;
         OdeSyst1[OdeSystAdaptInd(index_rec)].push_back({rate, OdeSystAdaptInd(index)});
     }
     
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     
     file >> OdeSyst2_size;
     for (short i=0; i!=OdeSyst2_size; ++i) {
         file >> index_rec >> rate >> index1 >> index2;
         OdeSyst2[OdeSystAdaptInd(index_rec)].push_back({rate, OdeSystAdaptInd(index1), OdeSystAdaptInd(index2)});
     }
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
     
     // read tab arrays
     Proteins.Rebuild_Instance(file, N_P_read, {0,1});
     GPComp.Rebuild_Instance(file, N_GPC_read, {0,7,10});
     PPComp.Rebuild_Instance(file, N_PPC_read, {0,8});
     Phosphorylates.Rebuild_Instance(file, N_Phosph_read, {0,4,6});
     Rebuild_PhosphorylationSwitches(file, N_Phosph_read);
     PPReact.Rebuild_Instance(file, N_PPR_read, {0});
     
     // finally recreate the addresses in OdeSyst correctly
     OdeSyst_ShiftAddresses();
}
     
     

inline void Genome::OdeSyst1_add(ode1 add_on, short OdeIndex)      // most important interface to OdeSyst1: Add one unit or one term to the dynamics of OdeSyst1 (specified by one ode1 struct 'add_on'). Add it to the row OdeIndex. Function adds this new ode1 struct at the end of the vector OdeSyst1[OdeIndex], and if an address pointer has been specified in add_on, then it writes the position (in the OdeSyst vector) into the variable that the pointer points to.
{
    OdeSyst1[OdeIndex].push_back(add_on);
    if (add_on.add != NULL)
        *add_on.add = OdeSyst1[OdeIndex].size() - 1;
}


inline void Genome::OdeSyst2_add(ode2 add_on, short OdeIndex)     // same as above, just for OdeSyst2
{
    OdeSyst2[OdeIndex].push_back(add_on);
    if (add_on.add != NULL)
        *add_on.add = OdeSyst2[OdeIndex].size() - 1;
}


inline void Genome::OdeSyst1_SetRate(const short& type, const short& ind, const short& ode_pos, float rate)   // Function that takes as input the type, index and ode_index of a arbitrary genetic instance (all as constant references) as well as a rate (float) that it sets at the corresponding position in OdeSyst1
{
    if (ode_pos>=0)             // check if the entry should still exists (ode_pos is set to -1 in a decaying complex when a decay product that has undergone phosphorylation switch is destroyed)
        OdeSyst1[OdeIndexSummands[type] + ind][ode_pos].rate = rate;      // use the vector OdeIndexSummands to get the corresponding summand (OdeSyst-offset) for the type. together with OdeInd access the rate.
}


inline void Genome::OdeSyst2_SetRate(const short& type, const short& ind, const short& ode_pos, float rate)   // analog to OdeSyst1_SetRate just for OdeSyst2
{
    if (ode_pos>=0)
        OdeSyst2[OdeIndexSummands[type] + ind][ode_pos].rate = rate;
}


inline void Genome::OdeSyst1_SetAddress(const short& type, const short& ind, short& ode_pos)   // Function that takes as input the type, index and ode_index of a arbitrary genetic instance and changes the address of the corresponding entry in OdeSyst1 to the address of ode_pos. This is needed in the function 'OdeSyst_ShiftAddresses' that must be called after each copying and assignment ('=') of Genome.
{
    if (ode_pos>=0)             // check if the entry should still exists (ode_pos is set to -1 in a decaying complex when a decay product that has undergone phosphorylation switch is destroyed)
        OdeSyst1[OdeIndexSummands[type] + ind][ode_pos].add = &ode_pos;      // use the vector OdeIndexSummands to get the corresponding summand (OdeSyst-offset) for the type. together with OdeInd access the address.
}


inline void Genome::OdeSyst2_SetAddress(const short& type, const short& ind, short& ode_pos)   // analog to OdeSyst1_SetAddress just for OdeSyst2
{
    if (ode_pos>=0)
        OdeSyst2[OdeIndexSummands[type] + ind][ode_pos].add = &ode_pos;
}


template<size_t T>
void Genome::OdeSyst_SetComplex(array<intfloat,T>& R, short complex_type, short complex_index)      // template function that takes a reference to a .tab-array row (either to GPComp.tab or PPComp.tab) as well as the type and index of the complex. Function adds contributions to dynamics (to OdeSyst) due to formation of complex, gene expression (in case of GPComp) and decay of complex and enters the references to OdeSyst in the designated fields in the respective .tab array (uses the functions OdeSyst1_add and OdeSyst2_add).
{   // determine indices
    short OdeIndex_C1 = OdeIndexSummands[R[1].i] + R[2].i;      // determine the OdeIndex for the first component with help of the array OdeIndexSummands
    short OdeIndex_C2 = OdeIndexSummands[R[4].i] + R[5].i;      // determine the OdeIndex for the first component with help of the array OdeIndexSummands
    short OdeIndex_Complex = OdeIndexSummands[complex_type] + complex_index;      // determine the OdeIndex for the product
    // set up complex formation
    OdeSyst2_add({-R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[3].i)}, OdeIndex_C1);    // update OdeSyst2 for Component 1
    OdeSyst2_add({-R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[6].i)}, OdeIndex_C2);    // update OdeSyst2 for Component 2
    OdeSyst2_add({R[0].f,OdeIndex_C1,OdeIndex_C2}, OdeIndex_Complex);   // update OdeSyst2 for Complex (assuming that the corresponding vector in OdeSyst2 is degenerate (has 0 elements) right now (the corresponding complex is assumed not to exist before and is being introduced via the function that calls OdeSyst_SetComplex) so that the entry for production appears in column 0 after push_back)
    
    // for GP Complexes
    if (complex_type == GPC) {
        short OdeIndex_exprProtein = OdeIndexSummands[P] + R[8].i;
        // set up gene expression
        OdeSyst1_add({R[7].f,OdeIndex_Complex,&(R[9].i)}, OdeIndex_exprProtein);  // update OdeSyst1 for expression of the protein
        // set up decay process
        OdeSyst1_add({R[10].f,OdeIndex_Complex,&(R[11].i)}, OdeIndex_C1);    // update OdeSyst1, Component1 for decay
        OdeSyst1_add({R[10].f,OdeIndex_Complex,&(R[12].i)}, OdeIndex_C2);    // update OdeSyst1, Component2 for decay
        OdeSyst1_add({-R[10].f,OdeIndex_Complex}, OdeIndex_Complex);         // update OdeSyst1, Complex for decay
    }
    // for PP Complexes: set up decay
    else if (complex_type == PPC) {
        
        OdeSyst1_add({-R[8].f,OdeIndex_Complex},OdeIndex_Complex);   // update OdeSyst1, Complex for decay
        
        if (R[9].i!=0)                                       // if decay mode is 1 or 2, i.e. there is at least one decay product
            OdeSyst1_add({R[8].f,OdeIndex_Complex,&(R[12].i)}, OdeIndex(R[10].i, R[11].i));  // update OdeSyst1, Component1 for decay (if necessary)
        if (R[9].i==2)                                                       // if decay mode == 2, i.e. there are two decay products that are identical to the two components
            OdeSyst1_add({R[8].f,OdeIndex_Complex,&(R[13].i)}, OdeIndex_C2);  // update OdeSyst1, Component2 for decay (if necessary)
    }
}


inline void Genome::OdeSyst1_delete(const short& type, const short& index, const short& ode_pos)    // function that deletes the dynamic contribution to element (type,index) at ode-position ode_pos from OdeSyst1 (after checking that the corresponding vector has not already been cleared (as it might be if the function is called during a deletion cascade triggered by a first random deletion of some other element) and that ode_pos has non-negative value (it can actually be negative (-1) for Phosphorylates as decay products, namely if a phosphorylated version of an element is chosen as a decay product via PhosphorylationSwitch and then this Phosphorylate gets destroyed while the unphosphorylated version of it stays alive, the reference to the decay product (from the PPComp or PPReact or Phosph) must be set to NULL which is indicated by setting the corresponding ode-position to -1; the check for non-negative ode-pos is therefor only important for deletions of PPComps and PPReacts and Phosphorylates where Phosphorylation Switch can happen))
{
    short ode_index = OdeIndexSummands[type] + index;
    if (Exists(type,index) && ode_pos>=0) {
        OdeSyst1[ode_index][ode_pos] = OdeSyst1[ode_index].back();          // exchange to-be-deleted-element with last element
        *(OdeSyst1[ode_index][ode_pos].add) = ode_pos;                      // for the formerly last element change the reference from .tab to it using the saved address
        OdeSyst1[ode_index].pop_back();                                     // destroy last element via pop_back
    }
}


inline void Genome::OdeSyst1_delete(const short& type, const short& index, const short& ode_pos, const short& builder_type, const short& builder_index)    // same as above for OdeSyst1 that also deletes the respective builder (builder_type, builder_index) from the associated builders vector of (type, index) using VecDel_element // could actually also be combined with the previous function.
{
    short ode_index = OdeIndexSummands[type] + index;
    if (Exists(type,index) && ode_pos>=0) {
        OdeSyst1[ode_index][ode_pos] = OdeSyst1[ode_index].back();
        *(OdeSyst1[ode_index][ode_pos].add) = ode_pos;
        OdeSyst1[ode_index].pop_back();
        VecDel_element(Instances_builders(type,index), {builder_type, builder_index});      // from the associated buidlers vector of (type, index) delete the respective builder (builder_type, builder_index) using VecDel_element
    }
}


inline void Genome::OdeSyst2_delete(const short& type, const short& index, const short& ode_pos)    // same as above just for OdeSyst2
{
    short ode_index = OdeIndexSummands[type] + index;
    if (Exists(type,index) && ode_pos>=0) {
        OdeSyst2[ode_index][ode_pos] = OdeSyst2[ode_index].back();
        *(OdeSyst2[ode_index][ode_pos].add) = ode_pos;
        OdeSyst2[ode_index].pop_back();
    }
}


inline void Genome::OdeSyst2_delete(const short& type, const short& index, const short& ode_pos, const short& builder_type, const short& builder_index)    // same as above just for OdeSyst2
{
    short ode_index = OdeIndexSummands[type] + index;
    if (Exists(type,index) && ode_pos>=0) {
        OdeSyst2[ode_index][ode_pos] = OdeSyst2[ode_index].back();
        *(OdeSyst2[ode_index][ode_pos].add) = ode_pos;
        OdeSyst2[ode_index].pop_back();
        VecDel_element(Instances_builders(type,index), {builder_type, builder_index});      // from the associated buidlers vector of (type, index) delete the respective builder (builder_type, builder_index) using VecDel_element
    }
}


void Genome::OdeSyst_ShiftAddresses()      // Function that sets the Addresses all components of OdeSyst1 and 2 to the tab-arrays correctly. This is necessary after copying or assigning the class Genome because all the pointers stored in OdeSyst would otherwise still point to the original class that is copied and not to the tab arrays of the new object. OdeSyst_ShiftAddresses could be used in a copy constructor of Genome but we rather use the default copy constructor and after each assignment or copy call OdeSyst_ShiftAddresses for the newly created object.
{
    short i;
    // Correct Addresses that point to GPComp
    {                       // block the statements so that the reference tab can be defined as the different tab arrays and the compiler does not complain about redefinition for other instantiations
        array< array<intfloat, 13>, N_GPC>& tab = GPComp.tab;     // create reference to the corresponding tab array
        for (i=0; i<N_GPC; ++i) {
            if (tab[i][0].f<0) continue;                            // check if the first element is set to -1 indicating that the element does not exist. In this case stop and continue with the next iteration.
            OdeSyst2_SetAddress(tab[i][1].i, tab[i][2].i, tab[i][3].i);        // production
            OdeSyst2_SetAddress(tab[i][4].i, tab[i][5].i, tab[i][6].i);
            OdeSyst1_SetAddress(P, tab[i][8].i, tab[i][9].i);                   // gene expression
            OdeSyst1_SetAddress(tab[i][1].i, tab[i][2].i, tab[i][11].i);        // decay
            OdeSyst1_SetAddress(tab[i][4].i, tab[i][5].i, tab[i][12].i);
        }
    }
    // Correct Addresses that point to PPComp
    {
        array< array<intfloat, 15>, N_PPC>& tab = PPComp.tab;
        for (i=0; i<N_PPC; ++i) {
            if (tab[i][0].f<0) continue;
            OdeSyst2_SetAddress(tab[i][1].i, tab[i][2].i, tab[i][3].i);         // production
            OdeSyst2_SetAddress(tab[i][4].i, tab[i][5].i, tab[i][6].i);
            if (tab[i][9].i!=0)                                                 // decay
                OdeSyst1_SetAddress(tab[i][10].i, tab[i][11].i, tab[i][12].i);
            if (tab[i][9].i==2)
                OdeSyst1_SetAddress(tab[i][4].i, tab[i][5].i, tab[i][13].i);
        }
    }
    // Correct Addresses that point to PPReact
    {
        array< array<intfloat, 15>, N_PPR>& tab = PPReact.tab;
        for (i=0; i<N_PPR; ++i) {
            if (tab[i][0].f<0) continue;
            OdeSyst2_SetAddress(tab[i][1].i, tab[i][2].i, tab[i][3].i);         // production
            OdeSyst2_SetAddress(tab[i][4].i, tab[i][5].i, tab[i][6].i);
            if (tab[i][8].i!=0)                                                 // decay
                OdeSyst2_SetAddress(tab[i][9].i, tab[i][10].i, tab[i][11].i);
            if (tab[i][8].i==2)
                OdeSyst2_SetAddress(tab[i][12].i, tab[i][13].i, tab[i][14].i);
        }
    }
    // Correct Addresses that point to Phosph
    {
        array< array<intfloat, 15>, N_Phosph>& tab = Phosphorylates.tab;
        for (i=0; i<N_Phosph; ++i) {
            if (tab[i][0].f<0) continue;
            OdeSyst1_SetAddress(tab[i][1].i, tab[i][2].i, tab[i][3].i);         // production
            OdeSyst1_SetAddress(tab[i][1].i, tab[i][2].i, tab[i][5].i);
            if (tab[i][7].i!=0)                                                  // decay
                OdeSyst1_SetAddress(tab[i][8].i, tab[i][9].i, tab[i][10].i);
            if (tab[i][7].i==2)
                OdeSyst1_SetAddress(tab[i][11].i, tab[i][12].i, tab[i][13].i);
        }
    }
    
}













/* old constructor (in the former implementation of the program) not to be used any more, just for reference


void Genome::Genome(int initial_prot_number) {
    for (int i=0; i<initial_prot_number; i++) {
        Proteins.push_back({i,unif_rand(),unif_rand(),i});
    }
    signature = initial_prot_number;
    Numbers[0] = Numbers[1] = initial_prot_number;
    
    genes1.reserve(5);
    proteins1.reserve(5);
    gpcomp1.reserve(5);
    ppcomp1.reserve(5);
    phosphorylates1.reserve(5);
    
    genes2.reserve(5);
    proteins2.reserve(5);
    gpcomp2.reserve(5);
    ppcomp2.reserve(5);
    phosphorylates2.reserve(5);
    
    genes1.resize(initial_prot_number);
    proteins1.resize(initial_prot_number);
    genes2.resize(initial_prot_number);
    proteins2.resize(initial_prot_number);
    for (int i=0; i<initial_prot_number; i++)   {
        genes1[i].reserve(5);
        proteins1[i].reserve(5);
        genes2[i].reserve(5);
        proteins2[i].reserve(5);
    }
    
    list<array<intfloat, 4>>::iterator it = Proteins.begin();
    for (int i=0; i<initial_prot_number; i++)   {
        proteins1[i][0] = {-i, -(*it)[2].f, P, i};
        proteins1[i][1] = {i, (*it)[1].f, G, i};
        it++;
    }
}

*/








/* alternative implementations to the above methods (not tested for compilability yet!!!)

template<short T>
void Genome::OdeSyst_SetComplex_without_using_OdeSyst_add(array<intfloat,T>& R, short complex_type, short complex_index)      // template function that takes a reference to a .tab-array row (either to GPComp.tab or PPComp.tab) as well as the type and index of the complex. Function adds contributions to dynamics (to OdeSyst) due to formation of complex, gene expression (in case of GPComp) and decay of complex and enters the references to OdeSyst in the designated fields in the respective .tab array.
{   // determine indices
    short OdeIndex_C1 = OdeIndexSummands[R[1].i] + R[2].i;      // determine the OdeIndex for the first component with help of the array OdeIndexSummands
    short OdeIndex_C2 = OdeIndexSummands[R[4].i] + R[5].i;      // determine the OdeIndex for the first component with help of the array OdeIndexSummands
    short OdeIndex_Complex = OdeIndexSummands[complex_type] + complex_index;      // determine the OdeIndex for the product
    
    // update OdeSyst2 for Component 1
    R[3] = OdeSyst2[OdeIndex_C1].size();                                 // fill the OdeSyst2 column index of the reaction contribution to the dynamics of the first component in GPComp.tab[3]
    OdeSyst2[OdeIndex_C1].push_back({-R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[3].i)});  // then fill OdeSyst2 for component 1 with the appropriate data
    
    // update OdeSyst2 for Component 2
    R[6] = OdeSyst2[OdeIndex_C2].size();                        // the same for the second component dynamics
    OdeSyst2[OdeIndex_C2].push_back({-R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[6].i)});
    
    // update OdeSyst2 for Complex
    OdeSyst2[OdeIndex_Complex].push_back({R[0].f,OdeIndex_C1,OdeIndex_C2});   // add the production reaction to the dynamics of the product in OdeSyst2 (assuming that the corresponding vector in OdeSyst2 is degenerate (has 0 elements) right now (the corresponding complex is assumed not to exist before and is being introduced via the function that calls OdeSyst2_SetReaction) so that the entry for production appears in column 0 after push_back)
    
    // for GP Complexes
    if (complex_type == GPC) {
        // update OdeSyst1 for expression of the protein
        short OdeIndex_exprProtein = OdeIndexSummands[P] + R[8].i;
        R[9] = OdeSyst1[OdeIndex_exprProtein].size();
        OdeSyst1[OdeIndex_Protein].push_back({R[7].f,OdeIndex_Complex,&(R[9].i)});
        
        // update OdeSyst1, Component1 for decay
        R[11] = OdeSyst1[OdeIndex_C1].size();                        // fill the OdeSyst1 column index of the decay contribution to the dynamics of the first decay product (=first component) in GPComp.tab[3]
        OdeSyst1[OdeIndex_C1].push_back({R[10].f,OdeIndex_Complex,&(R[11].i)});  // then fill OdeSyst1 for component 1 with the appropriate data
        
        // update OdeSyst1, Component2 for decay
        R[12] = OdeSyst1[OdeIndex_C2].size();                        // the same for the second decay product dynamics
        OdeSyst1[OdeIndex_C2].push_back({R[10].f,OdeIndex_Complex,&(R[12].i)});
        
        // update OdeSyst1, Complex for decay
        OdeSyst1[OdeIndex_Complex].push_back({-R[10].f,OdeIndex_Complex});
    }
    // for PP Complexes
    else if (complex_type == PPC) {         // in case of a PP Complex
        
        // update OdeSyst1, Complex for decay
        OdeSyst1[OdeIndex_Complex].push_back({-R[8].f,OdeIndex_Complex});   // add decay of complex to dynamics in OdeSyst1
        
        if (R[9].i!=0) {                                      // if decay mode is 1 or 2, i.e. there is at least one decay product
            // update OdeSyst1, Component1 for decay (if necessary)
            short OdeIndex_DecayProduct1 = OdeIndexSummands[R[10].i] + R[11].i;      // determine the OdeIndex for the first decay product with help of the array OdeIndexSummands
            R[12] = OdeSyst1[OdeIndex_DecayProduct1].size();                        // fill the OdeSyst1 column index of the decay contribution to the dynamics of the first decay product
            OdeSyst1[OdeIndex_DecayProduct1].push_back({R[8].f,OdeIndex_Complex,&(R[12].i)});  // then fill OdeSyst1 for decay product 1 with the appropriate data
        }
        if (R[9]==2) {                                                      // if decay mode == 2, i.e. there are two decay products that are identical to the two components
            // update OdeSyst1, Component2 for decay (if necessary)
            R[13] = OdeSyst1[OdeIndex_C2].size();                        // fill the OdeSyst1 column index of the decay contribution to the dynamics of the first decay product
            OdeSyst1[OdeIndex_C2].push_back({R[8].f,OdeIndex_Complex,&(R[13].i)});  // then fill OdeSyst1 for decay product 1 with the appropriate data
        }
    }
}


void Genome::OdeSyst_SetReaction_without_using_OdeSyst_add(array<intfloat,15>& R)      // function similar to OdeSyst2_SetCFReaction for Pure Reactions (the first part is indeed identical to OdeSyst2_SetCFReaction). Takes a reference to a PPReact.tab-array row from which it takes all relevant information to add the contributions to dynamics of component 1 and 2 and of the product(s) to OdeSyst2. For component 1 and 2 and the product(s) it fills the OdeSyst2 column Index into the corresponding positions in the PPReact.tab-array.
{
    short OdeIndex_C1 = OdeIndexSummands[R[1].i] + R[2].i;      // determine the OdeIndex for the first component with help of the array OdeIndexSummands
    short OdeIndex_C2 = OdeIndexSummands[R[4].i] + R[5].i;      // determine the OdeIndex for the first component with help of the array OdeIndexSummands
    
    R[3] = OdeSyst2[OdeIndex_C1].size();                        // fill the OdeSyst2 column index of the reaction contribution to the dynamics of the first component in GPComp.tab[3]
    OdeSyst2[OdeIndex_C1].push_back({-R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[3].i)});  // then fill OdeSyst2 for component 1 with the appropriate data
    
    R[6] = OdeSyst2[OdeIndex_C2].size();                        // the same for the second component dynamics
    OdeSyst2[OdeIndex_C2].push_back({-R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[6].i)});
    
    if (R[8].i!=0)                                       // if decay mode is 1 or 2, i.e. there is at least one product
    {
        short OdeIndex_P1 = OdeIndexSummands[R[9].i] + R[10].i;      // determine the OdeIndex for the first product with help of the array OdeIndexSummands
        R[11] = OdeSyst2[OdeIndex_P1].size();                        // fill the OdeSyst2 column index of the reaction contribution to the dynamics of the first product in PPReact.tab[11]
        OdeSyst2[OdeIndex_P1].push_back({R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[11].i)});  // then fill OdeSyst2 for product 1 with the appropriate data
    }
    if (R[8]==2)                                                       // if there is a second decay product
    {
        short OdeIndex_P2 = OdeIndexSummands[R[12].i] + R[13].i;      // determine the OdeIndex for the first product with help of the array OdeIndexSummands
        R[14] = OdeSyst2[OdeIndex_P2].size();                        // fill the OdeSyst2 column index of the reaction contribution to the dynamics of the first product in PPReact.tab[11]
        OdeSyst2[OdeIndex_P2].push_back({R[0].f,OdeIndex_C1,OdeIndex_C2,&(R[14].i)});  // then fill OdeSyst2 for product 1 with the appropriate data
    }
}



void Genome::ChooseDecayProduct_alternative(short Complex_type, short Complex_index, short d, short& DecProd_type, short& DecProd_index)
{
    if (d==1) {
        DecProd_type = Complex_type;
        DecProd_index = Complex_index;
    }
    else {
        array<intfloat,15> R = Instances[Complex_type]->tab[Complex_index];
        short size_component1 = (R[1]==PPC) ? Instances[R[1]]->tab[R[2]][7].i : 0;
        if (d  <= size_component1+2)
            ChooseDecayProduct(short R[1], short R[2], short d-1, short& DecProd_type, short& DecProd_index);
        else
            ChooseDecayProduct(short R[4], short R[5], short d-(size_component1+2), short& DecProd_type, short& DecProd_index);
    }
}


void Genome::PhosphorylationSwitch_alternative(short& type, short& index)
{
    short number_phosphorylates = Instances[type]->tab[index].back().i;
    vector<array<short,2>>& B = Instances[type]->builders[index];
    short i = 0;
    short d;
    
    if (type != Phosph) {
        d = ChooseInd(number_phosphorylates+1);
        while (d > 0)  {
            if(B[i][0] == Phosph)  {
                type = B[i][0];
                index = B[i][1];
                d--;
            }
            i++;
        }
        return;
    }
    
    if (type == Phosph) {
        d = ChooseInd(number_phosphorylates+2);
        if (d==0) {
            type = Phosphorylates.tab[index][!!!].i;
            index = Phosphorylates.tab[index][!!!].i;
            return;
        }
        while (d > 1)  {
            if(B[i][0] == Phosph)  {
                type = B[i][0];
                index = B[i][1];
                d--;
            }
            i++;
        }
        return;
    }
}

*/


