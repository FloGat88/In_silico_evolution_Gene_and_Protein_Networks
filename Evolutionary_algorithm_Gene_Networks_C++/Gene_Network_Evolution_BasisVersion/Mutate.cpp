//  Mutate.cpp

#include "RandomNumberGenerator.h"
#include "Genome.h"


void Genome::MutateRate()
{
    short d = ChooseInd_respective_weights(array< size_t, 5> ({Proteins.cons01.size()*2, GPComp.cons01.size()*3, PPComp.cons01.size()*2, Phosphorylates.cons01.size()*3, PPReact.cons01.size()}) );     // weights according to the number of rates in the different element types
    switch (d) {
        case 0:
            Proteins_MutateRate();
            break;
        case 1:
            GPComp_MutateRate();
            break;
        case 2:
            PPComp_MutateRate();
            break;
        case 3:
            Phosphorylates_MutateRate();
            break;
        case 4:
            PPReact_MutateRate();
            break;
        default:
            cerr << "ERROR in MutateRate" << '\n';
    }
    
}


void Genome::Add_Element()
{
    short d = ChooseInd_respective_weights(array< size_t, 5> ({Proteins.empty.size(), GPComp.empty.size(), PPComp.empty.size(), Phosphorylates.empty.size(), PPReact.empty.size() }) );         // equal weights
    
    switch (d) {
        case 0:
            Add_Protein();
            break;
        case 1:
            Add_GPComp();
            break;
        case 2:
            Add_PPComp();
            break;
        case 3:
            Add_Phosphorylate();
            break;
        case 4:
            Add_PPReact();
            break;
        case -1:
            cout << "Cannot add new element so delete one" << '\n';
            Delete_Element();
            break;
        default:
            cerr << "ERROR in Add_Element" << '\n';
    }
    
}


void Genome::Delete_Element()
{
    short d = ChooseInd_respective_weights(array< size_t, 5> ({Proteins.cons0.size(), GPComp.cons0.size()*3, PPComp.cons0.size()*3, Phosphorylates.cons0.size()*3, PPReact.cons0.size()*3 }) );         // weight for Proteins as the most fundamental type is chosen a factor 1/3 smaller than all other weights to favour the emergence of complex networks
    
    switch (d) {
        case 0:
            Delete_Protein();
            break;
        case 1:
            Delete_GPComp();
            break;
        case 2:
            Delete_PPComp();
            break;
        case 3:
            Delete_Phosphorylate();
            break;
        case 4:
            Delete_PPReact();
            break;
        case -1:
            cout << "cannot delete new element so add one" << '\n';
            Add_Element();
            break;
        default:
            cerr << "ERROR in Delete_Element" << '\n';
    }
    
}


void Genome::Mutate(int times)     // times gives the number of mutations that occur per function call
{
    for (int i = 1; i<=times; ++i)  {
        //cout << i << '\n';
        short d = ChooseInd_respective_weights(array< size_t, 3> ({3,2,1}) );   // favour creation of elements over their deletion and rate mutation over any topological change of the network
        switch (d) {
            case 0:              
                MutateRate();
                break;
            case 1:
                Add_Element();
                break;
            case 2:
                Delete_Element();
                break;
            default:
                cerr << "ERROR in Mutate" << '\n';
                break;
        }
    }
}
