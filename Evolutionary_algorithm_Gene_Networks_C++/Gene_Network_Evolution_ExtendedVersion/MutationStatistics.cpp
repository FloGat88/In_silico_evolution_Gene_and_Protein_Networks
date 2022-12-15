//
//  MutationsAnalyzer.cpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 08/09/16.
//  Copyright © 2016 Florian Gartner. All rights reserved.
//

#include <dirent.h>
#include "EvoAlg.h"



void EvoAlg::MutationStatistics(const vector<short> MutationTypes)
{
    constexpr size_t N_stat = 1000;                  // N_statistics: number of samples for each individual (i.e. for each file in the folder) and each number of mutations
    constexpr size_t N_mutations = 10;               // maximum number of mutations in a single individual up to which the statistics is to be made
    
    auto GetGeneration = [](string file_name)  {                                 // lambda expression that returns the corresponding generation number out ouf a file name
        return stoi(file_name.substr(10,file_name.find_first_of("_")-10));
    };
    auto CheckFile = [](string file_name)    {
        string substr = file_name.substr(0,10);
        return substr.compare("Generation")==0 ? true : false;
    };
    
    array<string, 4> OutputFileNames = {"MutationStatistics.txt","MutateRateStatistics.txt","AddElementStatistics.txt","DeleteElementStatistics.txt"};
    
    
    for(short d=0; d<MutationTypes.size(); ++d)   {
        short mtype=MutationTypes[d];
    
        vector<int> generations;                         // saves the number of generation for each individual in the finder
        vector<float> fitness;                             // saves the corresponding fitness for each individual in the folder
        array< vector< array<short,N_stat> >, N_mutations> mutation_kind;           // saves the kind of mutation that has happened in the individual previous in the 3rd dimension or in the original individual if index of third dim is 0, leading to the corresponding fitness difference that is saved in mutated_fitness_diff
        array< vector< array<float,N_stat> >, N_mutations> mutated_fitness_diff;    // saves in the first dimension the fitness differences in all the sample mutations to a certain individual (2nd dimension) and with a specified number of mutations (3rd dimension) (1... N_mutations)

        short m = 0;                            // index for the individuals in the folder
    
        generations.reserve(100);               // reserve space for vectors for better efficiency
        fitness.reserve(100);
        for (short i=0; i<N_mutations; ++i)    {
            mutation_kind[i].reserve(100);
            mutated_fitness_diff[i].reserve(100);
        }
    
        DIR *hdir;
        struct dirent *entry;
        hdir = opendir("MutationStatistics");      // open the directory with name "MutationStatistics" and make a pointer?? hdir into that directory to subsequently read out all files in it
        short chd = chdir("MutationStatistics");     // and also change into that directory
        if(chd != 0)
            cout << "could not change directory\n";
   
        do {
            entry = readdir(hdir);          // discard all entries that belong to the folder but are not the required text files like for example "." (current directory), ".." (directory one level up), Ds_store??. Check this by looking whether the filename starts with "Generation"
        }   while(CheckFile(entry->d_name)!=true);
        
        while (entry)   {                   // as long as still another file can be found and opend in the directory continue
            cout << entry->d_name << '\n';
            
            Individual individual, mutated_individual;           // individual is the original individual that is read out from the data (files in the directory) and mutated_individual denotes the mutated version thereof
            Rebuild_Individual_from_Data(entry->d_name, individual);           // recreate individual from text file in the directory and save it on individual
            
            generations.push_back(GetGeneration(entry->d_name));                // puch back a new entry line in generations and fitness and a new individual-plane in mutation_kind and mutated_fitness_diff
            fitness.push_back(individual.Fitness);
            
            for (short i=0; i<N_mutations; ++i)    {
                mutation_kind[i].push_back({});
                mutated_fitness_diff[i].push_back({});
            }

            for (short i=0; i<N_stat; ++i)   {                              // loop over entries in first dimension (samples)
                mutated_individual = individual;
                mutated_individual.OdeSyst_ShiftAddresses();              // shift the addresses in OdeSyst to the new object
                for (short t=0; t<N_mutations; ++t)   {                     // loop over third dimension (mutation numbers)
                    mutation_kind[t][m][i] = mutated_individual.Mutate(1,mtype);      // first make one further mutation and save its kind in mutation_kind
                    mutated_fitness_diff[t][m][i] =  FitnessEvaluator(mutated_individual) - fitness[m];    // calculate the fitness difference to the original individual and save it in mutated_fitness_diff
                }
            }
            do {
                entry = readdir(hdir);                  // read out the next file in the directory
            }   while(entry && CheckFile(entry->d_name)!=true);         // read next file if this file is not a real individual data file but the directory has not yet ended.
            ++m;
        }
    
        ofstream file (OutputFileNames[mtype]);         // create a file in the directory "MutationStatistics" with corresponding name specified in OutputFileNames in which all the data will be written
        if (file.is_open()) {
            file << "size\n";
            file << fitness.size() << "   " <<  N_stat << "   " << N_mutations << "\n\n";     // first write the sizes of the three dimensions (number of individuals, samples per individual per mutation number) and maximum mutation number
        
        
            file << "generation/fitness\n";                 // print generation and fitness
            for (short i=0; i<fitness.size(); ++i)
                file << left << setw(8) << generations[i] << setw(12) << fitness[i] << '\n';
            file << '\n';
        
            for (short t=0; t<N_mutations; ++t)  {          // for each mutation number (3rd dimension)
                file << "mutation_kind:  " << t+1 << " Mutations" << '\n';
                for (short i=0; i<mutation_kind[0].size(); ++i)  {       // and each individual (2nd dim)
                    for (short j=0; j<N_stat; ++j)                  // print respective mutation kind for each samples (1st dim)
                        file << left << setw(4) << mutation_kind[t][i][j];
                    file << '\n';
                }
                file << '\n';
        
                file << "mutated_fitness_diff:  " << t+1 << " Mutations" << '\n';
                for (short i=0; i<mutated_fitness_diff[0].size(); ++i)  {
                    for (short j=0; j<N_stat; ++j)
                        file << left << setprecision(6) << setw(15) << mutated_fitness_diff[t][i][j];       // print respective fitness difference for each samples (1st dim)
                    file << '\n';
                }
                file << '\n';
            }
        }
        else
            cout << "Could not open file\n";
        file.close();
    
        closedir(hdir);             // close the directory again
        chd = chdir("..");     // change back direcotory
        if(chd != 0)
            cout << "could not change directory\n";
    }
}



void EvoAlg::CompleteDeleteElementStatistics()
{
    constexpr short N_deletions = 2;        // number of deletions up to which all possible combinations of deletions are tried out. Choose either 1, 2 or 3.
    auto GetGeneration = [](string file_name)  {   // lambda expression that returns the corresponding generation number out ouf a file name
        return stoi(file_name.substr(10,file_name.find_first_of("_")-10));
    };
    auto CheckFile = [](string file_name)    {
        string substr = file_name.substr(0,10);
        return substr.compare("Generation")==0 ? true : false;
    };
    
    
    DIR *hdir;
    struct dirent *entry;
    hdir = opendir("MutationStatistics");      // open the directory with name "MutationStatistics" and make a pointer?? hdir into that directory to subsequently read out the files in it
    short chd = chdir("MutationStatistics");     // and also change into that directory
    if(chd != 0)
        cout << "could not change directory\n";
    
    do {
        entry = readdir(hdir);          // discard all entries that belong to the folder but are not the required text files like for example "." (current directory), ".." (directory one level up), Ds_store??. Check this by looking whether the filename starts with "Generation"
    }   while(CheckFile(entry->d_name)!=true);
    
    if(!entry)   {                   // if no file could be found in the directory
        cout << "Could not open file\n";
        exit(16);
    }
    
    Individual individual, mutant1, mutant2, mutant3;           // individual is the original individual that is read out from the data (file in the directory) and mutant1, mutant2 and mutant3 are the mutated version thereof with 1,2 or 3 mutations, respectively.
    
    Rebuild_Individual_from_Data(entry->d_name, individual);           // recreate individual from text file in the directory and save it on individual
    
    int generation = GetGeneration(entry->d_name);                // get the generation from the fle name
    float fitness_parent = individual.Fitness;
    short N_el = individual.GetNumberElements();
    
    vector<short> mutation_kind_1, mutation_kind_2, mutation_kind_3;   // in these matrices save the kind of element that is deleted (1=Protein, 2=GPComp, 3=PPComp, 4=Phosphorylate, 5=PPReact)
    vector<float> mutated_fitness_diff_1, mutated_fitness_diff_2, mutated_fitness_diff_3; // and here the corresponding fitenss differences with the parent
    mutation_kind_1.reserve(N_el);                 // reserve maximal sizes that the vecotrs can reach. In reality, the vectors for 2 and 3 deletions will be smaller because when one element is deleted also all its builders are deleted, therefore, e.g. in mutant1 often more elements than one will already have been deleted and the fewer combinations of two and three simulataneous deletions exist.
    mutation_kind_2.reserve(N_el*N_el);
    mutation_kind_3.reserve(N_el*N_el*N_el);
    mutated_fitness_diff_1.reserve(N_el);
    mutated_fitness_diff_2.reserve(N_el*N_el);
    mutated_fitness_diff_3.reserve(N_el*N_el*N_el);
    
    if(N_deletions>0)
        for(short t1=1; t1<=5; ++t1)                // go through all types
            for(short i1=0; i1<individual.GetNumberElements(t1); ++i1) {     // and for each type through all elements (i1 is the all-index, not the tab-index)
                mutant1 = individual;
                mutant1.OdeSyst_ShiftAddresses();
                mutant1.Delete_Element(t1,i1);          // delete the corresponding element
                mutation_kind_1.push_back(t1);                // save type of deleted element in mutation_kind_1
                mutated_fitness_diff_1.push_back(FitnessEvaluator(mutant1) - fitness_parent);         // and calculate the fitness difference to the parent and save the value in mutated_fitness_diff_1
                if(N_deletions>1)       // in case N_deletions is at least 2 continue and for each mutant from the first round perform each second possible deletion mutation
                    for(short t2=1; t2<=5; ++t2)
                        for(short i2=0; i2<mutant1.GetNumberElements(t2); ++i2) {
                            mutant2 = mutant1;
                            mutant2.OdeSyst_ShiftAddresses();
                            mutant2.Delete_Element(t2,i2);
                            mutation_kind_2.push_back(t2);
                            mutated_fitness_diff_2.push_back(FitnessEvaluator(mutant2) - fitness_parent);
                            if(N_deletions>2)       // if N_deletions = 3 continue again and perform each possible triple deletion mutation combination
                                for(short t3=1; t3<=5; ++t3)
                                    for(short i3=0; i3<mutant2.GetNumberElements(t3); ++i3) {
                                        mutant3 = mutant2;
                                        mutant3.OdeSyst_ShiftAddresses();
                                        mutant3.Delete_Element(t3,i3);
                                        mutation_kind_3.push_back(t3);
                                        mutated_fitness_diff_3.push_back(FitnessEvaluator(mutant3) - fitness_parent);
                                    }
                        }
            }

    
    ofstream file ("CompleteDeleteElementStatistics.txt");         // create a file in the directory "MutationStatistics" in which all the data will be written
    if (file.is_open()) {
        file << "size\n";
        file << 1 << "   " <<  N_el << "   " << N_deletions << "\n\n";     // first write number of individuals (always 1), number of elements and maximum mutation number (either 1, 2 or 3)
        
        
        file << "generation/fitness\n";                 // print generation and fitness
        file << left << setw(8) << generation << setw(12) << fitness_parent << '\n';
        file << '\n';
        
        auto PrintDeletionStatistics = [N_el, &file](short n_del, vector<short>& mutation_kind, vector<float>& mutated_fitness_diff)  {   // lambda expression to print the kind of the deleted element and resulting fitness difference between mutatnt and parent
            file << "mutation_kind:  " << n_del << endl;
            for(short i=0; i<mutation_kind.size(); ++i)
                file << left << setw(4) << mutation_kind[i];
            file << endl;
            file << "mutated_fitness_diff:   " << n_del << endl;
            for(short i=0; i<mutated_fitness_diff.size(); ++i)
                file << left << setprecision(6) << setw(15) << mutated_fitness_diff[i];
            file << endl << endl;
        };
        
        if(N_deletions>0)                   // call the above lambda function for the required numbers of deletions
            PrintDeletionStatistics(1, mutation_kind_1, mutated_fitness_diff_1);
        if(N_deletions>1)
            PrintDeletionStatistics(2, mutation_kind_2, mutated_fitness_diff_2);
        if(N_deletions>2)
            PrintDeletionStatistics(3, mutation_kind_3, mutated_fitness_diff_3);
        
    }
    else
        cout << "Could not open file\n";
    
    file.close();
    
    closedir(hdir);             // close the directory again
    chd = chdir("..");     // change back direcotory
    if(chd != 0)
        cout << "could not change directory\n";
}


