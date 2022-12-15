//  EvoAlg.h

#ifndef ____EvoAlg__
#define ____EvoAlg__

#include "Individual.h"
#include "Fitness_Evaluation.h"
#include "LamSimAnn_Energy.hpp"
#include "RandomNumberGenerator.h"
#include <dirent.h>
#include <thread>                       // for parallelisation
//#include <boost/filesystem.hpp>

/////////////////////////////////////////////// Additional Parameters ///////////////////////////////////////////////////////////////
constexpr short N = 8;             // number of Individuals, population size
constexpr short Print_Fittest_N_Individuals = 1;        // number of Individuals with highest rank (highest fitness) being printed at each check point

constexpr bool Print_Fittest_Lineage = false;     // if true print the complete lineage (only every ...th generation) of the fittest individual in the last generation (the output could then be used for analysis by MutationStatistics for example)
constexpr unsigned short Evolution_History_step_size = 10;    // only relevant if Print_Fittest_Lineage = true. In this case determines the number of generations that are left out in the record of Evolution_History between two successive records. (determines the crudeness of the record of the lineages)

constexpr bool rebuild_population_from_data = false;
constexpr int GenerationCountRestart = 0;
constexpr int MaxNumberGenerations = 1000; // for the other constants, Essential_Protein_Number and the maximal numbers of each element, N_P, N_GPC,..., see Genome.h

constexpr int FitnessTrajectory_size = MaxNumberGenerations-GenerationCountRestart+1;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct FitnessIndex_struct {float fit; short ind;};

struct FitnessTrajectory_struct {
    float mean;
    float max;
    FitnessTrajectory_struct() : mean(0), max(0)   {};
    FitnessTrajectory_struct(float mean,float max) : mean(mean), max(max)   {};
};


template<short N, class T>
class Evolve_FitterHalfSelection
{
private:
    T& FitnessEvaluator;
    const short EvolverID = 1;
    
public:
    Evolve_FitterHalfSelection(T& FitnessEvaluator) : FitnessEvaluator(FitnessEvaluator) {};
    short GetID() { return EvolverID;};
    
    void operator () (array<Individual, N>& Population, array< FitnessIndex_struct , N>& FitnessValues, float& Fitness_max, float& Fitness_tot)
    {
        nth_element(FitnessValues.begin(), FitnessValues.begin()+N/2, FitnessValues.end(), [](const FitnessIndex_struct& struct1, const FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});
#ifdef PARALLELIZE
        vector<thread> threads;
#endif
        for (short i=0; i<N/2; ++i) {
            short ind_reproduce = FitnessValues[i].ind;
            short ind_die = FitnessValues[i+N/2].ind;
            
            Population[ind_die] = Population[ind_reproduce];
            Population[ind_die].Mutate(2);
#ifndef NDEBUG
            Population[ind_die].Test_All();
#endif
#ifndef PARALLELIZE
            FitnessEvaluator(&Population[ind_die]);
#else
            threads.push_back(thread{FitnessEvaluator, &(Population[ind_die])});
#endif
        }
        
#ifdef PARALLELIZE
        for (std::thread& t:threads) {
            if (t.joinable())
                t.join();
        }
#endif
        for (short i=N/2; i<N; ++i) {
            short ind_die = FitnessValues[i].ind;
            float fit_old = FitnessValues[i].fit;
            float fit_new = Population[ind_die].Fitness;
            FitnessValues[i].fit = fit_new;
            Fitness_tot += (fit_new - fit_old);
            if (fit_new > Fitness_max)
                Fitness_max = fit_new;
        }
    };
};


       /*     if (Fitness_max < 23)   {         // Addendum: as long as the homogeneity bound is not exceeded mutate also the reproducing individual and hence all individuals in order to more easily escape from that local minimum
                Population[ind_reproduce].Mutate(2);
                float fit_new = FitnessValues[i].fit = FitnessEvaluator(Population[ind_reproduce]);
                Fitness_tot += (fit_new - fit_old);
                if (fit_new > Fitness_max)
                    Fitness_max = fit_new;
            }   */
    


template<short N, class T>
class Evolve_Stochastic
{
private:
    T& FitnessEvaluator;
    const short EvolverID = 2;
    
public:
    Evolve_Stochastic(T& FitnessEvaluator) : FitnessEvaluator(FitnessEvaluator) {};
    Individual FittestIndividual;                            // as good solutions can always be lost again in this evolution scheme, safe the fittest solution that has been encoutered
    short GetID() { return EvolverID;};
    
    void operator () (array<Individual, N>& Population, array< FitnessIndex_struct , N>& FitnessValues, float& Fitness_max, float& Fitness_tot)
    {
        short weights[N];
        for(short i=0; i<N; ++i)
            weights[i] = N-i;
        
        for (short i=0; i<N/2; ++i) {
            short p = ChooseInd_respective_weights(weights,N);

            nth_element(FitnessValues.begin(), FitnessValues.begin()+p, FitnessValues.end(), [](const FitnessIndex_struct& struct1, const FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});
        
            short ind_reproduce = FitnessValues[p].ind;
            short FitVal_ind_die = ChooseInd(N);
            short ind_die = FitnessValues[FitVal_ind_die].ind;
            float fit_old = Population[ind_die].Fitness;
            
            if(ind_die != ind_reproduce)
                Population[ind_die] = Population[ind_reproduce];
            
            float r=UnifRand();     // random number that determines if and how many mutations happen in the individual
            bool Mutation = false;          // indicates later if indeed a mutation has happened
            if(r > 0.5)  {
                Population[ind_die].Mutate(2);
                Mutation = true;
            }
 //           if(r > 0.75)
  //              Population[ind_die].Mutate(1);
            
            if(Mutation)              // only if a mutation has happend calculate new Fitness of the individual
                FitnessEvaluator(Population[ind_die]);
            
            float fit_new = FitnessValues[FitVal_ind_die].fit = Population[ind_die].Fitness;
                Fitness_tot += (fit_new - fit_old);
                if (fit_new > Fitness_max)  {           // If a new fitness record has been archieved, update Fitness_max and save this Individual as FittestIndividual
                    Fitness_max = fit_new;
                 //   FittestIndividual = Population[ind_die];
                 //   FittestIndividual.OdeSyst_ShiftAddresses();
                }
        }
    };
};


template<short N, class T>
class Evolve_Independent
{
private:
    T& FitnessEvaluator;
    const short EvolverID = 3;
    
public:
    Evolve_Independent(T& FitnessEvaluator) : FitnessEvaluator(FitnessEvaluator) {};
    short GetID() { return EvolverID;};
    
    void operator () (array<Individual, N>& Population, array< FitnessIndex_struct , N>& FitnessValues, float& Fitness_max, float& Fitness_tot)
    {
        for (short i=0; i<N; ++i)    {
            Individual individual_mutated = Population[i];
            individual_mutated.Mutate(1);
            float fit_new = FitnessValues[i].fit = FitnessEvaluator(individual_mutated);    // the indices in FitnessValues remain unchanged in this Evolver, therefore FitnessValues[i] indices simply the ith Individual (see initializztion of FitnessValues) and it is necessary only to change the corresponding Fitness entry in FitnessValues.
            float fit_old = Population[i].Fitness;
            
            if(fit_new >= fit_old)  {
                Population[i] = individual_mutated;
                Fitness_tot += (fit_new - fit_old);
                if (fit_new > Fitness_max)
                    Fitness_max = fit_new;
            }
        }
    }
};





class EvoAlg {
    private:
 
// choose FitnessEvaluator-and Evovler type
//        typedef FitnessEvaluation_DrosophilaGapgenes<N_discretisation>  FitnessEvaluator_Type;
        typedef FitnessEvaluation_Polarisation<N_discretisation>  FitnessEvaluator_Type;
        typedef Evolve_FitterHalfSelection< N, FitnessEvaluator_Type >  Evolver_Type;
    
        FitnessEvaluator_Type  FitnessEvaluator = FitnessEvaluator_Type(L_small);
        Evolver_Type  Evolver = Evolver_Type(FitnessEvaluator);

    
    void PrintIndividualData(array<Individual, N>& Population, array< FitnessIndex_struct , N>& FitnessValues, int generation_count, short rank)
    {
        if(Evolver.GetID() != 3)
            nth_element(FitnessValues.begin(), FitnessValues.begin()+rank-1, FitnessValues.end(), [](const FitnessIndex_struct& struct1, const FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});
        
        short index = FitnessValues[rank-1].ind;
        ofstream file ("Output_data_files/Generation"+to_string(generation_count)+"_Rank"+to_string(rank)+"_Fitness_"+to_string(Population[index].Fitness)+".txt");
        if (file.is_open()) {
            file << "Fitness:\n";
            file << Population[index].Fitness << endl;
            Population[index].PrintGenome(file);
        }
        else
            cout << "Could not open file\n";
        file.close();
    };

    
    void PrintIndividualData(Individual& individual, int generation_count, short rank)     // overloaded method to print directly a single Individual into the folder "MutationStatistics" (not search it in an array of Individulas according to the rank); to be used for example to print the Individuals of the lineage from Common_Ancesters in Record_Ancestry; only a reference to a single Individual needs to be provided and the corresponding generation.
    {
        ofstream file ("MutationStatistics/Generation"+to_string(generation_count)+"_Rank"+to_string(rank)+"_Fitness_"+to_string(individual.Fitness)+".txt");
        if (file.is_open()) {
            file << "Fitness:\n";
            file << individual.Fitness << endl;
            individual.PrintGenome(file);
        }
        else
            cout << "Could not open file\n";
        file.close();
    };
    
    
    class Record_Ancestry {
    private:
        EvoAlg& evoalg;          // save a reference on the containing class EvoAlg in this nested class to be able to use its member functions (especially "PrintIndividualData" in PrintLineage) within this nested class
        deque<array<Individual, N>> Evolution_History;     // stores the recent evolutionary history (at a certain discretisation (Evolution_History_step_size)) until coalescence to a single common ancester happened, the former common ancestry of the lineage is then transferred to Common_Lineage
        short EvoHist_size;             // index of Evolution_History at which ancestry_index_2 in the Individuals has been initialized the last time, to check for coalescence. When Coalescence has occurred the common ancester line up to index EvoHist_size in Evolution_History can be transferred to "Common_Ancesters" (to save storage) and the EvoHist_size is subsequently set to the actual size of Evolution_History;
        vector<Individual> Common_Lineage;   // saves the former common ancestry of the population (resp. the lineage of the fittest individual at the end) at a discretisation of Evolution_History_step_size
        vector<short> generation;    // saves the corresponding generation of the Individuals in 'Common_Lineage'
        vector<short> rank;         // saves the corresponding rank of the Individuals in 'Common_Lineage'
        
    public:
        Record_Ancestry(EvoAlg &evoalg) : evoalg(evoalg), EvoHist_size(0) {
            if(Print_Fittest_Lineage == true)  {           // only if Print_Fittest_Lineage = true and so the class will be really needed reserve space for the vectors
                Common_Lineage.reserve((MaxNumberGenerations-GenerationCountRestart)/Evolution_History_step_size + 1);   // if Print_Fittest_Lineage_N == true reserve enough space in Common_Lineage to save the whole lineage in Compressed form (the Population array is first saved as a whole in Evolution_History and can be transferred to Common_Ancesters as soon as coalescence has occured).
                generation.reserve((MaxNumberGenerations-GenerationCountRestart)/Evolution_History_step_size + 1);
                rank.reserve((MaxNumberGenerations-GenerationCountRestart)/Evolution_History_step_size + 1);
            }
        };
        
        bool Check_Coalescence(const array<Individual, N> &Population)   {      // Checks whether coalescence has apperaed in the population 'Population' (relative to the generation that is saved in 'Evolution_History[EvoHist_size-1]') by checking whether all ancestry_index_2 (that refer to that generation) in the population coincide.
            bool coalescence = true;
            short ancestry_index_2 = Population[0].ancestry_index_2;
            for(short i=1; i<N; ++i)
                if(Population[i].ancestry_index_2 != ancestry_index_2)
                    coalescence = false;
            return coalescence;
        };
    private:
        void Transfer_to_Common_Lineage(short EvoHist_size, short ancestry_index)   {     // Transfer data from Evolutionary_History to the compressed representation of the common Lineage in 'Common_Lineage', which can be done after Coalescence has appeared. Only the first 'EvoHist_size' columns (Population arrays) of Evolution_History will be transferred where 'ancestry_index' gives the index of the common ancester in the latest generation therein (at position ancestry_index-1)
            Common_Lineage.resize(Common_Lineage.size() + EvoHist_size);            // reserve enough space for the Transfer in Common_Lineage and in rank
            rank.resize(rank.size() + EvoHist_size);
            for(short i=EvoHist_size-1, j=Common_Lineage.size()-1; i>=0; --i,--j)   {
                Common_Lineage[j] = Evolution_History[i][ancestry_index];
                rank[j] = GetRank(Evolution_History[i], ancestry_index);      // get the rank of the individual and save it in 'rank'
                ancestry_index = Evolution_History[i][ancestry_index].ancestry_index_1;
            }
             Evolution_History.erase(Evolution_History.begin(), Evolution_History.begin()+EvoHist_size);
        };
        
        short GetRank(array<Individual, N> &Population_array, short index)   {    // get rank of the Individual with index 'index' within Population_array, (i.e. the number of individuals therein with higher fitness plus one). Is used in Transfer_to_Common_Lineage: Population_array is there a reference to a column (one time slice) of Evolution_History and as index the 'ancestry_index' from 'Transfer_to_Common_Lineage' is handed in.
            float test_fitness = Population_array[index].Fitness;
            short number_of_fitter = 0;
            for(short i=0; i<N; ++i)
                if(Population_array[index].Fitness > test_fitness)
                    ++number_of_fitter;
            return number_of_fitter + 1;
        };
        
    public:
        void Record_Generation(array<Individual, N> &Population, short generation_count)  {
            Evolution_History.push_back(Population);      // Copy of Indviduals to Evolution_History. Individuals in EvolutionHistory and Common_Lineage will only be printed to the data but not be Fitnessevaluated any more.
            generation.push_back(generation_count);
            for(short i=0; i<N; ++i)
                Population[i].ancestry_index_1 = i;
            if(Check_Coalescence(Population) == true)   {
                Transfer_to_Common_Lineage(EvoHist_size, Population[0].ancestry_index_2);
                EvoHist_size = Evolution_History.size();
                for(short i=0; i<N; ++i)
                    Population[i].ancestry_index_2 = i;
            }
        };
        
        void Print_Lineage(short index)  {
            Transfer_to_Common_Lineage(Evolution_History.size(), index);    // transfer all remaining elements of "Evolution_History" to "Common_Ancesters" (call the function with the index of the individual of the last generation (usually the fittest one))
            for(short i=0; i<Common_Lineage.size(); ++i)
                evoalg.PrintIndividualData(Common_Lineage[i], generation[i], rank[i]);   // call the routine from the containing EvoAlg class to which evoalg is a reference.
        };
        
    };

public:
/*
    void MutationStatistics (const vector<short> MutationTypes);     // for the implementatin see MutationStatistics.cpp
    void CompleteDeleteElementStatistics();          // for the implementatin see MutationStatistics.cpp. Make the complete MutationStatistics for all combinations of element deletions (up to 3 simulataneous deletions). Only supported for a single individual. As for MutationStatistics, make a folder called MutationStatistics that contains the data file of one individual who's mutation statistics is to be analyzed
*/
    void operator() ()
    {
        array< FitnessTrajectory_struct , (MaxNumberGenerations-GenerationCountRestart)+1> FitnessTrajectory;
        int generation_count = GenerationCountRestart;
    start_anew:
 //       initial_step = 1.0e-1;     // !!!!!!!!! remove again !!!!!!
        array<Individual, N> Population;
        array< FitnessIndex_struct , N> FitnessValues;
        float Fitness_max = -1e9, Fitness_tot = 0;
        Record_Ancestry AncestryRecorder(*this);        // instantiate an instance of Record_Ancestry (that will only come to use if Print_Fittest_Lineage == true) to record and finally print the Lineage that produces the fittest individual of the last generation
        
    // create initial Population either randomly or recreate from data and initialize FitnessValues, Fitness_max, Fitness_tot and write first entry into FitnessTrajectory
        if(rebuild_population_from_data)  {
            // Rebuild_Population_from_Data(Population);
            cout << "Rebuild from data not yet implemented/n";
            exit(6);
        }
        else
            for(short i=0; i<N; ++i)  {         // evaluate Fitness for all Individuals
                Population[i].Create_Random_Initial_Network();        // create random initial networks
                FitnessEvaluator(&Population[i]);                     // evaluate these random networks
            }
        
        for (short i=0; i<N; ++i) {
            FitnessValues[i].ind = i;
            FitnessValues[i].fit = Population[i].Fitness;
            Fitness_tot += FitnessValues[i].fit;
            if(FitnessValues[i].fit > Fitness_max)
                Fitness_max = FitnessValues[i].fit;
            Population[i].ancestry_index_1 = i;
            Population[i].ancestry_index_2 = i;      // initialize the ancester index according to the index in Population; this number will be inherited to an individual's offspring so that it can be seen who its ancestor was if the Population array of the ancester is saved somewhere (it is indeed saved in Evolution_History and Common_Lineage after all Evolution_History_step_size generations, afterwards ancester_index is initialized to the index in Population again )
        }
        FitnessTrajectory[0] = {Fitness_tot/N, Fitness_max};
        
        
    // let popultion evolve until either the maximum Fitness exceeds a certain value or the maximum generation count is reached
        while (Fitness_max < 95 && generation_count<MaxNumberGenerations) {
            ++generation_count;
            cout << "Calculating " << generation_count << "th generation\n";
  //          initial_step = generation_count <= 499  ?  1.0e-1 : 1.0e-4;       // !!!!!!!! remove this again !!!!!!!!
            Evolver(Population, FitnessValues, Fitness_max, Fitness_tot);
#ifndef NDEBUG
            for(short i=0; i<N; ++i) Population[i].Test_All();
#endif
            FitnessTrajectory[generation_count-GenerationCountRestart] = {Fitness_tot/N, Fitness_max};
            
            if ( (generation_count <= 100 && generation_count % 50 == 0) || (generation_count > 100 && generation_count % 50 == 0))  {
                for(short i=Print_Fittest_N_Individuals; i>=1; --i)
                    PrintIndividualData(Population, FitnessValues, generation_count, i);
                cout << "maximum Fitness found: " << Fitness_max << '\n';
            }
            
            if ( Print_Fittest_Lineage && (generation_count == 1 || generation_count % Evolution_History_step_size == 0) )
                AncestryRecorder.Record_Generation(Population, generation_count);
            
            /*        if (generation_count % 500 == 0 && Fitness_max < 22.5)
             goto start_anew;   */
        }
            //   const char dir_path[] = "Output_data_files\\TEST";
            //    boost::filesystem::path dir(dir_path);
            //   if(boost::filesystem::create_directory(dir)) {
            //      std::cout << "Success" << "\n";
            //   }
            //mkdir("last_generation");
        
        for(short i=N; i>=1; --i)        // save fittest 4 individuals in data
            PrintIndividualData(Population, FitnessValues, generation_count, i);
        
        if(Print_Fittest_Lineage)   {
            nth_element(FitnessValues.begin(), FitnessValues.begin(), FitnessValues.end(), [](const FitnessIndex_struct& struct1, const FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});    // find fittest Individual in the last generation (set it to front of 'FitnessValues')
            AncestryRecorder.Print_Lineage(FitnessValues[0].ind);   // print lineage of the fittest individual of the last generation from AncestryRecorder
        }

        // myFile.open( "Fitness_Evolution.txt", ios::out | ios::app );
        ofstream file ("Output_data_files/Fitness_Evolution.txt", ios::app);
        if (file.is_open())  {
            file << "Fitness Evolution" << endl;
            for (short i = GenerationCountRestart > 0 ? 1 : 0; i<=(generation_count - GenerationCountRestart); ++i)
                file << i+GenerationCountRestart << '\t' << setw(10) << FitnessTrajectory[i].mean << '\t' << setw(10) << FitnessTrajectory[i].max << '\n';
            file << endl;
            file << "Mutations" << endl;
            nth_element(FitnessValues.begin(), FitnessValues.begin(), FitnessValues.end(), [](const FitnessIndex_struct& struct1, const FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});       // find fittest individual to print its Mutations vector
            Population[FitnessValues[0].ind].Print_Mutations(file);
        }
        file.close();
    };
};



#endif /* defined(____EvoAlg__) */
