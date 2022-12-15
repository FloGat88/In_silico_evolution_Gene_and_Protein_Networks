//  EvoAlg.h

#ifndef ____EvoAlg__
#define ____EvoAlg__

#include "Individual.h"
#include "Dynamics_and_Fitness.h"
//#include <boost/filesystem.hpp>

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
    
public:
    Evolve_FitterHalfSelection(T& FitnessEvaluator) : FitnessEvaluator(FitnessEvaluator) {};
    
    void operator () (array<Individual, N>& Individuals, array< FitnessIndex_struct , N>& FitnessValues, float& Fitness_max, float& Fitness_tot)
    {
        nth_element(FitnessValues.begin(), FitnessValues.begin()+N/2, FitnessValues.end(), [](FitnessIndex_struct& struct1, FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});
        
        for (short i=0; i<N/2; ++i) {
            short ind_reproduce = FitnessValues[i].ind;
            short ind_die = FitnessValues[i+N/2].ind;
            float fit_old = Individuals[ind_die].Fitness;
            
            Individuals[ind_die] = Individuals[ind_reproduce];
            Individuals[ind_die].OdeSyst_ShiftAddresses();              // shift the addresses in OdeSyst to the new object
            Individuals[ind_die].Mutate(2);
            
            float fit_new = FitnessValues[i+N/2].fit = FitnessEvaluator(Individuals[ind_die]);
            Fitness_tot += (fit_new - fit_old);
            if (fit_new > Fitness_max)
                Fitness_max = fit_new;
            
       /*     if (Fitness_max < 23)   {         // Addendum: as long as the homogeneity bound is not exceeded mutate also the reproducing individual and hence all individuals in order to more easily escape from that local minimum
                Individuals[ind_reproduce].Mutate(2);
                float fit_new = FitnessValues[i].fit = FitnessEvaluator(Individuals[ind_reproduce]);
                Fitness_tot += (fit_new - fit_old);
                if (fit_new > Fitness_max)
                    Fitness_max = fit_new;
            }   */
                
        }
        
    }
    
};




class EvoAlg {
    private:
        constexpr static short N = 20;                     // number of Individuals
        constexpr static float L1 = 10.0;                  // length of the small system
        constexpr static float L2 = 30.0;                  // length of the large system (relevant only for FitnessEvaluation_2Gaussian_Peaks_scaling)
        constexpr static float D = 0.25;                   // diffusion constant
        constexpr static double T_max = 100.0;     // maximal time for the simulations of the dynamics to run: choose it such that T_max >> time scale for diffusion over the complete system (L^2/D)
        constexpr static short N_discretisation = 20;      // number of pieces with which the system sizes are discretized
        constexpr static int GenerationCountRestart = 5000;
        constexpr static int MaxNumberGenerations = 10000; // for the other constants, Essential_Protein_Number and the maximal numbers of each element, N_P, N_GPC,..., see Genome.h
    
        constexpr static int FitnessTrajectory_size = MaxNumberGenerations-GenerationCountRestart+1;
 
// FitnessEvaluator and Evolver for Complementary Gradients
        FitnessEvaluation_ComplementaryGradients<N_discretisation> FitnessEvaluator = FitnessEvaluation_ComplementaryGradients<N_discretisation>(L1, D, T_max);
        Evolve_FitterHalfSelection< N, FitnessEvaluation_ComplementaryGradients<N_discretisation> > Evolver = Evolve_FitterHalfSelection< N, FitnessEvaluation_ComplementaryGradients<N_discretisation> >(FitnessEvaluator);

    
    void PrintIndividualData(array<Individual, N>& Individuals, array< FitnessIndex_struct , N>& FitnessValues, int generation_count, short rank)
    {
        rank = rank - 1;
        nth_element(FitnessValues.begin(), FitnessValues.begin()+rank, FitnessValues.end(), [](FitnessIndex_struct& struct1, FitnessIndex_struct& struct2) -> bool {return (struct1.fit > struct2.fit);});
        short index = FitnessValues[rank].ind;
        ofstream file ("Generation"+to_string(generation_count)+"_Rank"+to_string(rank+1)+"_Fitness_"+to_string(FitnessValues[rank].fit)+".txt");
        if (file.is_open()) {
            file << "Fitness:\n";
            file << FitnessValues[rank].fit << "  (" << Individuals[index].Scaling_Fitness[0] << ' ' << Individuals[index].Scaling_Fitness[1] << ")\n";
            file << "Constants\n";
            file << N << ' ' << L1 << ' ' << L2 << ' ' << D << ' ' << T_max << ' ' << N_discretisation << ' ' << InitialRateMax << '\n';
            //    file << N_P << ' ' << N_GPC << ' ' << N_PPC << ' ' << N_Phosph << ' ' << N_PPR << '\n';
            Individuals[index].PrintGenome(file);
        }
        else
            cout << "Could not open file\n";
        file.close();
    };

    
    void Rebuild_Individual_from_Data(string file_name, Individual& individual)
    {
        ifstream file (file_name);
        if (file.is_open())  {
            short N_read, N_L_read;
            float L1_read, L2_read, D_read, InitialrateMax_read;
            double T_max_read;
            file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
            file >> individual.Fitness;
            file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
            file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
            file >> N_read >> L1_read >> L2_read >> D_read >> T_max_read >> N_L_read >> InitialrateMax_read;
            file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
            individual.Rebuild_Genome(file);
            if (L1_read != L1 || L2_read != L2 || D_read != D || T_max_read != T_max || N_L_read != N_discretisation)   {
                cout << "Recalculate Fitness because constants do not coincide\n";
                FitnessEvaluator(individual);
            }
        }
        else
            cout << "Could not open file\n";
        file.close();
    };
    
    
    void Rebuild_Population_from_Data(array<Individual, N>& Individuals)
    {
        
        Rebuild_Individual_from_Data("Generation5000_Rank1_Fitness_81.523689.txt" ,Individuals[0]);
        Rebuild_Individual_from_Data("Generation5000_Rank2_Fitness_81.521996.txt" ,Individuals[1]);
        Rebuild_Individual_from_Data("Generation5000_Rank3_Fitness_81.521996.txt" ,Individuals[2]);
        Rebuild_Individual_from_Data("Generation5000_Rank4_Fitness_81.521996.txt" ,Individuals[3]);
        for (short i=4; i<N; ++i)
        {
            Individuals[i] = Individuals[i-4];              // for the other individuals copy the rebuilt ones
            Individuals[i].OdeSyst_ShiftAddresses();        // shift the addresses in OdeSyst to the new objects
            Individuals[i].Mutate(20);                      // Mutation rush for those individuals
            FitnessEvaluator(Individuals[i]);
        }
    };

    
    void Create_Random_Population(array<Individual, N>& Individuals)
    {
        for (short i=0; i<N; ++i) {
            Individuals[i].Diversify();
            FitnessEvaluator(Individuals[i]);
        }
    };

    
public:
    void operator() ()
    {
        array< FitnessTrajectory_struct , (MaxNumberGenerations-GenerationCountRestart)+1> FitnessTrajectory;
        int generation_count = GenerationCountRestart;
start_anew:
        array<Individual, N> Individuals;
        array< FitnessIndex_struct , N> FitnessValues;
        float Fitness_max = -1e9, Fitness_tot = 0;
    
    // create initial Population either randomly or recreate from data and initialize FitnessValues, Fitness_max, Fitness_tot and write first entry into FitnessTrajectory
     //   Create_Random_Population(Individuals);
        Rebuild_Population_from_Data(Individuals);
        for (short i=0; i<N; ++i) {
            FitnessValues[i].ind = i;
            FitnessValues[i].fit = Individuals[i].Fitness;
            Fitness_tot += FitnessValues[i].fit;
            if(FitnessValues[i].fit > Fitness_max)
                Fitness_max = FitnessValues[i].fit;
        }
        FitnessTrajectory[0] = {Fitness_tot/N, Fitness_max};
        
    // let popultion evolve until either the maximum Fitness exceeds a certain value or the maximum generation count is reached
        while (Fitness_max < 95 && generation_count<MaxNumberGenerations) {
            ++generation_count;
            cout << "Calculating " << generation_count << "th generation\n";
            Evolver(Individuals, FitnessValues, Fitness_max, Fitness_tot);
            FitnessTrajectory[generation_count-GenerationCountRestart] = {Fitness_tot/N, Fitness_max};
            
            if (generation_count < 250 && generation_count % 5 == 0)  {
                PrintIndividualData(Individuals, FitnessValues, generation_count, 2);
                PrintIndividualData(Individuals, FitnessValues, generation_count, 1);
            }
            
            else if (generation_count % 50 == 0)  {
                PrintIndividualData(Individuals, FitnessValues, generation_count, 2);
                PrintIndividualData(Individuals, FitnessValues, generation_count, 1);
            }
            /*        if (generation_count % 500 == 0 && Fitness_max < 22.5)
             goto start_anew;   */
            }
            //   const char dir_path[] = "Output_data_files\\TEST";
            //    boost::filesystem::path dir(dir_path);
            //   if(boost::filesystem::create_directory(dir)) {
            //      std::cout << "Success" << "\n";
            //   }
            //mkdir("last_generation");
            for(short i=4; i>=1; --i)        // save fittest 4 individuals in data
                PrintIndividualData(Individuals, FitnessValues, generation_count, i);
            
            // myFile.open( "Fitness_Evolution.txt", ios::out | ios::app );
            // ofstream file ("Fitness_Evolution.txt", ios::app);
            ofstream file ("Fitness_Evolution.txt", ios::app);
            for (short i = GenerationCountRestart > 0 ? 1 : 0; i<=(generation_count - GenerationCountRestart); ++i)
            file << i+GenerationCountRestart << '\t' << setw(10) << FitnessTrajectory[i].mean << '\t' << setw(10) << FitnessTrajectory[i].max << '\n';
            file.close();
    };
};





    
    

    
    
    
        
    
    
    
    
    
    
    
    
    




#endif /* defined(____EvoAlg__) */
