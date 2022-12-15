function Evolution()
    N = 8;          % population size
    generation_count = int16(0);
    Fitness_max = 0;
    
    for i=N:-1:1
        Individuals(i) = Network;
        Individuals(i).FitnessPolarisation_pdepe();
    end
    
    while generation_count<1000
        generation_count = generation_count + 1
        [~,ind_repr] = maxk([Individuals.Fitness],N/2);
        [~,ind_die] = mink([Individuals.Fitness],N/2);
        Individuals(ind_die)=Individuals(ind_repr);
        %% sequential Fitness evaluation
        for i=ind_die
            Individuals(i).Mutate(2);
 %           Individuals(i).FitnessPolarisation_pdepe();
        end
 %       if max([Individuals.Fitness]) > Fitness_max
 %           Fitness_max = max([Individuals.Fitness]);
 %       end
 
 %% parallel fitness evaluation with parfor
        Individuals_EvalFitness = Individuals(ind_die);
        fitness = zeros(1,N/2);
        parfor i=1:N/2
            fitness(i) = Individuals_EvalFitness(i).FitnessPolarisation_pdepe()
        end
        for i=1:N\2
            Individuals(i).Fitness = fitness(i);
            if fitness(i) > Fitness_max
                Fitness_max = fitness(i);
            end
        end
    end
end

