#include "inc/common.h"
#include "inc/json.hpp"
#include "inc/xoshiro.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>

using namespace std;
using json = nlohmann::json;

// structure to hold simulation results
struct simulationStruct
{
    vector<int64_t> N;      // susceptible population
    vector<int64_t> M;      // mutant population
    int step;               // time step
};
simulationStruct sim_function(
    double lambda, double mu, double lambda1, double mu1, double beta,
    long long int Nf, int threshold, int max_steps)
{
    // Simulation setup --------------------------------------

    vector<int64_t> N(max_steps, 0);
    vector<int64_t> M(max_steps, 0);
    N[0] = 10LL;  // initial susceptible population
    M[0] = 0LL;   // initial mutant population
    int step = 0;

    // Initialise simulation --------------------------------------

    for (int i = 0; i < max_steps - 1; ++i)
    {
        int64_t new_N = 0; // susceptible divisions - deaths - mutants
        int64_t new_M = 0; // mutant divisions - deaths
        int64_t mutations = 0; // new mutants from susceptible population

        // stop simulation if population exceeds threshold or goes extinct
        if (N[i] >= Nf || M[i] >= Nf || N[i] <= 0 || M[i] < 0)
        {
            break;
        }
        
        // Susceptible population --------------------------------------
        
        // if population is small, explicit simulation
        if (N[i] < threshold)
        {
            for (int j = 0; j < N[i]; ++j)
            {
                
                double rand_div = get_rand_uniform();
                if (rand_div < lambda) // parent cell divides
                {
                    
                    double rand_mut = get_rand_uniform();
                    if (rand_mut < beta)
                    {
                        mutations++;  // offspring is mutant
                    }
                    else
                    {
                        new_N++;  // offspring is susceptible
                    }
                }
                
                double rand_death = get_rand_uniform();
                if (rand_death < mu)
                {
                    new_N--;  // parent cell dies
                }
            }
        }
        else // if population is large, use binom approximation
        {
            // sample number of division events
            int64_t divisions = get_rand_binom(N[i], lambda);
            
            // sample number of death events
            int64_t deaths_N = get_rand_binom(N[i], mu);

            // sample number of division events that produce mutants
            mutations = get_rand_binom(divisions, beta);
            
            // calculate net change
            new_N = (divisions - mutations) - deaths_N;
        }
        
        // Mutant population --------------------------------------
        
        // if population is small, explicit simulation
        if (M[i] < threshold && M[i] > 0)
        {
            for (int j = 0; j < M[i]; ++j)
            {
                double rand_div = get_rand_uniform();
                if (rand_div < lambda1)
                {
                    new_M++;  // cell divides
                }
                
                double rand_death = get_rand_uniform();
                if (rand_death < mu1)
                {
                    new_M--;  // cell dies
                }
            }
        }
        else if (M[i] > 0) // if population is large, use binom approximation
        {
            // sample number of division events
            int64_t divisions_M = get_rand_binom(M[i], lambda1);
            
            // sample number of death events
            int64_t deaths_M = get_rand_binom(M[i], mu1);
            
            // calculate net change
            new_M = divisions_M - deaths_M;
        }

        // update populations. Must be >=0
        N[i + 1] = max(N[i] + new_N, 0LL);
        M[i + 1] = max(M[i] + new_M + mutations, 0LL);
        
        step++;
    }
    
    return {N, M, step};
}


int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    
    // load parameters from json
    std::ifstream f("parameters.json");
    json input_parameters = json::parse(f);

    std::cout << "Input parameters: " << std::endl;
    for (const auto &[key, value] : input_parameters.items())
    {
        std::cout << key << ": " << value << std::endl;
    }

    // extract parameters
    double lambda = input_parameters["lambda"];      // susceptible division probability
    double mu = input_parameters["mu"];              // susceptible death probability
    double lambda1 = input_parameters["lambda1"];    // mutant division probability
    double mu1 = input_parameters["mu1"];            // mutant death probability
    double beta = input_parameters["beta"];          // mutation probability per division
    long long int Nf = input_parameters["Nf"];       // max population threshold
    int threshold = input_parameters["threshold"];   // simulation threshold for binom approximation
    int max_steps = input_parameters["max_steps"];   // maximum time steps
    int sims = input_parameters["sims"];             // number of simulations
    int test = input_parameters["test"];             // test mode flag (0=batch, 1=single)

    int complete_sim = 0;  // counter for valid simulations
   
    if (test == 1) // run single simulation for debugging
    {
        simulationStruct out = sim_function(lambda, mu, lambda1, mu1, beta, Nf, threshold, max_steps);
        
        cout << "Step, N, M, mutant_proportion" << endl;
        for (int i = 0; i <= out.step; ++i)
        {
            double total = out.N[i] + out.M[i];
            double prop = (total > 0) ? static_cast<double>(out.M[i]) / total : 0.0;
            
            cout << i << ", "
                 << out.N[i] << ", "
                 << out.M[i] << ", "
                 << round(prop * 1e8) / 1e8 << ", " <<endl;
        }

        cout << endl << "Total steps completed: " << out.step << endl;
    }
    else if (test == 0) // run batch simulations and write results to csv
    {
        ofstream file("out.csv");
        for (int i = 0; i < sims; ++i)
        {
            simulationStruct out = sim_function(lambda, mu, lambda1, mu1, beta, Nf, threshold, max_steps);

            // remove failed simulations (timestep < 10, population < 10 or >1e12)
            if (out.step < 10 || out.N[out.step] < 10 || out.N[out.step] > 1e12)
            {
                // silently skip
                continue;
            }
            
            // extract last 10 timesteps to capture point of exceeding threshold
            for (int j = out.step - 9; j <= out.step; j++)
            {
                file << i << "," << j << "," << out.N[j] << "," << out.M[j] << "\n";
            }
            complete_sim++;
        }
        file.close();
    }

    cout << endl << "Completed non-extinct simulations: " << complete_sim << " / " << sims << endl;
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    
    return 0;
}