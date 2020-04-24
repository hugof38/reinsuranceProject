#include "pch.h"
#include "framework.h"
#include "Dll2.h"

#include <math.h>
#include <iostream>
#include <random>
#include <chrono> 
//#include <future>

//static std::mutex s_collec;

static std::random_device r;
static std::default_random_engine generator{ r() };
static std::uniform_real_distribution<double> uDistribution(0.0, 1.0);
static std::normal_distribution<double> nDistribution(0.0, 1.0);

static auto lambda_tranche3 = [](double sample) {return std::max(sample - 20.0, 0.0); };

static auto lambda_AAL = [](double agr) {return std::max(agr - 5.0, 0.0); };

double tranche(double agr, double min_val, double max_val)
{
    return std::min(min_val, std::max(agr - max_val, 0.0));
}

DLL2_API  int simul_poiss(const double lambd)
{
    double p = exp(-lambd);
    int k = 0;
    double s = p;
    double u = uDistribution(generator);

    while (u > s)
    {
        ++k;
        p = p * lambd / k;
        s += p;
    }
    return(k);
}

static void simul_collec_loop(double* collec_simul, const double lambd, const double mu, const double sigma)
{
    int cmpt_poiss;
    double u;
    double sample;
    double AAL;

    double agregate_year[4] = {};

    cmpt_poiss = simul_poiss(lambd);

    for (int j = 0; j < cmpt_poiss; ++j)
    {
        u = nDistribution(generator);
        sample = exp(mu + sigma * u);

        agregate_year[0] += sample;
        agregate_year[1] += tranche(sample, 5.0, 5.0);
        agregate_year[2] += tranche(sample, 10.0, 10.0);
        agregate_year[3] += lambda_tranche3(sample);

    }

    //std::lock_guard<std::mutex> lock(s_collec);
    collec_simul[0] += agregate_year[0];

    AAL = lambda_AAL(agregate_year[1]);
    collec_simul[1] += AAL;
    collec_simul[2] += tranche(AAL, 5.0, 0.0);
    collec_simul[3] += tranche(AAL, 5.0, 5.0);
    collec_simul[4] += tranche(AAL, 5.0, 10.0);
    collec_simul[9] += tranche(AAL, 5.0, 15.0);

    collec_simul[5] += tranche(agregate_year[2], 10.0, 0.0);
    collec_simul[10] += tranche(agregate_year[2], 10.0, 10.0);
    collec_simul[6] += agregate_year[3];

    collec_simul[7] += agregate_year[2];

    collec_simul[8] += agregate_year[1];

}

DLL2_API void simul_collective(double* collec_simul, const int nbSimul, const double lambd, const double mu, const double sigma)
{

    for (int i = 0; i < nbSimul; ++i)
    {
        simul_collec_loop(collec_simul, lambd, mu, sigma);

        //std::async(std::launch::async, simul_collec_loop, collec_simul, lambd, mu, sigma);
    }

    for (int i = 0; i < 11; ++i)
    {
        collec_simul[i] = collec_simul[i] / nbSimul;
    }
}

DLL2_API int __stdcall wrap_simul_col(double* result_simul_time, const int nbSimul, const double lambd, const double mu, const double sigma)
{
    double collec_simul[11] = {};

    auto start = std::chrono::high_resolution_clock::now();
    simul_collective(collec_simul, nbSimul, lambd, mu, sigma);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    for (int i = 0; i < 11; ++i)
    {
        result_simul_time[i] = collec_simul[i];
    }
    result_simul_time[11] = duration.count();
    return(0);
}

