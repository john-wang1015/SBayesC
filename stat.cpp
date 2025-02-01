#include "stat.hpp"
#include <random>
using namespace std;

void Stat::seedEngine(const int seed){
    if (seed) {
        srand(seed);
        engine.seed(seed);
    } else {
        srand((int)time(NULL));
        engine.seed((int)time(NULL));
    }
}

float Stat::Normal::sample(const float mean, const float variance){
    static std::random_device rd;
    static std::mt19937 engine(rd());  // Fresh random seed
    std::normal_distribution<float> dist(mean, sqrtf(variance));
    return dist(engine);
}

float Stat::Uniform::sample(float a, float b) {
    return a + (b - a) * ranf();
}

float Stat::InvChiSq::sample(const float df, const float scale) {
    static std::random_device rd;
    static std::mt19937 engine(rd());
    std::gamma_distribution<float> gamma_dist(0.5f * df, 0.5f); 
    float gamma_sample = gamma_dist(engine);
    return scale * df / gamma_sample;  
}

float Stat::Gamma::sample(const float shape, const float scale) {
    static std::random_device rd;
    static std::mt19937 engine(rd());  
    std::gamma_distribution<float> gamma_dist(shape, scale);  
    return gamma_dist(engine); 
}

float Stat::Beta::sample(const float a, const float b) {
    static std::random_device rd;
    static std::mt19937 engine(rd());

    std::gamma_distribution<float> gamma_a(a, 1.0f);
    std::gamma_distribution<float> gamma_b(b, 1.0f);

    float x = gamma_a(engine);
    float y = gamma_b(engine);

    return x / (x + y);  
}

unsigned Stat::Bernoulli::sample(const float p){
    static std::random_device rd;
    static std::mt19937 engine(rd());
    std::bernoulli_distribution dist(p);
    return dist(engine);
}

unsigned Stat::Bernoulli::sample(const VectorXf &p){
    static std::random_device rd;
    static std::mt19937 engine(rd());
    float rnd = std::generate_canonical<float, 10>(engine); 
    float cum = 0;
    unsigned ret = 0;
    for (unsigned i = 0; i < p.size(); ++i) {
        cum += p[i];
        if (rnd < cum) {
            ret = i;
            break;
        }
    }
    return ret;
}

float Stat::NormalZeroMixture::sample(const float mean, const float variance, const float p){
    return bernoulli.sample(p) ? normal.sample(mean, variance) : 0;
}

// Sample Dirichlet

VectorXf Stat::Dirichlet::sample(const int n, const VectorXf &irx){
    VectorXf ps(n);
    double sx = 0.0;
    for (int i = 0; i < n; i++)
    {
        ps[i] = gamma.sample(irx(i), 1.0);
        sx += ps[i];
    }
    ps = ps / sx;
    return ps;
}