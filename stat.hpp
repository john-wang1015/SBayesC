#ifndef STAT_HPP
#define STAT_HPP

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

namespace Distributions{


    class Gamma {
        public:
            float sample(const float shape, const float scale);
    };

    class scaledInvChi{
        public:
            float sample(const float df, const float scale);
    };

    class beta{
        public:
            float sample(const float a, const float b);
    };

    class Dirichlet{
        public:
            Gamma gamma;
            VectorXf sample(const int n, const VectorXf &irx);
    };

}


#endif // STAT_HPP