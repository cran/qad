#include <Rcpp.h>
#include <cmath>
#include <cstdio>
using namespace Rcpp;

namespace std
{
    template<>
    struct hash<pair<double, double>>
    {
        size_t operator()(pair<double, double> const& p) const
        {
            auto hash1 = hash<double>{}(p.first); 
            auto hash2 = hash<double>{}(p.second); 
            return hash1*1000000000 + hash2; 
        }
    };
}


std::unordered_map<std::pair<double, double>, R_xlen_t>* pair_range (const NumericVector& x, const NumericVector& y, const NumericVector& x_range, const NumericVector& y_range) {
    std::unordered_map<std::pair<double, double>, R_xlen_t>* result = new std::unordered_map<std::pair<double, double>, R_xlen_t>();
    for (R_xlen_t i = 0; i < x.length(); i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();
        if (x_range[(R_xlen_t)x[i]-1] > 1 && y_range[(R_xlen_t)y[i]-1] > 1) {
            std::pair<double, double> key(x[i], y[i]);
            if ( result->find(key) == result->end() ) {
                (*result)[key] = 0;
            }
            (*result)[key] = (*result)[key] + 1;
        }
    }


    return result;
}


NumericVector range(const NumericVector& x) {
    NumericVector result(x.length());

    for (R_xlen_t i = 0; i < x.length(); i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();
        result[(R_xlen_t) x[i] - 1] = result[(R_xlen_t) x[i] - 1] + 1;
    }

    return result;
}


//' @export
// [[Rcpp::export]]
NumericMatrix build_checkerboard_weights (const NumericVector& X, const NumericVector& Y, R_xlen_t resolution) {
    R_xlen_t sample_size = std::min(X.length(), Y.length());
    R_xlen_t x_upper;
    R_xlen_t x_lower;
    R_xlen_t y_upper;
    R_xlen_t y_lower;
    NumericMatrix result(resolution, resolution);
    NumericVector x_range = range(X);
    NumericVector y_range = range(Y);
    R_xlen_t rx;
    R_xlen_t ry;
    R_xlen_t rp;
    double lambda_x;
    double lambda_y;
    std::unordered_map<std::pair<double, double>, R_xlen_t>* p_range = pair_range(X, Y, x_range, y_range);

    for (R_xlen_t i = 0; i < sample_size; i++) {
        if(i%100000 == 0) Rcpp::checkUserInterrupt();
        rx = x_range[(R_xlen_t)X[i] - 1];
        ry = y_range[(R_xlen_t)Y[i] - 1];
        if (rx > 1 && ry > 1) {
            std::pair<double, double> key(X[i], Y[i]);
            rp = (*p_range)[key];
            (*p_range)[key] = 0;
        } else {
            rp = 1;
        }

        if (rp != 0) {

            x_upper = std::ceil((X[i]) / (double)sample_size * resolution);
            x_lower = std::max(std::ceil((X[i] - (double)rx) / (double)sample_size * resolution), 1.);
            y_upper = std::ceil((Y[i]) / (double)sample_size * resolution);
            y_lower = std::max(std::ceil((Y[i] - (double)ry) / (double)sample_size * resolution), 1.);

            for (R_xlen_t x = x_lower; x <= x_upper; x++) {
                lambda_x = std::min(X[i], (double)x / (double)resolution * sample_size) - std::max(X[i] - (double)rx, (double)(x - 1) / (double)resolution * sample_size);
                for (R_xlen_t y = y_lower; y <= y_upper; y++) {
                    lambda_y = std::min(Y[i], (double)y / (double)resolution * sample_size) - std::max(Y[i] - (double)ry, (double)(y - 1) / (double)resolution * sample_size);
                    result(x-1, y-1) = result(x-1, y-1) + (lambda_x * lambda_y * rp)/(double)(sample_size * rx * ry);
                }
            }

        }
    }

    delete p_range;

    return result;
}

//' @export
// [[Rcpp::export]]
double local_kernel_integral(const NumericMatrix& A, R_xlen_t x, R_xlen_t y, R_xlen_t N, double y_sum) {
    double ak;
    double d;
    double k;
    double i_D;
    double k_D;
    double y0;
    double I1;
    double I2;

    ak = A(x-1, y-1);
    d = (y_sum + ak - y*ak)*N;
    k = (ak*N)*N;

    i_D = d;
    k_D = k - 1;

    if ((i_D + k_D * (double)(y-1) / (double)N)*(i_D + k_D * (double)y / (double)N) >= 0) {
        return std::abs(i_D / (double)N + k_D / 2. * ((double)(y*y)/(double)(N*N) - (double)((y-1)*(y-1))/(double)(N*N)));
    } else {
        y0 = i_D / -k_D;
        I1 = std::abs(i_D * (y0 - (double)(y-1) / (double)N) + k_D / 2. * (y0*y0 - (double)((y-1)*(y-1))/(double)(N*N)));
        I2 = std::abs(i_D * ((double)y / (double)N - y0) + k_D / 2. * ((double)(y*y)/(double)(N*N) - y0*y0));
        return I1+I2;
    }
}

//' @export
// [[Rcpp::export]]
double D1_Pi(const NumericMatrix& A, R_xlen_t resolution) {
    double result = 0;
    double y_sum;
    for (R_xlen_t x = 0; x < resolution; x++) {
        if(x%100000 == 0) Rcpp::checkUserInterrupt();
        y_sum = 0;
        for (R_xlen_t y = 0; y < resolution; y++) {
            result = result + local_kernel_integral(A, x+1, y+1, resolution, y_sum);
            y_sum = y_sum + A(x, y);
        }
    }

    return result/resolution;
}
