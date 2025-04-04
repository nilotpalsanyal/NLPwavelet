//######## NLPwavelet package C++ functions #########//

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <cmath>  // include the cmath library for 'fabs' if using
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double doublefact(int x) {
    int y = (x % 2 == 0) ? x - 1 : x;
    double numerator = std::tgamma(y + 2); // Gamma function: factorial(n) = Gamma(n+1)
    double denominator = std::pow(2, (y + 1) / 2) * std::tgamma((y + 3) / 2);
    return numerator / denominator;
}

// [[Rcpp::export]]
double M1_func(int r, double d, double tau1_l, double sigmasq) {
    double sqrt_term = std::sqrt(tau1_l / (1.0 + tau1_l));
    double scaled_d = d / std::sqrt(sigmasq);
    double factor = sqrt_term * scaled_d;
    
    arma::vec i = arma::regspace(0, r);

    // Apply std::lgamma element-wise using .transform()
    arma::vec lgamma_2i1 = i;
    lgamma_2i1.transform([](double val) { return std::lgamma(2 * val + 1); });

    arma::vec lgamma_ri = i;
    lgamma_ri.transform([r](double val) { return std::lgamma(r - val + 1); });

    double lgamma_2r1 = std::lgamma(2 * r + 1);

    arma::vec coeffs = arma::exp(lgamma_2r1 - lgamma_2i1 - lgamma_ri - (r - i) * std::log(2));
    arma::vec powers = arma::pow(arma::vec(i.n_elem, arma::fill::value(factor)), 2 * i);
    
    double M1 = arma::dot(coeffs, powers);
    M1 /= doublefact(2 * r - 1);
    
    return M1;
}

// [[Rcpp::export]]
double M2_func(int r, double d, double tau1_l, double sigmasq) {
    double sqrt_term = std::sqrt(tau1_l / (1.0 + tau1_l));
    double scaled_d = d / std::sqrt(sigmasq);
    double factor = sqrt_term * scaled_d;
    
    arma::vec i = arma::regspace(1, r + 1);

    // Apply std::lgamma element-wise using .transform()
    arma::vec lgamma_2i = i;
    lgamma_2i.transform([](double val) { return std::lgamma(2 * val); });

    arma::vec lgamma_r1i = i;
    lgamma_r1i.transform([r](double val) { return std::lgamma(r + 2 - val); });

    double lgamma_2r2 = std::lgamma(2 * r + 2);

    arma::vec coeffs = arma::exp(lgamma_2r2 - lgamma_2i - lgamma_r1i - (r + 1 - i) * std::log(2));
    arma::vec powers = arma::pow(arma::vec(i.n_elem, arma::fill::value(factor)), 2 * i - 1);
    
    double M2 = arma::dot(coeffs, powers);
    M2 /= doublefact(2 * r - 1);
    
    return M2;
}

// [[Rcpp::export]]
double h_func(double x, double d, double nu, double tau2_l, double sigmasq) {
    double log_term = std::log(std::pow(std::abs(x), -nu - 1));
    double exp_term = -(1.0 / (2 * sigmasq)) * (x * x - 2 * x * d) - (tau2_l * sigmasq) / (x * x);
    return std::exp(std::min(709.0, log_term + exp_term));
}

// [[Rcpp::export]]
double L_func(double x, double d, double nu, double tau2_l, double sigmasq) {
    return -(nu + 1) * std::log(std::abs(x)) - (1.0 / (2 * sigmasq)) * (x * x - 2 * x * d) - (tau2_l * sigmasq) / (x * x);
}

// [[Rcpp::export]]
double L_dd_func(double x, double nu, double tau2_l, double sigmasq) {
    return (nu + 1) / (x * x) - 1.0 / sigmasq - (6 * tau2_l * sigmasq) / (x * x * x * x);
}

arma::cx_vec find_roots(const arma::vec& coeffs) {
    int n = coeffs.n_elem - 1;
    arma::mat companion = arma::zeros(n, n);
    
    // Set subdiagonal elements to 1
    for (int i = 1; i < n; i++) {
        companion(i, i - 1) = 1.0;
    }
    
    // Set the last row using the polynomial coefficients
    companion.row(0) = -coeffs.subvec(1, n).t() / coeffs(0);
    
    // Compute the eigenvalues, which correspond to the roots
    return arma::cx_vec(arma::eig_gen(companion));
}

// [[Rcpp::export]]
Rcpp::List Lap_approx(double d, double nu, double tau2_l, double sigmasq) {
    // Define polynomial coefficients for: x^4 - d*x^3 + (nu+1)*sigmasq*x - 2*tau2_l*sigmasq^2 = 0
    arma::vec coeffs = {1, -d, (nu + 1) * sigmasq, 0, -2 * tau2_l * sigmasq * sigmasq};

    // Find polynomial roots
    arma::cx_vec roots = find_roots(coeffs);

    // Extract real roots
    arma::vec real_roots;
    for (size_t i = 0; i < roots.n_elem; i++) {
        if (std::abs(roots(i).imag()) < 1e-6) { // Ensure the root is real
            real_roots.insert_rows(real_roots.n_elem, arma::vec({roots(i).real()}));
        }
    }

    if (real_roots.n_elem == 0) {
        return Rcpp::List::create(Rcpp::Named("d.max") = 1e-50, Rcpp::Named("sigma.d.max") = NA_REAL);
    }

    // Compute L_func and L_dd_func values
    arma::vec L_values(real_roots.n_elem);
    arma::vec L_dd_values(real_roots.n_elem);

    for (size_t i = 0; i < real_roots.n_elem; i++) {
        double x = real_roots(i);
        L_values(i) = L_func(x, d, nu, tau2_l, sigmasq);
        L_dd_values(i) = L_dd_func(x, nu, tau2_l, sigmasq);
    }

    // Find indices where L_dd_values < 0
    arma::uvec negative_L_dd_indices = arma::find(L_dd_values < 0);

    if (negative_L_dd_indices.n_elem == 0) {
        return Rcpp::List::create(Rcpp::Named("d.max") = 1e-50, Rcpp::Named("sigma.d.max") = NA_REAL);
    }

    // Select the root with the maximum L_func value
    arma::uvec best_indices = negative_L_dd_indices.elem(arma::find(L_values(negative_L_dd_indices) == L_values(negative_L_dd_indices).max()));
    
    double d_max;
    if (best_indices.n_elem == 1) {
        d_max = real_roots(best_indices(0));
    } else {
        // If multiple candidates, choose the one with the highest L_dd value
        d_max = real_roots(best_indices(arma::index_max(L_dd_values(best_indices))));
    }

    // Compute sigma.d.max
    double L_dd_dmax = L_dd_func(d_max, nu, tau2_l, sigmasq);
    double sigma_d_max = std::sqrt(-1.0 / L_dd_dmax);

    return Rcpp::List::create(Rcpp::Named("d.max") = d_max, Rcpp::Named("sigma.d.max") = sigma_d_max);
}

// [[Rcpp::export]]
double lhood_contrib(double d, int r, double nu, double M1, double gamma1_l, double gamma2_l, double tau1_l, double tau2_l, double sigmasq, double d_star, double sigma_star) {
    double term1 = std::exp(std::max(-745.0, std::log(gamma1_l) + (-r - 0.5) * std::log(1 + tau1_l) + std::log(M1) - d * d / (2 * sigmasq * (1 + tau1_l))));
    double term2 = std::exp(std::max(-745.0, std::log(1 - gamma1_l) + std::log(gamma2_l) + nu / 2 * std::log(tau2_l * sigmasq) - std::lgamma(nu / 2) - d * d / (2 * sigmasq) + std::log(std::sqrt(2 * M_PI) * sigma_star) + std::log(h_func(d_star, d, nu, tau2_l, sigmasq))));
    double term3 = std::exp(std::max(-745.0, std::log(1 - gamma1_l) + std::log(1 - gamma2_l) - d * d / (2 * sigmasq)));
    return std::log(term1 + term2 + term3);
}

// [[Rcpp::export]]
Rcpp::List post_odds_func(double d, int r, double nu, double M1, double sigmasq, double gamma1_l, double gamma2_l, double tau1_l, double tau2_l, double d_star, double sigma_star) {
    double O1 = (gamma1_l != 0 && gamma2_l != 1) ? (gamma1_l / ((1 - gamma1_l) * (1 - gamma2_l)) * std::pow(1 + tau1_l, -r - 0.5) * M1 * std::exp(1 / (2 * sigmasq) * tau1_l / (1 + tau1_l) * d * d)) : 0;
    double O2 = gamma2_l / (1 - gamma2_l) * std::pow(tau2_l * sigmasq, nu / 2) * (1 / std::tgamma(nu / 2)) * std::sqrt(2 * M_PI) * sigma_star * h_func(d_star, d, nu, tau2_l, sigmasq);
    return Rcpp::List::create(Rcpp::Named("O1") = O1, Rcpp::Named("O2") = O2);
}

// [[Rcpp::export]]
Rcpp::List post_mixprobs_func(Rcpp::List post_odds) {
    double O1 = post_odds["O1"];
    double O2 = post_odds["O2"];
    double p1, p2;
    if(O1 != INFINITY && O2 != INFINITY) {
        p1 = O1 / (O1 + O2 + 1);
        p2 = O2 / (O1 + O2 + 1);
    } else if(O1 == INFINITY && O2 == INFINITY) {
        p1 = 0.5;
        p2 = 0.5;
    } else {
        p1 = (O1 == INFINITY) ? 1 : 0;
        p2 = (O2 == INFINITY) ? 1 : 0;
    }
    return Rcpp::List::create(Rcpp::Named("p1") = p1, Rcpp::Named("p2") = p2);
}

// [[Rcpp::export]]
double post_odds_func_indiv(std::string method, double d, double r, double sigmasq, 
                            double gamma_l, double tau_l, double M1 = 1, 
                            Rcpp::Nullable<double> d_star = R_NilValue, 
                            Rcpp::Nullable<double> sigma_star = R_NilValue) {
    double O;
    
    if (method == "mom") {
        O = (gamma_l / (1 - gamma_l)) * pow(1 + tau_l, -r - 0.5) * M1 * 
            exp((1 / (2 * sigmasq)) * (tau_l / (1 + tau_l)) * pow(d, 2));
        if (std::isnan(O)) O = INFINITY;
    } else {
        if (d_star.isNull() || sigma_star.isNull()) {
            Rcpp::stop("d_star and sigma_star must be provided if method is not mom.");
        }
        double d_star_val = Rcpp::as<double>(d_star);
        double sigma_star_val = Rcpp::as<double>(sigma_star);
        
        O = (gamma_l / (1 - gamma_l)) * pow(tau_l * sigmasq, r / 2) * 
            (1.0 / R::gammafn(r / 2)) * sqrt(2 * M_PI) * sigma_star_val * 
            h_func(d_star_val, d, r, tau_l, sigmasq);
    }
    
    return O;
}

// [[Rcpp::export]]
double lhood_contrib_indiv(std::string method, double d, double r, 
                                double gamma_l, double tau_l, double sigmasq, 
                                double M1 = 1, 
                                Rcpp::Nullable<double> d_star = R_NilValue, 
                                Rcpp::Nullable<double> sigma_star = R_NilValue) {    
    double out;
    double lower_bound = -745; // Prevent underflow in exponentiation

    if (method == "mom") {
        out = log( exp(std::max(lower_bound, log(gamma_l) + (-r - 0.5) * log(1 + tau_l) + 
                              log(M1) - d * d / (2 * sigmasq * (1 + tau_l)))) 
                 + exp(std::max(lower_bound, log(1 - gamma_l) - d * d / (2 * sigmasq))) );
    } else {
        if (d_star.isNull() || sigma_star.isNull()) {
            Rcpp::stop("d_star and sigma_star must be provided if method is not mom.");
        }
        double d_s = Rcpp::as<double>(d_star);
        double sigma_s = Rcpp::as<double>(sigma_star);

        // Assuming h_func is another Rcpp function
        // double h_val = Rcpp::as<double>(Rcpp::Function("h_func")(d_s, d, r, tau_l, sigmasq));
        double h_val = h_func(d_s, d, r, tau_l, sigmasq);

        out = log( exp(std::max(lower_bound, log(gamma_l) + (r / 2) * log(tau_l * sigmasq) - 
                              log(R::gammafn(r / 2)) - d * d / (2 * sigmasq) + 
                              log(sqrt(2 * M_PI) * sigma_s) + log(h_val))) 
                 + exp(std::max(lower_bound, log(1 - gamma_l) - d * d / (2 * sigmasq))) );
    }

    return out;
}



// library(microbenchmark)
// microbenchmark(doublefact(100),doublefact1(100),times=100,check="equal")
// microbenchmark(M1_func(1,2,0.4,2),M1_func1(1,2,0.4,2),times=100,check="equal")
// microbenchmark(M2_func(1,2,0.4,2),M2_func1(1,2,0.4,2),times=100,check="equal")
// microbenchmark(h_func(2,1,2,0.4,2),h_func1(2,1,2,0.4,2),times=100,check="equal")
// microbenchmark(L_func(2,1,2,0.4,2),L_func1(2,1,2,0.4,2),times=100,check="equal")
// microbenchmark(L_dd_func(2,2,0.4,2),L.dd_func1(2,2,0.4,2),times=100,check="equal")
// microbenchmark(Lap_approx(2,1,2,2),Lap_approx1(2,1,2,2),times=100,check="equal")
// microbenchmark(indiv_lhood_contrib(2,1,1,1.5,0.3,0.2,1.5,1.5,2,2,2),indiv_lhood_contrib1(2,1,1,1.5,0.3,0.2,1.5,1.5,2,2,2),times=100,check="equal")
// microbenchmark(post_odds_func(2,1,1,1.5,2,0.3,0.2,1.5,1.5,2,2),post_odds_func1(2,1,1,1.5,2,0.3,0.2,1.5,1.5,2,2),times=100,check="equal")
// microbenchmark(post_mixprobs_func(list(O1=0.5,O2=0.3)),post_mixprobs_func1(list(O1=0.5,O2=0.3)),times=100,check="equal")
// microbenchmark(post_odds_func_indiv("mom",2,1,1.5,2,0.3,1.5,2,2),post_odds_func_indiv1("mom",2,1,1.5,2,0.3,1.5,2,2),times=100,check="equal")




