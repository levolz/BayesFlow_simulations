#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double pexGAUS(double q,
               double mu = 5., double sigma = 1., double tau = 1.,
               bool lower_tail = true, bool log_p = false) {

	if (tau <= 0 || sigma <= 0) {
		double cdf = NA_REAL;
		return cdf;
	}

	double cdf;
	if (sigma < 1e-4){
		cdf = R::pexp(q - mu, tau, lower_tail, log_p);
		return cdf;
	}

	if (!traits::is_infinite<REALSXP>(q)){
		if (tau > .05 * sigma) {
			double z_i = q - mu - (sigma * sigma) / tau;
				cdf = R::pnorm((q - mu) / sigma, 0., 1., true, false) - std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) + (std::pow((mu + (sigma * sigma / tau)), 2) - mu * mu - 2. * q * (sigma * sigma / tau)) / (2. * sigma * sigma));
		} else {
			cdf = R::pnorm(q, mu, sigma, true, false);
		}
	} else {
		if (q < 0) {
			cdf = 0.;
		} else {
			cdf = 1.;
		}
	}

	if (!lower_tail) {
		cdf = 1. - cdf;
	}
	if (log_p) {
		cdf = std::log(cdf);
	}
	return cdf;
}


// [[Rcpp::export]]
NumericVector dEXG(NumericVector x,
  double mu = 5., double sigma = 1., double tau = 1, double a = 0, bool log_d = false)
{

	double Fa = pexGAUS(a,mu,sigma,tau);
  double cFa = 1 - Fa;

	int n = x.size();
	if (tau <= 0 || sigma <= 0) {
		NumericVector pdf(n, NA_REAL);
		return pdf;
	}

	NumericVector pdf(n);
	if (sigma < 1e-4){
		for (int i = 0; i < n; i++){
			pdf[i] = (R::dexp(x[i] - mu, tau, false))/cFa;
		}
		return pdf;
	}

	for (int i = 0; i < n; i++){
		if (tau > .05 * sigma){
			double z_i = x[i] - mu - (sigma * sigma) / tau;
			pdf[i] = std::exp(- std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau +
			         std::log(R::pnorm(z_i / sigma, 0., 1., true, false)))/cFa;
		} else {
			pdf[i] = std::exp(R::dnorm(x[i], mu, sigma, true))/cFa;
		}
	}
	if (log_d){
		for(int i = 0; i < n; i++){
			pdf[i] = std::log(pdf[i]);
		}
	}
	return pdf;
}



// [[Rcpp::export]]
NumericVector pEXG(NumericVector q,
  double mu = 5., double sigma = 1., double tau = 1., double a = 0)

{

  double Fa = pexGAUS(a,mu,sigma,tau);
  double cFa = 1 - Fa;


 	int n = q.size();
	if (tau <= 0 || sigma <= 0) {
		NumericVector cdf(n, NA_REAL);
		return cdf;
	}

	NumericVector cdf(n);
	if (sigma < 1e-4){
		for (int i = 0; i < n; i++){
			cdf[i] = (R::pexp(q[i] - mu, tau,true,false)-Fa)/cFa;
		}
		return cdf;
	}

	for (int i = 0; i < n; i++){
		if (!traits::is_infinite<REALSXP>(q[i])){
			if (q[i] <= a) {
			  cdf[i] = 0;
			} else {
			  if (tau > .05 * sigma){
				  double z_i = q[i] - mu - (sigma * sigma) / tau;
				  cdf[i] = ((R::pnorm((q[i] - mu) / sigma, 0., 1., true, false) -
				           std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) +
				           (std::pow((mu + (sigma * sigma / tau)), 2) - mu * mu -
				           2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma)))-Fa)/cFa;
			  } else {
				  cdf[i] = (R::pnorm(q[i], mu, sigma, true, false)-Fa)/cFa;
			  }
			}
		} else {
			if (q[i] < 0) {
				cdf[i] = 0.;
			} else {
				cdf[i] = 1;
			}
		}
	}

	return cdf;
}


const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)

// [[Rcpp::export]]
double pigt0(double t, double k = 1., double l = 1.){
	//if (t <= 0.){
	//  return 0.;
	//}
	double mu = k / l;
	double lambda = k * k;

	double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. + t/mu), 0., 1., true, false);
	double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. - t/mu), 0., 1., true, false);

	return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

// [[Rcpp::export]]
double digt0(double t, double k = 1., double l = 1.){
	//if (t <= 0.) {
	//  return 0.;
	//}
	double lambda = k * k;
	double e;
	if (l == 0.) {
		e = -.5 * lambda / t;
	} else {
		double mu = k / l;
		e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
	}
	return std::exp(e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI));
}

// [[Rcpp::export]]
double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
	if (t <= 0.){
		return 0.;
	}
	if (a < threshold){
		return pigt0(t, k, l);
	}

	double sqt = std::sqrt(t);
	double lgt = std::log(t);
	double cdf;

	if (l < threshold){
		double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
		double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1;

		double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
		double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

		cdf = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
	} else {
		double t1a = std::exp(- .5 * std::pow(k - a - t * l, 2) / t);
		double t1b = std::exp(- .5 * std::pow(a + k - t * l, 2) / t);
		double t1 = std::exp(.5* (lgt - M_LN2 - L_PI)) * (t1a - t1b);

		double t2a = std::exp(2. * l * (k - a) + R::pnorm(- (k - a + t * l) / sqt, 0., 1., true, true));
		double t2b = std::exp(2. * l * (k + a) + R::pnorm(- (k + a + t * l) / sqt, 0., 1., true, true));
		double t2 = a + (t2b - t2a) / (2. * l);

		double t4a = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
		double t4b = 2. * R::pnorm((k - a) / sqt - sqt * l, 0., 1., true, false) - 1.;
		double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;

		cdf = .5 * (t4 + t2 + t1) / a;
	}
	if (cdf < 0. || std::isnan(cdf)) {
		return 0.;
	}
	return cdf;
}

// [[Rcpp::export]]
double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10){
	if (t <= 0.){
		return 0.;
	}
	if (a < threshold){
		return digt0(t, k, l);
	}
	double pdf;
	if (l < threshold){
		double term = std::exp(- (k - a) * (k - a) / (2. * t)) - std::exp(- (k + a) * (k + a) / (2. * t));
		pdf = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(term) - M_LN2 - std::log(a));
	} else {
		double sqt = std::sqrt(t);

		double t1a = - std::pow(a - k + t * l, 2) / (2. * t);
		double t1b = - std::pow(a + k - t * l, 2) / (2. * t);
		double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

		double t2a = 2. * R::pnorm((- k + a) / sqt + sqt * l, 0., 1., true, false) - 1.;
		double t2b = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
		double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);

		pdf = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
	}
	if (pdf < 0. || std::isnan(pdf)) {
		return 0.;
	}
	return pdf;
}

// [[Rcpp::export]]
NumericVector pWald(NumericVector t, double v,
                    double B, double A, double t0){
	int n = t.size();
	NumericVector cdf(n);
	for (int i = 0; i < n; i++){
		t[i] = t[i] - t0;
		if (t[i] <= 0){
			cdf[i] = 0.;
		} else {
			cdf[i] = pigt(t[i], B + .5 * A, v, .5 * A);
		}
	}
	return cdf;
}

// [[Rcpp::export]]
NumericVector dWald(NumericVector t, double v,
                    double B, double A, double t0){
	int n = t.size();
	NumericVector pdf(n);
	for (int i = 0; i < n; i++){
		t[i] = t[i] - t0;
		if (t[i] <= 0){
			pdf[i] = 0.;
		} else {
			pdf[i] = digt(t[i], B + .5 * A, v, .5 * A);
		}
	}
	return pdf;
}


// [[Rcpp::export]]
NumericVector n1Wald(NumericMatrix dt,
  NumericVector v, NumericVector B, NumericVector A, double t0 = 0){

	int n = v.size();
	NumericVector out(dt.nrow());
	out = dWald(dt(0, _), v[0], B[0], A[0], t0);
	for (int i = 1; i < n; i++){
		out = out * (1-pWald(dt(i, _), v[i], B[i], A[i], t0));
	}
	return out;
}

// [[Rcpp::export]]
NumericVector n1Waldexg(NumericMatrix dt,
                        NumericVector v, NumericVector B, NumericVector A,
                        double t0, double minEXG = 0, int Si = 1){
	int n = v.size();
	NumericVector out(dt.nrow());
	Si = Si - 1;
	if (Si==0){
		out = dEXG(dt(0, _), v[0], B[0], A[0], minEXG);
	} else {
		out = dWald(dt(0, _), v[0], B[0], A[0], t0);
	}
	for (int i = 1; i < n; i++){
	  if (Si==i){
		  out = out = out * (1. - pEXG(dt(i, _), v[i], B[i], A[i], minEXG));
	  } else {
		  out = out * (1. - pWald(dt(i, _), v[i], B[i], A[i], t0));
	  }
	}
	return out;
}


// [[Rcpp::export]]
NumericVector stopfn_rdex(NumericVector t,
                          NumericVector vj, NumericVector Bj, NumericVector Aj,
                          double t0, double SSD, double minEXG, int Si = 1){

	int n = vj.size();
  NumericVector tmp( n * t.size());
  tmp = rep_each(t, n);
	NumericMatrix dt(n , t.size(), tmp.begin());
	Si = Si - 1;
  for (int i = 0; i < n; i++){
    if (Si == i ) {
      dt(i, _) = dt(i, _) + t0;
    } else {
      dt(i, _) = dt(i, _) + SSD;
    }
  }
	return n1Waldexg(dt, vj, Bj, Aj, t0 = 0, minEXG, Si+1);
}




