#include "gamma_fit.h"

double gamme_mle_update(double logmeanx, double meanlogx, double a) {
	double ia = 1 / a;
  boost::math::digamma(2);
	double z = ia + (meanlogx - logmeanx + log(a) - boost::math::digamma(a)) / (pow(abs(a), 2) * (ia - boost::math::trigamma(a)));
	return 1/z;
}



boost::math::gamma_distribution<> fit_gamma(int n, double sumx, double logsumx) {
	const int maxiter = 1000;
	const double tol = 1e-16;

	double meanx = sumx / n;
	double logmeanx = log(meanx);
	double meanlogx = logsumx / n;
	double a = (logmeanx - meanlogx) / 2;


	bool converged = false;
	int t = 0;
	while (!converged && t < maxiter) {
		t++;
		double a_old = a;
		a = gamme_mle_update(logmeanx, meanlogx, a);
		converged = abs(a - a_old) <= tol;
	}

	boost::math::gamma_distribution<> distr(a, meanx/a);
	return distr;
}

boost::math::gamma_distribution<> fit_gamma(const arma::mat &vals) {
	int n = vals.size();
	double sumx = 0;
	double logsumx = 0;

	for (int i = 0; i < n; i++) {
		sumx += vals[i];
		logsumx += log(vals[i]);
	}

	return fit_gamma(n, sumx, logsumx);
}
