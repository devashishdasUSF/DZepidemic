#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include<math.h>

class seric
{
public:
	double betaI;
	double betaE;

	double alpha;
	double gamma;
	
	double nuI;
	double nuE;
	double kappa;

	double N;
	double e0;

	arma::vec cases;

	void intialize0(const Rcpp::List & x)
	{
		betaI  = x["beta_I"];
		betaE  = x["beta_E"];

		alpha  = x["alpha"];
		gamma  = x["gamma"];
		
		nuI    = x["nu_I"];
		nuE    = x["nu_E"];
		kappa  = x["kappa"];

		N      = x["N"];
		e0     = x["e0"];

		cases = Rcpp::as<arma::vec>(x["cases"]);
	}

	Rcpp::List sim(const double Tmax = 200)
	{
		double s = N - 1.0, e = e0, ec = 0, ic = 0.0;
		double i = 0., r = 0.0, c = 0.0;
		double t = 0.0, dt = 0.01;
		double d_s_e, d_e_i, d_i_r;
		double d_ec_ic, d_ic_r;
		double d_e_ec, d_i_ic;
		int ndays = (int)Tmax;

		Rcpp::NumericVector S(ndays), E(ndays), I(ndays), R(ndays), C(ndays);
		Rcpp::NumericVector Ic(ndays), Ec(ndays);
		int day = 0, counter = 0;
		while(t <= Tmax)
		{
			counter++;
			t += dt;

			d_s_e    = betaI/N*s*i*dt + betaE/N*s*e*dt + kappa * (betaI/N*s*ic*dt + betaE/N*s*ec*dt);
			d_e_i    = alpha * e * dt;
			d_ec_ic  = alpha * ec * dt;
			d_i_r    = gamma * i * dt;
			d_ic_r   = gamma * ic * dt;
			d_e_ec   = nuE * e * dt;
			d_i_ic   = nuI * i * dt;

			s += -d_s_e;
			e += d_s_e - d_e_i - d_e_ec;
			i += d_e_i - d_i_r - d_i_ic;
			r += d_i_r + d_ic_r;

			ec += d_e_ec - d_ec_ic;
			ic += d_ec_ic + d_i_ic - d_ic_r;
			c += d_e_ec + d_i_ic;

			if(counter % (int)(1/dt) == 0)
			{
				S[day] = s;
				E[day] = e;
				I[day] = i;
				R[day] = r;
				C[day] = c;

				Ic[day] = ic;
				Ec[day] = ec;

				day++;
			}
		}
		Rcpp::List simResults;
		simResults["s"] = S;
		simResults["e"] = E;
		simResults["i"] = I;
		simResults["r"] = R;
		simResults["c"] = C;
		simResults["ic"] = Ic;
		simResults["ec"] = Ec;
		return simResults;
	}
};

//[[Rcpp::export]]
Rcpp::List sim(const Rcpp::List &params, const double Tmax = 200.0)
{

	seric model;
	model.intialize0(params);
	Rcpp::List simResults = model.sim(Tmax);

	return simResults;
}



