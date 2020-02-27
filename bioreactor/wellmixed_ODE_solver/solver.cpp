#include<iostream>
#include<vector>
#include<fstream>
#include <boost/math/special_functions/gamma.hpp>

const int X=0;
const int O2=1;
const int G=2;
const int Xy=3;
const int A=4;
const int B=5;

void get_rhs(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars)
{
  // JL -- not sure how you want to handle constants. I'm just going to dump them here
  // ------- MODEL CONSTANTS -------
  double y_xs = 0.009;
  double y_as = 1.01;
  double y_bs = 0.88;
  double y_os = 0.0467;

  double x_max = 11;
  double qs_max = 17;
  double o2_max = 0.214;

  double K_e = 0.214;
  double K_s = 31;

  double alpha_s = 3;
  double beta_s = 12;
  double alpha_e = 3;
  double beta_e = 12;

  double kLa = 5;

  // calculate chi_i
  // FIXME: we need a library for calculating the gamma incomplete function.
  // Using general "gammainc" for now. In my python code I'm using scipy.special.gammainc
  // (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainc.html)
  // There are a couple libraries that I'm seeing on google, but I don't feel comfortable
  // making that decision.
  //
  // Also, a word of explanation about the form of the gamma incomplete call. This function
  // is usually expressed as g(a, x), where a and x are obligate positive. Therefore, we need
  // to restrict x to be no lower than 0. Additionally, we need to make sure we don't hit a
  // divide by zero error, so we add a very small value to the denominator.
  double chi_s = boost::math::gamma_p(alpha_s, beta_s*std::max(solnvec[2]/(solnvec[3]+1e-8),0.0));
  double chi_e = boost::math::gamma_p(alpha_e, beta_e*std::max(solnvec[1]/(solnvec[4]+1e-8),0.0));
  double chi_p = 0.3;

  // calculate q_s
  double F_s = (solnvec[2] + solnvec[3])/(solnvec[2] + solnvec[3] + K_s);
  double F_e = (solnvec[1] + beta_e*solnvec[4])/(solnvec[1] + beta_e*solnvec[4] + K_e);
  double q_s = qs_max*F_s*F_e;

  // calculate intermediate rates
  double rar = chi_p*y_as*q_s*X;
  double rbr = (1-chi_p)*y_bs*q_s*X;

  double our = -chi_e*y_os*q_s*X;
  double otr = kLa *(o2_max - solnvec[1]);
  double rae = -(1-chi_e)*y_os*q_s*X;
  double rbe = -rae;


  // calculate final rates
  rhs[0] = y_xs*q_s*solnvec[0]*(1 - solnvec[0]/x_max);
  rhs[1] = our+otr;
  rhs[2] = -chi_s     *q_s*solnvec[0];
  rhs[3] = -(1-chi_s) *q_s*solnvec[0];
  rhs[4] = rar+rae;
  rhs[5] = rbr+rbe;

}

void advance(std::vector<double>& solnvec,int nvars,double t_now,double t_adv,double dt)
{
    double current_time=t_now;
    double final_time=t_now+t_adv;

    std::vector<double> rhs(nvars);
    std::vector<double> solnvec_n(nvars);

    while(current_time < final_time)
    {
        current_time += dt;

        //at current time level n
        solnvec_n=solnvec;

        //Doing RK23

        //stage 1
        get_rhs(rhs,solnvec,current_time,nvars);
        for(int i=0;i<nvars;i++)
        {
            solnvec[i] = solnvec_n[i] + 0.5*rhs[i]*dt;
        }
        
        //stage 2
        get_rhs(rhs,solnvec,current_time,nvars);
        for(int i=0;i<nvars;i++)
        {
            solnvec[i] = solnvec_n[i] + rhs[i]*dt;
        }
    }
}

int main()
{
    int nvars=6;
    std::vector<double> solnvec(nvars);


    double final_time=24.0;
    double advance_time=0.5;
    double current_time=0.0;
    double dt=advance_time/10.0;

    std::ofstream outfile("timehist.dat");

    //set initial conditions
    solnvec[X]  = 0.5;
    solnvec[O2] = 0.214;
    solnvec[G]  = 500.0;
    solnvec[Xy] = 250.0;
    solnvec[A]  = 0.0;
    solnvec[B]  = 0.0;

    //write initial condition
    outfile<<current_time<<"\t";
    for(int i=0;i<nvars;i++)
    {
        outfile<<solnvec[i]<<"\t";
    }
    outfile<<"\n";


    while(current_time < final_time)
    {
        current_time += advance_time;

        advance(solnvec,nvars,current_time,advance_time,dt);

        outfile<<current_time<<"\t";
        for(int i=0;i<nvars;i++)
        {
            outfile<<solnvec[i]<<"\t";
        }
        outfile<<"\n";
    }

    outfile.close();
    return(0);
}
