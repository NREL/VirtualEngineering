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
  // ------- MODEL CONSTANTS -------
  double y_xs = 0.009;
  double y_as = 1.01;
  double y_bs = 0.88;
  double y_os = 0.0467;

  double x_max = 11;
  double qs_max = 17;
  double o2_max = 0.214;

  double K_e = 0.0214;
  double K_s = 31;

  double alpha_s = 3;
  double beta_s = 12;
  double alpha_e = 1;
  double beta_e = 1e3;

  double kLa = 5;

  // calculate chi_i
  // A word of explanation about the form of the gamma incomplete call. This function
  // is usually expressed as g(a, x), where a and x are obligate positive. Therefore, we need
  // to restrict x to be no lower than 0. Additionally, we need to make sure we don't hit a
  // divide by zero error, so we condition use of the gamma function on the denominator being
  // non-zero.
  double chi_s = 1;
  double chi_e = 1;

  if (solnvec[Xy]>1e-8) {
    double sRatio = std::max(solnvec[G]/(solnvec[Xy]), 0.0);
    chi_s = boost::math::gamma_p(alpha_s, beta_s*sRatio);
  }
  
  if (solnvec[A]>1e-8) {
    double eRatio = std::max(solnvec[O2]/(solnvec[A]), 0.0);
    chi_e = boost::math::gamma_p(alpha_e, beta_e*eRatio);
  }
  
  double chi_p = 0.3;

  // calculate q_s
  double F_s = (solnvec[G] + solnvec[Xy])/(solnvec[G] + solnvec[Xy] + K_s);
  double F_e = (solnvec[O2] + solnvec[A]/beta_e)/(solnvec[O2] + solnvec[A]/beta_e + K_e);
  double q_s = qs_max*F_s*F_e;

  // calculate intermediate rates
  double rar = chi_p*y_as*q_s*solnvec[X];
  double rbr = (1-chi_p)*y_bs*q_s*solnvec[X];

  double our = -chi_e*y_os*q_s*solnvec[X];
  double otr = kLa *(o2_max - solnvec[O2]);
  double rae = -(1-chi_e)*y_as*q_s*solnvec[X];
  double rbe = -rae;


  // calculate final rates
  rhs[X] = y_xs*q_s*solnvec[X]*(1 - solnvec[X]/x_max);
  //rhs[O2] = our+otr;
  rhs[02] = 0;
  rhs[G] = -chi_s     *q_s*solnvec[X];
  rhs[Xy] = -(1-chi_s) *q_s*solnvec[X];
  rhs[A] = rar+rae;
  rhs[B] = rbr+rbe;

}

double get_our(std::vector<double> solnvec,double t,int nvars)
{
  // ------- MODEL CONSTANTS -------
  double y_os = 0.0467;

  double x_max = 11;
  double qs_max = 17;
  double o2_max = 0.214;

  double K_e = 0.0214;
  double K_s = 31;

  double alpha_s = 3;
  double beta_s = 12;
  double alpha_e = 1;
  double beta_e = 1e3;


  // calculate chi_i
  double chi_s = 1;
  double chi_e = 1;

  if (solnvec[Xy]>1e-8) {
    double sRatio = std::max(solnvec[G]/(solnvec[Xy]), 0.0);
    chi_s = boost::math::gamma_p(alpha_s, beta_s*sRatio);
  }
  
  if (solnvec[A]>1e-8) {
    double eRatio = std::max(solnvec[O2]/(solnvec[A]), 0.0);
    chi_e = boost::math::gamma_p(alpha_e, beta_e*eRatio);
  }
  
  double chi_p = 0.3;

  // calculate q_s
  double F_s = (solnvec[G] + solnvec[Xy])/(solnvec[G] + solnvec[Xy] + K_s);
  double F_e = (solnvec[O2] + solnvec[A]/beta_e)/(solnvec[O2] + solnvec[A]/beta_e + K_e);
  double q_s = qs_max*F_s*F_e;

  double our = -chi_e*y_os*q_s*solnvec[X];

  return our;
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
    double dt=advance_time/50.0;

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
