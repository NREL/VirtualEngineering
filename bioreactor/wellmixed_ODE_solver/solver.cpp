#include<iostream>
#include<vector>
#include<fstream>

const int X=0;
const int O2=1;
const int G=2;
const int Xy=3;
const int A=4;
const int B=5;

void get_rhs(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars)
{
    double lambda=1.0;

    for(int i=0;i<nvars;i++)
    {
        rhs[i] = -lambda*solnvec[i];
    }
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
