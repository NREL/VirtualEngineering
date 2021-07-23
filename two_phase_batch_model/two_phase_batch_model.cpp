#include<iostream>
#include<vector>
#include<fstream>


// parameters that are "fixed"
const double MwE = 65000.0;
const double MwG = 162.0;
const double Mwg = 180.0;
const double rhol = 1000.0;
const double rhoT = rhol;
const double rGg = MwG/Mwg;

// kinetics model parameters
const double kF = 5000.0;
const double kR = 1.0*kF;
const double KdF = 0.03773;
const double KdR = 0.5079;
const double KI = 0.01;


void get_adsorption(std::vector<double>& adsorption, double CGF, double CGRA, double CET, double cg, double epsl) {
    /*
    Adsorption equilibrium of enzyme between facile and recalcitrant glucan,
    and accounting for inhibition by glucose (and other sugars if present)
    */
    
    double CEGF = CET/(1.0 + KdF/KdR*CGRA/CGF + epsl*KdF/CGF*(1.0 + cg/KI));
    double CEGR = CET/(1.0 + KdR/KdF*CGF/CGRA + epsl*KdR/CGRA*(1.0 + cg/KI));
    
    adsorption[0] = CEGF;
    adsorption[1] = CEGR;
}


void get_rhs(std::vector<double>& rhs, std::vector<double> solnvec, double fET, double fX0, double fL0) {
    /*
    RHS of system of ODES
    */

    // unpack
    double fGF = solnvec[0];
    double fGR = solnvec[1];
    double fg = solnvec[2];
    
    // convenience term
    double mwGrho = MwG/rhoT;
    
    // molar concentrations
    double CGF = fGF/mwGrho;
    double CGR = fGR/mwGrho;
    double CET = rhoT/MwE*fET;
    
    // calculate liquid molar sugar concentration
    double fG = fGF + fGR;
    double fis = fG + fX0 + fL0;
    double fliq = 1.0 - fis;
    double epsl = rhoT/rhol*fliq;
    double cg = rhol/fliq/Mwg*fg;

    // enzyme adsorption -- with CGRA = CGR (no limited exposure)
    std::vector<double> adsorption(2);
    get_adsorption(adsorption, CGF, CGR, CET, cg, epsl);
    double rF = kF*adsorption[0];
    double rR = kR*adsorption[1];

    // ODEs
    double dfGF = -mwGrho*rF;
    double dfGR = -mwGrho*rR;
    double dfg = Mwg/rhoT*(rF + rR);
    
    rhs[0] = dfGF;
    rhs[1] = dfGR;
    rhs[2] = dfg;
}
    

void advance_soln(std::vector<double>& solnvec, double fET, double fX0, double fL0, double dt, int nvars) {

    std::vector<double> solnvec_0(nvars);
    std::vector<double> k1(nvars);
    std::vector<double> k2(nvars);
    std::vector<double> k3(nvars);
    std::vector<double> k4(nvars);

    // Store the initial value
    solnvec_0 = solnvec;

    // Calculate k1: the slope at the beginning of the interval
    get_rhs(k1, solnvec_0, fET, fX0, fL0);

    // Calculate k2: the slope at the midpoint of the interval following k1    
    for(int i=0; i<nvars; i++) {
        solnvec[i] = solnvec_0[i] + 0.5*dt*k1[i];            
    }

    get_rhs(k2, solnvec, fET, fX0, fL0);

    // Calculate k3: the slope at the midpoint of the interval following k2
    for(int i=0; i<nvars; i++) {
        solnvec[i] = solnvec_0[i] + 0.5*dt*k2[i];            
    }

    get_rhs(k3, solnvec, fET, fX0, fL0);

    // Calculate k4: the slope at the end of the interval following k3
    for(int i=0; i<nvars; i++) {
        solnvec[i] = solnvec_0[i] + dt*k3[i];            
    }

    get_rhs(k4, solnvec, fET, fX0, fL0);

    // Update the solution
    for(int i=0; i<nvars; i++) {
        solnvec[i] = solnvec_0[i] + 1.0/6.0*dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);          
    }
}


int main() {

    // Declare time-stepping parameters
    double t_final = 100.0;
    double dt = 0.1;
    double tt = 0.0;
    int t_steps = (int)(t_final/dt);

    // These values should be obtained from the notebook
    double lmbde = 0.03;
    double fis0 = 0.05;
    double dilution_factor = 0.5;
    double xG0 = 1.0;
    double xX0 = 0.0;
    double rhog0 = 0.0;
    double rhox0 = 1.0*dilution_factor;
    double rhof0 = 0.0*dilution_factor;

    double conversion_xylan = 0.5;
    double yF0 = 0.2 + 0.6*conversion_xylan;

    // initial conditions in model variables
    double fG0 = xG0*fis0;
    double fGF0 = yF0*fG0;
    double fGR0 = (1-yF0)*fG0;
    double fX0 = xX0*fis0;
    double fL0 = (1 - xG0 - xX0)*fis0;
    double fliq0 = 1 - fis0;
    double fg0 = rhog0*fliq0/rhol;
    double fET = lmbde*fG0;

    // Initialize all vectors
    int nvars = 3;
    std::vector<double> solnvec(nvars);

    // Store initial conditions
    solnvec[0] = fGF0;
    solnvec[1] = fGR0;
    solnvec[2] = fg0;

    // Initialize output file
    std::ofstream outfile("eh_results.csv");
    outfile.precision(9);

    // Write header line
    outfile<<"# time, fGF, fGR, fg, rhog\n";  

    // Save initial condition
    outfile<<0.0<<", ";
    for(int i=0; i<nvars; i++) {
        outfile<<solnvec[i]<<", ";
    }
    outfile<<0.0<<"\n";


    double fGF, fGR, fg;
    double fG, fis, fliq, rhog; 
    double convF, convR, conv, mG0, mG, mb;

    for(int i=0; i<t_steps; i++) {
        advance_soln(solnvec, fET, fX0, fL0, dt, nvars);
        tt += dt;

        // unpack and convert output quantities (maybe not all necessary in solver)
        fGF = solnvec[0];
        fGR = solnvec[1];
        fg = solnvec[2];

        fG = fGF + fGR;
        fis = fG + fX0 + fL0;
        fliq = 1.0 - fis;
        rhog = rhol/fliq*fg;

        convF = (fGF0 - fGF)/fG0;
        convR = (fGR0 - fGR)/fG0;
        conv = 1.0 - fG/fG0;
        mG0 = fG0 + fg0*rGg;
        mG = fG + fg*rGg;
        mb = 1.0 - mG/mG0;

        // Write data for this timestep
        outfile<<tt<<", ";
        for(int i=0; i<nvars; i++) {
            outfile<<solnvec[i]<<", ";
        }
        outfile<<rhog<<"\n";
    }


    outfile.close();

    return 0;
}
