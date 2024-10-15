#include"EHReactModel.H"

namespace reactmodel
{
    const int FGL=0;
    const int RGL=1;
    const int XLN=2;
    const int LGN=3;
    const int GLS=4;
    const int XLS=5;
    const int SLG=6;
    const int EZT=7;
    const int nspecies=8;

    const double MwE = 65000.0;
    const double MwG = 162.0;
    const double Mwg = 180.0;
    const double MwX = 166.129; // mol weight xylan
    const double Mwx = 150.13; // mol weight xylose
    const double MwL = 200.; // approx mol weight of lignin

    const double rhol = 1000.0;
    const double rhoT = rhol;
    const double rGg = MwG/Mwg;
    const double rXx = MwX/Mwx;

    // kinetics model parameters
    const double KdR = 0.05000000000000001;
    const double kL = 729.4513412205487;
    const double kR = 14712.525849904296;
    const double kF = 14712.525849904296;
    const double kX = 9999.999988133452;
    const double kappaRF = 9.33804072835234;
    const double kappaRL = 49.99999999999999;
    const double kappaRX = 11.281636078910811;
    const double kappaRs = 49.99999999971725;
    const double mGR1 = 1.;
    const double mX1 = 1.;

    void get_adsorption(std::vector<double>& Cenzymes,
            std::vector<double>& solnvec) 
    {
        /*
           Determine the equilibrium adsorption of the enzymes
           */

        // unpack
        double fGF = solnvec[0];
        double fGR = solnvec[1];
        double fX = solnvec[2];
        double fL = solnvec[3];
        double fg = solnvec[4];
        double fx = solnvec[5];
        double fsL = solnvec[6];
        double fET = solnvec[7];


        // Calculate Accessibility
        double fRis = fGR + fX + fL;
        double X_RGR = fGR/fRis;
        // this part is not in the paper -- did we neglect it? I think so,
        // essentially setting MGR1 = mX1 = 1; JJS 7/25/21
        double X_RGRA = 0.5*(3 - mGR1)*X_RGR + 0.5*(mGR1 - 1)*pow(X_RGR,3);
        double X_RX = fX/fRis;
        double X_RXA = 0.5*(3 - mX1)*X_RX + 0.5*(mX1 - 1)*pow(X_RX,3);

        double CGF = (rhoT/MwG) * fGF;
        double CGRA = (rhoT/MwG) * fRis * X_RGRA;
        double CXA = (rhoT/MwX) * fRis * X_RXA;

        double lamFR = CGF/CGRA;
        double lamXR = CXA/CGRA;

        // calculate denominator
        double fis = fRis + fGF;
        double eps_l = rhoT/rhol * (1 - fis);
        double csL = rhoT/MwL * fsL;
        double css = rhoT/Mwg*fg + rhoT/Mwx*fx;

        double denom = 1 + kappaRF*lamFR + kappaRX*lamXR + eps_l/CGRA * (KdR + kappaRL*csL + kappaRs*css);
        double CET = (rhoT/MwE) * fET;

        double CEGR = CET/denom;
        double CEGF = kappaRF * lamFR * CET/denom;
        double CEX = kappaRX * lamXR * CET/denom;

        Cenzymes[0] = CEGR;
        Cenzymes[1] = CEGF;
        Cenzymes[2] = CEX;
    }


    void get_rhs(std::vector<double>& rhs, std::vector<double> solnvec) 
    {
        /*
           RHS of system of ODES
           */

        // unpack
        double fGF = solnvec[0];
        double fGR = solnvec[1];
        double fX = solnvec[2];
        double fL = solnvec[3];
        double fg = solnvec[4];
        double fx = solnvec[5];
        double fsL = solnvec[6];
        double fET = solnvec[7];

        // enzyme adsorption 
        std::vector<double> Cenzymes(3);
        get_adsorption(Cenzymes, solnvec);
        double CEGR = Cenzymes[0];
        double CEGF = Cenzymes[1];
        double CEX = Cenzymes[2];

        // moles of lignin
        double CL = (rhoT/MwE) * fL;

        // calculate molar rates (rt = "r-tilda")
        double rtGF = kF*CEGF;
        double rtGR = kR*CEGR;
        double rtX = kX*CEX;
        double rtL = kL*CL*(rtGR + rtX);

        // calculate rates
        double R_GF = - MwG/rhoT * rtGF;
        double R_GR = - MwG/rhoT * rtGR;
        double R_X  = - MwX/rhoT * rtX;
        double R_g = Mwg/rhoT * (rtGF + rtGR);
        double R_x = Mwx/rhoT * rtX;
        double R_L = - MwL/rhoT * rtL;
        double R_sL = - R_L;
        double R_ET = 0.; // reaction term for enzymes is zero

        rhs[0] = R_GF;
        rhs[1] = R_GR;
        rhs[2] = R_X;
        rhs[3] = R_L;
        rhs[4] = R_g;
        rhs[5] = R_x;
        rhs[6] = R_sL;
        rhs[7] = R_ET;
    }


    void advance_soln(std::vector<double>& solnvec, double dt) 
    {

        std::vector<double> solnvec_0(nspecies);
        std::vector<double> k1(nspecies);
        std::vector<double> k2(nspecies);
        std::vector<double> k3(nspecies);
        std::vector<double> k4(nspecies);

        // Store the initial value
        solnvec_0 = solnvec;

        // Calculate k1: the slope at the beginning of the interval
        get_rhs(k1, solnvec_0);

        // Calculate k2: the slope at the midpoint of the interval following k1    
        for(int i=0; i<nspecies; i++) {
            solnvec[i] = solnvec_0[i] + 0.5*dt*k1[i];     
        }

        get_rhs(k2, solnvec);

        // Calculate k3: the slope at the midpoint of the interval following k2
        for(int i=0; i<nspecies; i++) {
            solnvec[i] = solnvec_0[i] + 0.5*dt*k2[i];            
        }

        get_rhs(k3, solnvec);

        // Calculate k4: the slope at the end of the interval following k3
        for(int i=0; i<nspecies; i++) {
            solnvec[i] = solnvec_0[i] + dt*k3[i];
        }

        get_rhs(k4, solnvec);

        // Update the solution
        for(int i=0; i<nspecies; i++) {
            solnvec[i] = solnvec_0[i] + 1.0/6.0*dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
    }

}
