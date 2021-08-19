#include<vector>
#include<scalar.H>
//#include <boost/math/special_functions/gamma.hpp>
//#include <gamma.hpp>

namespace reactmodel
{
    extern const int FGL ; //facile glucan
    extern const int RGL;  //recalcitrant glucan
    extern const int XLN;  //Xylan 
    extern const int LGN;  //Xylan 
    extern const int GLS;  //glucose
    extern const int XLS;  //Xylose 
    extern const int SLG;  //soluble lignin
    extern const int EZT;  //Enzyme total
    extern const int nspecies;
    
    extern const double MwE;
    extern const double MwG;
    extern const double Mwg;
    extern const double MwX; // mol weight xylan
    extern const double Mwx; // mol weight xylose
    extern const double MwL; // approx mol weight of lignin

    extern const double rhol;
    extern const double rhoT;
    extern const double rGg;
    extern const double rXx;


    // kinetics model parameters
    extern const double KdR;
    extern const double kL;
    extern const double kR;
    extern const double kF;
    extern const double kX;
    extern const double kappaRF;
    extern const double kappaRL;
    extern const double kappaRX;
    extern const double kappaRs;
    extern const double mGR1;
    extern const double mX1;


    void get_adsorption(std::vector<double>& Cenzymes,
        std::vector<double>& solnvec);

    void get_rhs(std::vector<double>& rhs, std::vector<double> solnvec);
    void advance_soln(std::vector<double>& solnvec,double dt);
}
