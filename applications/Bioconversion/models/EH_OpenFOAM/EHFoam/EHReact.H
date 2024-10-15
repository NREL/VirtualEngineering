{
    phixs==phis-phifs-phirs-phils;
    const vectorField& C = mesh.C();
    std::vector<double> solnvec(reactmodel::nspecies);
    std::vector<double> solnvec0(reactmodel::nspecies);
    std::vector<double> cenzymes(3);
    int outsteps=10;
    double solidsvfrac_cutoff=1e-10;
    double reactsolids_cutoff=1e-8;
    double conc_cutoff=1e-10;

    //some known issues
    /*1) when there are no solids and glucose in a cell with a 
    large enzyme concentrations, the reactions tend give some negative
    values for the products and bound enzymes
     
    2) I am assuming the concentrations of dissolved solids are with
    respect to the liquid. ie cg/cx/cl are all in kg per m3 of liquid
    Hence, in reactions i need a liquid volume fraction scaling to get the 
    correct mass fraction. That said, the transport equation needs 
    liquid volume fraction multiplier for true mass conservation. Its ok for 
    now as we are dealing with dilute systems

    3) Somehow the conversion when there is no lignin and xylan is not 
    matching our earlier cellulose model for the well mixed and intermediate 
    mixing case. Need to ask JJS about this.*/
    

    Info << "Starting reaction update..\n";
    std::ofstream os_timehist;
    if(Pstream::master() && singlecellrunflag)
    {
        os_timehist.open("timehist.dat");
    }

    #include "updateAvgVals.H"

    if(Pstream::master())
    {
        //write integrated quantities
        os_intquants<<reaction_time<<"\t";
        for(unsigned int sp=0;sp<avgvals.size();sp++)
        {
            os_intquants<<avgvals[sp]<<"\t";
        }
        os_intquants<<"\n";
    }

    //loop over cells, populate solution vector with average quantities and local O2
    //integrate reactions at each cell to prescribed time
    forAll(C,i)
    {
        double solidsvfrac=std::max(phis[i],solidsvfrac_cutoff);
        double liqvfrac=1.0-solidsvfrac;

        double totalmass=(solidsvfrac*rhos.value()+liqvfrac*rhol.value());
        double boundenz;

        solnvec[reactmodel::FGL]  = std::max(phifs[i],solidsvfrac_cutoff)*rhos.value()/totalmass;
        solnvec[reactmodel::RGL]  = std::max(phirs[i],solidsvfrac_cutoff)*rhos.value()/totalmass;
        solnvec[reactmodel::XLN]  = std::max(phixs[i],solidsvfrac_cutoff)*rhos.value()/totalmass;
        solnvec[reactmodel::LGN]  = std::max(phils[i],solidsvfrac_cutoff)*rhos.value()/totalmass;

        //cg is concentration of glucose in kg/(m3 of liquid)
        solnvec[reactmodel::GLS]  = std::max(cg[i],conc_cutoff)*liqvfrac/totalmass; 
        solnvec[reactmodel::XLS]  = std::max(cx[i],conc_cutoff)*liqvfrac/totalmass;
        solnvec[reactmodel::SLG]  = std::max(cl[i],conc_cutoff)*liqvfrac/totalmass;
        solnvec[reactmodel::EZT]  = (std::max(cef[i],conc_cutoff)+std::max(ceb[i],conc_cutoff))*liqvfrac/totalmass;

        for(int sp=0;sp<reactmodel::nspecies;sp++)
        {
                solnvec0[sp]=solnvec[sp];
        }

        double integ_time=0.0;
        int step=0;
        while(integ_time < reaction_update_time.value())
        {
            if(Pstream::master() && step%outsteps==0 && singlecellrunflag)
            {   
                os_timehist << integ_time <<"\t";
                for(unsigned int sp=0;sp<solnvec.size();sp++)
                {
                    os_timehist << solnvec[sp] <<"\t";
                }
                os_timehist<<"\n";
                os_timehist.flush();
            }
            reactmodel::advance_soln(solnvec,reaction_time_step.value());
            integ_time += reaction_time_step.value();
            step++;
        }

        reactmodel::get_adsorption(cenzymes,solnvec);
        boundenz=cenzymes[0]+cenzymes[1]+cenzymes[2]; //mol/m3

        /*if(boundenz<0.0)
        {
           Info << "\nbound enzymes:"<< cenzymes[0]<<"\t"<<cenzymes[1]<<"\t"<<cenzymes[2]<<"\n";
           Info << "solnvec0:";
           for(int sp=0;sp<reactmodel::nspecies;sp++)
           {
                Info << solnvec0[sp] <<"\t";
           }
           Info<<"\n";
           Info << "solnvec:";
           for(int sp=0;sp<reactmodel::nspecies;sp++)
           {
                Info << solnvec[sp] <<"\t";
           }
           Info<<"\n";

           Info<<"solids vfracs phis,phifs,phirs,phils,phixs:"<<phis[i]<<"\t"<<
               phifs[i]<<"\t"<<phirs[i]<<"\t"<<phils[i]<<"\t"<<phixs[i]<<"\n\n";

        }*/

        phifs[i] = solnvec[reactmodel::FGL]*totalmass/rhos.value(); 
        phirs[i] = solnvec[reactmodel::RGL]*totalmass/rhos.value(); 
        phixs[i] = solnvec[reactmodel::XLN]*totalmass/rhos.value(); 
        phils[i] = solnvec[reactmodel::LGN]*totalmass/rhos.value();

        phis[i]=phifs[i]+phirs[i]+phixs[i]+phils[i];
        double newliqvfrac=1.0-phis[i];


        cg[i]  = solnvec[reactmodel::GLS]*totalmass/newliqvfrac; 
        cx[i]  = solnvec[reactmodel::XLS]*totalmass/newliqvfrac; 
        cl[i]  = solnvec[reactmodel::SLG]*totalmass/newliqvfrac; 
        ceb[i] = boundenz*reactmodel::MwE;
        cef[i] = (solnvec[reactmodel::EZT]*totalmass/newliqvfrac)-ceb[i];
    }

    Info << "finished reaction update..\n";

    #include "updateAvgVals.H"
    /*if(Pstream::master())
    {
        //write integrated quantities
        os_intquants<<reaction_time<<"\t";
        for(unsigned int sp=0;sp<avgvals.size();sp++)
        {
            os_intquants<<avgvals[sp]<<"\t";
        }
        os_intquants<<"\n";
    }*/

    os_intquants.flush();
}
