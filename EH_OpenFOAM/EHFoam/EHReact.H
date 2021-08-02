{
    phixs==phis-phifs-phirs-phils;
    const vectorField& C = mesh.C();
    std::vector<double> solnvec(reactmodel::nspecies);
    std::vector<double> cenzymes(3);
    int outsteps=10;
    double solidsvfrac_cutoff=1e-10;
    double conc_cutoff=1e-10;
    std::vector<double> avgvals(10);

    Info << "Starting reaction update..\n";
    std::ofstream os_timehist;
    if(Pstream::master() && singlecellrunflag)
    {
        os_timehist.open("timehist.dat");
    }

    avgvals[0]=fvc::domainIntegrate(phis).value()/reactorvolume;
    avgvals[1]=fvc::domainIntegrate(phifs).value()/reactorvolume;
    avgvals[2]=fvc::domainIntegrate(phirs).value()/reactorvolume;
    avgvals[3]=fvc::domainIntegrate(phils).value()/reactorvolume;
    avgvals[4]=fvc::domainIntegrate(phixs).value()/reactorvolume;

    avgvals[5]=fvc::domainIntegrate(ceb*(1-phis)).value()/reactorvolume;
    avgvals[6]=fvc::domainIntegrate(cef*(1-phis)).value()/reactorvolume;
    avgvals[7]=fvc::domainIntegrate(cg*(1-phis)).value()/reactorvolume;
    avgvals[8]=fvc::domainIntegrate(cx*(1-phis)).value()/reactorvolume;
    avgvals[9]=fvc::domainIntegrate(cl*(1-phis)).value()/reactorvolume;

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
    /*if(Pstream::master())
    {
        //write integrated quantities
        os_intquants<<"after reaction:"<<
            fvc::domainIntegrate(phis).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(phifs).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(phirs).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(phils).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(phixs).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(ceb*(1-phis)).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(cef*(1-phis)).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(cg*(1-phis)).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(cx*(1-phis)).value()/reactorvolume<<"\t"<<
            fvc::domainIntegrate(cl*(1-phis)).value()/reactorvolume<<"\n";
    }*/
    os_intquants.flush();
}
