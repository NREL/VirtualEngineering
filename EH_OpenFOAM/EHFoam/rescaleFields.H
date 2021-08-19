{
    Info<<"Rescaling fields..\n";
    std::vector<double> currentavgvals(10);

    currentavgvals[0]=fvc::domainIntegrate(phis).value()/reactorvolume;
    currentavgvals[1]=fvc::domainIntegrate(phifs).value()/reactorvolume;
    currentavgvals[2]=fvc::domainIntegrate(phirs).value()/reactorvolume;
    currentavgvals[3]=fvc::domainIntegrate(phils).value()/reactorvolume;
    currentavgvals[4]=fvc::domainIntegrate(phixs).value()/reactorvolume;

    /*currentavgvals[5]=fvc::domainIntegrate(ceb*(1-phis)).value()/reactorvolume;
    currentavgvals[6]=fvc::domainIntegrate(cef*(1-phis)).value()/reactorvolume;
    currentavgvals[7]=fvc::domainIntegrate(cg*(1-phis)).value()/reactorvolume;
    currentavgvals[8]=fvc::domainIntegrate(cx*(1-phis)).value()/reactorvolume;
    currentavgvals[9]=fvc::domainIntegrate(cl*(1-phis)).value()/reactorvolume;*/
    
    currentavgvals[5]=fvc::domainIntegrate(ceb).value()/reactorvolume;
    currentavgvals[6]=fvc::domainIntegrate(cef).value()/reactorvolume;
    currentavgvals[7]=fvc::domainIntegrate(cg).value()/reactorvolume;
    currentavgvals[8]=fvc::domainIntegrate(cx).value()/reactorvolume;
    currentavgvals[9]=fvc::domainIntegrate(cl).value()/reactorvolume;

    phis  = phis* (avgvals[0]/currentavgvals[0]);
    phifs = phifs*(avgvals[1]/currentavgvals[1]);
    phirs = phirs*(avgvals[2]/currentavgvals[2]);
    phils = phils*(avgvals[3]/currentavgvals[3]);
    phixs = phixs*(avgvals[4]/currentavgvals[4]);

    ceb = ceb*(avgvals[5]/currentavgvals[5]);
    cef = cef*(avgvals[6]/currentavgvals[6]);
    cg  =  cg*(avgvals[7]/currentavgvals[7]);
    cx  =  cx*(avgvals[8]/currentavgvals[8]);
    cl  =  cl*(avgvals[9]/currentavgvals[9]);
}
