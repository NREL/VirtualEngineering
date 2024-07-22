avgvals[0]=fvc::domainIntegrate(phis).value()/reactorvolume;
avgvals[1]=fvc::domainIntegrate(phifs).value()/reactorvolume;
avgvals[2]=fvc::domainIntegrate(phirs).value()/reactorvolume;
avgvals[3]=fvc::domainIntegrate(phils).value()/reactorvolume;
avgvals[4]=fvc::domainIntegrate(phixs).value()/reactorvolume;

/*avgvals[5]=fvc::domainIntegrate(ceb*(1-phis)).value()/reactorvolume;
avgvals[6]=fvc::domainIntegrate(cef*(1-phis)).value()/reactorvolume;
avgvals[7]=fvc::domainIntegrate(cg*(1-phis)).value()/reactorvolume;
avgvals[8]=fvc::domainIntegrate(cx*(1-phis)).value()/reactorvolume;
avgvals[9]=fvc::domainIntegrate(cl*(1-phis)).value()/reactorvolume;*/

avgvals[5]=fvc::domainIntegrate(ceb).value()/reactorvolume;
avgvals[6]=fvc::domainIntegrate(cef).value()/reactorvolume;
avgvals[7]=fvc::domainIntegrate(cg).value()/reactorvolume;
avgvals[8]=fvc::domainIntegrate(cx).value()/reactorvolume;
avgvals[9]=fvc::domainIntegrate(cl).value()/reactorvolume;
