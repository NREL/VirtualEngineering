/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "Polynomial.H"
#include<vector>
#include"EHReactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    
    double prvs_update_time=0.0;
    double time=runTime.value();
    double reaction_time=0.0;
    double smallvalue=1e-10;
    dimensionedScalar smallconc("smallconc",dimDensity,scalar(smallvalue));
    dimensionedScalar smallsolids("smallsolids",dimless,scalar(smallvalue));
    
    std::ofstream os_intquants;
    if(Pstream::master())
    {
        os_intquants.open("integrated_quantities.dat");
    }
    
    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        
        fis==phis*rhos/rho;
    
        if(solveChemistryFlag)
        {   
            /*Info<<"time, prvs_update_time, fluid_update_time:"<<time<<"\t"<<
               prvs_update_time<<"\t"<<fluid_update_time<<"\n";*/

            if((time-prvs_update_time) >= fluid_update_time.value() && time>fluid_steadystate_time.value())
            {
                #include"EHReact.H"
                prvs_update_time=time;
                reaction_time += reaction_update_time.value(); 
            }
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        time=runTime.value();

        if(solveTransportFlag)
        {
            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    mesh.update();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & Uf();

                            #include "correctPhi.H"

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }


                setlvel==C1*(exp(-C2*phis/phi_max) - exp(-C2))/(1.0-exp(-C2))
                    *pos(phi_max-phis)*g/mag(g);

                #include "specEqns.H"

                //prevent negative values
                phis  = max(phis,  smallsolids);
                phifs = max(phifs, smallsolids);
                phirs = max(phirs, smallsolids);
                phils = max(phils, smallsolids);
                cg    = max(cg,    smallconc);
                cx    = max(cx,    smallconc);
                cl    = max(cl,    smallconc);
                cef   = max(cef,   smallconc);
                ceb   = max(ceb,   smallconc);

                phixs==phis-phifs-phirs-phils;
                if(rescaleflag)
                {
                    #include "rescaleFields.H"
                }

                visc == mu_l*pow((1.0-phis/phi_inf),-n_mu);
                visc.max(minvisc);
                visc.min(maxvisc);
                rho  == phis*rhos + (1.0-phis)*rhol;
                rhoPhi = fvc::interpolate(rho)*phi; 

                if(solveNavierStokesFlag)
                {
                    #include "UEqn.H"
                    // --- Pressure corrector loop
                    while (pimple.correct())
                    {
                        #include "pEqn.H"
                    }
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    os_intquants.close();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
