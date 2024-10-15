Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field visc\n" << endl;
volScalarField visc
(
    IOobject
    (
        "visc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
rho.oldTime();

Info<< "Reading field total solids\n" << endl;
volScalarField phis
(
    IOobject
    (
        "phis",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field facile cellulose in solids\n" << endl;
volScalarField phifs
(
    IOobject
    (
        "phifs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field recalcitrant cellulose in solids\n" << endl;
volScalarField phirs
(
    IOobject
    (
        "phirs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field lignin in solids\n" << endl;
volScalarField phils
(
    IOobject
    (
        "phils",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "setting xylan in solids\n" << endl;
volScalarField phixs
(
    IOobject
    (
        "phixs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phis-phifs-phirs-phils
);

Info<< "Reading field bound enzymes\n" << endl;
volScalarField ceb
(
    IOobject
    (
        "ceb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field free enzymes\n" << endl;
volScalarField cef
(
    IOobject
    (
        "cef",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field glucose\n" << endl;
volScalarField cg
(
    IOobject
    (
        "cg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field xylose\n" << endl;
volScalarField cx
(
    IOobject
    (
        "cx",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field soluble lignin\n" << endl;
volScalarField cl
(
    IOobject
    (
        "cl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField setlvel
(
    IOobject
    (
        "setlvel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("setvel", dimVelocity, vector(0.0,0.0,0.0))
);

#include "createPhi.H"

IOdictionary EHprops
(
 IOobject
 (
  "EHProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
 )
);

Info<< "Reading liquid viscosity\n" << endl;

dimensionedScalar mu_l
(
 EHprops.lookup("mu_l")
);

Info<< "Reading viscosity exponent\n" << endl;

dimensionedScalar n_mu
(
 EHprops.lookup("n_mu")
);

Info<< "Reading solids fraction scaling for viscosity\n" << endl;

dimensionedScalar phi_inf
(
 EHprops.lookup("phi_inf")
);

Info<< "Reading maximum solids fraction\n" << endl;

dimensionedScalar phi_max
(
 EHprops.lookup("phi_max")
);

Info<< "Reading solids diffusion coefficient \n" << endl;

dimensionedScalar Ks
(
 EHprops.lookup("Ks")
);

Info<< "Reading dissolved species diffusion coefficient \n" << endl;

dimensionedScalar Kd
(
 EHprops.lookup("Kd")
);

Info<< "Reading settling constant C1 \n" << endl;

dimensionedScalar C1
(
 EHprops.lookup("C1")
);

Info<< "Reading settling constant C2 \n" << endl;

dimensionedScalar C2
(
 EHprops.lookup("C2")
);

Info<< "Reading solids density \n" << endl;

dimensionedScalar rhos
(
 EHprops.lookup("rhos")
);

Info<< "Reading liquid density \n" << endl;

dimensionedScalar rhol
(
 EHprops.lookup("rhol")
);

Info<< "reading subcycling parameters \n" << endl;

dimensionedScalar reaction_update_time
(
 EHprops.lookup("reaction_update_time")
);

dimensionedScalar fluid_steadystate_time
(
 EHprops.lookup("fluid_steadystate_time")
);

dimensionedScalar fluid_update_time
(
 EHprops.lookup("fluid_update_time")
);

dimensionedScalar reaction_time_step
(
 EHprops.lookup("reaction_time_step")
);

rho==phis*rhos + (1.0-phis)*rhol;

Info<< "setting fis\n" << endl;
volScalarField fis
(
    IOobject
    (
        "fis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phis*rhos/rho
);

bool solveNavierStokesFlag = EHprops.lookupOrDefault("solveNavierStokes",  true);
bool solveTransportFlag = EHprops.lookupOrDefault("solveTransport",  true);
bool solveChemistryFlag = EHprops.lookupOrDefault("solveChemistry",  true);
bool singlecellrunflag  = EHprops.lookupOrDefault("singlecellrun", false);
bool rescaleflag  = EHprops.lookupOrDefault("rescalefields", false);

singlePhaseTransportModel laminarTransport(U, phi);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

visc==mu_l*pow((1.0-phis/phi_inf),-n_mu);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
dimensionedScalar maxvisc=mu_l*pow((1.0-phi_max/phi_inf),-n_mu);
dimensionedScalar minvisc=mu_l;

setlvel==C1*(exp(-C2*phis/phi_max) - exp(-C2))/(1.0-exp(-C2))*pos(phi_max-phis)*g/mag(g);

volScalarField philiq=(1.0-phis);
double reactorvolume=fvc::domainIntegrate(phis).value()
    +fvc::domainIntegrate(philiq).value();
Info<<"Volume of the domain:"<<reactorvolume<<"\n";

#include "createMRF.H"
#include "createFvOptions.H"

std::vector<double> avgvals(10);

#include "updateAvgVals.H"
