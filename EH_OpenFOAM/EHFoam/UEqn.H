MRF.correctBoundaryVelocity(U);

Info<<"max visc:"<<max(visc)<<"\n";
Info<<"min visc:"<<min(visc)<<"\n";

    fvVectorMatrix UEqn
    (
     fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
     + MRF.DDt(rho, U)
     - fvc::div(visc*dev2(T(fvc::grad(U))))
     - fvm::laplacian(visc, U)
     + turbulence->divDevRhoReff(rho, U)
     ==
     fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
        (
         UEqn
         ==
         fvc::reconstruct
         (
          (
           - ghf*fvc::snGrad(rho)
           - fvc::snGrad(p_rgh)
          ) * mesh.magSf()
         )
        );

    fvOptions.correct(U);
}
