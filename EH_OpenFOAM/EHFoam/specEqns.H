{
    volScalarField Kdeff("Kdeff",Kd*(1.0-phis));

    // Solve solids equation
    fvScalarMatrix phisEqn
        (
         fvm::ddt(phis) + fvm::div(phi,phis) + fvm::div(fvc::flux(setlvel), phis)
         ==
         fvm::laplacian(Ks,phis)
        );
    phisEqn.relax();
    phisEqn.solve();


    // Solve facile solids equation
    fvScalarMatrix phifsEqn
        (
         fvm::ddt(phifs) + fvm::div(phi,phifs) + fvm::div(fvc::flux(setlvel), phifs)
         ==
         fvm::laplacian(Ks,phifs)
        );
    phifsEqn.relax();
    phifsEqn.solve();
    
    // Solve recalcitrant solids equation
    fvScalarMatrix phirsEqn
        (
         fvm::ddt(phirs) + fvm::div(phi,phirs) + fvm::div(fvc::flux(setlvel), phirs)
         ==
         fvm::laplacian(Ks,phirs)
        );
    phirsEqn.relax();
    phirsEqn.solve();
    
    // Solve lignin equation
    fvScalarMatrix philsEqn
        (
         fvm::ddt(phils) + fvm::div(phi,phils) + fvm::div(fvc::flux(setlvel), phils)
         ==
         fvm::laplacian(Ks,phils)
        );
    philsEqn.relax();
    philsEqn.solve();

    // Solve bound enzyme equation
    fvScalarMatrix cebEqn
        (
         fvm::ddt(ceb) + fvm::div(phi,ceb) + fvm::div(fvc::flux(setlvel), ceb)
         ==
         fvm::laplacian(Kdeff,ceb)
        );

    cebEqn.relax();
    cebEqn.solve();

    // Solve free enzyme equation
    fvScalarMatrix cefEqn
        (
         fvm::ddt(cef) + fvm::div(phi, cef)
         ==
         fvm::laplacian(Kdeff,cef)
        );
    cefEqn.relax();
    cefEqn.solve();

    // Solve glucose equation
    fvScalarMatrix cgEqn
        (
         fvm::ddt(cg) + fvm::div(phi, cg)
         ==
         fvm::laplacian(Kdeff,cg)
        );
    cgEqn.relax();
    cgEqn.solve();
    
    // Solve xylose equation
    fvScalarMatrix cxEqn
        (
         fvm::ddt(cx) + fvm::div(phi, cx)
         ==
         fvm::laplacian(Kdeff,cx)
        );
    cxEqn.relax();
    cxEqn.solve();
    
    // Solve soluble lignin equation
    fvScalarMatrix clEqn
        (
         fvm::ddt(cl) + fvm::div(phi, cl)
         ==
         fvm::laplacian(Kdeff,cl)
        );
    clEqn.relax();
    clEqn.solve();


}
