{
    volScalarField Kdeff("Kdeff",Kd*(1.0-phis));

        //settling velocity
        volVectorField setlvel
        (
         "setlvel",
         C1*(exp(-C2*phis/phi_max) - exp(-C2))/(1.0-exp(-C2))*pos(phi_max-phis)*g/mag(g)
        );

        //std::cout<<"max settlevel:"<<max(mag(setlvel)).value()<<"\n";
        //std::cout<<"min settlevel:"<<min(mag(setlvel)).value()<<"\n";

    // Solve solids equation
    fvScalarMatrix phisEqn
        (
         fvm::ddt(phis) + fvm::div(phi,phis) + fvm::div(fvc::flux(setlvel), phis)
         //fvm::ddt(phis) + fvm::div(phi,phis) + fvc::div(setlvel*phis)
         //fvm::ddt(phis) + fvm::div(phi,phis) + phis*fvc::div(setlvel) + (setlvel & fvc::grad(phis))
         ==
         fvm::laplacian(Ks,phis)
        );
    phisEqn.relax();
    phisEqn.solve();

     
    // Solve facile solids equation
    fvScalarMatrix phifsEqn
        (
         fvm::ddt(phifs) + fvm::div(phi,phifs) + fvm::div(fvc::flux(setlvel), phifs)
         //fvm::ddt(phifs) + fvm::div(phi,phifs) + fvc::div(setlvel*phifs)
         //fvm::ddt(phifs) + fvm::div(phi,phifs) + phifs*fvc::div(setlvel) + (setlvel & fvc::grad(phifs))
         ==
         fvm::laplacian(Ks,phifs)
        );
    phifsEqn.relax();
    phifsEqn.solve();

    // Solve bound enzyme equation
    fvScalarMatrix cebEqn
        (
         fvm::ddt(ceb) + fvm::div(phi,ceb) + fvm::div(fvc::flux(setlvel), ceb)
         //fvm::ddt(ceb) + fvm::div(phi,ceb) + fvc::div(setlvel*ceb)
         //fvm::ddt(ceb) + fvm::div(phi,ceb) + fvc::div(setlvel*ceb) + ceb*fvc::div(setlvel) + (setlvel & fvc::grad(ceb))
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

}
