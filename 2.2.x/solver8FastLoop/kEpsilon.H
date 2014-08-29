if (turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }

    tmp<volTensorField> tgradU2 = fvc::grad(U2);
    volScalarField G(2*nut2*(tgradU2() && dev(symm(tgradU2()))));
    tgradU2.clear();

    #include "wallFunctions.H"
    #include "kineticTheoryModel.H"

    // const dimensionedScalar Vave;
    //const dimensionedScalar R;
    //const dimensionedScalar da;
    // volScalarField Theta_;
    //scalar da_ = da.value();
    // const phaseModel& phase1,

    dictionary phase1Dictionary
    (
        transportProperties.subDict("phase1")
    );
    
    dimensionedScalar da(phase1Dictionary.lookup("d"));

    
    dimensionedScalar tausgSmall
    (
        "small",
        dimTime,
        SMALL
    );


    volScalarField Ur(mag(U1 - U2)); 

    volScalarField Rep = da*Ur/nu2; 

    const volScalarField ones
    (
        IOobject
        (
            "ones",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        scalar(1)
    );

    tmp<volScalarField> Ew
    (
        new volScalarField
        (
            IOobject
            (
                "Ew",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimMass/dimLength/pow(dimTime, 3.0),
                0.0
            )
        )
    );

    tmp <volScalarField> Cw =
        (
            pos(sign(Rep - 150.0*ones))*(1.0 - pos(sign(Rep - 310.0*ones)))
        )*(10.0/3.0)
      + pos(sign(Rep - 310.0*ones))*8.0;

    // Must assign dimension of mug to modelling constant "1.2" for
    // 310 <= Rep < 610 in Table 7 of Rao et al. (2011), AIChE Journal Vol 0,
    // Nr. 0, after Lun, C. (2000), "Numerical simulation of dilute turbulent
    // gas-solid flows". Int J Multiphase Flow 26, 1707-1736. Hence, constant
    // is named "CLun", after Lun.
    dimensionedScalar CLun
    (
        "CLun",
        dimMass/dimLength/dimTime,
        scalar(1.2)
    );

    tmp <volScalarField> mut =
        (
            pos(sign(Rep - 150.0*ones))*(1.0 - pos(sign(Rep - 310.0*ones)))
        )*0.017*Rep*rho2*nu2
      + (
            pos(sign(Rep - 310.0*ones))*(1.0 - pos(sign(Rep - 610.0*ones)))
        )*(CLun + pow(Rep, 2.0)*rho2*nu2)
      + pos(sign(Rep - 610.0*ones))*0.029*Rep*rho2*nu2;

    Ew = 12.0*Cw*mut*k/pow(da, 2.0);

    volScalarField Ik =
      - (2.0*k - kineticTheory.ksg())/(kineticTheory.tausg() + tausgSmall)
      + (Ew/rho2);
     

    // Dissipation equation
    
    //Info<<"dist: "<<scalar(0.2)*Ik<<"\n";
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon)
      + fvm::div(phi2, epsilon)
      - fvm::Sp(fvc::div(phi2), epsilon)
      - fvm::laplacian
        (
            alpha1Eps*nuEff2, epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
      ==
        (C1*G*epsilon/k)
      - fvm::Sp(C2*epsilon/k, epsilon)
      + fvm::Sp(scalar(2.5)*Ik/(k*(1-alpha1)), epsilon) 
    );

    // Info<< "Step 16" << endl;
    #include "wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), 1.0e-15));

    // Info<< "Step 17" << endl;

    //  Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(k)
      + fvm::div(phi2, k)
      - fvm::Sp(fvc::div(phi2), k)
      - fvm::laplacian
        (
            alpha1k*nuEff2, k,
            "laplacian(DkEff,k)"
        )
      ==
        G
      - fvm::Sp(epsilon/k, k)
      + fvm::Sp(Ik*alpha1/(k*(1 - alpha1)), k)
    );

    kEqn.relax();
    kEqn.solve();

    // Info<< "Step 18" << endl;
    k.max(dimensionedScalar("zero", k.dimensions(), 1.0e-8));

    //- Re-calculate turbulence viscosity
    nut2 = Cmu*sqr(k)/epsilon;

    #include "wallViscosity.H"
}

nuEff2 = nut2 + nu2;
