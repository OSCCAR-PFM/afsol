//  *** SETUP THE GRAVITY TERMS FOR PERIODIC AND NON-PERIODIC CASE

    Info << "\nReading g and phiMean and creating gPhase1 & gPhase2" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    uniformDimensionedVectorField gPhase1 = g;
    uniformDimensionedVectorField gPhase2 = g;

    dimensionedScalar phiMean = sum(mesh.V()*alpha1) / sum(mesh.V());
    dimensionedScalar rhoMean = phiMean.value()      * rho1
                              + (1.-phiMean.value()) * rho2;

    if (interfacialProperties.found("useDynamicPressure"))
    {
        Info << "\nPressure contains only dynamic pressure \n" << endl;
        Info << "\nCalculating rhoMean and setting gPhase1 and gPhase2 \n" << endl;

        gPhase1  =  (rho1.value()-rhoMean.value()) / rho1.value() * g ;
        gPhase2  =  (rho2.value()-rhoMean.value()) / rho2.value() * g ;
        Info << "rhoMean         : " << rhoMean.value() << endl; 
        Info << "phiMean         : " << phiMean.value() << endl;
        Info << "gPhase1(solid)  : " << gPhase1 << endl;
        Info << "gPhase2(fluid)  : " << gPhase2 << endl;
    }
    else
    {
         Info << "\nPressure includes hydrostatic and dynamic pressure \n" << endl;
    }

// *** SETUP THE CORRECTION FACTOR TO ENFORCE CONSTANT MOMENTUM OF DOMAIN ***


    uniformDimensionedVectorField gAdjust = g; //adjustable gravity-like term to ensure constant solid flow if necessary
    gAdjust = scalar(0.) * gAdjust;            //initialize with zero

