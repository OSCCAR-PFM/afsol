/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "kineticTheoryModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const Foam::phaseModel& phase1,
    const Foam::phaseModel& phase2,
    const Foam::volVectorField& U2,
    const Foam::volScalarField& k,
    const Foam::volScalarField& alpha1,
    const Foam::dragModel& drag1
)
:
    phase1_(phase1),
    phase2_(phase2),
    U1_(phase1.U()),
    U2_(U2),
    k_(k),
    alpha1_(alpha1),
    phi1_(phase1.phi()),
    drag1_(drag1),

    rho1_(phase1.rho()),
    da_(phase1.d()),
    nu1_(phase1.nu()),
    rho2_(phase2.rho()),
    nu2_(phase2.nu()),
    

    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            U1_.time().constant(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            U1_.time().constant(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    kineticTheory_(kineticTheoryProperties_.lookup("kineticTheory")),
    equilibrium_(kineticTheoryProperties_.lookup("equilibrium")),

    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    conductivityModel_
    (
        conductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            kineticTheoryProperties_
        )
    ),
    granularPressureModel_
    (
        granularPressureModel::New
        (
            kineticTheoryProperties_
        )
    ),
    frictionalStressModel_
    (
        frictionalStressModel::New
        (
            kineticTheoryProperties_
        )
    ),
   // rho2_(transportProperties_.lookup("rho2")),
   // nu2_(transportProperties_.lookup("nu2")),
    e_(kineticTheoryProperties_.lookup("e")),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    alphaMinFriction_(kineticTheoryProperties_.lookup("alphaMinFriction")),
    Fr_(kineticTheoryProperties_.lookup("Fr")),
    eta_(kineticTheoryProperties_.lookup("eta")),
    p_(kineticTheoryProperties_.lookup("p")),
    phi_(dimensionedScalar(kineticTheoryProperties_.lookup("phi"))*M_PI/180.0),
    R_(kineticTheoryProperties_.lookup("R")),
    Vave_(kineticTheoryProperties_.lookup("vave")),
    Theta_
    (
        IOobject
        (
            "Theta",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh()
    ),
    mu1_
    (
        IOobject
        (
            "mu1",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    
    tausg_
    (
        IOobject
        (
            "tausg",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
    ),
    
    ksg_
    (
        IOobject
        (
            "ksg",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
    ),
    pa_
    (
        IOobject
        (
            "pa",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    
    
    ppMagf_
    (
        IOobject
        (
            "ppMagf",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    
    kappa_
    (
        IOobject
        (
            "kappa",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    gs0_
    (
        IOobject
        (
            "gs0",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    gs0Prime_
    (
        IOobject
        (
            "gs0prime",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    
pfH_
    (
        IOobject
        (
            "pfH",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    muf_
    (
        IOobject
        (
            "muf",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModel::solve(const volTensorField& gradU1t)
{
    if (!kineticTheory_)
    {
        return;
    }

    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    surfaceScalarField phi(1.5*rho1_*phi1_*fvc::interpolate(alpha1_));

    volTensorField dU(gradU1t.T());    //fvc::grad(U1_);
    volSymmTensorField D(symm(dU));

    // NB, drag = K*alpha1*alpha2,
    // (the alpha1 and alpha2 has been extracted from the drag function for
    // numerical reasons)
    volScalarField Ur(mag(U1_ - U2_));
    volScalarField alpha2Prim(alpha1_*(1.0 - alpha1_)*drag1_.K(Ur));

    // Calculating the radial distribution function (solid volume fraction is
    //  limited close to the packing limit, but this needs improvements)
    //  The solution is higly unstable close to the packing limit.
    gs0_ = radialModel_->g0
    (
        min(max(alpha1_, scalar(1e-6)), alphaMax_ - 0.01),
        alphaMax_
    );

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha1_,
            gs0_,
            rho1_,
            e_
        )
    );

    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha1_, Theta_, gs0_, rho1_, da_, e_);

    // particle viscosity (Table 3.2, p.47)
    mu1_ = viscosityModel_->mu1(alpha1_, Theta_, gs0_, rho1_, da_, e_);

    dimensionedScalar Tsmall
    (
        "small",
        dimensionSet(0 , 2 ,-2 ,0 , 0, 0, 0),
        1.0e-6
    );

    dimensionedScalar TsmallSqrt = sqrt(Tsmall);
    volScalarField ThetaSqrt(sqrt(Theta_));

    // dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        12.0*(1.0 - sqr(e_))*sqr(alpha1_)*rho1_*gs0_*(1.0/da_)*ThetaSqrt/sqrtPi
    );

    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1(3.0*alpha2Prim);
    volScalarField J2
    (
        0.25*sqr(alpha2Prim)*da_*sqr(Ur)
       /(max(alpha1_, scalar(1e-6))*rho1_*sqrtPi*(ThetaSqrt + TsmallSqrt))
    );
    
   /* dictionary phase2Dictionary
    (
        transportProperties.subDict("phase2")
    );

    dimensionedScalar rho2(phase2Dictionary.lookup("rho"));
    dimensionedScalar nu2(phase2Dictionary.lookup("nu"));*/
    
    
      Info<< "Step 1" << endl;
    dimensionedScalar ST = da_*da_*Vave_*rho1_/(18*rho2_*nu2_*2*R_);
    Info<< "Step 2" << endl;
    volScalarField tauc_ = da_*sqrt(M_PI)/(24*(alpha1_+1.0e-6)*gs0_*ThetaSqrt);
    Info<< "Step 3" << endl;
    volScalarField tauD_ = rho1_*(alpha1_+1.0e-6)/((drag1_.K(Ur));
    Info<< "Step 4" << endl;
    
     tausg_ = tauc_;
    
   if (ST.value()<100)
    
    {
    
      tausg_ = tauD_;
      
     }
     
     else 
     
     {
     
       tausg_ = tauc_;
       
     }
     
   //  Info<<"dist: "<<tausg_<<"\n";
     dimensionedScalar tausgsmall
    (
        "small",
        dimensionSet(0, 0, 1, 0, 0),
        1.0e-6
    );
    
    //Info<<"dist: "<<k_<<"\n";
    ksg_ = sqrt(6*k_*Theta_);
    
    Info<< "Step 5" << endl;
    
    volScalarField IT = rho2_*(ksg_-3.0*Theta_)/(tausg_+tausgsmall);
    
  // Info<<"dist: "<<J2-J1<<"\n";
    Info<< "Step 6" << endl;
    // bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ = (4.0/3.0)*sqr(alpha1_)*rho1_*da_*gs0_*(1.0+e_)*ThetaSqrt/sqrtPi;

    // stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau(2.0*mu1_*D + (lambda_ - (2.0/3.0)*mu1_)*tr(D)*I);

Info<< "Step 7" << endl;
    if (!equilibrium_)
    {
        // construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20
        // no grad infront of Ps
        // wrong sign infront of laplacian
       // Info<<"dist: "<<IT<<"\n";
        fvScalarMatrix ThetaEqn
        (
            fvm::ddt(1.5*alpha1_*rho1_, Theta_)
          + fvm::div(phi, Theta_, "div(phi,Theta)")
         ==
            fvm::SuSp(-((PsCoeff*I) && dU), Theta_)
          + (tau && dU)
          + fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
          + fvm::Sp(-gammaCoeff, Theta_)
        //  + fvm::Sp(-J1, Theta_)
        //  + fvm::Sp(J2/(Theta_ + Tsmall), Theta_)
          + fvm::Sp(IT/((Theta_+ Tsmall)), Theta_)
          
        );
Info<< "Step 8" << endl;
        ThetaEqn.relax();
        ThetaEqn.solve();
    }
    
    else
    {
        // equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1(2.0*(1.0 + e_)*rho1_*gs0_);
        volScalarField K3
        (
            0.5*da_*rho1_*
            (
                (sqrtPi/(3.0*(3.0-e_)))
               *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha1_*gs0_)
               +1.6*alpha1_*gs0_*(1.0 + e_)/sqrtPi
            )
        );

        volScalarField K2
        (
            4.0*da_*rho1_*(1.0 + e_)*alpha1_*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0
        );

        volScalarField K4(12.0*(1.0 - sqr(e_))*rho1_*gs0_/(da_*sqrtPi));

        volScalarField trD(tr(D));
        volScalarField tr2D(sqr(trD));
        volScalarField trD2(tr(D & D));

        volScalarField t1(K1*alpha1_ + rho1_);
        volScalarField l1(-t1*trD);
        volScalarField l2(sqr(t1)*tr2D);
        volScalarField l3
        (
            4.0
           *K4
           *max(alpha1_, scalar(1e-6))
           *(2.0*K3*trD2 + K2*tr2D)
        );

        Theta_ = sqr((l1 + sqrt(l2 + l3))/(2.0*(alpha1_ + 1.0e-4)*K4));
    }
// Info<<"dist: "<<tausg<<"\n";
    Theta_.max(1.0e-15);
    Theta_.min(1.0e+3);

    volScalarField pf
    (
        frictionalStressModel_->frictionalPressure
        (
            alpha1_,
            alphaMinFriction_,
            alphaMax_,
            Fr_,
            eta_,
            p_
           
        )
    );

    pfH_ = pf;

    PsCoeff += pf/(Theta_+Tsmall);

    PsCoeff.min(1.0e+10);
    PsCoeff.max(-1.0e+10);

    // update particle pressure
    pa_ = PsCoeff*Theta_;

    // frictional shear stress, Eq. 3.30, p. 52
    volScalarField muf
    (
        frictionalStressModel_->muf
        (
            alpha1_,
            alphaMax_,
            pf,
            D,
            phi_
            
        )
    );

    // add frictional stress
    mu1_ += muf;
    mu1_.min(1.0e+2);
    mu1_.max(0.0);

    Info<< "kinTheory: max(Theta) = " << max(Theta_).value() << endl;

    volScalarField ktn(mu1_/rho1_);

    Info<< "kinTheory: min(nu1) = " << min(ktn).value()
        << ", max(nu1) = " << max(ktn).value() << endl;

    Info<< "kinTheory: min(pa) = " << min(pa_).value()
        << ", max(pa) = " << max(pa_).value() << endl;
}

volScalarField& Foam::kineticTheoryModel::ppMagf(const volScalarField& alphaUpdate)
{
    volScalarField alpha1 = alphaUpdate;

    gs0_ = radialModel_->g0(min(alpha1, alphaMinFriction_), alphaMax_);
    gs0Prime_ = radialModel_->g0prime(min(alpha1, alphaMinFriction_), alphaMax_);
    
    // Computing ppMagf
    ppMagf_ = Theta_*granularPressureModel_->granularPressureCoeffPrime
    (
        alpha1,
        gs0_,
        gs0Prime_,
        rho1_,
        e_
    );
    
    volScalarField ppMagfFriction = frictionalStressModel_->frictionalPressurePrime
    (
        alpha1,
        alphaMinFriction_,
        alphaMax_,
        Fr_,
        eta_,
        p_
    );
    
    // NOTE: this might not be appropriate if J&J model is used (verify)
    forAll(alpha1, cellI)
    {
        if(alpha1[cellI] >= alphaMinFriction_.value())
        {
         ppMagf_[cellI] = ppMagfFriction[cellI];
        }
    }
    
    ppMagf_.correctBoundaryConditions();
    
    return ppMagf_;
}



// ************************************************************************* //
