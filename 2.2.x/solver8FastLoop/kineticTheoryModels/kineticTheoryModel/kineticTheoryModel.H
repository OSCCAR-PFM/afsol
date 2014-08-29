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

Class
    Foam::kineticTheoryModel

Description

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
{
    // Private data

        const phaseModel& phase1_;
        const phaseModel& phase2_;
        const volVectorField& U1_;
        const volVectorField& U2_;
        const volScalarField& k_;
        const volScalarField& alpha1_;
        const surfaceScalarField& phi1_;

        const dragModel& drag1_;

        const dimensionedScalar& rho1_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nu1_;
        const dimensionedScalar& rho2_;
        const dimensionedScalar& nu2_;

        //- dictionary holding the modeling info
        IOdictionary kineticTheoryProperties_;
        IOdictionary transportProperties_;

        //- use kinetic theory or not.
        Switch kineticTheory_;

        //- use generation == dissipation
        Switch equilibrium_;

        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        autoPtr<conductivityModel> conductivityModel_;

        autoPtr<kineticTheoryModels::radialModel> radialModel_;

        autoPtr<granularPressureModel> granularPressureModel_;

        autoPtr<frictionalStressModel> frictionalStressModel_;

         //- coefficient of restitution
      //  const dimensionedScalar rho2_;
        
        //- coefficient of restitution
       // const dimensionedScalar nu2_;
                      
        //- coefficient of restitution
        const dimensionedScalar e_;

        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- min value for which the frictional stresses are zero
        const dimensionedScalar alphaMinFriction_;

        //- material constant for frictional normal stress
        const dimensionedScalar Fr_;

        //- material constant for frictional normal stress
        const dimensionedScalar eta_;

        //- material constant for frictional normal stress
        const dimensionedScalar p_;

        //- angle of internal friction
        const dimensionedScalar phi_;
        
         //- coefficient of restitution
        const dimensionedScalar R_;
        
        //- coefficient of restitution
        const dimensionedScalar Vave_;

        //- The granular energy/temperature
        volScalarField Theta_;

        //- The granular viscosity
        volScalarField mu1_;

        //- The granular bulk viscosity
        volScalarField lambda_;
        volScalarField tausg_;
        volScalarField ksg_;

        //- The granular pressure
        volScalarField pa_;
        
        volScalarField ppMagf_;


        //- The granular temperature conductivity
        volScalarField kappa_;

        //- The radial distribution function
        volScalarField gs0_;
        volScalarField gs0Prime_;

	volScalarField pfH_;
        volScalarField muf_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


public:

    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volVectorField& U2,
            const volScalarField& alpha1,
            const volScalarField& k,
            const dragModel& drag1
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        void solve(const volTensorField& gradU1t);

        bool on() const
        {
            return kineticTheory_;
        }

        const volScalarField& mu1() const
        {
            return mu1_;
        }

        const volScalarField& pa() const
        {
            return pa_;
        }
        
        volScalarField& ppMagf(const volScalarField& alphaUpdate);

        const volScalarField& lambda() const
        {
            return lambda_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }
        
        const volScalarField& tausg() const
        {
            return tausg_;
        }
        
        const volScalarField& ksg() const
        {
            return ksg_;
        }
        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
