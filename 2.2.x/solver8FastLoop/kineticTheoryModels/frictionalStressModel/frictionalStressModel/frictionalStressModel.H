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
    Foam::frictionalStressModel

SourceFiles
    frictionalStressModel.C

\*---------------------------------------------------------------------------*/

#ifndef frictionalStressModel_H
#define frictionalStressModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class frictionalStressModel Declaration
\*---------------------------------------------------------------------------*/

class frictionalStressModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        frictionalStressModel(const frictionalStressModel&);

        //- Disallow default bitwise assignment
        void operator=(const frictionalStressModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("frictionalStressModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        frictionalStressModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        frictionalStressModel(const dictionary& dict);


    // Selectors

        static autoPtr<frictionalStressModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~frictionalStressModel();


    // Member Functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alpha1,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const = 0;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alpha1f,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const = 0;

        virtual tmp<volScalarField> muf
        (
            const volScalarField& alpha1,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volSymmTensorField& D,
            const dimensionedScalar& phi
        ) const = 0;


	virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alpha1,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p,
            const dimensionedScalar& rho1,
            const dimensionedScalar& da,
    	    const volSymmTensorField& D
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alpha1,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& n,
            const dimensionedScalar& p,
            const dimensionedScalar& rho1,
            const dimensionedScalar& da,
    	    const volSymmTensorField& D
        ) const;

        virtual tmp<volScalarField> muf
        (
            const volScalarField& alpha1,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volSymmTensorField& D,
            const dimensionedScalar& phi,
            const dimensionedScalar& rho1,
            const dimensionedScalar& da
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
