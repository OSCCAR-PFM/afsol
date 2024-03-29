/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(frictionalStressModel, 0);

    defineRunTimeSelectionTable(frictionalStressModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frictionalStressModel::frictionalStressModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::frictionalStressModel::~frictionalStressModel()
{}





Foam::tmp<Foam::volScalarField> Foam::frictionalStressModel::
frictionalPressure
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
) const
{
	return frictionalPressure(alpha1, alphaMinFriction, alphaMax, Fr, eta, p);
}

Foam::tmp<Foam::volScalarField> Foam::frictionalStressModel::
frictionalPressurePrime
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
) const
{
	return frictionalPressurePrime( alpha1, alphaMinFriction, alphaMax, Fr, eta,p);
}

Foam::tmp<Foam::volScalarField> Foam::frictionalStressModel::muf
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const dimensionedScalar& phi,
    const dimensionedScalar& rho1,
    const dimensionedScalar& da
) const
{
	return muf( alpha1, alphaMax, pf, D,phi);
}



// ************************************************************************* //
