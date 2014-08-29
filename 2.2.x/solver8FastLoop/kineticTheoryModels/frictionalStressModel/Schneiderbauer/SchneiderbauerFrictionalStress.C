/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

// Schellander David
// 20.03.2012
// david.schellander@jku.at

#include "SchneiderbauerFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SchneiderbauerFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SchneiderbauerFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SchneiderbauerFrictionalStress::SchneiderbauerFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SchneiderbauerFrictionalStress::~SchneiderbauerFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SchneiderbauerFrictionalStress::
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


    tmp<volScalarField> p_f
    (
        new volScalarField
        (
            IOobject
            (
                "p_f",
                alpha1.mesh().time().timeName(),
                alpha1.mesh()
            ),
            alpha1.mesh(),
            dimensionedScalar("p_f", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );

   // volScalarField& muff = mag(D);
    Info<<"da: "<<da.value()<<"\n";
    Info<<"rho1: "<<rho1.value()<<"\n";

    scalar dNormCell = scalar(0.0);

    scalar b = scalar(0.2); //Forterre and Pouliquen 2008


    forAll (D, celli)
    {

           if(alphaMinFriction.value() - alpha1[celli] <= 0.0)
        {

           /* dNormCell = sqrt(2.0/9.0*(sqr(D[celli].xx()) + sqr(D[celli].yy()) + sqr(D[celli].zz()))
                                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                                  + sqr(D[celli].yz()));*/
            dNormCell = sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                              + sqr(D[celli].yy() - D[celli].zz())
                              + sqr(D[celli].zz() - D[celli].xx()))
                              + sqr(D[celli].xy()) + sqr(D[celli].xz())
                              + sqr(D[celli].yz()));

            dNormCell = max(dNormCell,scalar(1.0e-6));

            p_f()[celli] = 4.0*rho1.value()*b*b*da.value()*da.value()*dNormCell*dNormCell/(max(alphaMax.value()-alpha1[celli],1e-6)*max(alphaMax.value()-alpha1[celli],1e-6));

            if(p_f()[celli] > 10000)
            {
                p_f()[celli] = 10000;

            }
        }
           else
           {
               p_f()[celli]=0.0;
           }
    }

    return p_f;
}


Foam::tmp<Foam::volScalarField> Foam::SchneiderbauerFrictionalStress::
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

    tmp<volScalarField> p_fPrim
    (
        new volScalarField
        (
            IOobject
            (
                "p_fPrim",
                alpha1.mesh().time().timeName(),
                alpha1.mesh()
            ),
            alpha1.mesh(),
            dimensionedScalar("p_fPrim", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );

   // volScalarField& muff = mag(D);


    scalar dNormCell = scalar(0.0);

    scalar b = scalar(0.2); //Forterre and Pouliquen 2008


    forAll (D, celli)
    {

        if(alphaMinFriction.value() - alpha1[celli] <= 0)
        {
            /*dNormCell = sqrt(2.0/9.0*(sqr(D[celli].xx()) + sqr(D[celli].yy()) + sqr(D[celli].zz()))
                                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                                  + sqr(D[celli].yz()));*/

            dNormCell = sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                              + sqr(D[celli].yy() - D[celli].zz())
                              + sqr(D[celli].zz() - D[celli].xx()))
                              + sqr(D[celli].xy()) + sqr(D[celli].xz())
                              + sqr(D[celli].yz()));

            dNormCell = max(dNormCell,scalar(1.0e-6));

            p_fPrim()[celli] = -2.0*4.0*rho1.value()*b*b*da.value()*da.value()*dNormCell*dNormCell/(max(alphaMax.value()-alpha1[celli],1e-6)*max(alphaMax.value()-alpha1[celli],1e-6)*max(alphaMax.value()-alpha1[celli],1e-6));

            if(p_fPrim()[celli] < -10000)
            {
                p_fPrim()[celli] = -10000;

            }
        }
        else
        {
            p_fPrim()[celli]=0.0;
        }
    }

    return p_fPrim;

}


Foam::tmp<Foam::volScalarField> Foam::SchneiderbauerFrictionalStress::muf
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

    // Creating muf assuming it should be 0 on the boundary which may not be
    // true
    tmp<volScalarField> tmuf
    (
        new volScalarField
        (
            IOobject
            (
                "muf",
                alpha1.mesh().time().timeName(),
                alpha1.mesh()
            ),
            alpha1.mesh(),
            dimensionedScalar("muf", dimensionSet(1, -1, -1, 0, 0), 0.0)
        )
    );

   // volScalarField& muff = mag(D);


    scalar dNormCell = scalar(0.0);

    //Calculation of internal friction coefficient (Jop et. al, Forterre and Pouliquen)
    scalar muIMin = scalar(0.382); //tan(20.9°)
    scalar muIMax = scalar(0.644); //tan(32.7°)
    scalar Is = scalar(10.0);
    scalar I0 = scalar(0.279); // value from Schneiderbauer et. al.


    forAll (D, celli)
    {

            /*dNormCell = sqrt(2.0/9.0*(sqr(D[celli].xx()) + sqr(D[celli].yy()) + sqr(D[celli].zz()))
                                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                                  + sqr(D[celli].yz()))+1.0e-15;*/

            dNormCell = sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                              + sqr(D[celli].yy() - D[celli].zz())
                              + sqr(D[celli].zz() - D[celli].xx()))
                              + sqr(D[celli].xy()) + sqr(D[celli].xz())
                              + sqr(D[celli].yz()))+1.0e-15;

            dNormCell = max(dNormCell,scalar(1.0e-6));

            if(pf[celli] > 1e-6)
            {
                Is = 2.0*dNormCell*da.value()/(sqrt(pf[celli]/rho1.value()));
            }
            else
            {
                Is = 10.0;
            }

            tmuf()[celli] = 0.5*pf[celli]*(muIMin+(muIMax-muIMin)/(I0/Is+1.0))/dNormCell;

            if(tmuf()[celli] > 10000)
            {
                tmuf()[celli] = 10000;

            }
        }

    return tmuf;
}


// ************************************************************************* //
