/*---------------------------------------------------------------------------*\
particleThetaSchneiderbauerSchellanderFvPatchScalarField

Copyright Information
    Copyright (C) 2012 David Schellander and Simon Schneiderbauer
    
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::particleThetaSchneiderbauerSchellanderFvPatchScalarField

SourceFiles
    particleThetaSchneiderbauerSchellanderFvPatchScalarField.C
    
Description
    Schneiderbauer & Schellander boundary condition for the granular temperature.
    
Authors
    Simon Schneiderbauer <simon.schneiderbauer@jku.at>
    David Schellander <david.schellander@jku.at>

Notes
    The code is known to work with OpenFOAM 1.7.x

\*---------------------------------------------------------------------------*/

#ifndef particleThetaSchneiderbauerSchellanderFvPatchScalarField_H
#define particleThetaSchneiderbauerSchellanderFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class particleThetaSchneiderbauerSchellanderFvPatch Declaration
\*---------------------------------------------------------------------------*/

class particleThetaSchneiderbauerSchellanderFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data
    
       	scalar restitutionCoefficient_;	// normal restitution coefficent
	scalar alpha2_;	// tangential restitution coefficient
	scalar muF_;	// wall friction for very small values
	scalar muS_;	// wall friction for sliding
    scalar sigma_;

public:

    //- Runtime type information
    TypeName("particleThetaSchneiderbauerSchellander");


    // Constructors

        //- Construct from patch and internal field
        particleThetaSchneiderbauerSchellanderFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<double, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        particleThetaSchneiderbauerSchellanderFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<double, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  particleThetaSchneiderbauerSchellanderFvPatchScalarField
        //  onto a new patch
        particleThetaSchneiderbauerSchellanderFvPatchScalarField
        (
            const particleThetaSchneiderbauerSchellanderFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<double, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new particleThetaSchneiderbauerSchellanderFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        particleThetaSchneiderbauerSchellanderFvPatchScalarField
        (
            const particleThetaSchneiderbauerSchellanderFvPatchScalarField&,
            const DimensionedField<double, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const DimensionedField<double, volMesh>& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new particleThetaSchneiderbauerSchellanderFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
