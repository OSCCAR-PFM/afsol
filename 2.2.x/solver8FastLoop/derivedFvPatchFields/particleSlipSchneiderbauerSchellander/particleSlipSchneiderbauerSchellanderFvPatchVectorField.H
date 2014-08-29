/*---------------------------------------------------------------------------*\
particleSlipSchneiderbauerSchellanderFvPatchVectorField

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
    Foam::particleSlipSchneiderbauerSchellanderFvPatchVectorField

SourceFiles
    particleSlipSchneiderbauerSchellanderFvPatchVectorField.C
    
Description
    Schneiderbauer & Schellander boundary condition ffor the velocity of the particle
    
Authors
    Simon Schneiderbauer <simon.schneiderbauer@jku.at>
    David Schellander <david.schellander@jku.at>

Notes
    The code is known to work with OpenFOAM 1.7.x

\*---------------------------------------------------------------------------*/

#ifndef particleSlipSchneiderbauerSchellanderFvPatchVectorField_H
#define particleSlipSchneiderbauerSchellanderFvPatchVectorField_H

#include "fvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class particleSlipSchneiderbauerSchellanderFvPatch Declaration
\*---------------------------------------------------------------------------*/


class particleSlipSchneiderbauerSchellanderFvPatchVectorField
:
    public partialSlipFvPatchVectorField
{
    // Private data
	scalar restitutionCoefficient_;	// normal restitution coefficent
	scalar alpha2_;	// tangential restitution coefficient
	scalar muF_;	// wall friction for very small velocities (sticky paricles)
	scalar muS_;	// wall friction sliding
    scalar sigma_;

public:

    //- Runtime type information
    TypeName("particleSlipSchneiderbauerSchellander");


    // Constructors

        //- Construct from patch and internal field
        particleSlipSchneiderbauerSchellanderFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        particleSlipSchneiderbauerSchellanderFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  particleSlipSchneiderbauerSchellanderFvPatchVectorField onto a new patch
        particleSlipSchneiderbauerSchellanderFvPatchVectorField
        (
            const particleSlipSchneiderbauerSchellanderFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new particleSlipSchneiderbauerSchellanderFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        particleSlipSchneiderbauerSchellanderFvPatchVectorField
        (
            const particleSlipSchneiderbauerSchellanderFvPatchVectorField&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<Vector<double>, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new particleSlipSchneiderbauerSchellanderFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
