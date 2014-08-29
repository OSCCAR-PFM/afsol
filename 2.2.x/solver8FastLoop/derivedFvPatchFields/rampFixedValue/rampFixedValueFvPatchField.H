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

Class
    Foam::rampFixedValueFvPatchField

Description
    Foam::rampFixedValueFvPatchField

SourceFiles
    rampFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef rampFixedValueFvPatchField_H
#define rampFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rampFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rampFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Reference value
        Field<Type> refValueHigh_;

	//- Reference value
        Field<Type> refValueLow_;

        //- startRamp
        scalar startRamp_;

        //- endRamp
        scalar endRamp_;

        //- Current time index
        label curTimeIndex_;


    // Private member functions

        //- Return current scale
        scalar currentScale() const;


public:

    //- Runtime type information
    TypeName("rampFixedValue");


    // Constructors

        //- Construct from patch and internal field
        rampFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rampFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rampFixedValueFvPatchField
        //  onto a new patch
        rampFixedValueFvPatchField
        (
            const rampFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rampFixedValueFvPatchField
        (
            const rampFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new rampFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rampFixedValueFvPatchField
        (
            const rampFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new rampFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const Field<Type>& refValueHigh() const
            {
                return refValueHigh_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& refValueHigh()
            {
                return refValueHigh_;
            }

	    //- Return the ref value
            const Field<Type>& refValueLow() const
            {
                return refValueLow_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& refValueLow()
            {
                return refValueLow_;
            }

            //- Return startRamp
            scalar startRamp() const
            {
                return startRamp_;
            }

            scalar& startRamp()
            {
                return startRamp_;
            }

            //- Return endRamp
            scalar endRamp() const
            {
                return endRamp_;
            }

            scalar& endRamp()
            {
                return endRamp_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "rampFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
