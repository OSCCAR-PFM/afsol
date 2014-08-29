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

#include "rampFixedValueFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar rampFixedValueFvPatchField<Type>::currentScale() const
{
    return
       min(1.0,max((this->db().time().value()-startRamp_)/(endRamp_-startRamp_),0.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
rampFixedValueFvPatchField<Type>::rampFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    refValueHigh_(p.size()),
    refValueLow_(p.size()),
    startRamp_(0.0),
    endRamp_(0.0),
    curTimeIndex_(-1)
{}


template<class Type>
rampFixedValueFvPatchField<Type>::rampFixedValueFvPatchField
(
    const rampFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    refValueHigh_(ptf.refValueHigh_, mapper),
    refValueLow_(ptf.refValueLow_, mapper),
    startRamp_(ptf.startRamp_),
    endRamp_(ptf.endRamp_),
    curTimeIndex_(-1)
{}


template<class Type>
rampFixedValueFvPatchField<Type>::rampFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    refValueHigh_("refValueHigh", dict, p.size()),
    refValueLow_("refValueLow", dict, p.size()),
    startRamp_(readScalar(dict.lookup("startRamp"))),
    endRamp_(readScalar(dict.lookup("endRamp"))),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(refValueLow_+(refValueHigh_-refValueLow_)*currentScale());
    }
}


template<class Type>
rampFixedValueFvPatchField<Type>::rampFixedValueFvPatchField
(
    const rampFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    refValueHigh_(ptf.refValueHigh_),
    refValueLow_(ptf.refValueLow_),
    startRamp_(ptf.startRamp_),
    endRamp_(ptf.endRamp_),
    curTimeIndex_(-1)
{}


template<class Type>
rampFixedValueFvPatchField<Type>::rampFixedValueFvPatchField
(
    const rampFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    refValueHigh_(ptf.refValueHigh_),
    refValueLow_(ptf.refValueLow_),
    startRamp_(ptf.startRamp_),
    endRamp_(ptf.endRamp_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void rampFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    refValueHigh_.autoMap(m);
    refValueLow_.autoMap(m);
}


template<class Type>
void rampFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const rampFixedValueFvPatchField<Type>& tiptf =
        refCast<const rampFixedValueFvPatchField<Type> >(ptf);

    refValueHigh_.rmap(tiptf.refValueHigh_, addr);
    refValueLow_.rmap(tiptf.refValueLow_, addr);
}


template<class Type>
void rampFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        patchField = refValueLow_+(refValueHigh_-refValueLow_)*currentScale();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void rampFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    refValueHigh_.writeEntry("refValueHigh", os);
    refValueLow_.writeEntry("refValueLow", os);
    os.writeKeyword("startRamp")
        << startRamp_ << token::END_STATEMENT << nl;
    os.writeKeyword("endRamp")
        << endRamp_ << token::END_STATEMENT << nl;
    this->writeEntry("value",os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
