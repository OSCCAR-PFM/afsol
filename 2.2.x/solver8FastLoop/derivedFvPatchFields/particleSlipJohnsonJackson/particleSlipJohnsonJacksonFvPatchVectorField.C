/*---------------------------------------------------------------------------*\
particleSlipJohnsonJacksonFvPatchVectorField

Copyright Information
    Copyright (C) 2008 Alberto Passalacqua 
    Copyright (C) 2008-2010 Juho Peltola
    
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

\*---------------------------------------------------------------------------*/

#include "particleSlipJohnsonJacksonFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
   // specularityCoefficient_(p.size()),
    restitutionCoefficient_(p.size()),
    muS_(p.size()),
    da_(p.size()),
    U1_(p.size())
{}


particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const particleSlipJohnsonJacksonFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(tdpvf, p, iF, mapper),
    //specularityCoefficient_(tdpvf.specularityCoefficient_),
    restitutionCoefficient_(tdpvf.restitutionCoefficient_),
    muS_(tdpvf.muS_)
{}


particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
   // specularityCoefficient_(readScalar(dict.lookup("specularityCoefficient"))),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    muS_(readScalar(dict.lookup("muS")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        partialSlipFvPatchVectorField::evaluate();
    }
}


particleSlipJohnsonJacksonFvPatchVectorField::particleSlipJohnsonJacksonFvPatchVectorField
(
    const particleSlipJohnsonJacksonFvPatchVectorField& tdpvf,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(tdpvf, iF),
  //  specularityCoefficient_(tdpvf.specularityCoefficient_),
    restitutionCoefficient_(tdpvf.restitutionCoefficient_),
    muS_(tdpvf.muS_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void particleSlipJohnsonJacksonFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    
    
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    const dictionary& kineticTheoryProperties = db().lookupObject<IOdictionary>
    (
        "kineticTheoryProperties"
    );

    dictionary phase1Dictionary
    (
        transportProperties.subDict("phase1")
    );

    dimensionedScalar rho1(phase1Dictionary.lookup("rho"));
    dimensionedScalar alphaMax(kineticTheoryProperties.lookup("alphaMax"));

    const fvPatchScalarField& alpha1 = 
	patch().lookupPatchField<volScalarField, scalar>("alpha1");

    const fvPatchScalarField& g0 = 
	patch().lookupPatchField<volScalarField, scalar>("gs0");

    const fvPatchScalarField& mu1 = 
	patch().lookupPatchField<volScalarField, scalar>("mu1");

    scalarField alpha1Patch = alpha1.patchInternalField() + 1.0e-6;
    scalarField ThetaPatch = max(alpha1Patch, 1.0e-6);
    scalarField g0Patch = g0.patchInternalField();
    scalarField mu1Patch = mu1.patchInternalField()+1.0e-6;

    const fvPatchVectorField& U1 =
        patch().lookupPatchField<volVectorField, vector>("U1");
    scalarField magU1Patch = mag(U1.patchInternalField()) + 1.0e-6;

    const fvPatchScalarField& muf =
        patch().lookupPatchField<volScalarField, scalar>("muf");

    const fvPatchScalarField& pf =
        patch().lookupPatchField<volScalarField, scalar>("pfH");

    if (db().foundObject<volScalarField>("Theta"))
    {
	const fvPatchScalarField& Theta = 
	  patch().lookupPatchField<volScalarField, scalar>("Theta");
	
	ThetaPatch = Theta.patchInternalField();
    }
    
    // The partial slip BC in OpenFOAM is implemented as
    //
    // valueFraction*U + (1-valueFraction)*grad(U) = 0
    //
    // To find valueFraction, we re-write Johnson and Jackson BC as
    //
    // 1/c*U + grad(U) = 0
    //
    // where
    //
    // c = valueFraction/(1 - valueFraction)
    //
    // As a consequence
    //
    // valueFraction = 1/(c + 1)


    // Revisiting Johnson-Jackson boundary condition

     scalarField mu_ = magU1Patch; //just to initialize

   forAll(magU1Patch,celli)
   {
      
            mu_[celli] = muS_;
      
   }

     scalarField k = scalar(7.0)/scalar(2.0)*mu_*(scalar(1.0)+restitutionCoefficient_);

     scalarField specularityCoefficientprim0 = -0.0012596 + 0.1064551*k - 0.0428147*pow(k,2) + 0.0097594*pow(k,3) - 0.0012508258*pow
     (k,4) + 0.0000836983*pow(k,5) - 0.00000226955*pow(k,6);
     
     scalarField r = magU1Patch/sqrt(3*ThetaPatch);
     scalarField rprim = 4*k/(7*sqrt(6*M_PI*specularityCoefficientprim0));

  /* scalarField specularityCoefficientprim = k; //just to initialize

   forAll(r,celli)
   {
     
     if (r[celli] < rprim[celli])
    {
	 specularityCoefficientprim = -7*sqrt(6*M_PI)*specularityCoefficientprim0*specularityCoefficientprim0*r/(8*k)+specularityCoefficientprim0;
    }

     else
    {
         specularityCoefficientprim = 2*k/(7*r*sqrt(6*M_PI));
    }
   
    }*/


      scalarField specularityCoefficientprim = neg(r - rprim)*(-7*sqrt(6*M_PI)*specularityCoefficientprim0*specularityCoefficientprim0*r/(8*k)+specularityCoefficientprim0)
      + pos(r - rprim)*(2*k/(7*r*sqrt(6*M_PI)));


      
     scalarField specularityCoefficient = specularityCoefficientprim/(1+da_*sqrt(3.0*ThetaPatch)
     *M_PI*specularityCoefficientprim*rho1.value()*alpha1Patch*g0Patch/(24.0*mu1Patch*alphaMax.value()));
     
  forAll(specularityCoefficient,celli)
   {

      if ((specularityCoefficient[celli] < 0) || (specularityCoefficient[celli] > 1))
    {
	FatalErrorIn
        (
            "particleSlipJohnsonJacksonFvPatchScalarField::"
            "updateCoeffs()"
        )   << "The value of the specularity coefficient has to be between 0 and 1."
            << abort(FatalError);
    }
   }
//   scalarField tau = (M_PI*rho1.value()*alpha1Patch*g0Patch*specularityCoefficient*sqrt(3.0*ThetaPatch))/(6.0*alphaMax.value());
//   +specularityCoefficient*pf/magU1Patch;
    
   scalarField tau = (M_PI*rho1.value()*alpha1Patch*g0Patch*specularityCoefficient*sqrt(3.0*ThetaPatch))/(6.0*alphaMax.value());
   
    scalarField c =tau/(mu1Patch)*scalar(2.0)/this->patch().deltaCoeffs();

    this->valueFraction() = c/(c + scalar(1));  
    
   /* scalarField c = (6.0*mu1Patch*alphaMax.value())/
(M_PI*rho1.value()*alpha1Patch*g0Patch*
specularityCoefficient*sqrt(3.0*ThetaPatch))*scalar(2.0)/this->patch().deltaCoeffs();

    this->valueFraction() = scalar(1)/(c + scalar(1));*/
    
 
   Info<<"dist: "<<specularityCoefficient<<"\n";


    partialSlipFvPatchVectorField::updateCoeffs();
}

// Write
void particleSlipJohnsonJacksonFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
 //   os.writeKeyword("specularityCoefficient")
 //      << specularityCoefficient_ << token::END_STATEMENT << nl;
  os.writeKeyword("restitutionCoefficient")
      << restitutionCoefficient_ << token::END_STATEMENT << nl;
  os.writeKeyword("muS")
      << muS_ << token::END_STATEMENT << nl;
    writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    particleSlipJohnsonJacksonFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
