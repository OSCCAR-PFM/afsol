/*---------------------------------------------------------------------------*\
particleThetaJohnsonJacksonFvPatchScalarField

Copyright Information
    Copyright (C) 2008 Alberto Passalacqua 
    Copyright (C) 2008-2010 Juho Peltola
    Copyright (C) 2010-2011 Alberto Passalacqua 
    
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

#include "particleThetaJohnsonJacksonFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(p.size()),
    muS_(p.size()),
    da_(p.size()),
    U1_(p.size())
    //specularityCoefficient_(p.size())
{
    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const particleThetaJohnsonJacksonFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    muS_(ptf.muS_)
   // specularityCoefficient_(ptf.specularityCoefficient_)
{}


particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    muS_(readScalar(dict.lookup("muS")))
   // specularityCoefficient_(readScalar(dict.lookup("specularityCoefficient")))
{
    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        particleThetaJohnsonJacksonFvPatchScalarField::evaluate();
    }
}


particleThetaJohnsonJacksonFvPatchScalarField::particleThetaJohnsonJacksonFvPatchScalarField
(
    const particleThetaJohnsonJacksonFvPatchScalarField& pivpvf,
    const DimensionedField<double, volMesh>& iF
)
:
    mixedFvPatchScalarField(pivpvf, iF),
    restitutionCoefficient_(pivpvf.restitutionCoefficient_),
    muS_(pivpvf.muS_)
    
   // specularityCoefficient_(pivpvf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleThetaJohnsonJacksonFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void particleThetaJohnsonJacksonFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void particleThetaJohnsonJacksonFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if ((restitutionCoefficient_ < 0) || (restitutionCoefficient_ > 1))
    {
	FatalErrorIn
        (
            "particleThetaJohnsonJacksonFvPatchScalarField::"
            "updateCoeffs()"
        )   << "The value of the restitution coefficient has to be between 0 and 1."
            << abort(FatalError);
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

    const fvPatchVectorField& U1 = 
	patch().lookupPatchField<volVectorField, vector>("U1");

    const fvPatchScalarField& alpha1 = 
	patch().lookupPatchField<volScalarField, scalar>("alpha1");

    const fvPatchScalarField& g0 = 
	patch().lookupPatchField<volScalarField, scalar>("gs0");
    
    const fvPatchScalarField& kappa = 
	patch().lookupPatchField<volScalarField, scalar>("kappa");

    const fvPatchScalarField& mu1 = 
	patch().lookupPatchField<volScalarField, scalar>("mu1");
      
    scalarField magU1Patch = mag(U1.patchInternalField());

    scalarField ThetaPatch = max(patchInternalField(), 1.0e-6);
    
    scalarField alpha1Patch = max(alpha1.patchInternalField(), 1.0e-6);
    
    scalarField g0Patch = g0.patchInternalField();
    
    scalarField kappaPatch = max(kappa.patchInternalField(), 1.0e-15);

    scalarField mu1Patch = mu1.patchInternalField()+1.0e-6;


 
    // The mixed BC in OpenFOAM is implemented as
    //
    // valueFraction*(Theta - ThetaRef) + (1-valueFraction)*(grad(Theta) - grad(Theta)_Ref) = 0
    //
    // To find valueFraction, we re-write Johnson and Jackson BC as
    //
    // c*(Theta - ThetaRef) + (grad(Theta) - grad(Theta)_Ref = 0
    //
    // where
    //
    // c = valueFraction/(1 - valueFraction)
    //
    // As a consequence
    //
    // valueFraction = 1/(c+1)
    //
    // We distinguish two cases, according to the value of the restituition
    // coefficient (See below)
    
    //scalarField delta = 1.0/patch().deltaCoeffs();
    
    scalarField zeroAlpha = pos(alpha1Patch - 1.0e-6);

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
  

   /* scalarField specularityCoefficientprim = magU1Patch; //just to initialize

   forAll(r,celli)
   {
     
     if (r <= 4*k/(7*sqrt(6*M_PI*specularityCoefficientprim0)))
    {
	specularityCoefficientprim = -7*sqrt(6*M_PI)*specularityCoefficientprim0*specularityCoefficientprim0*r/(8*k)+specularityCoefficientprim0;
    }

     else
    {
        specularityCoefficientprim = 2*k/(7*r*sqrt(6*M_PI));
    }
//   
    }*/

	scalarField specularityCoefficientprim = neg(r - rprim)*(-7*sqrt(6*M_PI)*specularityCoefficientprim0*specularityCoefficientprim0*r/(8*k)+specularityCoefficientprim0)
      + pos(r - rprim)*(2*k/(7*r*sqrt(6*M_PI)));



      
     scalarField specularityCoefficient = specularityCoefficientprim/(1+da_*sqrt(3.0*ThetaPatch)
     *M_PI*specularityCoefficientprim*rho1.value()*alpha1Patch*g0Patch/(24.0*mu1Patch*alphaMax.value()));


Info<< "specularityCoefficient = "
            << "  Min(specularityCoefficient) = " << min(specularityCoefficient)
            << "  Max(specularityCoefficient) = " << max(specularityCoefficient)
            << endl;

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

    if (restitutionCoefficient_ != 1.0)
    {
      	// If the restitution coefficient is < 1, Johnson and Jackson BC
	// can be written in the form
	//
	// c*(Theta - ThetaRef) + grad(Theta) = 0
	//
	// with
	//
	// c = valueFraction/(1 - valueFraction)
	//
	// grad(Theta)_Ref = 0
	//
	// and
	//
	// valueFraction = c/(c + 1)

        //this->refValue() = 2.0*specularityCoefficient_*sqr(magU1Patch)
        //    /(3.0*(scalar(1) - sqr(restitutionCoefficient_)));

        //Info<<"e: "<<restitutionCoefficient_<<"\n";



        this->refValue() = 2.0*specularityCoefficient*magU1Patch*magU1Patch
                    /(3.0*(scalar(1) - restitutionCoefficient_*restitutionCoefficient_));

	this->refGrad() = 0.0;
	
           scalarField c = M_PI*alpha1Patch*rho1.value()*g0Patch
            *(scalar(1) - restitutionCoefficient_*restitutionCoefficient_)*sqrt(3.0*ThetaPatch)
            /(4.0*kappaPatch*alphaMax.value())*scalar(2.0)/this->patch().deltaCoeffs();
        this->valueFraction() = c/(c + scalar(1));

        //test
        this->refGrad() = -c*(ThetaPatch-2.0*specularityCoefficient*magU1Patch*magU1Patch
                             /(3.0*(scalar(1) - restitutionCoefficient_*restitutionCoefficient_)));


        this->valueFraction() = 0.0;
        this->refValue() = 0.0;  
        
        
        /*	this->refValue() = 2.0*specularityCoefficient*sqr(magU1Patch)
/(3.0*(scalar(1) - sqr(restitutionCoefficient_)));

this->refGrad() = 0.0;

scalarField c = -M_PI*alpha1Patch*rho1.value()*g0Patch
*(scalar(1) - sqr(restitutionCoefficient_))*sqrt(3.0*ThetaPatch)
/(4.0*kappaPatch*alphaMax.value())*scalar(2.0)/this->patch().deltaCoeffs();

this->valueFraction() = c/(c + scalar(1)); */

    }
    else
    {
	// If the restitution coefficient is 1, the BC degenerates in the form
	//
	// grad(Theta) - grad(Theta)_Ref = 0
	//
	// with 
	//
	// ThetaRef = 0
	//
	// and
	//
	// valueFraction = 0
	
/*	this->refValue() = 0.0;

        //this->refGrad() = zeroAlpha*M_PI*specularityCoefficient_*alpha1Patch*rho1.value()
            //*g0Patch*sqrt(3.0*ThetaPatch)*sqr(magU1Patch)/(6.0*alphaMax.value()*kappaPatch);

        this->refGrad() = zeroAlpha*M_PI*specularityCoefficient*alpha1Patch*rho1.value()
            *g0Patch*sqrt(3.0*ThetaPatch)*magU1Patch*magU1Patch/(6.0*alphaMax.value()*kappaPatch);

	this->valueFraction() = 0.0; */
	
	this->refValue() = 0.0;

        this->refGrad() = M_PI*specularityCoefficient*alpha1Patch*rho1.value()
         *g0Patch*sqrt(3.0*ThetaPatch)*sqr(magU1Patch)/(6.0*kappaPatch);

        this->valueFraction() = 0.0;
	
    } 

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void particleThetaJohnsonJacksonFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("muS")
      << muS_ << token::END_STATEMENT << nl;
   // os.writeKeyword("specularityCoefficient")
      //  << specularityCoefficient_ << token::END_STATEMENT << nl;
    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    particleThetaJohnsonJacksonFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
