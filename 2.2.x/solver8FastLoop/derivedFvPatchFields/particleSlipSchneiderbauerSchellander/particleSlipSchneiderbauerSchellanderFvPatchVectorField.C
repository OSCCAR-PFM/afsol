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

    
Description
    Schneiderbauer & Schellander boundary condition ffor the velocity of the particle
    
Authors
    Simon Schneiderbauer <simon.schneiderbauer@jku.at>
    David Schellander <david.schellander@jku.at, www.schellcom.at>

\*---------------------------------------------------------------------------*/

#include "particleSlipSchneiderbauerSchellanderFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleSlipSchneiderbauerSchellanderFvPatchVectorField::particleSlipSchneiderbauerSchellanderFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    restitutionCoefficient_(p.size()),
    alpha2_(p.size()),
    muF_(p.size()),
    sigma_(p.size())
 //   muS_(p.size())
{}


particleSlipSchneiderbauerSchellanderFvPatchVectorField::particleSlipSchneiderbauerSchellanderFvPatchVectorField
(
    const particleSlipSchneiderbauerSchellanderFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(tdpvf, p, iF, mapper),
    restitutionCoefficient_(tdpvf.restitutionCoefficient_),
    alpha2_(tdpvf.alpha2_),
    muF_(tdpvf.muF_),
    sigma_(tdpvf.sigma_)
   // muS_(tdpvf.muS_)
{}


particleSlipSchneiderbauerSchellanderFvPatchVectorField::particleSlipSchneiderbauerSchellanderFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    alpha2_(readScalar(dict.lookup("alpha2"))),
    muF_(readScalar(dict.lookup("muF"))),
    sigma_(readScalar(dict.lookup("sigma")))
   // muS_(readScalar(dict.lookup("muS")))
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


particleSlipSchneiderbauerSchellanderFvPatchVectorField::particleSlipSchneiderbauerSchellanderFvPatchVectorField
(
    const particleSlipSchneiderbauerSchellanderFvPatchVectorField& tdpvf,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(tdpvf, iF),
    restitutionCoefficient_(tdpvf.restitutionCoefficient_),
    alpha2_(tdpvf.alpha2_),
    muF_(tdpvf.muF_),
    sigma_(tdpvf.sigma_)
  // muS_(tdpvf.muS_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void particleSlipSchneiderbauerSchellanderFvPatchVectorField::updateCoeffs()
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

    const fvPatchVectorField& U1 =
        patch().lookupPatchField<volVectorField, vector>("U1");

    const fvPatchScalarField& alpha1 = 
	patch().lookupPatchField<volScalarField, scalar>("alpha1");

    const fvPatchScalarField& g0 = 
	patch().lookupPatchField<volScalarField, scalar>("gs0");

    const fvPatchScalarField& mu1 = 
	patch().lookupPatchField<volScalarField, scalar>("mu1");

    const fvPatchScalarField& pf =
        patch().lookupPatchField<volScalarField, scalar>("pfH");

    const fvPatchScalarField& muf =
        patch().lookupPatchField<volScalarField, scalar>("muf");

    vectorField U1Patch = U1.patchInternalField();
    scalarField magU1Patch = mag(U1.patchInternalField());
    scalarField alpha1Patch = alpha1.patchInternalField() + 1.0e-7;
    scalarField ThetaPatch = max(alpha1Patch, 1.0e-6);
    scalarField g0Patch = g0.patchInternalField();
    scalarField mu1Patch = mu1.patchInternalField()+1.0e-6;
    
    if (db().foundObject<volScalarField>("Theta"))
    {
	const fvPatchScalarField& Theta = 
	  patch().lookupPatchField<volScalarField, scalar>("Theta");
	
        ThetaPatch = Theta.patchInternalField()+1.0e-9;
    }
    
    // The partial slip BC in OpenFOAM is implemented as
    //
    // valueFraction*U + (1-valueFraction)*grad(U) = 0
    //
    // To find valueFraction, we re-write Schneiderbauer BC as
    //
    // c*U + grad(U)*delta n = 0
    //
    // where
    //
    // c = valueFraction/(1 - valueFraction)
    //
    // As a consequence
    //
    // valueFraction = c/(c + 1)
    //
    // tau = - mu gradU
    //
    // - grad U = tau/mu
    //
    // and so
    //
    // c = tau* delta n /(mu*U)
    //
    // and to find delta n you should look into the partialSlip BC file
    //
    // delta n = 2/this->patch().deltaCoeffs();
    

   if (db().foundObject<volScalarField>("Theta"))
   {
       const fvPatchScalarField& Theta =
         patch().lookupPatchField<volScalarField, scalar>("Theta");

       ThetaPatch = Theta+1.0e-9;
   }

   scalarField mu_ = magU1Patch; //just to initialize

   forAll(magU1Patch,celli)
   {
   

       //if(magU1Patch[celli] < 1.0e-4)
      // {
            mu_[celli] = muF_;
      // }
      // else
       //{
         //   mu_[celli] = muS_;
       //}

   }
   
   
      
    
    scalarField u = magU1Patch;
    scalarField v = magU1Patch;
    scalarField Un = (this->patch().nf() & U1Patch);
    scalarField alpha = magU1Patch;
    scalarField gama = magU1Patch;
    forAll(U1Patch,celli)
    
   {
   
  /* if(sigma_ = 0.0)
   {
   
    gama[celli] = sigma_;
   }
   
   else if(sigma_>0.0) 
   {*/
   
  // alpha[celli] = (atan(mag(U1Patch[celli].z())/mag(U1Patch[celli].x())))*180/M_PI;
   // alpha[celli] = scalar(4.0);
   alpha[celli] = atan(((sqrt(scalar(1.5)*ThetaPatch[celli]))/2+mag(U1Patch[celli].z()))/(mag(U1Patch[celli].x())))*180/M_PI;
   
      gama[celli] = (-0.398942*exp(-4050/(sigma_*sigma_))*sigma_+0.398942*exp(-0.5*alpha[celli]*alpha[celli]/(sigma_*sigma_))*sigma_)/(0.5*erf(63.6396/sigma_)+0.5*erf(0.707107*alpha[celli]/sigma_));
    
   // }
   
       
     u[celli]=(U1Patch[celli].x()*cos(gama[celli]*M_PI/180)+Un[celli]*sin(gama[celli]*M_PI/180));

     v[celli]=(-U1Patch[celli].x()*sin(gama[celli]*M_PI/180)+Un[celli]*cos(gama[celli]*M_PI/180));
     
    // Info<<"dist: "<<cos(gama_*M_PI/180)<<"\n";
   //  Info<<"dist: "<<v[celli]<<"\n";

    }
    
 //   forAll(v,celli)
 //  {

   //    if(v[celli] < 0)
       
  //    {


    scalarField mu0 = 7.0/2.0*(1.0+restitutionCoefficient_)/(1.0+alpha2_)*mu_;

    scalarField us = (u+mu0*v)/(sqrt(2.0*(ThetaPatch+scalar(1e-6)))*mu0);

    scalarField erfU1 = erf(us);

    scalar eta = 1.0/2.0*(1.0+restitutionCoefficient_);

    scalarField tau = -rho1.value()*alpha1Patch*g0Patch*eta*mu_/mu0*(u*v*(1-erfU1)+sqrt(2*ThetaPatch/
    M_PI)*(exp(-v*v/(2*(ThetaPatch+scalar(1e-6))))-exp(-us*us))*mu0*(v)+(v*v+(ThetaPatch+scalar(1e-6)))*mu0*(erf(v/sqrt(2*(ThetaPatch+scalar(1e-6))))-erfU1));


    scalarField c = tau/((u*cos(gama*M_PI/180)-v*sin(gama*M_PI/180))*mu1Patch)*scalar(2.0)/this->patch().deltaCoeffs();

    this->valueFraction() = c/(c+scalar(1));  
    
        
   // Info<<"dist: "<<u*cos(gama*M_PI/180)+v*sin(gama*M_PI/180)<<"\n";
  // Info<<"dist: "<<tau<<"\n";

    partialSlipFvPatchVectorField::updateCoeffs();


    //Info<<"dist: "<<scalar(1)/this->patch().deltaCoeffs()<<"\n";

// } 
// }
   
}

// Write
void particleSlipSchneiderbauerSchellanderFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha2")
        << alpha2_ << token::END_STATEMENT << nl;
    os.writeKeyword("muF")
        << muF_ << token::END_STATEMENT << nl;
   // os.writeKeyword("muS")
    //    << muS_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigma")
        << sigma_ << token::END_STATEMENT << nl;
    writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    particleSlipSchneiderbauerSchellanderFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
