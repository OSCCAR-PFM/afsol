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

Description
    Schneiderbauer & Schellander boundary condition for the granular temperature.
    
Authors
    Simon Schneiderbauer <simon.schneiderbauer@jku.at>
    David Schellander <david.schellander@jku.at, www.schellcom.at>

\*---------------------------------------------------------------------------*/

#include "particleThetaSchneiderbauerSchellanderFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleThetaSchneiderbauerSchellanderFvPatchScalarField::particleThetaSchneiderbauerSchellanderFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(p.size()),
    alpha2_(p.size()),
    muF_(p.size()),
    sigma_(p.size())
  //  muS_(p.size())
{
    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


particleThetaSchneiderbauerSchellanderFvPatchScalarField::particleThetaSchneiderbauerSchellanderFvPatchScalarField
(
    const particleThetaSchneiderbauerSchellanderFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    alpha2_(ptf.alpha2_),
    muF_(ptf.muF_),
    sigma_(ptf.sigma_)
    
   // muS_(ptf.muS_)
{}


particleThetaSchneiderbauerSchellanderFvPatchScalarField::particleThetaSchneiderbauerSchellanderFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    alpha2_(readScalar(dict.lookup("alpha2"))),
    muF_(readScalar(dict.lookup("muF"))),
    sigma_(readScalar(dict.lookup("sigma")))
   // muS_(readScalar(dict.lookup("muS")))
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
        particleThetaSchneiderbauerSchellanderFvPatchScalarField::evaluate();
    }
}


particleThetaSchneiderbauerSchellanderFvPatchScalarField::particleThetaSchneiderbauerSchellanderFvPatchScalarField
(
    const particleThetaSchneiderbauerSchellanderFvPatchScalarField& pivpvf,
    const DimensionedField<double, volMesh>& iF
)
:
    mixedFvPatchScalarField(pivpvf, iF),
    restitutionCoefficient_(pivpvf.restitutionCoefficient_),
    alpha2_(pivpvf.alpha2_),
    muF_(pivpvf.muF_),
    sigma_(pivpvf.sigma_)
   // muS_(pivpvf.muS_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleThetaSchneiderbauerSchellanderFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void particleThetaSchneiderbauerSchellanderFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void particleThetaSchneiderbauerSchellanderFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if ((restitutionCoefficient_ < 0) || (restitutionCoefficient_ > 1))
    {
	FatalErrorIn
        (
            "particleThetaSchneiderbauerSchellanderFvPatchScalarField::"
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
      
    scalarField magU1Patch = mag(U1.patchInternalField());
    vectorField U1Patch = U1.patchInternalField();
      

    //scalarField ThetaPatch = max(patchInternalField(), 1.0e-6);
    scalarField ThetaPatch = patchInternalField()+1.0e-9;
    
    scalarField alpha1Patch = max(alpha1.patchInternalField(), 1.0e-6);
    
    scalarField g0Patch = g0.patchInternalField();
    
    scalarField kappaPatch = max(kappa.patchInternalField(), 1.0e-15);

    const fvPatchScalarField& pf =
        patch().lookupPatchField<volScalarField, scalar>("pfH");
 
    // The mixed BC in OpenFOAM is implemented as
    //
    // valueFraction*(Theta - ThetaRef) + (1-valueFraction)*(grad(Theta) - grad(Theta)_Ref) = 0
    //
    // To find valueFraction, we re-write thes BC as
    //
    // c*(Theta - ThetaRef) + (grad(Theta) - grad(Theta)_Ref) = 0
    //
    // where
    //
    // c = valueFraction/(1 - valueFraction)
    //
    // As a consequence
    //
    // valueFraction = c/(c+1)
    //
    // grad(Theta)_Ref = 0
    //

    scalarField mu_ = magU1Patch; //just to initialize

    forAll(magU1Patch,celli)
    {

       // if(magU1Patch[celli] < 1.0e-4)
      //  {
             mu_[celli] = muF_;
       // }
       // else
       // {
        //     mu_[celli] = muS_;
       // }

    }
    
    scalarField u = magU1Patch;
    scalarField v = magU1Patch;
    scalarField Un = (this->patch().nf() & U1Patch);
    scalarField alpha =  magU1Patch;
    scalarField gama =  magU1Patch;
    scalarField nz =  magU1Patch;
    forAll(U1Patch,celli)
   {
   
   /*  if(sigma_ = 0)
   {
   
    gama[celli] = sigma_;
   }
   
   else if(sigma_>0.0)
   {*/
   
  // alpha[celli] = (atan(mag(U1Patch[celli].z())/mag(U1Patch[celli].x())))*180/M_PI;
   // alpha[celli] = scalar(4.0);
    alpha[celli] = atan(((sqrt(scalar(1.5)*ThetaPatch[celli]))/2+mag(U1Patch[celli].z()))/(mag(U1Patch[celli].x())))*180/M_PI;
   
      gama[celli] =(-0.398942*exp(-4050/(sigma_*sigma_))*sigma_+0.398942*exp(-0.5*alpha[celli]*alpha[celli]/(sigma_*sigma_))*sigma_)/(0.5*erf(63.6396/sigma_)+0.5*erf(0.707107*alpha[celli]/sigma_));
      
    
  //  }
    //  Info<<"dist: "<<gama[celli]<<"\n";
     
      
     u[celli]=(U1Patch[celli].x()*cos(gama[celli]*M_PI/180)+Un[celli]*sin(gama[celli]*M_PI/180));

     v[celli]=(-U1Patch[celli].x()*sin(gama[celli]*M_PI/180)+Un[celli]*cos(gama[celli]*M_PI/180));
     
     
      vectorField NORMAL = this->patch().nf();
       nz[celli] = NORMAL[celli].z();
     
     }
     
   //  forAll(v,celli)
  // {

    //   if(v[celli] < 0)
       
    //   {

  /*  scalar eta = 1.0/2.0*(1.0+restitutionCoefficient_);
    scalarField mu0 = 7.0/2.0*(1.0+restitutionCoefficient_)/(1.0+alpha2_)*mu_;
    scalarField us = magU1Patch/(sqrt(2.0*ThetaPatch)*mu0);
    scalarField erfU1 = erf(us);

    scalarField tau = -mu_*(alpha1Patch*eta*rho1.value()*g0Patch*ThetaPatch*erfU1); */
    
 
    
    scalar eta = 1.0/2.0*(1.0+restitutionCoefficient_);
    scalarField mu0 = 7.0/2.0*(1.0+restitutionCoefficient_)/(1.0+alpha2_)*mu_;
    scalarField us = (u+mu0*v)/(sqrt(2.0*(ThetaPatch+scalar(1e-6)))*mu0);
    scalarField erfU1 = erf(us);
    
   // scalarField tau = -mu_*(alpha1Patch*eta*rho1.value()*g0Patch*ThetaPatch*erfU1);

    scalarField tau = -rho1.value()*alpha1Patch*g0Patch*eta*mu_/mu0*(u*v*(1-erfU1)+sqrt(2*ThetaPatch/
    M_PI)*(exp(-v*v/(2*(ThetaPatch+scalar(1e-6))))-exp(-us*us))*mu0*(v)+(v*v+(ThetaPatch+scalar(1e-6)))*mu0*(erf(v/sqrt(2*(ThetaPatch+scalar(1e-6))))-erfU1)); 

    // Info<<"eta: "<<eta<<"\n";
  //  scalarField tauSU1 = tau*magU1Patch;
   
   scalarField N = eta*rho1.value()*alpha1Patch*g0Patch*((v*v+(ThetaPatch+scalar(1e-6)))*(1-erf(v/sqrt
     (2*ThetaPatch)))-sqrt(2*ThetaPatch/M_PI)*(v)*exp(-v*v/(2*(ThetaPatch+scalar(1e-6)))));
   
   scalarField tauSU1 = (tau*(u*cos(gama*M_PI/180)-v*sin(gama*M_PI/180))-N*(u*sin(gama*M_PI/180)+v*cos(gama*M_PI/180))); 
   
   scalarField ul = u/(sqrt(2.0*ThetaPatch)*mu0);
    scalarField erfU2 = erf(ul);
    
    scalarField fact1 = scalar(2.0)*mu_*u*u*(scalar(2.0)*eta-mu0);
    scalarField fact2 = ThetaPatch*(scalar(14.0)*mu_*eta-4*mu0*(scalar(1.0)+mu_)-scalar(6.0)*mu_*mu0*mu0*eta);
    scalarField fact3 = sqrt(ThetaPatch)*(scalar(4.0)*(eta-scalar(1.0))+scalar(6.0)*mu_*mu_*eta);
    scalarField fact4 = sqrt(scalar(2.0)*M_PI)*mu_*u*erfU2;
     
    // scalarField D0 = -(alpha1Patch*rho1.value()*g0Patch*eta*sqrt(ThetaPatch)/(mu0*mu0*sqrt(2.0*M_PI))*(exp(-(ul*ul))*(mu_*(fact1+fact2))+mu0*mu0*sqrt(ThetaPatch)*(fact3-fact4))); 
     
     scalarField D0 = -(alpha1Patch*rho1.value()*g0Patch*eta*sqrt(ThetaPatch)/(mu0*mu0*sqrt(2.0*M_PI))*(exp(-(ul*ul))*(mu_*(fact1+fact2))+mu0*mu0*sqrt(ThetaPatch)*(fact3-fact4)));
     
     scalarField fact5 =  ThetaPatch*(exp(-(ul*ul))*scalar(8.0)*mu0+scalar(8.0)*mu0*mu0*mu0*(scalar(1.0)-exp(-(ul*ul)))+exp(-(ul*ul))*mu_*(scalar(8.0)*mu0+exp(-(ul*ul))*scalar(18.0)*mu0*mu0*eta-scalar(28.0)*eta));
    scalarField fact6 = exp(-(ul*ul))*scalar(4.0)*mu_*u*u*(scalar(2.0)*eta-mu0);
    scalarField fact7 = scalar(3.0)*mu0*mu0*ThetaPatch*(scalar(4.0)*eta-scalar(4.0)+scalar(6.0)*mu_*mu_*eta);
    scalarField fact8 = mu_*(scalar(7.0)*mu_*(u*u+scalar(2.0)*ThetaPatch)*scalar(2.0)*eta-scalar
    (4.0)*mu0*(u*u+scalar(2.0)*ThetaPatch)*(scalar(1.0)+mu_)-scalar(18.0)*mu_*mu0*mu0*ThetaPatch*eta)*
    (scalar(1.0)-erfU2);
     
     scalarField D1 = -(alpha1Patch*rho1.value()*g0Patch*eta*(v)/(4*mu0*mu0*mu0*sqrt(ThetaPatch*M_PI))
     *((sqrt(2.0)*mu_*u*(fact5-fact6))-mu0*sqrt(ThetaPatch*M_PI)*(fact7+fact8)));
     
     
     scalarField fact9 = exp(-ul*ul)*scalar(2.0)*mu_*mu_*u*u*u*u*(scalar(2.0)*eta-mu0);
  scalarField fact10 = exp(-ul*ul)*(scalar(4.0)*mu0*(scalar(1.0)+mu_)-scalar(7.0)*mu_*scalar(2.0)*eta+mu_*mu0*mu0*scalar(2.0)*eta+scalar(2.0)*mu0*mu0*mu0*(mu_-scalar(1.0)));
  scalarField fact11 = exp(-ul*ul)*(scalar(4.0)*mu0*(scalar(1.0)+mu_)-scalar(7.0)*mu_*scalar(2.0)*eta+scalar(9.0)*mu_*mu0*mu0*scalar(2.0)*eta);
  scalarField fact12 = scalar(4.0)*(eta-scalar(1.0))+scalar(3.0)*mu_*mu_*scalar(2.0)*eta;
  scalarField fact13 = scalar(4.0)*sqrt(M_PI)*mu_*pow(mu0,4)*pow(ThetaPatch,1.5)*u*erfU2;
  
  scalarField D2 = -eta*rho1.value()*alpha1Patch*g0Patch*v*v/(scalar(4.0)*mu0*mu0*mu0*mu0*sqrt
  (M_PI*ThetaPatch*ThetaPatch*ThetaPatch))*(sqrt(2.0)*(fact9-mu_*ThetaPatch*u*u*(fact10)-
  mu0*mu0*ThetaPatch*ThetaPatch*(mu_*(fact11)-scalar(3.0)*mu0*mu0*(fact12)))-fact13);
  
  
  scalarField fact14 = exp(-ul*ul)*scalar(9.0)*mu_*scalar(2.0)*eta+scalar(4.0)*mu0*(scalar(1.0)-exp(-ul*ul));
   scalarField fact15 = exp(-ul*ul)*scalar(2.0)*mu_*pow(u,4)*(scalar(1.0)+restitutionCoefficient_-mu0);
   scalarField fact16 = exp(-ul*ul)*(scalar(7.0)*mu_*scalar(2.0)*eta-scalar(4.0)*mu0*(scalar(1.0)+mu_)-scalar
   (3.0)*mu_*mu0*mu0*scalar(2.0)*eta+scalar(2.0)*pow(mu0,3));
   scalarField fact17 = scalar(4.0)*(eta-scalar(1.0))+scalar(3.0)*mu_*mu_*scalar(2.0)*eta*erfU2;
   
   scalarField D3 = -eta*rho1.value()*alpha1Patch*g0Patch*(v)*v*v/(scalar(12.0)*pow(mu0,5)*sqrt
   (M_PI)*pow(ThetaPatch,2.5))*(sqrt(2.0)*mu_*u*(pow
   (mu0,4)*ThetaPatch*ThetaPatch*fact14-fact15-ThetaPatch*u*u*(fact16))-scalar(3.0)*sqrt(M_PI)*pow
   (mu0,5)*pow(ThetaPatch,2.5)*fact17);
    
    
    scalarField N0 = eta*rho1.value()*alpha1Patch*g0Patch*ThetaPatch;
    
  scalarField D = D0+D1+D2+D3;
 // scalarField D = D0/N0*N;
   
    
   //scalarField c = -(tauSU1 + term2)/(kappaPatch*ThetaPatch);

   this->valueFraction() =0.0;
   this->refValue() = 0.0;
   this->refGrad() = -(-tauSU1 + D)*nz/(kappaPatch);

   //Info<<"term: "<<term2<<"\n";
      /* Info<<"dist: "<<fact1<<"\n";
       Info<<"dist: "<<fact2<<"\n";
       Info<<"dist: "<<fact3<<"\n";
       Info<<"dist: "<<fact4<<"\n";
       Info<<"dist: "<<term2<<"\n";
       Info<<"dist: "<<kappaPatch<<"\n";*/
       //Info<<"dist: "<<D<<"\n";
       

    mixedFvPatchScalarField::updateCoeffs();
//}
//}
}


// Write
void particleThetaSchneiderbauerSchellanderFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha2")
        << alpha2_ << token::END_STATEMENT << nl;
    os.writeKeyword("muF")
        << muF_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigma")
        << sigma_ << token::END_STATEMENT << nl;
  //  os.writeKeyword("muS")
     //   << muS_ << token::END_STATEMENT << nl;
    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    particleThetaSchneiderbauerSchellanderFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
