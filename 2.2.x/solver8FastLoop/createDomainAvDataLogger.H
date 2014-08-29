    fileName probeDir = "domainAvData";
    probeDir = probeDir/mesh.time().timeName();

    mkDir(probeDir);

    OFstream* sPtrUSlip 		= new OFstream(probeDir/"uSlip");
    OFstream* sPtrIntegralMomentum 	= new OFstream(probeDir/"integralMomentum");
    OFstream* sPtrMeanFCorrSlipVel 	= new OFstream(probeDir/"meanFCorrSlipVel");
    OFstream* sPtrFluidMaxVel 		= new OFstream(probeDir/"fluidMaxVel");
    OFstream* sPtrFluidMinVel 		= new OFstream(probeDir/"fluidMinVel");
    OFstream* sPtrParticleMaxVel 	= new OFstream(probeDir/"particleMaxVel");
    OFstream* sPtrParticleMinVel 	= new OFstream(probeDir/"particleMinVel");
