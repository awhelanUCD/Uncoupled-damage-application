//initialise total, positive and negative energy release rates
scalar Y,YPos,YNeg; 

const scalarField& DEpsilonPEqI = DEpsilonPEq.internalField();
scalarField& DDamageI = DDamage.internalField();
const scalarField& EI = E.internalField();
const scalarField& nuI = nu.internalField();

forAll(DDamageI, cellI)
{
    scalar traceSigma = Foam::tr(sigmaI[cellI]);

    //denominator associated with the positive stress tensor
    scalar denomPositive = 1/(2*EI[cellI]);

    //denominator associated with the negative stress tensor
    scalar denomNegative = h.value()/(2*EI[cellI]);

    YPos = denomPositive*(
                          (1.0 + nuI[cellI])*
                          (sigmaPositiveI[cellI]&&sigmaPositiveI[cellI])-nuI[cellI]*Foam::pow(macauley(traceSigma), 2.0)
                         );
    YNeg = denomNegative*(
                          (1.0 + nuI[cellI])*
                          (sigmaNegativeI[cellI]&&sigmaNegativeI[cellI])-nuI[cellI]*Foam::pow(macauley(-traceSigma), 2.0)
                         );

    Y = -YPos - YNeg;

    //update damage increment
    DDamageI[cellI] = Foam::pow(-Y/s0.value(), b.value())*DEpsilonPEqI[cellI];

}

forAll(DDamage.boundaryField(), patchI)
{

    // Take references to boundary fields
    const scalarField& EP = E.internalField();
    const scalarField& nuP = nu.internalField();
    const scalarField& DEpsilonPEqP = DEpsilonPEq.boundaryField()[patchI];
    scalarField& DDamageP = DDamage.boundaryField()[patchI];
    const symmTensorField& sigmaP = sigma.boundaryField()[patchI];
    const symmTensorField& sigmaPositiveP = sigmaPositive.boundaryField()[patchI];
    const symmTensorField& sigmaNegativeP = sigmaNegative.boundaryField()[patchI];

    forAll(DDamageP, faceI)
    {
        scalar traceSigma = Foam::tr(sigmaP[faceI]);

        //denominator associated with the positive stress tensor
        scalar denomPositive = 1/(2*EP[faceI]);

        //denominator associated with the negative stress tensor
        scalar denomNegative = h.value()/(2*EP[faceI]);

        YPos = denomPositive*(
                              (1.0 + nuP[faceI])*
                              (sigmaPositiveP[faceI]&&sigmaPositiveP[faceI])-nuP[faceI]*Foam::pow(macauley(traceSigma), 2.0)
                             );
        YNeg = denomNegative*(
                              (1.0 + nuP[faceI])*
                              (sigmaNegativeP[faceI]&&sigmaNegativeP[faceI])-nuP[faceI]*Foam::pow(macauley(-traceSigma), 2.0)
                             );

        Y = -YPos - YNeg;

        //update damage increment
        DDamageP[faceI] = Foam::pow(-Y/s0.value(), b.value())*DEpsilonPEqP[faceI];

    }
}


