/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    uncoupledDamage

Description
   This application calculates uncoupled damage fields. 
   
   The traixiality, decomposed stress tensor, the max principal stress and 
   lode angle parameter fields are calculated. The particular uncoupled 
   damage law can then be easily modified to a users preference.

   The user has the option to print out the traixiality, decomposed stress tensor,
   the principal stress and the lode angle parameter by running this application 
   with the option -verbose.

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "mechanicalModel.H"
#include "logVolFields.H"
#include "Eigen/Dense"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Function to caluclate principal stresses. Currently set to only
// assign value for max prinicpal stress. 
void calculateMaxPrincipalStress
(
  const symmTensor& sigma,
  scalar& sigmaPMax
)
{
    

    // Convert tau to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        sigma.xx(), sigma.xy(), sigma.xz(),
        sigma.xy(), sigma.yy(), sigma.yz(),
        sigma.xz(), sigma.yz(), sigma.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();
    label iMax = -1;
    label iMid = -1;
    label iMin = -1;
    const scalar a = EVals[0];
    const scalar b = EVals[1];
    const scalar c = EVals[2];


    if (a < b)
    {
        if (a < c)
        {
            if (b < c)
            {
                // a < b
                // a < c
                // b < c
                // a < b < c
                iMin = 0;
                iMid = 1;
                iMax = 2;
            }
            else
            {
                // a < b
                // a < c
                // b > c
                // a < c < b
                iMin = 0;
                iMid = 2;
                iMax = 1;
            }
        }
        else
        {
            // a < b
            // a > c
            // c < a < b
            iMin = 2;
            iMid = 0;
            iMax = 1;
        }
    }
    else
    {
        if (b < c)
        {
            if (a < c)
            {
                // a > b
                // b < c
                // a < c
                // b < a < c
                iMin = 1;
                iMid = 0;
                iMax = 2;
            }
            else
            {
                // a > b
                // b < c
                // a > c
                // b < c < a
                iMin = 1;
                iMid = 2;
                iMax = 0;
            }
        }
        else
        {
            // a > b
            // b > c
            // c < b < a
            iMin = 2;
            iMid = 1;
            iMax = 0;
        }
    }

    sigmaPMax = EVals[iMax];
}

// Decompose stress tensor into it's positive and negative compinents
void decompose
(
  const symmTensor& sigma,
  symmTensor& sigmaPositive,
  symmTensor& sigmaNegative
)
{
    

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        sigma.xx(), sigma.xy(), sigma.xz(),
        sigma.xy(), sigma.yy(), sigma.yz(),
        sigma.xz(), sigma.yz(), sigma.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();
    //Initialise Vectors which will hold positive and negative Eigenvalues respectively
    Eigen::Vector3d EPosVals=EVals;
    Eigen::Vector3d ENegVals=EVals;

    //Set positive and negative Eignvalue vectors
    for (int i=0;i<3;i++)
    {
        if (EVals(i)>0)
        {
            ENegVals(i)=0.0;
        }
        else
        {
            EPosVals(i)=0.0;
        }
    }
    
    Eigen::Matrix3d P = EPosVals.asDiagonal();
    Eigen::Matrix3d N = ENegVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();
    Eigen::Matrix3d resultPos = EVecs*P*EVecs.inverse();
    Eigen::Matrix3d resultNeg = EVecs*N*EVecs.inverse();

    sigmaPositive.xx()=resultPos(0,0);
    sigmaPositive.xy()=resultPos(0,1);
    sigmaPositive.xz()=resultPos(0,2);
    sigmaPositive.yy()=resultPos(1,1);
    sigmaPositive.yz()=resultPos(1,2);
    sigmaPositive.zz()=resultPos(2,2);

    sigmaNegative.xx()=resultNeg(0,0);
    sigmaNegative.xy()=resultNeg(0,1);
    sigmaNegative.xz()=resultNeg(0,2);
    sigmaNegative.yy()=resultNeg(1,1);
    sigmaNegative.yz()=resultNeg(1,2);
    sigmaNegative.zz()=resultNeg(2,2);

}

// Calulate the lode angle parameter
void calcLodeAngle(const symmTensor& S, scalar& lodeAngle)
{
    const scalar q = Foam::sqrt(3.0/2.0)*Foam::mag(S);
    scalar zheta = (27.0/2.0)*(Foam::det(S)/(Foam::pow(q, 3.0)));

    // This loop ensures that there is no floating point exception 
    // when numerical error leads to a value for zheta of e.g. 1.000001

    if (zheta > 1.0)
    {
        zheta = 1;
    }
    else if (zheta < -1.0)
    {
        zheta = -1;
    }

    lodeAngle = 1 - (2.0/3.14)*(Foam::acos(zheta));

    // This loop ensures that there is no floating point exception 
    // when numerical error leads to a value for the lodeAngle of e.g. 1.000001 

    if (lodeAngle > 1.0)
    {
        lodeAngle = 1.0;
    }
    else if (lodeAngle < -1.0)
    {
        lodeAngle = -1.0;
    }
}

// Macauley function
scalar macauley(scalar x)
{

    if (x < 0) 
    {
        x = 0;
    }
    return x;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("verbose", "");
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();

    bool verbose = false;
    if (args.optionFound("verbose"))
    {
        verbose = true;
    }

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    // Create accumulated damage field
    volScalarField Damage
    (
        IOobject
        (
            "damage",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
        //dimensionedScalar("zero", dimensionSet(2, -2, -4, 0, 0, 0, 0), 0.0)
    );

    volScalarField DDamage
    (
        IOobject
        (
            "deltaDamage",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );


    // Create mechanical model

    mechanicalModel mechanical(mesh);

    // Lame's first parameter
    const volScalarField lambda = mechanical.lambda();

    // Lame's second parameter -- shear modulus
    const volScalarField mu = mechanical.mu();

    // Young's modulus
    const volScalarField E = mu*(3.0*lambda + 2.0*mu)/(lambda + mu);

    // Poisson's ration
    const volScalarField nu = lambda/(2*(lambda + mu));

    // Read dictionary to define parameters
    IOdictionary damDict
    (
        IOobject
        (
            "ductileDamageDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Lookup parameters for damage law
#   include "damageLaws/uncoupledLemaitreCrackClosure/damageParameters.H"

    // Field to hold max principal stresses
    volScalarField maxPrincipalStress
    (
        IOobject
        (
            "maxPrincipalStress",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );

    // Lode angle field
    volScalarField lodeParameter
    (
        IOobject
        (
            "lodeParameter",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );

   // Positive and negative stress tensors

    volSymmTensorField sigmaPositive
    (
        IOobject
        (
            "sigmaPositive",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    volSymmTensorField sigmaNegative
    (
        IOobject
        (
            "sigmaNegative",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    // Create a field to store the triaxiality
    volScalarField triaxiality
    (
        IOobject
        (
            "triaxiality",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );


    // Create a field to store the old-time epsilonPEq field
    volScalarField epsilonPEq_0
    (
        IOobject
        (
            "epsilonPEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );

    // Calculate damage for all time-steps

    for (label i=startTime; i<endTime; i++)
    {

        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        Info << "    Reading mesh" << endl;
        mesh.readUpdate();

        IOobject sigmaheader
        (
            "sigmaCauchy",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject epsilonPEqheader
        (
            "epsilonPEq",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check sigma exists
        if (sigmaheader.headerOk() && epsilonPEqheader.headerOk())
        {
            Info<< "    Reading sigmaCauchy" << endl;
            const volSymmTensorField sigma(sigmaheader, mesh);

            Info<< "    Reading epsilonPEq" << endl;
            const volScalarField epsilonPEq(epsilonPEqheader, mesh);

            // Calculate DEpsilonPEq
            const volScalarField DEpsilonPEq = epsilonPEq - epsilonPEq_0;
    
            // Update the old-time epsilonPEq
            epsilonPEq_0 = 1.0*epsilonPEq;

            // For the first time-step (i.e. 0 time), we are not able to
            // calculate the increment as we do not haev the old-time
            if (startTime == 0)
            {
                Info<< "    Skipping initial time" << endl;
                continue;
            }

            // Calculate equivalent stress
            const volScalarField sigmaEq =
                sqrt((3.0/2.0)*magSqr(dev(sigma)))
              + dimensionedScalar("SMALL", dimPressure, SMALL);

            // Calculate hydrostatic stress
            const volScalarField sigmaHyd = (1.0/3.0)*tr(sigma);

            // Calculate constraint/trixiality
            triaxiality = sigmaHyd/sigmaEq;

            // Read displacement increment
            const volVectorField DU
            (
                IOobject
                (
                    "DU",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            // Calculate gradient of DU
            const volTensorField gradDU = fvc::grad(DU);

            // Calculate relative deformation gradient
            const volTensorField relF = I + gradDU.T();

            // Calculate true strain increment
            // Is this correct for logarithic strain? or should it be:
            // DEpsilon = 0.5*log(F.T() & F) - 0.5*log(F_old.T() & F_old) ...?
            const volSymmTensorField DEpsilon = 0.5*log(symm(relF.T() & relF));

            // Calculate equivalent true strain increment
            const volScalarField DEpsilonEq =
                sqrt((2.0/3.0)*magSqr(dev(DEpsilon)));

            // Calculate lode angle, postive and negative compnenets of the stress tensor
            // and the maximun principal stress for each cell 

            // Take references to internal fields
            scalarField& lodeParameterI = lodeParameter.internalField();
            scalarField& maxPrincipalStressI = maxPrincipalStress.internalField();
            symmTensorField& sigmaPositiveI = sigmaPositive.internalField();
            symmTensorField& sigmaNegativeI = sigmaNegative.internalField();
            const symmTensorField& sigmaI = sigma.internalField();

            forAll(lodeParameterI, cellI)
            {
                // Calculate the maximun principal stress for a given cell
                calculateMaxPrincipalStress
                ( 
                    sigmaI[cellI],
                    maxPrincipalStressI[cellI] 
                );

                // Decompose the stress tensor into it's positive and negative components for a given cell
                decompose
                (
                    sigmaI[cellI],
                    sigmaPositiveI[cellI],
                    sigmaNegativeI[cellI]
                );
   
                // Caculate the lode angle for each cell               
                calcLodeAngle
                (
                    dev(sigmaI[cellI]),
                    lodeParameterI[cellI]
                );

       
            }

            forAll(lodeParameter.boundaryField(), patchI)
            {

                // Take references to boundary fields
                scalarField& maxPrincipalStressP = maxPrincipalStress.boundaryField()[patchI];
                const symmTensorField& sigmaP = sigma.boundaryField()[patchI];
                symmTensorField& sigmaPositiveP = sigmaPositive.boundaryField()[patchI];
                symmTensorField& sigmaNegativeP = sigmaNegative.boundaryField()[patchI];
                scalarField& lodeParameterP = lodeParameter.boundaryField()[patchI];

                forAll(sigmaP, faceI)
                {
                    // Calculate the maximun principal stress for a given cell
                    calculateMaxPrincipalStress
                    (
                        sigmaP[faceI],
                        maxPrincipalStressP[faceI] 
                    );

                    // Decompose the stress tensor into it's positive and negative components for a given cell
                    decompose
                    (
                        sigmaP[faceI],
                        sigmaPositiveP[faceI],
                        sigmaNegativeP[faceI]
                    );
 
                    // Caculate the lode angle for each cell
                    calcLodeAngle
                    (
                        dev(sigmaP[faceI]),
                        lodeParameterP[faceI]
                    );
                }
             }

            // Calculate damage increment
           #   include "damageLaws/uncoupledLemaitreCrackClosure/damageLaw.H"

            // Increment total Damage
            Damage += DDamage;

            // Write Damage field
            Damage.write();

            Info<< nl << "Max Damage: " << max(Damage).value() << endl;

            // Option to print other fields
            if (verbose == true)
            {
                maxPrincipalStress.write();
                lodeParameter.write();
                sigmaPositive.write();
                sigmaNegative.write();
                triaxiality.write();
            }
        }
        else
        {
            Info<< "    Ensure both sigmaCauchy and epsilonPEq files exist to "
                << "obtain damage" << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
