const volScalarField term1 = D1*pow(2.71828,-D2*triaxiality);
const volScalarField term2 = D5*pow(2.71828,-D6*triaxiality);
const volScalarField term3 = D3*pow(2.71828,-D4*triaxiality);
const volScalarField epsilonF = (0.5*(term1+term2)-term3)*pow(lodeParameter, 2.0)
                               + 0.5*(term1-term2)*lodeParameter+term3;
DDamage = DEpsilonPEq/epsilonF;

