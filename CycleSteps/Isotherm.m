function q=Isotherm(y, P, T, isotherm_parameters)
%hypotheticalisotherm- Calculate the molar loadings of a two component mixture
%   Calculate the total loading of a two component mixture loading. It is
%   assumed that the isotherm is able to be represented by a dual site
%   langmuir competitive isotherm with temperature dependent parameters.
%   The isotherm can be in terms of concentration [mol/m^3] or Partial
%   pressure [Pa]. 
%   Inputs:
%   y: mole fraction of component one [-]. It is assumed the mole fraction
%   of component 2 is 1-y
%   P: Total pressure of the gas [Pa]
%   T: Temperature of the gas [K]
%   isothermparams: Parameters for the isotherm. View script Input_PSA for
%   structure
%   input_units: specify whether the isotherm is in terms of concentration
%   of partial pressure

    R=8.314;

    q_s_b_1=isotherm_parameters(1);
    q_s_d_1=isotherm_parameters(3);
    q_s_b_2=isotherm_parameters(2);
    q_s_d_2=isotherm_parameters(4);
    b_1=isotherm_parameters(5);
    d_1=isotherm_parameters(7);
    b_2=isotherm_parameters(6);
    d_2=isotherm_parameters(8);
    deltaU_b_1=isotherm_parameters(9);
    deltaU_d_1=isotherm_parameters(11);
    deltaU_b_2=isotherm_parameters(10);
    deltaU_d_2=isotherm_parameters(12);


    B_1=b_1*exp(-deltaU_b_1/R./T);
    D_1=d_1*exp(-deltaU_d_1/R./T);
    B_2=b_2*exp(-deltaU_b_2/R./T);
    D_2=d_2*exp(-deltaU_d_2/R./T);

    if isotherm_parameters(13) == 0

        P_1=y.*P;
        P_2=(1-y).*P;
        input_1=P_1;
        input_2=P_2;

    elseif isotherm_parameters(13) ==1
        C_1=y.*P./R./T;
        C_2=(1-y).*P./R./T;
        input_1=C_1;
        input_2=C_2;
    else
        error('Please specify whether the isotherms are in terms of Concentration or Partial Pressure')

    end


    q1_b=q_s_b_1.*B_1.*input_1./(1+B_1.*input_1+B_2.*input_2);
    q1_d=q_s_d_1.*D_1.*input_1./(1+D_1.*input_1+D_2.*input_2);

    q1=q1_b+q1_d;

    q2_b=q_s_b_2.*B_2.*input_2./(1+B_1.*input_1+B_2.*input_2);
    q2_d=q_s_d_2.*D_2.*input_2./(1+D_1.*input_1+D_2.*input_2);

    q2=q2_b+q2_d;
	
	q = [q1, q2] ;

end