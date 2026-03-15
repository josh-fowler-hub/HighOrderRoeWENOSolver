function [t,rhoR,rhoL,uL,uR,PR,PL,R,g,middle] = test_case(test)
% test cases for shock tube, based on Toro's book

switch test
    case 1 % modified sod shock tube
        % test 1
        fprintf('Case 1: Toro Test #1 (Modified Shock Tube) \n');
        t = 0.2;
        rhoR = 0.125;
        rhoL = 1;
        uL = 0.75;
        uR = 0;
        PR = 0.1;
        PL = 1;
        R = 8314/1;
        g = 1.4;
        middle = 0.3;
    case 2
        % test 2
        fprintf('Case 2: Toro Test #2 (123 Test) \n');
        t = 0.15;
        rhoR = 1;
        rhoL = 1;
        uL = -2;
        uR = 2;
        PR = 0.4;
        PL = 0.4;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
    case 3
        % test 3
        fprintf('Case 3: Toro Test #3 (Left Blast Wave) \n');
        t = 0.012;
        rhoR = 1;
        rhoL = 1;
        uL = 0;
        uR = 0;
        PR = 0.01;
        PL = 1000;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
    case 4
        % test 4
        fprintf('Case 4: Toro Test #4 (Right Blast Wave) \n');
        t = 0.035;
        rhoR = 5.99242;
        rhoL = 5.99924;
        uL = 19.5975;
        uR = -6.19633;
        PR = 46.0950;
        PL = 460.894;
        R = 8314/1;
        g = 1.4;
        middle = 0.4;
    case 5
        % test 5
        fprintf('Case 5: Toro Test #5 (Resulting Shock from Toro Test 4 & 5) \n');
        t = 0.012;
        rhoR = 1;
        rhoL = 1;
        uL = -19.59745;
        uR = -19.59745;
        PR = 0.01;
        PL = 1000;
        R = 8314/1;
        g = 1.4;
        middle = 0.8;
    case 6 % Configuration 1, Sod's Problem
        fprintf('Case 6: Sods problem \n');
        PL = 1;
        PR = 0.1;
        uL = 0;
        uR = 0;
        rhoL = 1;
        rhoR = 0.125;
        t = 0.1;
        R = 8314/1;
        g = 1.4;
        middle = 0.3;
        
    case 7 % Configuration 2, Left Expansion and right strong shock
        fprintf('Case 7: Left Expansion and right strong shock \n');
        PL = 1000;
        PR =  0.1;
        uL = 0;
        uR = 0;
        rhoL = 3;
        rhoR = 2;
        t = 0.02;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 8 % Configuration 3, Right Expansion and left strong shock
        fprintf('Case 8: Right Expansion and left strong shock \n');
        PL = 7;
        PR = 10;
        uL = 0;
        uR = 0;
        rhoL = 1;
        rhoR = 1;
        t = 0.1;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 9 % Configuration 4, Double Shock
        fprintf('Case 9: Double Shock \n');
        PL = 450;
        PR = 45;
        uL = 20;
        uR = -6;
        rhoL = 6;
        rhoR = 6;
        t = 0.02;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 10 % Configuration 5, Double Expansion
        fprintf('Case 10: Double Expansion \n');
        PL = 40;
        PR = 40;
        uL = -2;
        uR = 2;
        rhoL = 1;
        rhoR = 2.5;
        t = 0.03;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;

    case 11 % Configuration 6, Cavitation
        fprintf('Case 11: Cavitation \n');
        PL = 0.4;
        PR = 0.4;
        uL = -2;
        uR = 2;
        rhoL = 1;
        rhoR = 1;
        t = 0.1;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 12 % Shocktube problem of G.A. Sod, JCP 27:1, 1978
        fprintf('Case 12: Shocktube problem of G.A. Sod, JCP 27:1, 1978\n');
        PL = 1.0;
        PR = 0.1;
        uL = 0.75;
        uR = 0;
        rhoL = 1;
        rhoR = 0.125;
        t = 0.17;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 13 % Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Case 13: Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997\n');
        PL = 3.528;
        PR =  0.571;
        uL = 0.698;
        uR =  0;
        rhoL = 0.445;
        rhoR = 0.5;
        t = 0.15;
        R = 8314/1;
        g = 1.4;
        middle = 0.5; 
      
    case 14 % Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Case 14: Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997\n');
        PL = 10.333;
        PR = 1;
        uL = 0.92;
        uR = 3.55;
        rhoL = 3.857;
        rhoR = 1;
        t = 0.09;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 15 % Shocktube problem with supersonic zone
        fprintf('Case 15: Shocktube problem with supersonic zone\n');
        PL = 1;
        PR = 0.02;
        uL = 0;
        uR = 0.00;
        rhoL = 1;
        rhoR = 0.02;
        t = 0.162;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;

    case 16 % Stationary shock
        fprintf('Case 16: Stationary shock\n');
        PL = 1.0;
        PR = 0.1;
        uL = -2.0;
        uR = -2.0;
        rhoL = 1.0;
        rhoR = 0.125;
        t = 0.1;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
        
    case 17 % right side of 2-d riemman case 12
        fprintf('Case 17: right side of 2-d riemman case 14\n');
        PL = 1.0;
        PR = 0.4;
        uL = 0.7276;
        uR = 0.0;
        rhoL = 1.0;
        rhoR = 0.5313;
        t = 0.1;
        R = 8314/1;
        g = 1.4;
        middle = 0.5; 
        
    case 18 % Stationary Shock
        fprintf('Case 18: right side of 2-d riemman case 15\n');
        PL = 0.1;
        PR = 0.676;
        uL = 1.2;
        uR = 0.723966942148760;
        rhoL = 1.0;
        rhoR = 1.657534246575342;
        t = 0.1;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
    case 19 % Stationary Shock
        fprintf('Case 19: User Specified Test\n');
        PL = 1;
        PR = 1.01;
        uL = 1.2;
        uR = 1.2;
        rhoL = 1.0;
        rhoR = 1.01;
        t = 0.2;
        R = 8314/1;
        g = 1.4; % Helium 1.66
        middle = 0.75;
    case 20 % Manufactured Solution
        PL = 1;
        PR = 1;
        uL = 1;
        uR = 1;
        rhoL = 1;
        rhoR = 1;
        t = 0.2;
        R = 8314/1;
        g = 1.4;
        middle = 0.5;
end
end