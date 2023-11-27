clear

cd '/Users/aishiguha/Documents/MATLAB/Code'



for rval = 0.8:0.1:0.8

     for g2val = 2.1:0.1:2.3

        for trl = 11:20

            tic
  
    
            %Variables to find Bell States
    
            r = rval; %squeezing strength
    
            gamma2 = g2val; %BSM Detector Threshold
    
            trial = trl;
    
            
            disp(r)
            disp(gamma2)
    
    
    
            sigma0 = 1/sqrt(2); % scale of vacuum fluctuations
    
            sigma = 1/sqrt(2);
    
            
    
            cd '/Users/aishiguha/Documents/MATLAB/BellStates'
    
            filename_load = 'BellStates_r-' + string(r) + '_g2-' + string(gamma2) + '_trial-' + string(trial) + '.mat';
    
            load(filename_load, 'z1H_m', 'z1V_m', 'z2H_m', 'z2V_m','z3H_m','z3V_m', 'z4H_m','z4V_m','z1p_m', 'z2p_m', 'z3p_m', 'z4p_m', 'z5p_m', 'z6p_m', 'z7p_m', 'z8p_m', 'BSM_Count' );
    
            %additonal
           
                    
            if BSM_Count < 300000
    
                limit = BSM_Count;
    
            else
    
                limit = 300000;
    
            end
    
                z1H_m =  z1H_m(1:limit);
    
                z1V_m =  z1V_m(1:limit);
    
                z2H_m =  z2H_m(1:limit);
    
                z2V_m =  z2V_m(1:limit);
    
                z3H_m =  z3H_m(1:limit);
    
                z3V_m =  z3V_m(1:limit);
    
                z4H_m =  z4H_m(1:limit);
    
                z4V_m =  z4V_m(1:limit);
    
                z1p_m =  z1p_m(1:limit);
    
                z2p_m =  z2p_m(1:limit);
    
                z3p_m =  z3p_m(1:limit);
    
                z4p_m =  z4p_m(1:limit);
    
                z5p_m =  z5p_m(1:limit);
    
                z6p_m =  z6p_m(1:limit);
    
                z7p_m =  z7p_m(1:limit);
    
                z8p_m =  z8p_m(1:limit);
    
            
    
            %setting QST Detector Threshold
    
            gamma1 = (0:0.1:3)'; 
    
            
    
            cd '/Users/aishiguha/Documents/MATLAB/Data'
    
            filename = 'Double_r-'+string(r) +'_g2-' +string(gamma2) + '_g1-'+string(gamma1(length(gamma1)))+','+string(gamma1(2)-gamma1(1))+','+string(gamma1(1)) + '_trial-' + string(trial) +'.mat';
    
            
    
            Fidelity = zeros(1,length(gamma1));
    
            BellStat = zeros(1,length(gamma1));
    
            M = length(gamma1);
    
            
    
            TotalSampled = 0;
    
            ArrayofEfficiencies = [];
    
            
    
            %Psi- from SPDC 1
    
            
    
            a1H = sigma0*( cosh(r)*z1H_m + sinh(r)*conj(z2V_m) ); %matrix mult, column x row
    
            a1V = sigma0*( cosh(r)*z1V_m - sinh(r)*conj(z2H_m) );
    
            a2H = sigma0*( cosh(r)*z2H_m - sinh(r)*conj(z1V_m) );
    
            a2V = sigma0*( cosh(r)*z2V_m + sinh(r)*conj(z1H_m) );
    
            
    
            %Psi- from SPDC 2
    
            a3H = sigma*( cosh(r)*z3H_m + sinh(r)*conj(z4V_m) ); %row x column
    
            a3V = sigma*( cosh(r)*z3V_m - sinh(r)*conj(z4H_m) );
    
            a4H = sigma*( cosh(r)*z4H_m - sinh(r)*conj(z3V_m) );
    
            a4V = sigma*( cosh(r)*z4V_m + sinh(r)*conj(z3H_m) );
    
            
    
            %Beamsplitter for heralding
    
            
    
            b5H = (1/sqrt(2))*(a2H + a4H); %replace sigma
    
            b5V = (1/sqrt(2))*(a2V + a4V);
    
            b6H = (1/sqrt(2))*(a2H - a4H);
    
            b6V = (1/sqrt(2))*(a2V - a4V);
    
            
    
            theta = 0;
    
            
    
            %polarizers for heralding
    
            c5H = b5H*cos(theta)^2 + b5V*cos(theta)*sin(theta) + sigma*repmat(z5p_m,M,1)*(1-cos(theta)^2) - sigma*repmat(z6p_m,M,1)*cos(theta)*sin(theta); %column vector of 1's of length r, repmat
    
            c5V = b5V*cos(theta)*sin(theta) + b5V*sin(theta)^2 - sigma*repmat(z5p_m,M,1)*cos(theta)*sin(theta) + sigma*repmat(z6p_m,M,1)*(1-sin(theta)^2);
    
            
    
            eta = pi/2;
    
            
    
            c6H = b6H*cos(eta)^2 + b6V*cos(eta)*sin(eta) + sigma*repmat(z7p_m,M,1)*(1-cos(eta)^2) - sigma*repmat(z8p_m,M,1)*cos(eta)*sin(eta);
    
            c6V = b6V*cos(eta)*sin(eta) + b6V*sin(eta)^2 - sigma*repmat(z7p_m,M,1)*cos(eta)*sin(eta) + sigma*repmat(z8p_m,M,1)*(1-sin(eta)^2);
    
            
    
            
    
            %Polarizing beamsplitters and QST
    
            
    
            % choose PBS basis for right side (detectors 1 & 2)
    
            %cycle through PBS's using for loops
    
            
    
            for x = 1:4
    
            for y = 1:4 % does putting the for loop here cover all realizations? How exactly do the realizations work? 
    
            switch x
    
               case 1 % H/V basis
    
                  b1H = repmat(a1H,M,1);
    
                  b1V = sigma0*repmat(z2p_m,M,1); % vacuum mode from unused input port
    
                  b2H = sigma0*repmat(z1p_m,M,1); % vacuum mode from unused input port
    
                  b2V = repmat(a1V,M,1);
    
               case 2 % D/A basis
    
                  % detectors 1 & 3
    
                  b1H = (a1H + a1V)/2 + sigma0*(repmat(z1p_m,M,1) - repmat(z2p_m,M,1))/2;
    
                  b1V = (a1H + a1V)/2 - sigma0*(repmat(z1p_m,M,1) - repmat(z2p_m,M,1))/2;
    
                  b2H = sigma0*(repmat(z1p_m,M,1) + repmat(z2p_m,M,1))/2 + (a1H - a1V)/2;
    
                  b2V = sigma0*(repmat(z1p_m,M,1) + repmat(z2p_m,M,1))/2 - (a1H - a1V)/2;
    
               case 3 % R/L basis
    
                  b1H =    (a1H - 1i*a1V)/2 +    sigma0*(repmat(z1p_m,M,1) + 1i*repmat(z2p_m,M,1))/2;
    
                  b1V = 1i*(a1H - 1i*a1V)/2 - 1i*sigma0*(repmat(z1p_m,M,1) + 1i*repmat(z2p_m,M,1))/2;
    
                  b2H =    sigma0*(repmat(z1p_m,M,1) - 1i*repmat(z2p_m,M,1))/2 +    (a1H + 1i*a1V)/2;
    
                  b2V = 1i*sigma0*(repmat(z1p_m,M,1) - 1i*repmat(z2p_m,M,1))/2 - 1i*(a1H + 1i*a1V)/2;
    
                case 4 % H/V basis for Z, do I change this
    
                  b1H = repmat(a1H,M,1);
    
                  b1V = sigma0*repmat(z2p_m,M,1); % vacuum mode from unused input port
    
                  b2H = sigma0*repmat(z1p_m,M,1); % vacuum mode from unused input port
    
                  b2V = repmat(a1V,M,1);
    
            end
    
            
    
            % choose PBS basis for left side (detectors 3 & 4)
    
            switch y
    
               case 1 % H/V basis (I)
    
                  b3H = repmat(a3H,M,1);
    
                  b3V = sigma0*repmat(z4p_m,M,1); % vacuum mode from unused input port, repmat is okay?
    
                  b4H = sigma0*repmat(z3p_m,M,1); % vacuum mode from unused input port
    
                  b4V = repmat(a3V,M,1);
    
               case 2 % D/A basis (X)
    
                  b3H = (a3H + a3V)/2 + sigma0*(repmat(z3p_m,M,1) - repmat(z4p_m,M,1))/2;
    
                  b3V = (a3H + a3V)/2 - sigma0*(repmat(z3p_m,M,1) - repmat(z4p_m,M,1))/2;
    
                  b4H = sigma0*(repmat(z3p_m,M,1) + repmat(z4p_m,M,1))/2 + (a3H - a3V)/2;
    
                  b4V = sigma0*(repmat(z3p_m,M,1) + repmat(z4p_m,M,1))/2 - (a3H - a3V)/2;
    
               case 3 % R/L basis (Y)
    
                  b3H =    (a3H - 1i*a3V)/2 +    sigma0*(repmat(z3p_m,M,1) + 1i*repmat(z4p_m,M,1))/2;
    
                  b3V = 1i*(a3H - 1i*a3V)/2 - 1i*sigma0*(repmat(z3p_m,M,1) + 1i*repmat(z4p_m,M,1))/2;
    
                  b4H =    sigma0*(repmat(z3p_m,M,1) - 1i*repmat(z4p_m,M,1))/2 +    (a3H + 1i*a3V)/2;
    
                  b4V = 1i*sigma0*(repmat(z3p_m,M,1) - 1i*repmat(z4p_m,M,1))/2 - 1i*(a3H + 1i*a3V)/2;
    
                case 4 % H/V basis for Z gate (can it be the exact same?)
    
                  b3H =  repmat(a3H,M,1);
    
                  b3V = sigma0*repmat(z4p_m,M,1); % vacuum mode from unused input port
    
                  b4H = sigma0*repmat(z3p_m,M,1); % vacuum mode from unused input port
    
                  b4V = repmat(a3V,M,1);
    
            end
    
            
    
            % Boolean detection outcomes
    
            D1 = abs(b1H) > gamma1 | abs(b1V) > gamma1; %highlight
    
    D2 = abs(b2H) > gamma1 | abs(b2V) > gamma1;
    
            D3 = abs(b3H) > gamma1 | abs(b3V) > gamma1;
    
            D4 = abs(b4H) > gamma1 | abs(b4V) > gamma1;
    
            
    
            % quadruples
    
            C1356 = sum(D1 & D3 & ~D2 & ~D4, 2); % Detectors 1 & 3 & 5 & 6, want to sum over columns
    
            C2356 = sum(D2 & D3 & ~D1 & ~D4, 2); % Detectors 2 & 3 & 5 & 6
    
            C1456 = sum(D1 & D4 & ~D2 & ~D3, 2); % Detectors 1 & 4 & 5 & 6
    
            C2456 = sum(D2 & D4 & ~D1 & ~D3, 2); % Detectors 2 & 4 & 5 & 6
    
            
    
    
    
            Sampled = C1356 + C2356 + C1456 + C2456; %This is just for the last coincidence count (So for bases ZZ)
    
            
    
            %fidelity and Bell Stat calculations
            disp(string(x) + ''+ string(y))
            if (x == 1) && (y == 1)
    
                E11 = (C1356 + C2456 + C1456 + C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 1) && (y == 2)
    
                E12 = (C1356 - C2456 - C1456 + C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 1) && (y == 3)
    
                E13 = (C1356 - C2456 - C1456 + C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 1) && (y == 4)
    
                E14 = (C1356 - C2456 - C1456 + C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 2) && (y == 1)
    
                E21 = (C1356 - C2456 + C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 2) && (y == 2)
    
                E22 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 2) && (y == 3)
    
                E23 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 2) && (y == 4)
    
                E24 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 3) && (y == 1)
    
                E31 = (C1356 - C2456 + C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 3) && (y == 2)
    
                E32 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 3) && (y == 3)
    
                E33 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 3) && (y == 4)
    
                E34 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 4) && (y == 1)
    
                E41 = (C1356 - C2456 + C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 4) && (y == 2)
    
                E42 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 4) && (y == 3)
    
                E43 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            if (x == 4) && (y == 4)
    
                E44 = (C1356 + C2456 - C1456 - C2356)./(C1356 + C2456 + C1456 + C2356);
    
            end
    
            
    
            TotalSampled = TotalSampled + Sampled; %Used for efficiency
    
            
    
            end
    
            end 
    
            
    
            I = [1 0; 0 1];
    
            X = [0 1; 1 0];
    
            Y = [0 -1i; 1i 0];
    
            Z = [1 0; 0 -1];
    
            braPsi = [0 1 -1 0];
    
            ketPsi = [0;1;-1;0];
    
            
    
            A = Z;
    
            A_p = X;
    
            B = -(Z + X)/(sqrt(2));
    
            B_p = (Z - X)/(sqrt(2));
    
            rhomatrix = zeros(4,4,length(gamma1));
    
            for i = 1:length(gamma1) %remember to change across parameters
    
                densityrho = 0.25*(kron(I,I)*E11(i) + kron(I,X)*E12(i) + kron(I, Y)*E13(i) + kron(I, Z)*E14(i) + kron(X,I)*E21(i) + kron(X,X)*E22(i) + kron(X,Y)*E23(i) + kron(X,Z)*E24(i) + kron(Y,I)*E31(i) + kron(Y,X)*E32(i) + kron(Y,Y)*E33(i) + kron(Y,Z)*E34(i) + kron(Z,I)*E41(i) + kron(Z,X)*E42(i) + kron(Z,Y)*E43(i) + kron(Z,Z)*E44(i));
    
                Fidelity(i) = 0.5*braPsi*densityrho*ketPsi;
    
                BellStat(i) = abs(trace(densityrho*(kron(A,B))) + trace(densityrho*(kron(A_p,B)))) + abs(trace(densityrho*(kron(A_p,B_p))) - trace(densityrho*(kron(A,B_p))));
                
                rhomatrix(:,:,i) = densityrho;
            end
    
            
    
            Efficiency = TotalSampled/(16*(BSM_Count)); %Do we divde by just N, or 16*N? Bc, I have N realizations every basis change
    
            
    
            
    
            %----------------------------------------------------------------------
    
            % Coincidence Detection Efficiency
    
            
    
            for i = 1:2
    
            for j = 1:2
    
            switch i
    
                case 1
    
                    b1H_b = a1H;
    
                    b1V_b = z2p_m;
    
                    b2H_b = z1p_m;
    
                    b2V_b = a1V;
    
                case 2
    
                    b1H_b = (a1H + a1V)*(1/sqrt(2));
    
                    b1V_b = z2p_m;
    
                    b2H_b = z1p_m;
    
                    b2V_b = (a1H - a1V)*(1/sqrt(2));
    
            end
    
            
    
            switch j
    
                case 1
    
                    b3H_b = a3H*(cos(pi/8)) + a3V*(sin(pi/8));
    
                    b3V_b = z4p_m;
    
                    b4H_b = z3p_m;
    
                    b4V_b = -a3H*(sin(pi/8)) + a3V*(cos(pi/8));
    
                case 2
    
                    b3H_b = a3H*(sin(pi/8)) + a3V*(cos(pi/8));
    
                    b3V_b = z4p_m;
    
                    b4H_b = z4p_m;
    
                    b4V_b = -a3H*(cos(pi/8)) + a3V*(sin(pi/8));
    
            end
    
            % Boolean detection outcomes
    
            D1_b = abs(b1H_b) > gamma1 | abs(b1V_b) > gamma1; %highlight
    
            D2_b = abs(b2H_b) > gamma1 | abs(b2V_b) > gamma1;
    
            D3_b = abs(b3H_b) > gamma1 | abs(b3V_b) > gamma1;
    
            D4_b = abs(b4H_b) > gamma1 | abs(b4V_b) > gamma1;
    
            
    
            
    
            % quadruples
    
            C1356_b = sum(D1_b & D3_b & ~D2_b & ~D4_b, 2); % Detectors 1 & 3 & 5 & 6, want to sum over columns
    
            C2356_b = sum(D2_b & D3_b & ~D1_b & ~D4_b, 2); % Detectors 2 & 3 & 5 & 6
    
            C1456_b = sum(D1_b & D4_b & ~D2_b & ~D3_b, 2); % Detectors 1 & 4 & 5 & 6
    
            C2456_b = sum(D2_b & D4_b & ~D1_b & ~D3_b, 2); % Detectors 2 & 4 & 5 & 6
    
            
    
            %Bob part in Alice/Bob Efficiency
    
            B1 = C1356_b./sum(D3_b&~D4_b,2);
    
            B2 = C1456_b./sum(~D3_b&D4_b,2);
    
            B3 = C2356_b./sum(D3_b&~D4_b,2);
    
            B4 = C2456_b./sum(~D3_b&D4_b,2);
    
            
    
            %Alice part in Alice/Bob Efficiency
    
            A1 = C1356_b./sum(D1_b&~D2_b,2);
    
            A2 = C1456_b./sum(D1_b&~D2_b,2);
    
            A3 = C2356_b./sum(~D1_b&D2_b,2);
    
            A4 = C2456_b./sum(~D1_b&D2_b,2);
    
            
    
            ABArray = [B1,B2,B3,B4,A1,A2,A3,A4];
    
            ABEfficiency_new = min(ABArray,[],2);
    
            ArrayofEfficiencies = [ArrayofEfficiencies , ABEfficiency_new];
    
            ABEfficiency = min(ArrayofEfficiencies,[],2);
    
            end
    
            end
    
            
            
            
    
            
    
            
    
            %Check on coincidence detection efficiency, put end statements in correct
    
            %place
    
            
    
            %What the
    
            %Take minimum of all of the above, and miminum over all bases..
    
            %[x y min(x,y)]
    
            
    
            %--------------------------------------------------------------------
    
            % Saving variables
    
            
    
            
    
            disp(BSM_Count)
    
            max_Fidelity = max(Fidelity);
    
            
    
             save(filename, 'gamma1', 'r', 'gamma2','Fidelity', 'BellStat', 'Efficiency', 'ABEfficiency','max_Fidelity', 'BSM_Count', 'trial', 'rhomatrix')
    
            

            
    
          
    
            toc
    
            
        end
     end
end


