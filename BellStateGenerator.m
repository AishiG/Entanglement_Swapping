clear
tic

cd '/Users/aishiguha/Documents/MATLAB/BellStates'


for rval = 2.0:0.1:2.1
    for g2val = 1.3:0.1:1.4
        for trl = 11:12

            N = 1000*1e3; % number of realizations (1000 = 1 ms)
            sigma0 = 1/sqrt(2); % scale of vacuum fluctuations
            sigma = 1/sqrt(2);
            r = rval; %squeezing strength
            gamma2 = g2val; %BSM detector threshold
            trial = trl;
            
            filename = 'test_BellStates_r-' + string(r) + '_g2-' + string(gamma2) + '_trial-' + string(trial) + '.mat';
            
            

            if isfile(filename)
                load(filename, 'z1H_m', 'z1V_m', 'z2H_m', 'z2V_m','z3H_m','z3V_m', 'z4H_m','z4V_m','z1p_m', 'z2p_m', 'z3p_m', 'z4p_m', 'z5p_m', 'z6p_m', 'z7p_m', 'z8p_m', 'BSM_Count');
                z1Hb = z1H_m;
                z1Vb = z1V_m;
                z2Hb = z2H_m;
                z2Vb = z2V_m;
                z3Hb = z3H_m;
                z3Vb = z3V_m;
                z4Hb = z4H_m;
                z4Vb = z4V_m;
                z1pb = z1p_m;
                z2pb = z2p_m;
                z3pb = z3p_m;
                z4pb = z4p_m;
                z5pb = z5p_m;
                z6pb = z6p_m;
                z7pb = z7p_m;
                z8pb = z8p_m;
            
                BSM_skip = BSM_Count;
            end
            
            

            %Initializing Arrays
            z1H_m = [];
            z1V_m =  [];
            z2H_m =  [];
            z2V_m =  [];
            z3H_m =  [];
            z3V_m =  [];
            z4H_m =  [];
            z4V_m =  [];
            
            %PBS Noise
            z1p_m =  [];
            z2p_m =  [];
            z3p_m =  [];
            z4p_m = [];
            %Polarizer Noise, include in BSM selection (how data work)
            z5p_m = [];
            z6p_m =  [];
            z7p_m =  [];
            z8p_m =  [];
            
            BSM = 0;
            
            while BSM <= 10000
            %Generating Bell States---------------------------------------
            % input random variables for the entanglement source
            z1H= complex(randn(1,N),randn(1,N))/sqrt(2);
            z1V = complex(randn(1,N),randn(1,N))/sqrt(2);
            z2H = complex(randn(1,N),randn(1,N))/sqrt(2);
            z2V = complex(randn(1,N),randn(1,N))/sqrt(2);
            z3H = complex(randn(1,N),randn(1,N))/sqrt(2);
            z3V = complex(randn(1,N),randn(1,N))/sqrt(2);
            z4H = complex(randn(1,N),randn(1,N))/sqrt(2);
            z4V = complex(randn(1,N),randn(1,N))/sqrt(2);
            
            %PBS Noise
            z1p = complex(randn(1,N),randn(1,N))/sqrt(2);
            z2p = complex(randn(1,N),randn(1,N))/sqrt(2);
            z3p = complex(randn(1,N),randn(1,N))/sqrt(2);
            z4p = complex(randn(1,N),randn(1,N))/sqrt(2);
            %Polarizer Noise, include in BSM selection (how data work)
            z5p = complex(randn(1,N),randn(1,N))/sqrt(2);
            z6p = complex(randn(1,N),randn(1,N))/sqrt(2);
            z7p = complex(randn(1,N),randn(1,N))/sqrt(2);
            z8p = complex(randn(1,N),randn(1,N))/sqrt(2);

            
            %Psi- from SPDC 1
            a1H = sigma0*( cosh(r)*z1H + sinh(r)*conj(z2V) );
            a1V = sigma0*( cosh(r)*z1V - sinh(r)*conj(z2H) );
            a2H = sigma0*( cosh(r)*z2H - sinh(r)*conj(z1V) );
            a2V = sigma0*( cosh(r)*z2V + sinh(r)*conj(z1H) );
            
            %Psi- from SPDC 2
            a3H = sigma*( z3H*cosh(r) + conj(z4V)*sinh(r) );
            a3V = sigma*( z3V*cosh(r) - conj(z4H)*sinh(r) );
            a4H = sigma*( z4H*cosh(r) - conj(z3V)*sinh(r) );
            a4V = sigma*( z4V*cosh(r) + conj(z3H)*sinh(r) );
            
            %Beamsplitter for heralding
            
            b5H = (1/sqrt(2))*(a2H + a4H); %replace sigma
            b5V = (1/sqrt(2))*(a2V + a4V);
            b6H = (1/sqrt(2))*(a2H - a4H);
            b6V = (1/sqrt(2))*(a2V - a4V);
            
            theta = 0;
            %polarizers for heralding
            c5H = b5H*cos(theta)^2 + b5V*cos(theta)*sin(theta) + sigma*z5p*(1-cos(theta)^2) - sigma*z6p*cos(theta)*sin(theta);
            c5V = b5V*cos(theta)*sin(theta) + b5V*sin(theta)^2 - sigma*z5p*cos(theta)*sin(theta) + sigma*z6p*(1-sin(theta)^2);
            
            eta = pi/2;
            c6H = b6H*cos(eta)^2 + b6V*cos(eta)*sin(eta) + sigma*z7p*(1-cos(eta)^2) - sigma*z8p*cos(eta)*sin(eta);
            c6V = b6V*cos(eta)*sin(eta) + b6V*sin(eta)^2 - sigma*z7p*cos(eta)*sin(eta) + sigma*z8p*(1-sin(eta)^2);
            
            %Detection on 5 and 6
            D56 = (abs(c5H) > gamma2 | abs(c5V) > gamma2) & (abs(c6H) > gamma2 | abs(c6V) > gamma2);
            
            % Picking out only the elements that gave Bell State Measurement
            z1H_m = [z1H_m, (nonzeros(D56.*z1H))'];
            z1V_m = [z1V_m, (nonzeros(D56.*z1V))'];
            z2H_m = [z2H_m, (nonzeros(D56.*z2H))'];
            z2V_m = [z2V_m, (nonzeros(D56.*z2V))'];
            z3H_m = [z3H_m, (nonzeros(D56.*z3H))'];
            z3V_m = [z3V_m, (nonzeros(D56.*z3V))'];
            z4H_m = [z4H_m, (nonzeros(D56.*z4H))'];
            z4V_m = [z4V_m, (nonzeros(D56.*z4V))'];
            
            %PBS Noise.
            z1p_m = [z1p_m, (nonzeros(D56.*z1p))'];
            z2p_m = [z2p_m, (nonzeros(D56.*z2p))'];
            z3p_m = [z3p_m, (nonzeros(D56.*z3p))'];
            z4p_m = [z4p_m, (nonzeros(D56.*z4p))'];
            
            %Polarizer Noise, include in BSM selection 
            z5p_m = [z5p_m, (nonzeros(D56.*z5p))'];
            z6p_m = [z6p_m, (nonzeros(D56.*z6p))'];
            z7p_m = [z7p_m, (nonzeros(D56.*z7p))'];
            z8p_m = [z8p_m, (nonzeros(D56.*z8p))'];
           
            
            BSM = length(z1H_m);
            disp(BSM)

            end
            
            
            if isfile('BellStates_r-' + string(r) + '_g2-' + string(gamma2) + '_trial-' + string(trial) + '.mat')
                z1H_m =  [z1Hb z1H_m];
                z1V_m =  [z1Vb z1V_m];
                z2H_m =  [z2Hb z2H_m];
                z2V_m =  [z2Vb z2V_m];
                z3H_m =  [z3Hb z3H_m];
                z3V_m =  [z3Vb z3V_m];
                z4H_m =  [z4Hb z4H_m];
                z4V_m =  [z4Vb z4V_m];
                z1p_m =  [z1pb z1p_m];
                z2p_m =  [z2pb z2p_m];
                z3p_m =  [z3pb z3p_m];
                z4p_m =  [z4pb z4p_m];
                z5p_m =  [z5pb z5p_m];
                z6p_m =  [z6pb z6p_m];
                z7p_m =  [z7pb z7p_m];
                z8p_m =  [z8pb z8p_m];
                BSM_Count = length(z1H_m);
                disp(BSM_Count)
            save(filename, 'z1H_m', 'z1V_m', 'z2H_m', 'z2V_m','z3H_m','z3V_m', 'z4H_m','z4V_m','z1p_m', 'z2p_m', 'z3p_m', 'z4p_m', 'z5p_m', 'z6p_m', 'z7p_m', 'z8p_m', 'BSM_Count'  );
            else
            BSM_Count = length(z1H_m);
            save(filename, 'z1H_m', 'z1V_m', 'z2H_m', 'z2V_m','z3H_m','z3V_m', 'z4H_m','z4V_m','z1p_m', 'z2p_m', 'z3p_m', 'z4p_m', 'z5p_m', 'z6p_m', 'z7p_m', 'z8p_m', 'BSM_Count'  );
            disp(BSM_Count)
            disp("New trial")
            end


        end
    end
end

disp(filename)
cd '/Users/aishiguha/Documents/MATLAB/Code'
toc