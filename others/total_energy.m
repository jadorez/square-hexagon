set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);
% Jad Zahar,   jad.zahar@epfl.ch
%%%%%%%% INTERNSHIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define variables and matrices-------------------------------------------
%--------------------------------------------------------------------------
n = [10:10:80]; % #of cells = N = n * n
%please input even values of n

%main parameters
J = 1;
Jp = 1;
S = .5;

%lattice geometry
lp = 1; % length of J1 bonds
l = 1/sqrt(2); %1/sqrt(2) length of J1' bonds
L = 2*l + lp; %length of primitive vectors(do not edit)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%% build and diagonalize
H = zeros(1, length(n));
%will contant total energy values
H_norm = H;
%will contant energy value per cell


for j = 1:1:length(H)
%we want to check for convergence so iterate for multiple grid-sizes

    N = n(j)*n(j);%  #of cells = N = n * n 

    % N values of k in the first BZ
    delta = (2*pi/L)/n(j);

    kx = [-(n(j)/2)*delta + delta/2: delta: (n(j)/2)*delta - delta/2];
    ky = [-(n(j)/2)*delta + delta/2: delta: (n(j)/2)*delta - delta/2];
    
    
    
    W = zeros([5, length(kx), length(ky)]); %will contain eigenvalues for each *k*
    
    
    
    %build A
    A = diag(S*0.5*[2*J+Jp, 2*J+Jp, 2*J+Jp, 2*J+Jp, 4*Jp]);
    
    for iterX = 1:1:length(kx)
        for iterY = 1:1:length(ky)
            x = kx(iterX);
            y = ky(iterY);
            %build B
            B =   S*0.5*[[  0          f(-x*lp)*J     f(y*lp)*J          0         f(x*l - y*l)*Jp]
                        [f(x*lp)*J         0                0        f(y*lp)*J     f(-x*l - y*l)*Jp]
                        [f(-y*lp)*J        0                0        f(-x*lp)*J    f(x*l + y*l)*Jp]
                        [    0          f(-y*lp)*J      f(x*lp)*J        0         f(-x*l + y*l)*Jp]
                        [f(-1*(x*l - y*l))*Jp  f(-1*(-x*l - y*l))*Jp  f(-1*(x*l + y*l))*Jp  f(-1*(-x*l + y*l))*Jp 0]];
        
            %find energy levels
            W(:,iterX, iterY) = eig((A+B)*(A-B));
        
        end
    end

    
    
    %%% compute total energy value
    H(j) = N*(-4*S*(S+1)*(J+Jp)) + N*(-2*S*(J+2*Jp)) + 0.5*sum(sum(sum(sqrt(W))));
    H_norm(j) = H(j)/N;
end


%% plot result
figure
plot(real(H_norm))



function result = f(a)
    result = -1*exp(1i*a);
end
