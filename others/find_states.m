set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);
% Jad Zahar,   jad.zahar@epfl.ch
%%%%%%%% INTERNSHIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jp = .7;

%fixed k value
kx = (pi/L);
ky = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define variables and matrices


method = 2;% 1.diagonalize I-*M
           % 2.diagonalise (A+B)*(A-B)

%bonds
J = 1;
S = .5;

%lattice geometry
lp = 1; % length of J1 bonds
l = 1/sqrt(2); %1/sqrt(2) length of J1' bonds
L = 2*l + lp; %length of primitive vectors

%% build and diagonalize

if method == 1
    E = zeros([10, length(kx)]); %will contain energy levels for each *k*
elseif method ==2
    E = zeros([5, length(kx)]); %will contain energy levels for each *k*
end


%build A
A = diag(S*0.5*[2*J+Jp, 2*J+Jp, 2*J+Jp, 2*J+Jp, 4*Jp]);

for n = 1:1:length(kx) 
    x = kx(n);
    y = ky(n);
    %build B
    B =   S*0.5*[[  0          f(-x*lp)*J     f(y*lp)*J          0         f(x*l - y*l)*Jp]
                [f(x*lp)*J         0                0        f(y*lp)*J     f(-x*l - y*l)*Jp]
                [f(-y*lp)*J        0                0        f(-x*lp)*J    f(x*l + y*l)*Jp]
                [    0          f(-y*lp)*J      f(x*lp)*J        0         f(-x*l + y*l)*Jp]
                [f(-1*(x*l - y*l))*Jp  f(-1*(-x*l - y*l))*Jp  f(-1*(x*l + y*l))*Jp  f(-1*(-x*l + y*l))*Jp 0]];
    assert(ishermitian(B));%check that B is hermitian

    %find energy levels
    if method ==1
        
        M = [A B
           -B -A];
        E(:,n) = eig(M);
    elseif method ==2
        E(:,n) = eig((A+B)*(A-B));
    end
    
    [V,~] = eig((A+B)*(A-B));

end



%in order to analyse the data more easily when using method 2, create
if method == 2
    trueE = sqrt(E);
end

disp("Energy values at $\vec{k}=(\pi/L,0)$ are :")
disp(real(trueE))
disp("The contribution of each original HP boson to the corresponding bogoliubov boson is in matrix V. (qualitative, this is not really the composition)")
disp(V)
disp("You can find which band is which column by looking at the order of the energies above.")


function result = f(a)
    result = -1*exp(1i*a);
end
