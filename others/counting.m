set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);
% Jad Zahar,   jad.zahar@epfl.ch
%This will give the quantum energy (E - Eclass) as well as the number of
%each magnon type for the ground state. All for multiple descending values
%of J' starting from 1 =J.
%%%%%%%% INTERNSHIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define variables and matrices-------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%% PARAMETER VALUES
Jp = [.99:-.005:.2];
%%%%%%%%%%%


n = 60; % #of cells = N = n * n
%main parameters
J = 1;
S = 20;

%lattice geometry
lp = 1; % length of J1 bonds
l = 1/sqrt(2); %1/sqrt(2) length of J1' bonds
L = 2*l + lp; %length of primitive vectors(do not edit)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%% build and diagonalize
H = zeros(1, length(Jp));
%will contant total energy values
H_norm = H;
%will contant energy value per cell
H_norm_q = H;
%will contain normalised (E - Eclassical)

N = n*n;%  #of cells = N = n * n 

% N values of k in the first BZ
delta = (2*pi/L)/n;
% sample points in the 1st BZ
kx = [-(n/2)*delta + delta/2: delta: (n/2)*delta - delta/2];
ky = [-(n/2)*delta + delta/2: delta: (n/2)*delta - delta/2];
W = zeros([5, length(kx), length(ky)]); %will contain eigenvalues for each case
magnons = zeros(5, length(Jp)); %will contain the number of each type of magnon for each J' value



for j = 1:1:length(H)
%we will solve for different values of Jp
    
    
    %build A
    A = diag(S*0.5*[2*J+Jp(j), 2*J+Jp(j), 2*J+Jp(j), 2*J+Jp(j), 4*Jp(j)]);
    
    for iterX = 1:1:length(kx)
        for iterY = 1:1:length(ky)
            x = kx(iterX);
            y = ky(iterY);
            %build B
            B =   S*0.5*[[  0          f(-x*lp)*J     f(y*lp)*J          0         f(x*l - y*l)*Jp(j)]
                        [f(x*lp)*J         0                0        f(y*lp)*J     f(-x*l - y*l)*Jp(j)]
                        [f(-y*lp)*J        0                0        f(-x*lp)*J    f(x*l + y*l)*Jp(j)]
                        [    0          f(-y*lp)*J      f(x*lp)*J        0         f(-x*l + y*l)*Jp(j)]
                        [f(-1*(x*l - y*l))*Jp(j)  f(-1*(-x*l - y*l))*Jp(j)  f(-1*(x*l + y*l))*Jp(j)  f(-1*(-x*l + y*l))*Jp(j) 0]];
        
            %find energy levels
            W(:,iterX, iterY) = eig((A+B)*(A-B));



            %find bogoliubov transformation matrix
            [V,~] = eig([A B; -B -A]);
            V_inv = inv(V);
            %now find the coefficients Bk
            % where a'a --bogo--> (Ai)a'a + (Bi)aa'+... (e.c. i=1..5)
            for X = 1:1:5
                new_M = V_inv*accumarray([X X], 1,[10, 10])*V;
                    %for this apply the bogo transf to a matrix that only holds
                    %an a'a term
                magnons(X,j) = (magnons(X,j)+ new_M(6,6)+new_M(7,7)+new_M(8,8)+new_M(9,9)+new_M(10,10))/N;

            end
       
        
        end
    end

    
    
    %%% compute total energy value
    H(j) = N*(-4*S*(S+1)*(J+Jp(j))) + N*(-2*S*(J+2*Jp(j))) + 0.5*sum(sum(sum(sqrt(W))));
    H_norm(j) = H(j)/N;
    H_norm_q(j) = (H(j) - N*(-4*S*(S+1)*(J+Jp(j))))/N;
end

%% plot results
figure
plot(Jp, real(magnons(1,:)), 'k-','HandleVisibility','off');
hold on
plot(Jp, real(magnons(2,:)),'k-', 'DisplayName',"$A,B,D,E$ bosons");
plot(Jp, real(magnons(3,:)),'k-','HandleVisibility','off');
plot(Jp, real(magnons(4,:)),'k-','HandleVisibility','off');
plot(Jp, real(magnons(5,:)),'r-', "DisplayName","$C$ boson");
grid on
xlabel("$J_1'/J_1$")
ylabel("n° of bosons")
legend("Location","best");





%________________check that W is a real spectrum
% hold off
% for a = [1 2 3 4 5]
%     for b = 1:1:60
%         plot(b, real(squeeze(W(a, 25, b))), 'r.');
%         hold on
%     end
% end
%_________________________________________



function result = f(a)
    result = -1*exp(1i*a);
end
