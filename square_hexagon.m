set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
% Jad Zahar,   jzahar@student.ethz.ch
%%%%%%%% INTERNSHIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function builds the visual panel
app

function app 
%build panel
panel = uifigure("Position",[400 300 500 300]);

%build J' slider
sld = uislider(panel,'slider','Limits',[0 1.2]);
    sld_label =uilabel(panel);
    sld_label.Interpreter = 'latex';
    sld_label.Text =  "$J_1'/J_1$";
    sld_label.FontSize = 20;
    sld_label.Position = [35 50 100 100];

%build method switch
swt = uiswitch(panel,"Position",[100 170 300 25]);
swt.Items = ["LSWT","PWT"];

%biuld ""keep figure" checkbox
cbx = uicheckbox(panel, "Text", "Keep figures", "Position",[350 100 300 25]);

%declare callback functions
sld.ValueChangedFcn = {@run, sld, swt, cbx};
swt.ValueChangedFcn = {@run, sld, swt, cbx};

end



function run(~, ~, slider, my_switch, checkbox)
%runs the program, takes UI elements values as inputs

% Initialize---------------------------------------------------------------

if checkbox.Value %if box is checked, keep figures
    fig = figure;
else
    clf;
end

%look at method input (PWT or LSWT)
method = my_switch.Value;


%bond values
J = 1;
Jp = slider.Value; %slider.Value
S = .5; % spin number
M = .5;%value for parameter M, this will only appear in PWT and as factor of S


%lattice geometry
lp = 1; % length of J1 bonds
l = 1/sqrt(2); %1/sqrt(2) length of J1' bonds
L = 2*l + lp; %length of primitive vectors
d = 2*(sqrt(2)*l + lp/sqrt(2));


%Path on the 1st BZ, TODO : make fineness of path a parameter, 300pt fornow
kx = (pi/L)*[0:0.01:1, ones([1 100]), 1:-0.01:0];
ky = (pi/L)*[0:0.01:1, 1:-0.01:0, zeros([1,100])];


% Failure
fail = false; % set to true if PWT method fails



% Calculations ------------------------------------------------------------

if method == "LSWT"
% with linear spin wave theory

    % Energies, will contain energy levels for each *k*
    E = zeros([5, length(kx)]);
    
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
    
        E(:,n) = eig((A+B)*(A-B));
    
    end

%  -   -   -   -   -   -   -  -   -   -   -   -   -   -   -   -   -   -   -

elseif method == "PWT"
% with plaquette wave theory (+lswt for free spin)

    % set C
    C = 0.5*Jp*sqrt(2*S)*sqrt(M)/sqrt(3);
    

    E = zeros([7, length(kx)]); %will contain energy levels for each *k*


    for n = 1:1:length(kx) 
        x = kx(n);
        y = ky(n);

        % build A
            A_pw = zeros(7);
            A_pw(5,1) = C*(-1 - f(d*x) - f(d*y)+f(d*(x + y)));
            A_pw(6,1) = C*(+f(d*x) - f(d*y))/sqrt(2);
            A_pw(7,1) = C*(-1 - f(d*(x+y)))/sqrt(2);
            A_pw = A_pw + A_pw';
            A_pw = A_pw + diag(-1*[0 J 0 0 J 0 0])+2*diag([0 J J J J J J]);
            A_pw = A_pw*0.5;
         % build B
            B_pw =   zeros(7);
            B_pw(1,2)=-1 - f(d*x)-f(d*y)+f(d*(x+y));
            B_pw(1,3)=(f(d*x)-f(d*y))/sqrt(2);
            B_pw(1,4)=(-1-f(d*(x+y)))/sqrt(2);
            B_pw = B_pw+B_pw';
            B_pw = B_pw*C*0.5;
    
        % find energy levels
        E(:,n) = eig((A_pw+B_pw)*(A_pw-B_pw));

    end

    % Check results
    fail = (any(imag(E) > 20*sqrt(eps), 'all') || any(real(E) < -100*sqrt(eps), 'all'));
    %failed if energy is complex or negative

end

% Figure-------------------------------------------------------------------
for m = [1:1:size(E,1)]
    plot(sqrt(abs(real(E(m,:)))), 'b.');
    hold on
end

%figure visuals
ylabel("$E$", "Rotation",0)
xline(0, '--')
xline(100, '--')
xline(200, '--')
xline(302, '--')

xticks([0 100 200 300])
xticklabels({'$\Gamma$','K','M','$\Gamma$'})

ylim([0 1.15])

title(append('$J_1'' = $' , num2str(Jp), '$J_1$'), 'FontSize', 20);



% if failed : make figure red
if fail
    title(append('$J_1'' = $' , num2str(Jp)) + "$J_1$ - Invalid energy levels obtained !", "FontSize",20, 'Color','r')
end

end


function result = f(a)
    result = -1*exp(1i*a);
end
