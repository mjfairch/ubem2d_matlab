% Builds a NACA 4-digit series airfoil with unit chord length.
%
% Input:
% -----------
% code_str : 4-digit string (default = '2315') encoding the foil's shape.
%   The first digit is M, the second digit is P, and the last two digits 
%   are TH, defined by:
%
%   M = maximum camber as a percentage of the chord
%   P = distance of max camber from leading edge in tens of percent of chord
%   TH = thickenss of the airfoil as a percentage of chord
%
%   For example, 2315 means a cambered airfoil whose max camber is 2% of the
%   chord, located 30% chordwise from the leading edge, and whose thickness
%   is 15% of the chord.
%
% npan : integer (default = 50); number of panels along the surface.
%
% te_clamp : boolean, default = true, trailing edge clamp.  When false,
%   the trailing edge has a slight nonzero thickness.  When true, the foil
%   adjusted so that the trailing edge has exactly zero thickness.
% 
% uniform : boolean, default = false, distribution of points along chord.
%   When false, points are uniformly distributed along the chord.
%   When true, points are distributed more densly near leading edge.
%
% Output:
% --------
% foil : Airfoil2d object representing the 4-digit series NACA airfoil
%
% References:
% -----------
% https://en.wikipedia.org/wiki/NACA_airfoil
%
% Example:
% --------
% foil = foil_naca4('2315',50);
% plot(foil.x,foil.y,'ko-'); axis('equal'); grid on;
% fprintf('Leading edge at x=%.2g, y=%.2g\n',foil.x(foil.le),foil.y(foil.le));
%
% Author:   Michael J. Fairchild (Princeton University)
% Date:     2016-09-22
function foil = foil_naca4(code_str, npan, te_clamp, uniform)
if (nargin < 1) 
    code_str = '2315'; 
end
if (nargin < 2)
    npan = 50; 
end
if (nargin < 3)
    te_clamp = true;
end
if (nargin < 4) 
    uniform = false; 
end

if (length(code_str) ~= 4)
    error('NACA code must be 4 digits'); 
end
c = 1;                                  % chord length
m = .01*str2double(code_str(1));        % camber parameter
p = .1*str2double(code_str(2));         % location of max camber
th = 0.01*str2double(code_str(3:4));    % thickness

if (npan < 4)
    error('Must specify at least four panels');
end
if (m == 0 && p ~= 0)
    error('m == 0 and p ~= 0 are inconsistent');
end
if (m ~= 0 && (p == 0 || p == 1))
    error('when m ~= 0, then p must lie in the interval 0 < p < 1');
end

coefs = [-.1015, .2843, -.3516, -.1260, 0];
if (te_clamp)
    coefs(1) = -.1036;
end
sqrt_coef = .2969;

% Determine abcissa
n = npan + 1;
N = floor(n/2)+mod(n,2);
% Abcissa along upper surface
if (uniform)
    xu = linspace(c,0,N);
else
    t = linspace(.5*pi,0,N);
    xu = c*(1-cos(t));
end
lep = N-1;
% Abcissa along lower surface defined in terms of those along upper surface
xl = zeros(1,n-N);
if (mod(n,2) == 0) % n even
    dxu = diff(xu);
    for i=1:N-1
        xl(i) = xu(end-(i-1))-.5*dxu(end-(i-1));
    end
else
    for i=1:N-1
        xl(i) = xu(end-i);
    end
end
xl(end) = xu(1);

% Construct the thickness distribution along upper and lower surfaces
ytu = 5*th*c*(sqrt_coef*sqrt(xu/c) + polyval(coefs,xu/c));
ytl = 5*th*c*(sqrt_coef*sqrt(xl/c) + polyval(coefs,xl/c));

% Compute the upper and lower surfaces
if (m == 0) % symmetric (i.e. noncambered) foil
    x = cat(2,xu,xl);
    y = cat(2,ytu,-ytl);
else
    % Compute camber line and its derivative along upper and lower surfaces
    ycu = zeros(size(xu));
    dycudx = zeros(size(xu));
    ycl = zeros(size(xl));
    dycldx = zeros(size(xl));
    for i=1:length(ycu)
        if (xu(i) <= p*c)
            ycu(i) = m*xu(i)/(p*p)*(2*p-xu(i)/c);
            dycudx(i) = 2*m/(p*p)*(p-xu(i)/c);
        else
            ycu(i) = m*(c-xu(i))/((1-p)*(1-p))*(1+xu(i)/c-2*p);
            dycudx(i) = 2*m/((1-p)*(1-p))*(p-xu(i)/c);
        end
    end
    for i=1:length(ycl)
        if (xl(i) <= p*c)
            ycl(i) = m*xl(i)/(p*p)*(2*p-xl(i)/c);
            dycldx(i) = 2*m/(p*p)*(p-xl(i)/c);
        else
            ycl(i) = m*(c-xl(i))/((1-p)*(1-p))*(1+xl(i)/c-2*p);
            dycldx(i) = 2*m/((1-p)*(1-p))*(p-xl(i)/c);
        end
    end
    % Upper and lower thickness distributions
    thu = atan(dycudx);
    thl = atan(dycldx);
    % x and y coordinates along airfoil boundary
    xU = xu - ytu.*sin(thu);
    yU = ycu + ytu.*cos(thu);
    xL = xl + ytl.*sin(thl);
    yL = ycl - ytl.*cos(thl);
    x = cat(2,xU,xL);
    y = cat(2,yU,yL);
end

foil = Airfoil2d(x,y,lep+1);