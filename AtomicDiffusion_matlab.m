% Computational solutions for atomic diffusion through plates, hollow cylinders, and spheres.
% Written by Gavin Latham and optimized by Caleb Pennell on 5/22/23.
% Written for use at ANCORP - this was a side project with no peer review process. Use at own discretion.

% The material science used in this code is based off of "The Mathematics of Diffusion" by J. Crank.

clc; clear; hold off; clf('reset');
fprintf('This code is to be used for modeling hydrogen concentration in various geometries made of 316 S.S. after a brakeout process\n');
D0 = 4.7e-7; % Temperature-independent preexponential diffusion coefficient (m^2/s)
k = 8.617e-5; % Boltzmann constant (eV/K)
celsius = input('What is the bakeout temperature in celcius? ');
T = celsius + 273.15; % Temperature (K)
t = input('How many minutes are you baking for? ') * 60; % time (s)
Co = .04; % Initial concentration of hydrogen in 316 S.S. (mol/m^3)
C1 = 5e-7; % Constant surface concentration (moles/m^3)

Q = .48; % Activation energy of diffusion (eV)
D = D0 * exp(-Q / (k * T)); % Diffusion coefficient (m^2/s) -> citation below

%E Hashimoto, T Kino,
%Hydrogen permeation through type 316 stainless steels and ferritic steel for a fusion reactor,
%Journal of Nuclear Materials,Volumes 133–134,1985,Pages 289-291,ISSN 0022-3115,https://doi.org/10.1016/0022-3115(85)90153-9.

check = 0;
while check == 0
  fprintf('\n');
  shape = input('What general shape is your part? Plate, tube, or sphere? ' , 's');
  shape = lower(shape);
  if strcmp(shape, "plate") == 1 || strcmp(shape, "sphere") == 1 || strcmp(shape, "tube") == 1
    check = 1;
  else
    fprintf('This is not a valid shape, please choose from the specified shapes.\n');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(shape, "plate") == 1

  h = .00001; % Step size (m)
  dt = ((h^2) / (2 * D)); % Time step dependent on physical step size, h
  fprintf('\n');
  x = input('What is the thickness of the plate in inches? ') / 39.37;
  x = [0:h:x];
  C = zeros(1,length(x));
  length = input('What is the length of the plate in inches? This is only necessary for the total moles of H at the end. ')/39.37;
  height = input('What is the height of the plate in inches? This is only necessary for the total moles of H at the end. ')/39.37;
  C(2:end - 1) = Co;  % Concentration profile (mol/m^3)
                      % Each element represents the concentration at a step, so boundary is index 1, thus 0
  tic;

  N = t / dt; % Total steps taken
  ti = 0; % Time at step i
  D2C = 0; % Fick's second law
  for i = 1:1:N
    ti = ti + dt;
    D2C = (C(3:end) - (2*C(2:end-1)) + C(1:end-2))/h^2; % Numerical solution to Fick's second law
    C(2:end - 1) = C(2:end - 1) + dt * D * D2C;
  end

  fprintf("\n");
  time_elapsed = toc

  plot(x,C)
  xlabel("Distance in (m)");
  ylabel("Hydrogen Concentration (mol/m^3)");
  str = sprintf("Hydrogen Concentration After %g Minute Bake at %gC", t/60, celsius);
  str2 = sprintf('%g minute(s)',t / 60);
  title(str);

  moles_of_H = sum((length.*height).*C(2:end - 1).*h); % (mol)
  fprintf('\nThere are %g moles of Hydrogen left in the steel plate\n',moles_of_H);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(shape, "tube") == 1

  b = input('What is the outer diameter in inches? ')/(2 * 39.37); % Outer radius (m) - this is changeable
  a = b - (input('What is the thickness of the tubing wall in inches? ') / 39.37); % Inner diameter (m) - this is changeable
  h = input('What is the height of your tube in inches? ') / 39.37;

  resolution = .01; % Resolution of the zero, so down to hundredreths, rounded
  max_alpha = 100000; % Max value of zeros (alphas) willing to look to
  max_i = max_alpha / resolution; % Amount of steps in the function to look for zeros
  base_i = 1; % Multiple used to efficiently find roots, eg. if first alpha is 1963.17
              % Then each alpha will be almost exactly a multiple of this, like 1960
              % So skip forward just under this amount

  x = 0:resolution:max_alpha + 1; % Physical profile through thickness (m)

  fprintf("\nStarting... please wait, this may take a little: \n");
  tic;

  ax = a.*x; % Precalculated for Crank's equation
  bx = b.*x;

  Uo = (besselj(0,ax).*bessely(0,bx)) - (besselj(0,bx).*bessely(0,ax)); % Used to determine alpha values, Crank's

  alpha = [];

  i = base_i / resolution; % i step size finding alphas (zeros of the function)
  while i < max_i
    if ((Uo(i) > 0 && Uo(i + 1) < 0 ) || (Uo(i) < 0 && Uo(i + 1) > 0)) %finding roots of Uo which is alpha

      if length(alpha) == 0 % First aplpha determines following solution pattern
          base_i = (i + 1) * resolution - 5;
      end

      alpha(end + 1) = (i + 1) * resolution; % Add zero to solution vector
      i = i + base_i / resolution;

    else
      i = i + 1;
    end

    if length(alpha) == 100 % 100 zeros is sufficient for accurate profile plot
      break
    end
  end


  step = ((b-a) / (length(alpha) - 1)); % Map alphas onto r
  r = [a:step:b]; % Thickness profile of tube, radial (m)

  summ = zeros(1, length(alpha)); % Empty summation vector
  for j = 1:1:length(alpha)

    aj = a.*alpha(j);
    bj = b.*alpha(j);
    rj = r.*alpha(j);

    besj_aj = besselj(0,aj);
    besj_bj = besselj(0,bj);

    U2 = besselj(0,rj).*bessely(0,bj) - besj_bj.*bessely(0,rj); % Equation pulled from Crank's text
    func = (besj_aj.*U2)./(besj_aj + besj_bj).*exp(-D * (alpha(j)^2) * t); % Same equation, two lines
    summ = summ + func;
  end

  C = (C1-Co).*(1 - (pi.*summ)) + Co; % Final piece of Crank's equation from text

  plot(r,C);
  hold on
  xlabel("Distance into wall (m)");
  ylabel("Hydrogen Concentration (mol/m^3)");
  str = sprintf("Hydrogen Concentration After Bake at %gC", celsius);
  str2 = sprintf('%g minute(s)',t / 60);
  title(str);
  legend(str2,'Location','northeastoutside')
  hold off

  fprintf('Done!\n\n');
  time_elapsed = toc

  dx = 1 / length(C);
  integral = sum((C.*dx),2); % (mol/m^2)

  % This is the integral multiplied by the cross-section area cutting
  % parallel to the tube face. It is also multiplied by the height.
  moles_of_H = integral * pi * (b^2 - a^2) * h; % (mol)
  fprintf('\nThere are %g moles of Hydrogen left in the steel tube\n',moles_of_H);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(shape, "sphere") == 1

  b = input('What is the outer radius in inches? ') / 39.37; % Outer radius (m)
  a = b - (input('What is the wall thickness in inches? ') / 39.37); % Inner radius (m)

  n_limit = 1000; % This is the amount of summations to do. The greater, the more precise the function.
  r_step = .00001; % Profile step (m)
  r = a:r_step:b;

  summ = zeros(1, length(r));
  for n = 1:1:n_limit

    % Equation pulled from Crank's text
    func = ((b.*cos(n.*pi) - a)./n).*(sin((n.*pi.*(r - a))./(b - a))).*exp((-D * (n^2) * (pi^2) * t) / ((b - a)^2));
    summ = summ + func;
  end

  C = (C1 - Co).*(1 + ((2./(pi.*r)).*summ)) + Co; % Final equation fromo Crank's

  plot(r,C)

  hold on

  xlim([a, b])
  ylim([0, 1.1 * max(C)])
  xlabel("Distance into wall (m)");
  ylabel("Hydrogen Concentration (mol/m^3)");
  str = sprintf("Hydrogen Concentration After Bake at %gC", celsius);
  str2 = sprintf('%g minute(s)',t / 60);
  title(str);
  legend(str2,'Location','northeastoutside')
  hold off

  SA = 4.*pi.*(r.^2); % (m^2)
  moles_of_H = sum(SA.*C.*r_step); % (mol)
  fprintf('\n');
  fprintf('There are %g moles of Hydrogen left in the steel sphere\n',moles_of_H);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


