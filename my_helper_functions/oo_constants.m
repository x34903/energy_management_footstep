FAKE_HEBI = 1 ;

PROP_METHOD = 1 ;% ode45
% PROP_METHOD = 2 ;% Runge Kutta
% PROP_METHOD = 3 ;% Trapezoidal Rule 
% PROP_METHOD = 4 ;% Euler Integration
params.PROP_METHOD = PROP_METHOD ;

params.min_dtn = 0.01 ;% necessary within integration methods

qconfig = [ pi/2 pi 0 0 ]';
grav = 10 ;

%%

% LEG COM properties
m_knee     = 0.25 ;
l_knee     = 0    ;% from hip
J_knee_COM = 0    ;

% Wheel properties
m_wheel     = 0.25 ;
r_wheel     = 0.1 ;
J_wheel_COM = m_wheel*r_wheel^2 ;

% hip properties
m_actuator = 0.25 ;
l_leg      = 1    ;% from toe
m_hip      = m_actuator + m_wheel ;

m_cm = 2*m_knee + m_wheel + m_actuator ;
l_cm = l_leg ;

% total mass = 1.0, all located at hip


%% ENERGY CONTROLLER


K_control = 1e0 ;



