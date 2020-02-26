% format compact; close all; clear; clc; 
% Build Rotations and Jacobian from Denavit Hartenberg Parameters
clear; clc; format compact;

oo_constants ;

%% DH Params (given)

syms q1_ q1dot_ t_

q_    = [ q1_   ].' ;
qdot_ = [ q1dot_].' ;

qconfig = [ pi/2]';
qr_ = q_ + qconfig ;

q     = [ q_  ; qdot_ ]   ;
n = length(q_)  ;% number of links
Nq = length(q) ;% number of states

params.grav = grav  ;
params.DH_n = 1     ;
params.Nq   = Nq    ;
params.q_   = q_    ;
params.qdot_= qdot_ ;
params.q    = q     ;
params.DH_q = sym('DH_q',[1,n])  ;
params.DH_l = sym('DH_l',[1,n])  ;
params.lc   = sym('lc',[1,n])    ;
params.m    = sym('m',[1,n])     ;
params.moi  = sym('moi',[3,3,n]) ;

%% LINK 1 % stance leg, toe to knee
link  = 1 ;
qi    = qr_(link);% q1 ;%pi / 2 + q1              ;% symbolic rotation

% individual link properties
mi  = m_knee + m_hip + m_knee; %************** ALTERED WITH MASS OF OTHER KNEE

% total link length
li  = l_leg ; 

% center of mass of link
lci = ( m_hip*l_leg + m_knee*(l_leg-l_knee) ) / ( m_hip + m_knee ) ;% leg COM offset from toe

% moment of inertia about link center of mass   
ji  = m_hip*( l_leg - lci )^2 ...
        + m_knee*( (l_leg-lci) - l_knee )^2 ...
      + J_wheel_COM ...
      + J_knee_COM ;% j about leg com

moi_i = diag( [0 0 ji] ) ;% link moment of inertia about its com

% store data in proper locations
params.DH_q(link)    = qi    ;
params.DH_l(link)    = li    ;
params.lc(link)      = lci   ;
params.m(link)       = mi    ;
params.moi(:,:,link) = moi_i ;


%% Rotation Matrix
% state = [q1].';

tic
[ B_MATRIX , C_MATRIX , G_MATRIX , OTHER ] = symbolic_dynamic_rotation( params ) ;
t_create_dynamics = toc

% create a bunch of useful functions based on the output
matlabFunction( OTHER.TRANSFORMATIONS , 'File', 'symbolic_functions/transformations1', 'Vars', {q}) ;
matlabFunction( OTHER.TRANSFORMATIONS , 'File', 'symbolic_functions/transformations_short1', 'Vars', {q_}) ;
matlabFunction( B_MATRIX , 'File', 'symbolic_functions/B_matrix1', 'Vars', {q}) ;
matlabFunction( C_MATRIX , 'File', 'symbolic_functions/C_matrix1', 'Vars', {q}) ;
matlabFunction( G_MATRIX , 'File', 'symbolic_functions/G_matrix1', 'Vars', {q}) ;
matlabFunction( G_MATRIX , 'File', 'symbolic_functions/G_matrix_short1', 'Vars', {q_}) ;


%% INTERNAL CONTROL

syms E_d % desired energy level (for control action) 

Ek = 1/2*params.m*params.lc*q1dot_.^2 ;
Ep = params.m*grav*params.lc*cos(q1_) ;
E = Ek + Ep ;
Ebar = E - E_d ;


% m_knee*lci^2 + (m_actuator+m_wheel+m_knee)*l_leg^2 + J_knee_COM+J_wheel_COM - B_MATRIX
syms flywheel_control_on cp_control_on
controllers = [ flywheel_control_on cp_control_on ]

% tau = -control_on*K_control*q1dot_.^2*Ebar ;
tau = -sign(q1dot_)*flywheel_control_on*K_control*q1dot_.^2*Ebar ;% negative tau adds energy when spinning negative

% syms tau

qddot_ = simplify( (B_MATRIX - J_wheel_COM) \ (  tau  - C_MATRIX*qdot_ - G_MATRIX ) ) 

% Edot = simplify( params.m*params.lc*q1dot_*qddot_ - params.m*grav*params.lc*sin(q1_)*q1dot_ )
% qddot_ = B_MATRIX \ ( [ 0 ; tau_sym ] - C_MATRIX*qdot_ - G_MATRIX ) ;

% tau_sym = sym('tau_sym', [ params.DH_n-1 , 1 ] ) ;
% qddot_ = simplify( B_MATRIX \ ( [ 0 ; tau_sym ] - C_MATRIX*qdot_ - G_MATRIX ) ) ;

% t_create_dynamicsa = toc;
% qddot_ = simplify(qddot_) ;% 27 minutes
% t_simplify_dynamics = toc - t_create_dynamicsa ;
% fprintf('time to simplify dynamics: %1.3f sec\n',toc) 

t_create_dynamicsa = toc ;
qddot_ = simplify( qddot_ ) ;% 27 minutes
t_simplify_dynamics = toc - t_create_dynamicsa ;
fprintf('time to simplify dynamics: %1.3f sec\n',toc) 

qdot = [qdot_;qddot_] ;

t_create_dynamics = toc;
fprintf('time to create dynamics: %1.3f sec\n',t_create_dynamics) 

% matlabFunction( qdot , 'File', 'symbolic_functions/dynamics1', 'Vars', {t_,q,tau_sym}) ;% 14 min w/ qddot simplify, 24 min w/out
% t_create_d = toc ;
% fprintf('time to create dynamics function: %1.3f sec\n',t_create_d) 
    
matlabFunction( qdot , 'File', 'symbolic_functions/dynamics1', 'Vars', {t_,q,E_d,controllers}) ;% 14 min w/ qddot simplify, 24 min w/out
t_create_d = toc ;
fprintf('time to create dynamics function: %1.3f sec\n',t_create_d) 
    
%% Steady State

syms slope_m 

q1dot_ss = sqrt( 4*grav/params.lc * sin( slope_m ) + 2*E_d / ( params.m * params.lc.^2 ) ) ; 

q2_ss = acos( sqrt( 2*E_d / params.m ) ./ ( q1dot_ss.*params.lc ) ) ;

q1_ss = -q2_ss / 2 ;

q_ss = [ q1_ss ; q2_ss ; q1dot_ss ; 0 ] ;

matlabFunction( q_ss , 'File', 'symbolic_functions/create_ss', 'Vars', { slope_m , E_d }) ;% 

% q_ss = create_ss( slope_m , E_d )



%% 

syms q1dot_recreate q2_ q2dot_ 'real'

q_config = [ q1_ ; q2_ ] ;

q1dot_recreate = sqrt( ( E_d - Ep ) / ( 1/2*params.m*params.lc ) ) ;

matlabFunction( q1dot_recreate , 'File', 'symbolic_functions/recreate_q1dot', 'Vars', { q_ , E_d }) ;% 


%%
syms xn hxn l_swing 'real'

[ q1_IK , q2_IK ] = get_IK( xn, hxn , l_swing ) ; 

q1dot_IK = recreate_q1dot( q1_IK , E_d ) ;% 

q_minus = [ q1_IK ; q2_IK ; q1dot_IK ; 0 ] ;

q_plus = reset_map2( q_minus ) ;


[Em , PEm , KEm ] = get_energy( q_minus , 0 ) ;
[Ep , PEp , KEp ] = get_energy( q_plus  , 0 ) ;

% PEp = PEp + (m_knee+m_hip+m_knee)*grav*hxn ;

dE  = simplify(  Ep -  Em ) ;
dPE = simplify( PEp - PEm ) ;
dKE = simplify( KEp - KEm ) ;

matlabFunction( dE , dPE , dKE , 'File', 'symbolic_functions/get_energy_change', 'Vars', { xn , hxn , l_swing , E_d } ) ;


%% IMPACT FUNCTION

l_toe = sqrt( l_swing.^2 + l_leg.^2 - 2*l_swing*l_leg*cos(q2_) ) ;

cos_phi = ( l_swing.^2 - l_leg^2 - l_toe.^2 ) / ( -2*l_leg*l_toe ) ;

sin_phi = sin(q2_) * l_swing / l_toe ; 

psi_real = atan2( sin_phi , cos_phi ) - pi/2 - q1_ ;

matlabFunction( psi_real , 'File', 'symbolic_functions/get_psi_real', 'Vars', { q_config ,  l_swing } ) ;









