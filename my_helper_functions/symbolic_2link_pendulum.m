% format compact; close all; clear; clc; 
% Build Rotations and Jacobian from Denavit Hartenberg Parameters
clear; clc; format compact;

%
oo_constants ;

%% DH Params (given)

syms q1_ q2_ q1dot_ q2dot_ t_

q_    = [ q1_    q2_    ].' ;
qdot_ = [ q1dot_ q2dot_ ].' ;

qconfig = [ pi/2 pi ]';
qr_ = q_ + qconfig ;

q     = [ q_  ; qdot_ ]   ;
n = length(q_)  ;% number of links
Nq = length(q) ;% number of states

params.grav = grav  ;
params.DH_n = 2     ;
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
mi  = m_knee + m_hip ; 

% total link length
li  = l_leg ; 

% center of mass of link
lci = ( m_hip*l_leg + m_knee*(l_leg-l_knee) ) / ( m_hip + m_knee ) ;% leg COM offset from toe

% moment of inertia about link center of mass   
ji  = m_hip*( l_leg - lci )^2 ...
      + m_knee*( (l_leg-lci) - l_knee )^2 ...
      + J_wheel_COM*0 ...
      + J_knee_COM ;% j about leg com

moi_i = diag( [0 0 ji] ) ;% link moment of inertia about its com

% store data in proper locations
params.DH_q(link)    = qi    ;
params.DH_l(link)    = li    ;
params.lc(link)      = lci   ;
params.m(link)       = mi    ;
params.moi(:,:,link) = moi_i ;


%% LINK 2 % stance leg, knee to hip
link  = link+1 ;
qi    = qr_(link);% %pi / 2 + q1              ;% symbolic rotation

% individual link properties
mi  = m_knee ;% m_thigh ; 

% total link length
li  = l_leg ;% l_thigh ;% length of thigh

% center of mass of link
lci = l_knee ;%  ;% knee COM offset from hip

% moment of inertia about link center of mass   
ji  = mi*lci^2 + J_knee_COM ;%  J_COM_thigh ;         % j about thigh com
 
moi_i = diag( [0 0 ji] ) ;% link moment of inertia about its com


% store data in proper locations
params.DH_q(link)    = qi ;
params.DH_l(link)    = li ;
params.lc(link)      = lci  ;
params.m(link)       = mi    ;
params.moi(:,:,link) = moi_i ;


%% Dynamics

tic
[ B_MATRIX , C_MATRIX , G_MATRIX , OTHER ] = symbolic_dynamic_rotation( params ) ;
t_create_dynamics = toc

% create a bunch of useful functions based on the output
matlabFunction( OTHER.TRANSFORMATIONS , 'File', 'symbolic_functions/transformations2', 'Vars', {q}) ;
matlabFunction( OTHER.TRANSFORMATIONS , 'File', 'symbolic_functions/transformations_short2', 'Vars', {q_}) ;
matlabFunction( B_MATRIX , 'File', 'symbolic_functions/B_matrix2', 'Vars', {q}) ;
matlabFunction( C_MATRIX , 'File', 'symbolic_functions/C_matrix2', 'Vars', {q}) ;
matlabFunction( G_MATRIX , 'File', 'symbolic_functions/G_matrix2', 'Vars', {q}) ;
matlabFunction( G_MATRIX , 'File', 'symbolic_functions/G_matrix_short2', 'Vars', {q_}) ;



tau_sym = sym('tau_sym', [ params.DH_n-1 , 1 ] ) ;
qddot_ = B_MATRIX \ ( [ 0 ; tau_sym ] - C_MATRIX*qdot_ - G_MATRIX ) ;

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

matlabFunction( qdot , 'File', 'symbolic_functions/dynamics2', 'Vars', {t_,q,tau_sym}) ;% 14 min w/ qddot simplify, 24 min w/out
t_create_d = toc ;
fprintf('time to create dynamics function: %1.3f sec\n',t_create_d) 
    
    
%% LINEARIZE DYNAMICS

%{
qdot_short = [ qdot(2) ; qdot(4) ] ;
q_short = [ q2_ q2dot_ ].';
dfdq = jacobian(qdot_short,q_short) ;% 3 sec
dfdu = jacobian(qdot_short,tau_sym) ;% 2 sec

% tic
% % XX min w/out dfdq,dfdu simplify!!   % XX min w/ dfdq,dfdu simplify!!
matlabFunction( dfdq , dfdu , 'File', 'symbolic_functions/linearize2', 'Vars', {q,tau_sym}) ;
% t_create_lin = toc ;
% fprintf('time to create linearization function: %1.3f sec\n',t_create_lin)     
    
save symbolic2.mat  

%% CREATE IMPACT FUNCTION

clear; clc;
load symbolic2.mat
%}
% set post-impact coordinate translation
syms q1dotp q2dotp 
qdotp = [ q1dotp ; q2dotp ] ; 

q1p =  q1_ + q2_ ; 
q2p = -q2_; 
qp_ = [ q1p ; q2p ] ;

matlabFunction( [qp_;0;0] , 'File', 'symbolic_functions/reset_map2_basic_only', 'Vars', { q } ) ;

%% MINUS - BEFORE COLLISION

Jvc = OTHER.Jvc ;
Jwc = OTHER.Jwc ;
H   = OTHER.TRANSFORMATIONS ;
H_COM = OTHER.H_COM ;
MoI = OTHER.MoI ;

L = sym('L',[n,1]) ;% store angular momentum
L(:) = 0 ;
for TD = n:-1:1 % so far, only accounts for planar vectors. adjust L/L_ to save 3D
    
    for i = 1:TD
        m_i = params.m(i) ;

        J_P  = Jvc(:,:,i) ;
        J_O = Jwc(:,:,i) ;  
        
        r_COM = -H_COM(1:3,end,i+1) ;% vector, stance foot to COM 
        r_TD = -H(1:3,end,TD+1) ;% vector, stance foot to swing foot touchdown
        r_i = r_TD - r_COM  ;% vector from touchdown point to current mass 
        
        v_i = simplify( J_P*qdot_ ) ;% velocity vector
        
        L_  =  mycross_noconj( r_i , m_i*v_i ) ...
                   + MoI(:,:,i)*J_O*qdot_ ;

        L(TD) = L(TD) + L_(3) ;
    end

end
Lm_ = simplify(L) ;% Angular momentum before impact! 

Qm = jacobian( Lm_ , qdot_ ) ;

%% PLUS - AFTER COLLISION
% get transformation information backawards

qr_ = qp_ + qconfig ;

params2 = params ;
params2.q_ = qp_ ;
params2.qdot_= qdotp ;
params2.q    = [ qp_ ; qdotp ]  ;
params2.runshort = 1 ;% rerun but don't compute B,C,G

for link = 1:n 

    params2.DH_q   = subs( params.DH_q , q_ , qp_ ) ;
    params2.DH_l(link)    = params.DH_l(n+1-link) ;
    params2.lc(link)      = params.DH_l(n+1-link) - params.lc(n+1-link)   ;
    params2.m(link)       = params.m(n+1-link) ;
    params2.moi(:,:,link) = params.moi(:,:,n+1-link) ;

end

tic
[ ~ , ~ , ~ , OTHERp ] = symbolic_dynamic_rotation( params2 ) ;

Jvc_p = OTHERp.Jvc ;
Jwc_p = OTHERp.Jwc ;
H_p  = OTHERp.TRANSFORMATIONS ;
H_COM_p = OTHERp.H_COM ;
MoI_p = OTHERp.MoI ;

L(:) = 0 ;
for TD = n:-1:1 % so far, only accounts for planar vectors. adjust L/L_ to save 3D
    
    for i = n:-1:TD %1:TD
        m_i = params2.m(i) ;

        J_P = Jvc_p(:,:,i) ;
        J_O = Jwc_p(:,:,i) ;  
        
        r_COM = - H_COM_p( 1:3 , end , i+1 )  ;
        r_TD = - H_p( 1:3 , end , TD ) ;% vector, stance foot to swing foot touchdown
        r_i = simplify( r_TD - r_COM ) ;% vector from touchdown point to current mass 
        
        v_i = simplify( J_P*qdotp ) ;% velocity vector         
        
        L_ = mycross_noconj( r_i , m_i*v_i ) ... 
                   + MoI_p(:,:,i) * J_O * qdotp ; 
        L(n+1-TD) = L(n+1-TD) + L_(3) ; 
    end

end

Lp_ = simplify( L ) ;% Angular momentum after impact!

Qp = jacobian(Lp_,qdotp) ;

qdot_plus = Qp \ Qm * qdot_ 


matlabFunction( [ qp_ ; qdot_plus ] , 'File', 'symbolic_functions/reset_map2', 'Vars', { q } ) ;


%% BACKWARDS RESET FUNCTION

q1m = q1_ + q2_ ; 
q2m = -q2_ ; 
qm_ = [ q1m ; q2m ] ;

Qp_ = subs( Qp , q_ , qp_ ) ;
Qm_ = subs( Qm , q_ , qp_ ) ;

qdot_minus = Qm_ \ Qp_ * qdot_ 

matlabFunction( [ qm_ ; qdot_minus ] , 'File', 'symbolic_functions/reset_map_bw2', 'Vars', { q } ) ;

%% ENERGY

syms hxn ;
KE = 0 ;
PE = 0 ;
for i = 1:n
%     KE = KE + 1/2*( params.m(i)*params.lc(i)^2 + params.moi(3,3,i) )*q(n+i)^2 ;   
%     PE = PE + params.m(i)*grav*( hxn + OTHER.TRANSFORMATIONS(2,4,i+1) ) ;
%     KE = KE + 1/2*( params.m(i)*params.lc(i)^2 + params.moi(3,3,i) )*q(n+i)^2 ;
    PE = PE + params.m(i)*grav*( hxn + OTHER.H_COM(2,4,i+1) ) ;    
end
KE = 1/2*qdot_.'*B_MATRIX*qdot_ ;
E_total = KE + PE ;

matlabFunction( E_total , PE , KE , 'File', 'symbolic_functions/get_energy', 'Vars', { q , hxn } ) ;

%% 

syms E0 hxn % Ediff

qd = [ 0*pi/180  0 ].' ;
Ef = get_energy( [qd(1) 0 qd(2) 0 ]' , hxn )  ;%%%**************** set Ef = Ed instead. 
% For CP, Ef = mgl. 
% Therefore, it solves for the point where Ediff = E0 - Ed. Goal is to
% regulate Ed, not to regulate E0. We need to make sure we have energy Ed
% post collision, which could be some value mgl+KE(at top - i.e. walking
% speed). Shouldn't change things much, just requires us to pass in Ed.
%
% I wonder if this will work to help increase velocity when we know a hill
% is coming up...
%

Ediff = E0 - Ef ;

syms q1_sym q1dot_sym 'real' 
assume( q1_sym < 0 & q1_sym > -pi ) ;

q2_sym = acos( cos(q1_sym) - hxn ) + 2*pi - q1_sym ;
assume( q2_sym < pi & q2_sym > 0 ) ;
assume( q2_sym , 'real' )

qm_sym = [ q1_sym q2_sym q1dot_sym 0 ].' ;

[Em_sym,PEm_sym,KEm_sym] = get_energy( qm_sym , 0) ;

qp_sym = reset_map2(qm_sym) ;
[Ep_sym,PEp_sym,KEp_sym] = get_energy( qp_sym , hxn ) ;

eqn = (KEm_sym - KEp_sym ) - Ediff 

q1dotm_2 = 2/(1*1^2) *(E0 - PEm_sym ) ;
assume( q1dotm_2 , 'real') ;
assume( q1dotm_2 , 'positive') ;

find_zero = expand( subs( eqn, q1dot_sym^2 , q1dotm_2 ) ) ;
% find_zero = @(q1_sym,hxn,E0) (981*hxn)/100 - (981*cos(q1_sym))/100 - (cos(q1_sym - acos(cos(q1_sym) - hxn))^2*(2*E0 - (981*cos(q1_sym))/50))/2 + 981/100 ; 
% EQN = @(q1_sym,hxn,E0) (981*hxn)/100 - (981*cos(q1_sym))/100 + (981*cos(q1_sym)^5)/100 + (981*hxn^2*cos(q1_sym)^3)/100 + (981*cos(q1_sym)*sin(q1_sym)^2)/100 - (981*cos(q1_sym)^3*sin(q1_sym)^2)/100 - E0*cos(q1_sym)^4 - E0*sin(q1_sym)^2 - (981*hxn*cos(q1_sym)^4)/50 + E0*cos(q1_sym)^2*sin(q1_sym)^2 + (981*hxn*cos(q1_sym)^2*sin(q1_sym)^2)/50 - (981*hxn^2*cos(q1_sym)*sin(q1_sym)^2)/100 + 2*E0*hxn*cos(q1_sym)^3 + (981*cos(q1_sym)^3*sin(q1_sym)*(2*hxn*cos(q1_sym) - cos(q1_sym)^2 - hxn^2 + 1)^(1/2))/50 - E0*hxn^2*cos(q1_sym)^2 + E0*hxn^2*sin(q1_sym)^2 - 2*E0*hxn*cos(q1_sym)*sin(q1_sym)^2 - 2*E0*cos(q1_sym)^2*sin(q1_sym)*(2*hxn*cos(q1_sym) - cos(q1_sym)^2 - hxn^2 + 1)^(1/2) - (981*hxn*cos(q1_sym)^2*sin(q1_sym)*(2*hxn*cos(q1_sym) - cos(q1_sym)^2 - hxn^2 + 1)^(1/2))/50 + 2*E0*hxn*cos(q1_sym)*sin(q1_sym)*(2*hxn*cos(q1_sym) - cos(q1_sym)^2 - hxn^2 + 1)^(1/2) + 981/100 ;
% matlabFunction( find_zero , 'File', 'symbolic_functions/find_zero', 'Vars', { q1_sym, hxn , Ediff , E0 } ) ;
matlabFunction( find_zero , 'File', 'symbolic_functions/find_zero', 'Vars', { q1_sym, hxn , E0 } ) ;


%% INVERSE KINEMATICS

syms xn l_swing 'real'


l_stance = l_leg

alpha = atan2( hxn, xn ) ;
q2_IK = acos( ( ( xn.^2 + hxn.^2 ) - l_stance.^2 - l_swing.^2 ) / ( -2*l_stance*l_swing) ) ;
beta  = acos( ( l_swing.^2 - l_stance.^2 - ( xn.^2 + hxn.^2 ) ) / ( -2*l_stance*sqrt( xn.^2 + hxn.^2 )) ) ;
q1_IK = beta - pi/2 + alpha ;

matlabFunction( q1_IK , q2_IK , 'File', 'symbolic_functions/get_IK', 'Vars', { xn, hxn , l_swing } ) ;



%% INVERSE KINEMATICS EQUAL LEGS

alpha = atan2( hxn, xn ) ;
q2_IK = acos( ( ( xn.^2 + hxn.^2 ) - 2*l_stance.^2 ) / ( -2*l_stance.^2) ) ;
q2_IK = piecewise( q2_IK < 0 , -q2_IK , q2_IK ) ;

beta = pi/2 - q2_IK ./ 2 ;
% beta  = acos( sqrt( xn.^2 + hxn.^2 ) / ( 2*l_stance ) ) ;
q1_IK = beta + alpha - pi/2;

matlabFunction( q1_IK , q2_IK , 'File', 'symbolic_functions/get_IK2', 'Vars', { xn, hxn } ) ;




