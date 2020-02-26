clear; clc; format compact ; close all ;
%
% this sets up proper file paths
localDir = fileparts(mfilename('fullpath')) ;      
restoredefaultpath ;% clear paths before adding
addpath(fullfile(localDir, 'symbolic_functions')) ;
addpath(fullfile(localDir, 'my_helper_functions')) ;
addpath(fullfile(localDir, 'helper_functions_from_others')) ;
%}

%% CREATE SYMBOLIC FUNCTIONS
% Necessary to run if symbolic functions do not already exist or parameters/constants chagne
%
% these three function create the symbolic functions used throughout the
% code.
symbolic_basic_functions ;
symbolic_2link_pendulum ;
symbolic_simple_pendulum ;
%}

%%
clear;
%
oo_constants ;

% notice discretization is finer in region that we will be walking in to
% collect more data points in that region. 
gamma_ARRAY = [ [-50:0.1:-30.1] [-30:0.01:0] [0.1:0.1:89] ]*pi/180 ;% [rad] possible terrain slopes (negative is downhill)
q2_ARRAY = [ [ 1:0.1:9.9] [ 10:0.01:89] ]*pi/180 ;% [rad] possible inner angles 


Vdi = 0 ;% desired velocity after impact

N = 10* length( gamma_ARRAY ) * length( q2_ARRAY ) ;


DATA_SAVE = zeros( floor(N/10) , 3 ) ;% preallocate for speed

tic
count = 0 ;

    
Ed_Kpi = 1/2 * m_cm * Vdi.^2 % [J] kinetic energy at top of other side
Ed_pi  = m_cm*grav*l_cm + Ed_Kpi ;% [J] total energy at top of other side
    

for i_gamma = 1:length( gamma_ARRAY )
	gammai = gamma_ARRAY( i_gamma ) ;% [slope of terrain/step]
        
    
    q2mi = 1 *pi/180 ;
    d_q2 = 0.1 * pi/180 ;
    check_q2 = [ 1 1 1 ];

        
    for i_q2 = 1:length( q2_ARRAY )
        q2mi = q2_ARRAY( i_q2 ) ;

        ri = 2*l_cm*sin(q2mi/2) ;% [m] length of footstep
        dhi = ri*sin(gammai) ;% height of foot touchdown h direction
        dxi = ri*cos(gammai) ;% distance of footstep in x direction


        % POST COLLISION INFORMATION
        q1dot_pi = -sqrt( 2/(m_cm*l_cm^2)*( Ed_pi - m_cm*grav*l_cm*cos( q2mi/2 + gammai ) ) ) ;% [rad/s] angular rate after collision

        qpi = [ gammai+q2mi/2 , -q2mi , q1dot_pi , 0 ].' ;% post-collision state

        %             [ Edi , PEdi , KEdi ] = get_energy( qpi , 0 ) ;% post-collision local energy     


        % PRE COLLISION INFORMATION
        qmi = reset_map_bw2( qpi ) ;% pre-collision state            

        [ E0i , PE0i , KE0i ] = get_energy( qmi , 0 ) ;% pre-collision energy      

        
        
        % only save if initial energy is within the bounds we expect to use
        % for walking
        if E0i < 25 && E0i > 10 
            count = count + 1 ;
           

            DATA_SAVE(count,:) = [ E0i , dhi , dxi ] ; %[ E0i , Ed_pi , dhi , dxi ] ; 
        end
      
        
    end


    
end
    
N_CP = count ;

toc
count 

DATA_SAVE = DATA_SAVE(1:count,:) ;

save('symbolic_functions/curve_of_capture_data.mat','DATA_SAVE','N_CP','-v7.3') 

%}
%% REGRESSION AND CREATE FUNCTION FOR CAPTURE POINT POLYNOMIAL
% above can be commented out if you want to just run regression below

clear; 
load curve_of_capture_data.mat

order = 6 ;% order of polynomial regression

tic
INPUTSx = DATA_SAVE(1:N_CP,[1,2]) ;%  [ E0i , dhi ]  
OUTPUTx = DATA_SAVE(1:N_CP,3) ;% [ dxi ] 
regx = my_regression_function( [ INPUTSx , OUTPUTx ] , order , 'regression_create_x_cp' ) 

tocx = toc 


save('symbolic_functions/curve_of_capture_regression.mat','INPUTSx','OUTPUTx','regx','order','-v7.3') ;
%}
%% PLOT REGRESSION RESULT
% comment above out if already run and nothing changed

clear; 
load curve_of_capture_regression.mat
%

x_train = INPUTSx ;
y_train = OUTPUTx ;

n_mesh = 100;

xmin  = min(x_train);
xmax = max(x_train);
X_E0 = linspace(xmin(1),xmax(1),n_mesh)'; 
X_dh = linspace(xmin(2),xmax(2),n_mesh)';
[ X_E0_mesh , X_dh_mesh ] = meshgrid( X_E0 , X_dh ) ; 

R2 = regx.R2 ;

for i1 = 1:n_mesh
for i2 = 1:n_mesh
%     Z1(i1,i2) = regression_create_x2( [ X_E0_mesh(i1,i2) , 10 , X_dh_mesh(i1,i2) ] )  ;
    Z1(i1,i2) = regression_create_x_cp( [ X_E0_mesh(i1,i2) , X_dh_mesh(i1,i2) ] )  ;    
end, end
% Z2 = recreate_x2( X_E0_mesh , X_dh_mesh ) ;


myfig = figure; 
myfig.Position = [493 492 748 505] ;
hold on;
surf(X_E0_mesh,Z1,X_dh_mesh); 
% surf(X_E0_mesh,Z2,X_dh_mesh); 
plot3(x_train(:,1),y_train,x_train(:,2),'bo')% comment this line to stop plotting data on top of mesh
% plot3(x_train(:,1),y_train,x_train(:,2),'bo') 
XX_train = x_train ;
YY_train = y_train ;
xlabel('energy before impact E^- [J]','fontsize',22,'position',[16.9,-0.1,-1.0],'rotation',-20);
ylabel('step size \Deltax [m]','fontsize',22,'position',[25.40,0.58,-1.0],'rotation',14);
zlabel('step height \Deltah [m]','fontsize',22);
% title(sprintf('Degree Polynomial Fit: %1.0f',order),'fontsize',22)
% title({'Curve of Equal Energy';sprintf('Degree Polynomial Fit: %1.0f',order)},'fontsize',26)
title('Curve of Capture','fontsize',26,'position',[14.22 , 0.81 , 0.84 ])
xlim([10 20])
ylim([0 max(y_train) ]) %max(max(Z1))]) ;
leg = legend({['curve of capture' newline...
                    '\Deltax_{cp} = R_{cp}( E^- , \Deltah )']},...
                    'location','east');
% leg = legend('Function Evaluation','Test Data','Training Data');
leg.FontSize = 20;
leg.Position = [ 0.56 , 0.68 , 0.33 , 0.09 ] ;
text(12,0.2,-0.75,sprintf('Degree Polynomial Fit: %1.0f',order),'fontsize',20)
text(12,0.2,-0.95,sprintf('R^2 = %1.4f',R2),'fontsize',20)
% text(15,0.2,-0.8,sprintf('R2_{adj}=%1.5f',R2_adj),'fontsize',20)
view(50 , 25)
grid on


coeff_table = regx.coeff_table
mu = regx.mu
sigma = regx.sigma 
R2 = regx.R2
order = order



%% FUNCTIONS

function [value, isterminal, direction] = impact_event( t , q , params )

    q1 = q(1) ;
    q2 = params.q2 ;
    hxn = params.hxn ;

    % point where foot hits the ground
    value = -hxn + cos(q1) - cos(q1+q2) ;
    
    isterminal = 1 ;
    direction = -1 ;
    
end




