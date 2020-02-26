clear; clc; format compact ; close all ;
%
% this sets up proper file paths
localDir = fileparts(mfilename('fullpath')) ;      
restoredefaultpath ;% clear paths before adding
addpath(fullfile(localDir, 'symbolic_functions')) ;
addpath(fullfile(localDir, 'my_helper_functions')) ;
addpath(fullfile(localDir, 'helper_functions_from_others')) ;


%% CREATE SYMBOLIC FUNCTIONS
% Necessary to run if symbolic functions do not already exist or parameters/constants chagne
%
% these three function create the symbolic functions used throughout the
% code.
symbolic_basic_functions ;
symbolic_2link_pendulum ;
symbolic_simple_pendulum ;
%}

%% ASSEMPBLE DATA HERE
% can hide all above once run once
%

clear; 
oo_constants ;

q1_ARRAY = [ [-89:0.1:-70.1] [-70:0.01:-10] [-10.1:0.1:89] ]*pi/180 ;
q2_ARRAY = [ [ 1:0.1:9.9] [ 10:0.01:60] [ 60.1:0.1:89] ]*pi/180 ;% [rad] possible inner angles 

N = 10*length( q1_ARRAY ) * length( q2_ARRAY ) ;

DATA_SAVE = zeros( N , 3 ) ;

tic
count = 0 ;
for i_q1 = 1:length( q1_ARRAY )
    q1mi = q1_ARRAY( i_q1 ) ; 

    PEm = m_cm*grav*l_cm*cos( q1mi ) ;
    
    for i_q2 = 1:length( q2_ARRAY )
        q2mi = q2_ARRAY( i_q2 ) ;
        
        PEp = m_cm*grav*l_cm*cos( q1mi + q2mi ) ;

        dPE = PEm - PEp ;
        
        dKE = -dPE ;
        
        q1dotmi_2 = dKE / (1/2*m_cm*l_cm^2*sin( q2mi )^2) ;
        
        KEm = 1/2*m_cm*( l_cm^2 * q1dotmi_2 ) ;
        KEp = 1/2*m_cm*( l_cm^2 * q1dotmi_2 * cos( q2mi )^2 ) ;
        KEm - KEp ;
        
        Em = KEm + PEm ;
        Ep = KEp + PEp ;
        
        if abs( Em - Ep ) > 1E-10
           Em
           Ep
           pause()
        end
        
        Ti = transformations_short2( [ q1mi ; q2mi ] ) ;
        XHi = Ti(1:2,end,end) ;
           
        dxi = XHi(1) ;
        dhi = XHi(2) ;           

        if ( Em >= m_cm*grav*l_cm ) && ( Em < 50 ) && ( dxi > 0 ) && ( dhi > -1 ) && ( dhi < 0.8 )  
           count = count + 1 ;

           DATA_SAVE(count,:) = [ Em , dhi , dxi ] ;
            
        end
        
        
    end
    
end

toc
count 

DATA_SAVE = DATA_SAVE(1:count,:) ;

save('symbolic_functions/curve_of_equal_energy_data.mat','DATA_SAVE','-v7.3') 

%}
%% REGRESSION AND CREATE FUNCTION FOR EQUAL ENERGY POLYNOMIAL
% above can be commented out if you want to just run regression below
%
clear; 
load curve_of_equal_energy_data.mat

order = 4 ;% order of polynomial regression

tic
INPUTSx = DATA_SAVE(:,[1:2]) ;%  [ E0 , dh ]  
OUTPUTx = DATA_SAVE(:,3) ;% [ dxi ] 
regx = my_regression_function( [ INPUTSx , OUTPUTx ] , order , 'regression_create_x_equal_energy' ) 

tocx = toc 

save('symbolic_functions/curve_of_equal_energy_regression.mat','INPUTSx','OUTPUTx','regx','order','-v7.3') ;
%}
%% PLOT REGRESSION RESULT
% comment above out if already run and nothing changed

clear; 
load curve_of_equal_energy_regression.mat

x_train = INPUTSx ;
y_train = OUTPUTx ;
R2 = regx.R2 ;

n_mesh = 100;

xmin = min(x_train);
xmax = max(x_train);
X_E0 = linspace(xmin(1),xmax(1),n_mesh)'; 
X_dh = linspace(xmin(2),xmax(2),n_mesh)';
[ X_E0_mesh , X_dh_mesh ] = meshgrid( X_E0 , X_dh ) ; 

for i1 = 1:n_mesh
for i2 = 1:n_mesh
    Z1(i1,i2) = regression_create_x_equal_energy( [ X_E0_mesh(i1,i2) , X_dh_mesh(i1,i2) ] )  ;    
end, end

myfig = figure; 
myfig.Position = [493 492 748 505] ;
hold on;
surf(X_E0_mesh,Z1,X_dh_mesh); 
plot3(x_train(:,1),y_train,x_train(:,2),'bo') ;% comment out to show mesh without data
XX_train = x_train ;
YY_train = y_train ;
xlabel('energy before impact E^- [J]','fontsize',22,'position',[16.9,-0.1,-1.0],'rotation',-22);
ylabel('step size \Deltax [m]','fontsize',22,'position',[25.9,0.6,-1.1],'rotation',14);
zlabel('step height \Deltah [m]','fontsize',22);
title('Curve of Equal Energy','fontsize',26,'position',[13.56,0.72,0.12])
xlim([10 25])
ylim([0 max(y_train) ]) %max(max(Z1))]) ;
leg = legend({['curve of equal energy' newline...
                    '\Deltax_{\DeltaE=0} = R_{\DeltaE=0}( E^- , \Deltah )']},...
                    'location','east');
leg.FontSize = 20;
leg.Position = [0.60,0.70,0.28,0.09 ] ;
text(20,0.6,-0.9,sprintf('R^2 = %1.4f',R2),'fontsize',20)
text(20,0.6,-0.8,sprintf('Degree Polynomial Fit: %1.0f',order),'fontsize',20)
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
       
    value = -hxn + cos(q1) - cos(q1+q2) ;
    
    isterminal = 1 ;
    direction = -1 ;
    
end




