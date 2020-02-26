%% FUNCTION : : DYNAMICS 
function [ B_MATRIX , C_MATRIX , G_MATRIX , OTHER ] = oo_Symbolic_Dynamic_Rotation( params ) 
    
    a.b.c = 1 ;
    RUNSHORT = 0 ;
    if isfield( params , 'runshort')
        if params.runshort
            RUNSHORT = 1 
        end
    end

    % get DH parameters
    DH    = get_DH( params ) ;
    DH_n  = size(DH,1)       ;
    Nq    = params.Nq     ;% number of full states
    Nq_   = params.Nq / 2 ;% number of position states
    q_    = params.q_     ;
    qdot_ = params.qdot_  ;
    
    
    % preallocate for link transformations (from previous link)
    H_DH_i = sym('H_DH_i',[4,4,DH_n]) ;    
    H_DH_i(:,:,:) = 0 ;
    
    % preallocate for link transformations (from world frame)
    H_i = sym('H_i',[4,4,DH_n+1]) ;    
    H_i(:,:,:) = 0 ;
    H_i(:,:,1) = eye(4) ;
    
    % preallocate for link rotation axis
    Z = zeros(3,DH_n+1); 

    % preallocate for link locations (world frame)
    O = sym('O',[3,DH_n+1]) ;
    O(:,:) = 0 ;

    % preallocate for link COM transformations (from previous link)
    H_DH_COM_i = sym('H_DH_COM_i',[4,4,DH_n]) ;    
    H_DH_COM_i(:,:,:) = 0 ;
    
    % preallocate for link COM transformations (from world frame)
    H_COM_i = sym('H_COM_i',[4,4,DH_n+1]) ;    
    H_COM_i(:,:,:) = 0 ;
    H_COM_i(:,:,1) = eye(4) ;

    % preallocate for COM locations (world frame)
    O_COM = sym('O_COM',[3,DH_n+1]) ;
    O_COM(:,:) = 0 ;

    for i = 1:DH_n
        
        % get DH parameters
        theta_i = DH(i,1) ;% rotation axis
        d_i     = DH(i,2) ;% *used for 3d only
        a_i     = DH(i,3) ;% length until end effector
        alpha_i = DH(i,4) ;% *used for 3d only
        
        % transformation of next link (link frame)
        H_DH_i(:,:,i+1) = single_H(theta_i,d_i,a_i,alpha_i) ;

        % transformation of next link (world frame)
        H_i(:,:,i+1) = simplify( H_i(:,:,i) * H_DH_i(:,:,i+1) );

        % location and rotation axis of joint
        O(:,i+1) = H_i(1:3,4,i+1) ; % location
        Z(:,i+1) = H_i(1:3,3,i+1) ; % rotation axis

        % get location of link COM for transformation matrix (link frame)
        lc_i = params.lc(i) ;
        H_DH_COM_i(:,:,i) = single_H(theta_i,d_i,lc_i,alpha_i) ;
        
        % link COM transformation matrix (world frame)
        H_COM_i(:,:,i+1) = H_i(:,:,i) * H_DH_COM_i(:,:,i) ;  

        % location of link center of mass (world frame)
        O_COM(:,i+1)  = H_COM_i(1:3,4,i+1) ; 

    end

    Z(:,1)     = H_i(1:3,3,1) ;     
    O_COM(:,1) = zeros(3,1) ; 

    % preallocate for COM jacobians (world frame)
    Jvc_COM_i = sym('Jvc_COM_i',[3,DH_n,DH_n]) ;
    Jvc_COM_i(:,:,:) = 0 ;
   
    Jwc_COM_i = zeros(3,DH_n,DH_n) ;
    
    % prallocate for moment of inertia (world frame)
    MoI_ = sym('MoI_',[3,3,DH_n]) ;  
    MoI_(:,:,:) = 0 ;    
    
    % NOTE: The following Jvc/Jwc setup is for rotational links only 
    % [see Siciliano p252 for more info]
    parfor i = 1:DH_n
%     for i = 1:DH_n
        
        % position of next "current link" com
        Oip1 = O_COM(:,i+1) ;

        % zero out the jth jacobian (works with parfor)
        Jvc_COM_j = sym('Jvc_COM_j',[3,DH_n]) ; 
        Jvc_COM_j(:,:) = 0 ;
        Jwc_COM_j = sym('Jwc_COM_j',[3,DH_n]) ; 
        Jwc_COM_j(:,:) = 0 ;   
        
        for j = 1:DH_n
            
            
            if j <= i % for each mass up to "current link"
                
                % position of current link
                Oj = O(:,j) ; 
                Z_ = Z(:,j) ;

                % jth element of ith jacobian
                Jvc_COM_j(:,j) = mycross_noconj( Z_ , Oip1-Oj ) ;%[Siciliano Eqn 7.18]
                Jwc_COM_j(:,j) = Z_ ;
               
            end
                
        end

        % place the j loop into the ith element
        Jvc_ = simplify( Jvc_COM_j ) ;
        Jvc_COM_i(:,:,i) = Jvc_ ;
        Jwc_COM_i(:,:,i) = Jwc_COM_j ;
        
        
        R_      = H_i(1:3,1:3,i+1) ; % rotation matrices         

        % moment of inertia (link frame)
        MoI_i   = params.moi(:,:,i) ;

        % moment of inertia (world frame)
        MoI_(:,:,i) = R_*MoI_i*(R_.') ;
    
    end

   
    t_setup = toc ;
    fprintf('%1.1f sec to setup (%1.1f sec elapsed total)\n',t_setup ,toc)       
    
    
    B_MATRIX = sym('B',[ DH_n , DH_n ]) ;
    B_MATRIX(:,:) = 0 ;
    C_MATRIX = ones( Nq_,Nq_ )*q_ ;% symbolic placeholder    
    G_MATRIX = sym('G_MATRIX',[Nq_,1]) ;
    G_MATRIX(:) = 0 ;    
    
   %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%  
   % CREATE MASS/INERTIA MATRIX ( B-matrix per Siciliano ) 
    
    if ~RUNSHORT
   
    %     parfor i = 1:DH_n
        for i = 1:DH_n

            MoI_i = MoI_(:,:,i)*0 ;

            
            m_i = params.m(i) ;

            J_P = Jvc_COM_i(:,:,i) ;
            J_O = Jwc_COM_i(:,:,i) ;  

            B_MATRIX = B_MATRIX + m_i*J_P.'*J_P ...
                  + J_O.'*MoI_i*J_O ;     

        end

        t_B_create = toc ;
        fprintf('%1.1f sec to create B matrix (%1.1f sec elapsed total)\n',t_B_create-t_setup,toc)   

        B_MATRIX = simplify(B_MATRIX)  ;

        t_B_simplify = toc ;
        fprintf('%1.1f sec to simplify B matrix (%1.1f sec elapsed total)\n',t_B_simplify-t_B_create,toc)    



        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%      
        % CREATE CORIOLIS MATRIX ( C-matrix per Siciliano )    

        parfor i = 1:Nq_ % ( 0.29 sec for Nq_=2 )
    %      for i = 1:Nq_ %  ( 0.20 sec for Nq_ = 2 ) 
            ck = sym('ck',[1,Nq_]) ;
            ck(:) = 0 ;
            for j = 1:Nq_
                for k = 1:Nq_

                    % Siciliano
                    dMij_dqk = jacobian( B_MATRIX(i,j) , q_(k) ) ;
                    dMik_dqj = jacobian( B_MATRIX(i,k) , q_(j) ) ;
                    dMjk_dqi = jacobian( B_MATRIX(j,k) , q_(i) ) ;

                    ck(k) = ( dMij_dqk + dMik_dqj - dMjk_dqi )* qdot_(k) ;
                end

                C_MATRIX(i,j) = 1/2*sum(ck) ;% Coriolis matrix {Eqn 3.10} 
            end
        end

        t_C_create = toc ;
        fprintf('%1.1f sec to create C matrix (%1.1f sec elapsed total)\n',t_C_create-t_B_simplify,toc)  

        C_MATRIX = simplify(C_MATRIX) ;

        t_C_simplify = toc ;
        fprintf('%1.1f sec to simplify C matrix (%1.1f sec elapsed total)\n',t_C_simplify - t_C_create , toc)    



        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%      
        % CREATE GRAVITY MATRIX 

        g0 = [ 0 -params.grav 0].' ;

        for i = 1:Nq_

           m_i   = params.m(i)  ;
           J_P_i = Jvc_COM_i(:,:,i) 
           gi = - m_i*g0.'*J_P_i ;

           G_MATRIX = ( G_MATRIX.' + gi ).'; 

        end


        t_G_create = toc ;
        fprintf('%1.1f sec to create G matrix (%1.1f sec elapsed total)\n',t_G_create-t_C_simplify,toc)    

        G_MATRIX = simplify( G_MATRIX ) ;

        t_G_simplify = toc ;
        fprintf('%1.1f sec to simplify G matrix (%1.1f sec elapsed total)\n',t_G_simplify-t_G_create,toc)
    
    end
    
    % 
    OTHER.TRANSFORMATIONS = simplify( H_i );
    OTHER.MoI = simplify( MoI_ ) ;
    OTHER.Jvc = simplify( Jvc_COM_i ) ;
    OTHER.Jwc = Jwc_COM_i ;
    OTHER.H_COM = simplify( H_COM_i ) ;

    
    fprintf('%1.2f sec to create %1.0i-link dynamics\n',toc,Nq_)
    

    
end

%% FUNCTION : : Transformation Matrix (single line) from Denavit Hartenberg Parameters
function H_DH_i = single_H(theta_i,d_i,a_i,alpha_i)

    H_DH_i =...
        [cos(theta_i), -sin(theta_i)*cos(alpha_i), sin(theta_i)*sin(alpha_i), a_i*cos(theta_i);
         sin(theta_i), cos(theta_i)*cos(alpha_i), -cos(theta_i)*sin(alpha_i), a_i*sin(theta_i);
         0, sin(alpha_i), cos(alpha_i), d_i;
         0, 0, 0, 1] ;
         
end

%% EDIT HERE : : FUNCTION : : Denavit Hartenberg
function DH = get_DH( params )
% Configuration from Denavit Hartenberg Parameters
    
    DH = sym('DH',[params.DH_n,4]) ;

    DH(:,1) = params.DH_q ;% rotation column 
    DH(:,3) = params.DH_l ;% translation column
    DH(:,[2,4]) = zeros(params.DH_n,2) ;% not used in 2d rotation only 
end

% function C = mycross_noconj(a,b) 
% 
%     C = [ a(2)*b(3) - a(3)*b(2) ;
%           a(3)*b(1) - a(1)*b(3) ;
%           a(1)*b(2) - a(2)*b(1) ] ;
%       
% end