%% Artificial Intelligence and Machine Learning for Engineering Design
% Prof: Levent Burak Kara
% Homework 2
% Question 1

% Steve Crews
% 13 Oct 18

%% 
function regress_out = my_regression( data , order , regression_function_name  )

    reg_factor = 1e-3 ;


    Nx = size(data,2)-1;

    coeff_base = create_coefficients( order , Nx ); 

    %% Part A No-standardization Eigenvalues
%     disp('Part A ')

    x_train = data(:,1:Nx);
    y_train = data(:,end);
    n = size(x_train,1);

    A_train = create_polys( x_train , coeff_base ) ;

    [ ~ , S , ~ ] = svd(A_train'*A_train);
    lambda = diag(S);
%     disp('eigenvalues=')
%     fprintf('     %3.2e\n\n',lambda)

    ratio = max(abs(lambda))/min(abs(lambda)) ;
    fprintf('Unscaled Eigenvalue Ratio: %3.2e\n',ratio) 


    %% Part D Standardization Eigenvalues 
%     disp(' '), disp('Part D ')

    x_train_scaled = x_train;
    sigmax = zeros(1,Nx) ;
    mux = zeros(1,Nx) ;
    for i = 1:Nx
        mux(i) = mean(x_train(:,i));
        sigmax(i) = std(x_train(:,i));
        x_train_scaled(:,i) = (x_train(:,i) - mux(i))/sigmax(i);
        % now mean = 0 and std dev = 1
    end

    A_train_scaled = create_polys( x_train_scaled , coeff_base ) ;

    Nreg = size(A_train_scaled,2) ;

    [ ~ , S_scaled, ~ ] = svd(A_train_scaled'*A_train_scaled);
    lambda_scaled = diag(S_scaled);
%     disp('eigenvalues=')
%     fprintf('     %3.2e\n\n',lambda_scaled)

    ratio_scaled = max(abs(lambda_scaled))/min(abs(lambda_scaled)) ;
    fprintf('Scaled Eigenvalue Ratio: %3.2e\n',ratio_scaled)

    
    c_train_scaled = (A_train_scaled'*A_train_scaled + eye( Nreg ) * reg_factor )\A_train_scaled'*y_train;
    

    
    %% R2
    % BETTER TO DO THIS WITH TEST DATA, COULD SPLIT UP DATA 70/30, 60/40
    % etc and try out
    
    my = mean(y_train); %mean of y data
    myvec = my*ones(n,1); %mean y as a vector
    y_meanshift = y_train - myvec;% mean shifted (zero mean) y

    SStot = y_meanshift'*y_meanshift ;

    f_predict = A_train_scaled*c_train_scaled ;    
    
    e = y_train - f_predict ;
    SSres = e'*e; %Just an alternative to SSres = e'e;

    R2 = 1 - SSres / SStot ;

    % adjusted R2 only interesting if using test data
    k = order*Nx;
    R2_adj = 1 - ( (1-R2)*(n-1) / (n-k-1) ) ;
    
    
    %% Create Function

    x_sym = sym('x_',[1,Nx]) ;

    x_sym_scaled = x_sym ;
    for i = 1:Nx
        x_sym_scaled(:,i) = (x_sym(:,i) - mux(i))/sigmax(i);
        % now mean = 0 and std dev = 1
    end

    A_sym = create_polys( x_sym_scaled , coeff_base ) ;

    f_sym = A_sym*c_train_scaled ;

    func_location_name = strcat( 'symbolic_functions/' , regression_function_name ) ;
    matlabFunction( f_sym , 'File', func_location_name , 'Vars', { x_sym }) ; 

    coeff_table = zeros(order+1) ;
    for cnt = 1:length(coeff_base)
        i = order - coeff_base( cnt,1 ) + 1 ;
        j = coeff_base( cnt,2 ) + 1 ;
        coeff_table( i,j ) = c_train_scaled(cnt) ;%Note: this table is not scaled with mean/std
    end

    regress_out.c = c_train_scaled ;
    regress_out.function_name = regression_function_name ;
    regress_out.mu = mux ;
    regress_out.sigma = sigmax ;
    regress_out.R2 = R2 ;
    regress_out.R2_adj = R2_adj ;
    regress_out.coeff_table = coeff_table ;
    
end


%% FUNCTIONS

function A = create_polys( x , coeff_base )

    Nx = size(x,2) ;
    n  = size(x,1) ;
    N_coeff = size(coeff_base,1) ; 
    
    if isa(x,'double')
        A = ones(n,N_coeff) ;
    else
        A = sym('A',[n,N_coeff]);
        A(:) = 1 ;
    end
    
        
    for i1 = 1:N_coeff
        x_base = ones(n,1) ;
        for i2 = 1:Nx 
            x_base = x_base .* x(:,i2) .^ ( coeff_base(i1,i2) ) ;
        end
        A( :,i1) = x_base ;
    end

end

function coeff_base = create_coefficients( order , xn ) 

    % Initialize
    coeff_base = zeros(order^xn,xn);

    % Create Colums Corresponding to Mathematical Base
    for ii=1:xn
        coeff_base(:,ii)=mod(floor((1:order^xn)/order^(ii-1)),order);
    end

    % Flip - Reduce - Augment
    coeff_base=fliplr(coeff_base); 
    coeff_base=coeff_base(sum(coeff_base,2)<=order,:); 
    Ab=diag(repmat(order,[1,xn])); 
    coeff_base=[coeff_base;Ab];

    % Degree Conditionals
    for ii=1:xn
        coeff_base=coeff_base(coeff_base(:,ii)<=order,:);
    end

end