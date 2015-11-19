function EKF4(filename, mode)
	% performe extended Kalman filtering with the 2D musculoskeletal model (gait2de)
	
	% Inputs:
	%	filename.........(string) Name of a file with simulated state trajectories (e.g. result21)
    %   mode.............(scalar) 1: measure all states, 2: measure states 1-18 and GRF

	global model

	% store some information about the model
	model.Nstates = 50;		% number of state variables
    if (mode == 1)
        model.Nmeas = 50;       % measure all states
    elseif (mode == 2)
        model.Nmeas = 24;		% 9 coordinates, 9 velocities, and 6 ground reactions
    elseif (mode == 3)
        model.Nmeas = 18;		% 9 coordinates, 9 velocities
    else
        error('mode value not recognized')
    end
    model.Nmus = 16;		% number of muscles
	model.dt = 0.0016;		% time step for the Kalman filter

	% covariances (these need more thinking to get good values!!)
% 	Qu = 0.00001*eye(model.Nmus);			% covariance in unknown muscle controls Decreasing Qu makes the results worst
% 	Qz = 0.01*eye(model.Nstates);		% process noise covariance
% 	R  = 0.01*eye(model.Nmeas);			% measurement noise covariance
    
    Qu =0.1*eye(model.Nmus);			% covariance in unknown muscle controls Decreasing Qu makes the results worst
	Qz = diag([ zeros(1,9) 1 1 1 zeros(1,38)]);		% process noise covariance
	R  = 0.0000001*eye(model.Nmeas);			% measurement noise covariance

    

	% w = [u,z], so Q = [Qu zeros ; zeros Qv]
	Q = [Qu zeros(model.Nmus, model.Nstates) ; zeros(model.Nstates, model.Nmus) Qz];

	% load the data file and generate simulated measurements
	% use a file that has many gait cycles!
	[t, ydata, xtrue] = measure(filename, R, model.dt);
    N = numel(t);           % number of samples

	% initialization of the Kalman filter
	xhat = xtrue(:,1); 				% initial state estimate, from the true data (later we should use something else)
	P =0.001*eye(model.Nstates); 	% estimation error covariance
	
	% determine observability at xhat, using eq. 1.154 in Simon's book
	[~, F] = process(xhat);
	[~, H] = meas(xhat);
	QQ = [];
	HFn = H;		% H*F^n, starting at n=0
	for i = 1:model.Nstates
		QQ = [QQ ; HFn];
		HFn = HFn * F;				% H*F^n = H*F^(n-1) * F
        disp([i rank(QQ)])
	end
	fprintf('Rank of observability matrix is %d.\n', rank(QQ));
    disp('Hit ENTER to continue');pause

	% create an array to store the results
	xhatArray = zeros(size(xtrue));

    xhatArray(:,1) = xhat;

	% run the Kalman filter through all time steps
	Pminus = P;        % make a copy because we will overwrite P during the iterations
    for k = 1:N-1
	
		% do equations 13.47
		[xhat, F, L] = process(xhat);	% the process model
		P = F*P*F' + L*Q*L';

%         % do equations 13.64
%         Pminus = P;        % make a copy because we will overwrite P during the iterations
%         Niterations=1;
%         for i = 1:Niterations
%             [h, H, M] = meas(xhat);            % the measurement model
%             K = P * H' / (H * P * H' + M * R * M'); % Kalman gain                                
%             xhat = xhat + K * (ydata(:,k) - h);     % state estimate measurement update         
%             P = (eye(model.Nstates) - K * H) * Pminus ;  % estimation error covariance measurement update    
%         end
        
        		% do equations 13.49
		[y, H, M] = meas(xhat);			% the measurement model
		K = P * H' / (H * P * H' + M * R * M'); % Kalman gain                                Eq 13.49(1)
        xhat = xhat + K * (ydata(:,k) - y); 	% state estimate measurement update          Eq 13.49(2)
		P = (eye(model.Nstates) - K * H) * Pminus ;  % estimation error covariance measurement update     Eq 13.49(3)

        
		% Save state estimate for plotting
		xhatArray(:, k+1) = xhat;   % Save data in arrays.		
	end

	% plot the results
	close all
	musclenames={'RIliopsoas', 'RGlutei', 'RHamstrings', 'RRectus', 'RVasti', 'RGastroc', 'RSoleus', 'RTibialisAnterior',
		'LIliopsoas', 'LGlutei', 'LHamstrings', 'LRectus', 'LVasti', 'LGastroc', 'LSoleus', 'LTibialisAnterior'};
		for x1=1:9
	figure(1)
	subplot(5,2,x1); hold on;
	plot(t,xtrue(x1,:),'r-', t,xhatArray(x1,:),'b--')
	set(gca,'FontSize',8);
	xlabel('time (s)'), ylabel(['Position #',num2str(x1)])
	title(['State #',num2str(x1)],'Color','R')
	RMSError = sqrt(sum((xtrue(x1,:)-xhatArray(x1,:)).^2)/N);
	disp(['RMS position estimation error = ',num2str(RMSError), ' rad'])

	grid
	box on
		end
		
		   for x1=10:18
	figure(2)
	subplot(5,2,x1-9); hold on;
	plot(t,xtrue(x1,:),'r-', t,xhatArray(x1,:),'b--')
	set(gca,'FontSize',8);
	xlabel('time (s)'), ylabel(['Velocity #',num2str(x1)])
	title(['State # ',num2str(x1)],'Color','R')
	RMSError = sqrt(sum((xtrue(x1,:)-xhatArray(x1,:)).^2)/N);
	disp(['RMS velocity estimation error = ',num2str(RMSError), ' rad/sec'])
	grid
	box on
		   end
		   for x1=19:34
	figure(3)
	subplot(4,4,x1-18); hold on;
	plot(t,xtrue(x1,:),'r-', t,xhatArray(x1,:),'b--')
	set(gca,'FontSize',8);
	xlabel('time (s)'), ylabel('Lce(cm)')
	title(['State # ',num2str(x1),musclenames(x1-18)],'Color','R')
	RMSError = sqrt(sum((xtrue(x1,:)-xhatArray(x1,:)).^2)/N);
	disp(['RMS muscle length estimation error = ',num2str(RMSError)])
	grid
	box on
		  end
		  for x1=35:50
	figure(4)
	subplot(4,4,x1-34); hold on;
	plot(t,xtrue(x1,:),'r-', t,xhatArray(x1,:),'b--')
	set(gca,'FontSize',8);
	% set(gca,'FontSize',12);
	xlabel('time (s)'), ylabel('Active State')
	title(['State #',num2str(x1),musclenames(x1-34)],'Color','R')
	RMSError = sqrt(sum((xtrue(x1,:)-xhatArray(x1,:)).^2)/N);
	disp(['RMS muscle activation estimation error = ',num2str(RMSError)])
	grid
	box on
		  end
	%

end
%================================================================
function [f, F, L] = process(x);
	global model
	% System process model in the form of eq. 13.44 and 13.46 in Simon's book:
	% 		x_k+1 = f(x_k, w);
	%
	% w consists of muscle controls (u) and other unknown errors and disturbances (z)
	% so w = [u;z] and covariance matrix of w is Q = [Qu 0 ; 0 Qz];
	%
	%
	% Inputs:
	%	x................(50 x 1) state of the system at time k
	%   
	% Outputs:
	%   f................(50 x 1) state of the system at time k+1
	%   F................(50 x 50) Jacobian df/dx
	%   L................(50 x 66) jacobian df/dw
	
	% continuous time dynamics: dx/dt = gait2de(x,u) + z
	% we evaluate this at u=0 and z=0
	u = zeros(model.Nmus,1);
	g = gait2de(x,u);

	% estimate dg/dx and dg/du by finite differences
	h = 1e-7;
	dg_dx = zeros(model.Nstates, model.Nstates);
	for i=1:model.Nstates
        xsave=x(i);
        x(i)=x(i)+h;
        gh = gait2de(x,u);
        dg_dx(:,i)=(gh-g)/h;
        x(i)=xsave;    
    end
    for i=1:model.Nmus
        usave=u(i);
        u(i)=u(i)+h;
        gh = gait2de(x,u);
        dg_du(:,i)=(gh-g)/h;
        u(i)=usave;
    end
	
	% dg/dz is identity matrix
	dg_dz = eye(model.Nstates);
	
	% discrete time dynamics by forward Euler method
	f = x + g * model.dt;
	
	% Jacobians F=df/dx and L=df/dw for the linearized discrete time model
	F = eye(model.Nstates) + dg_dx*model.dt;
	L = [dg_du  dg_dz] * model.dt;

end
%================================================================
function [y, H, M] = meas(x)
% Measurement model y = h(x) + v (assuming additive noise v)
% Using notation from Simon's book eq. 13.44 and 13.48
%
% Inputs:
%	x................(Nstates x 1) State of the system
%
% Outputs:
%	y................(Nmeas x 1) Measured variables
%	H................(Nmeas x Nstates) Jacobian dh/dx
%   M................(Nmeas x Nmeas) Jacobian dh/dv

	global model
	
	y = zeros(model.Nmeas,1);
	u = zeros(model.Nmus,1);
    
    if (model.Nmeas == model.Nstates)
        y = x;
        H = eye(model.Nstates);
    elseif (model.Nmeas == 24)
        % compute the ground reactions and their derivatives
        [~,grf] = gait2de(x,u);
        dgrf_dx = zeros(size(grf,1), model.Nstates);
        h = 1e-7;
        for i=1:model.Nstates
            xsave=x(i);
            x(i)=x(i)+h;
            [~,grfh] = gait2de(x,u);
            dgrf_dx(:,i)=(grfh-grf)/h;
            x(i)=xsave;
        end

        % y consists of states 1-18 and the grf
        y = [x(1:18) ; grf];

        % the Jacobian dy/dx
        H = [eye(18) zeros(18,32) ; dgrf_dx];
        
    elseif (model.Nmeas == 18)
        % y consists of states 1-18 
        y = x(1:18);

        % the Jacobian dy/dx
        H = [eye(18) zeros(18,32)];
    else
        error('meas: no code yet for this number of measurements');   
    end
        
	% the Jacobian dh/dv
	M  = eye(model.Nmeas);      % Eq 13.48(2)
		
end
%================================================================
function [t, ydata, x] = measure(file, R, dt);
% Simulate measurements using state histories from file
%
% Inputs:
%	file.............(string) Name of the mat file with state-time history
%	R................(Nmeas x Nmeas) Covariance matrix of measurement noise
%	dt...............(scalar) Sampling time interval
%
% Outputs:
%	t.................(1 x N) time stamps of the measurements
%	data..............(Nmeas x N) the measurements
%   x.................(Nstates x N) resampled state time history

	% do the Cholesky decomposition of R
	A = chol(R);
%     A = sqrt(R);
	Nmeas = size(R,1);

	load(file);
	
	% resample the states to sampling time dt
	t1 = linspace(0,Result.dur,size(Result.x,2));		% original time stamps
	t = 0:dt:Result.dur;								% the new samples
	x = interp1(t1, Result.x', t)';						% resampled states
	Nsamples = numel(t);
	
	% simulate the measurement system and add noise
	% noise is simulated using https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
	ydata = zeros(Nmeas, Nsamples);
	for i = 1:numel(t)
		ydata(:,i) = meas(x(:,i)) + A*randn(Nmeas,1);
	end

end
%================================================================