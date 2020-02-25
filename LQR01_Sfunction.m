function [sys,x0,str,ts,simStateCompliance] = LQR01_Sfunction(t,x,u,flag)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 19;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%

ts  = [0.05 0];% ts_为ts的第一位，修改采样周期时下面的也要同时修改
global ts_1;
ts_1=ts(1,1);

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

sys = [];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
%%
%设置输入的期望路径
%     %直线
%     target_path_x=25*sin(0.2*t);
%     target_path_y=0;
%     target_path_heading=0.2*t;
%     target_path_curvature=0;
% %圆
% target_path_x=50*sin(0.2*t);
% target_path_y=50-50*cos(0.2*t);
% target_path_heading=0.2*t;
% target_path_curvature=1/50;
%RoadXY107

% time_gain=1;
% target_path_x=time_gain*t;
% target_path_y=road_y_function(target_path_x);
% target_path_heading=road_heading_function(target_path_x);
% target_path_curvature=road_curvature_function(target_path_x);

% %直线
% target_path_x_end=1000;
% target_path_x_pointDistance=0.1;
% target_path_x=0:target_path_x_pointDistance:target_path_x_end;
% target_path_size=target_path_x_end/target_path_x_pointDistance+1;
% target_path_y=zeros(1,target_path_size);
% target_path_heading=zeros(1,target_path_size);
% target_path_curvature=zeros(1,target_path_size);
% 
% %阶跃换道
% target_path_x_end=100;
% target_path_x_pointDistance=0.1;
% target_path_x=0:target_path_x_pointDistance:target_path_x_end;
% target_path_size=target_path_x_end/target_path_x_pointDistance+1;
% target_path_y=zeros(1,target_path_size);
% target_path_y(200:end)=2;
% target_path_heading=zeros(1,target_path_size);
% target_path_curvature=zeros(1,target_path_size);

%RX107
target_path_end=700;
target_path_space=0.1;
target_path_size=target_path_end/target_path_space;
target_path_x=0:target_path_space:target_path_end;
target_path_y=zeros(1,target_path_size+1);
target_path_heading=zeros(1,target_path_size+1);
target_path_curvature=zeros(1,target_path_size+1);
for n=1:target_path_size+1
    target_path_y(n)=road_y_function(target_path_x(n));
    target_path_heading(n)=road_heading_function(target_path_x(n));
    target_path_curvature(n)=road_curvature_function(target_path_x(n));
end
% plot(target_path_x,target_path_y);
% plot(target_path_x,target_path_heading);
% plot(target_path_x,target_path_curvature);
%%
%车俩参数，carsim内部的参数
mass_=1592;%整车质量(kg)(1111/60/60)
cf_=40207*2;%前轮侧偏刚度的2倍(N/rad)
cr_=40207*2;%后轮侧偏刚度的2倍(N/rad)
lf_=1.180;%前悬(m)
lr_=1.770;%后悬(m)
wheelbase_=lf_+lr_;%轴距(m)
iz_=2488;%转动惯量(kg*m^2)
tire_single_direction_max_degree=40.4*pi/180;
%设置参数
global ts_1;
ts_=ts_1;
Ts_=0.01;%离散化时间常数
filter_tire_angle=5*ts_;
%%
%建立状态空间方程
basis_state_size_=4;
preview_window=0;%预估暂时不考虑
matrix_size=basis_state_size_+preview_window;

%初始化A矩阵
matrix_A=zeros(basis_state_size_,basis_state_size_);
%matrix_Ad=zeros(basis_state_size_,basis_state_size_);
%matrix_Adc=zeros(matrix_size,matrix_size);
matrix_A_coeff=zeros(matrix_size,matrix_size);

matrix_A(1,2)=1;
matrix_A(2,3)=(cf_+cr_)/mass_;
matrix_A(3,4)=1;
matrix_A(4,3)=(cf_*lf_-cr_*lr_)/iz_;

matrix_A_coeff(2,2)=-(cf_+cr_)/mass_;
matrix_A_coeff(2,4)=(cr_*lr_-cf_*lf_)/mass_;
matrix_A_coeff(3,4) = 1.0; 
matrix_A_coeff(4,2)=(lr_*cr_-lf_*cf_)/iz_;
matrix_A_coeff(4,4)=-(lf_*lf_*cf_+lr_*lr_*cr_)/iz_;

%初始化B矩阵
matrix_B=zeros(basis_state_size_,1);
%matrix_Bd=zeros(basis_state_size_,1);
%matrix_Bdc=zeros(matrix_size,1);
matrix_B(2,1)=cf_/mass_;
matrix_B(4,1)=cf_*lf_/iz_;
matrix_Bd=Ts_.*matrix_B;%B矩阵离散化

%%
%读取当前信息
linear_v=u(5);
% linear_v=10/3.6;
if linear_v==0
    linear_v=0.001;
end
if t==0
    current_path_x=0;
    current_path_y=0;
    current_path_heading=0;
    angle_v=0;
    previous_tire_angle=0;
% elseif t==ts_
%     current_path_x=0;
%     current_path_y=-3.11;
%     current_path_heading=-0.7621989;
%     angle_v=0;
%     previous_tire_angle=0;
else
    current_path_x=u(1);
    current_path_y=u(2);
    current_path_heading=u(3);   
    angle_v=u(4);                
    
    previous_tire_angle=u(6);
end

%%
% if t>=ts_
%     
%     path_back1_x=time_gain*(t-ts_);
%     path_back1_y=road_y_function(path_back1_x);
%     distanceQ_back1=PointDistanceSquare(path_back1_x,path_back1_y,current_path_x,current_path_y);
%     
%     path_center_x=time_gain*(t);
%     path_center_y=road_y_function(path_center_x);
%     distanceQ_center=PointDistanceSquare(path_center_x,path_center_y,current_path_x,current_path_y);
% 
%     distanceQ_array=[distanceQ_back1,distanceQ_center,distanceQ_front10,...
%         distanceQ_front60,distanceQ_front110,distanceQ_front200,...
%         distanceQ_front300,distanceQ_front400,distanceQ_front500];
%     [distanceQ_min,id]=min(distanceQ_array);
%  
%     switch id
%         case 1
%             min_id=-1;
%             target_path_nearest_x=path_back1_x;
%             target_path_nearest_y=road_y_function(target_path_nearest_x);
%             target_path_nearest_heading=road_heading_function(target_path_nearest_x);
%             target_path_nearest_curvature=road_curvature_function(target_path_nearest_x);
%         case 2
%             min_id=0;
%             target_path_nearest_x=path_center_x;
%             target_path_nearest_y=road_y_function(target_path_nearest_x);
%             target_path_nearest_heading=road_heading_function(target_path_nearest_x);
%             target_path_nearest_curvature=road_curvature_function(target_path_nearest_x);
%     end
% else
%     min_id=0;
%     target_path_nearest_x=target_path_x;
%     target_path_nearest_y=target_path_y;
%     target_path_nearest_heading=target_path_heading;
%     target_path_nearest_curvature=target_path_curvature;
% end
%%
% point_range=800;
% target_path_choose_x=zeros(point_range,1);
% target_path_choose_y=zeros(point_range,1);
% distance_choose=zeros(point_range,1);
% 
% for n=0:1:point_range-1
%     target_path_choose_x(n+1)=time_gain*(t+n*ts_);
%     target_path_choose_y(n+1)=road_y_function(target_path_choose_x(n+1));
%     distance_choose(n+1)=PointDistanceSquare(target_path_choose_x(n+1),target_path_choose_y(n+1),...
%         current_path_x,current_path_y);
% end
% [distance_choose_min,id]=min(distance_choose);
% min_id=id-1;
% target_path_nearest_x=time_gain*(t+min_id*ts_);
% target_path_nearest_y=road_y_function(target_path_nearest_x);
% target_path_nearest_heading=road_heading_function(target_path_nearest_x);
% target_path_nearest_curvature=road_curvature_function(target_path_nearest_x);

distance_choose=11*ones(1,target_path_size); 
for n=1:target_path_size
    if(dot([current_path_x-target_path_x(n),current_path_y-target_path_y(n)],[1,tan(target_path_heading(n))])<=0)
        distance_choose(n)=PointDistanceSquare(target_path_x(n),target_path_y(n),...
            current_path_x,current_path_y);
    end
end
% for n=1:target_path_size
%     if(dot([current_path_x-target_path_x(n),current_path_y-target_path_y(n)],[1,0])<=0)
%         distance_choose(n)=PointDistanceSquare(target_path_x(n),target_path_y(n),...
%             current_path_x,current_path_y);
%     end
% end
[~,index_min]=min(distance_choose);

target_path_nearest_x=target_path_x(index_min);
target_path_nearest_y=target_path_y(index_min);
target_path_nearest_heading=target_path_heading(index_min);
target_path_nearest_curvature=target_path_curvature(index_min);

%%
%初始化状态State，反馈增益K
matrix_State=zeros(matrix_size,1);
matrix_K=zeros(1,matrix_size);
matrix_u=0;

%设计二次型规划参数Q，R
%600_100_600_100_1;1000_100_600_200
matrix_Q=eye(matrix_size);
matrix_Q(1,1)=30;
matrix_Q(2,2)=8;
matrix_Q(3,3)=13;
matrix_Q(4,4)=20;
matrix_R=2;
max_num_iteration=40;

%%
%计算状态量
dx=current_path_x-target_path_nearest_x;
dy=current_path_y-target_path_nearest_y;
%横向误差
lateral_error=cos(target_path_nearest_heading)*dy-sin(target_path_nearest_heading)*dx;
%航向误差
heading_error=current_path_heading-target_path_nearest_heading;
%横向误差变化率
lateral_error_rate=linear_v*sin(heading_error);
%航向变化率
heading_rate=angle_v;
target_heading_rate=target_path_nearest_curvature*linear_v;
%航向误差变化率
heading_error_rate=heading_rate-target_heading_rate;

%%
%更新状态矩阵State
matrix_State(1,1)=lateral_error;
matrix_State(2,1)=lateral_error_rate;
matrix_State(3,1)=heading_error;
matrix_State(4,1)=heading_error_rate;

%更新A的矩阵
matrix_A(2,2)=matrix_A_coeff(2,2)/linear_v;
matrix_A(2,4)=matrix_A_coeff(2,4)/linear_v;
matrix_A(4,2)=matrix_A_coeff(4,2)/linear_v;
matrix_A(4,4)=matrix_A_coeff(4,4)/linear_v;

%将A矩阵离散化
matrix_i=eye(size(matrix_A));
matrix_Ad=(inv(matrix_i-Ts_*0.5*matrix_A))*(matrix_i+Ts_*0.5*matrix_A);

%%
%求解反馈量
matrix_K(1,:)=SolveLQRProblem(matrix_Ad,matrix_Bd,matrix_Q,matrix_R,max_num_iteration);

matrix_u=-(matrix_K(1,:)*matrix_State(:,1));
tire_angle_feedback=matrix_u;

%求解前馈量
kv=lr_*mass_/2/cf_/wheelbase_-lf_*mass_/2/cr_/wheelbase_;
tire_angle_feedforward_1=wheelbase_*target_path_nearest_curvature;
tire_angle_feedforward_2=kv*linear_v*linear_v*target_path_nearest_curvature;
tire_angle_feedforward_3=lr_*target_path_nearest_curvature-lf_*mass_*linear_v*linear_v*target_path_nearest_curvature/2/cr_/wheelbase_;
tire_angle_feedforward_4=tire_angle_feedforward_1+tire_angle_feedforward_2-matrix_K(1,3)*tire_angle_feedforward_3;
tire_angle_feedforward=tire_angle_feedforward_4;
%合并
tire_angle_unlimited=tire_angle_feedback+tire_angle_feedforward;
% tire_angle_unlimited=rem(tire_angle_unlimited,2*pi);
tire_angle=tire_angle_unlimited;

if tire_angle<-tire_single_direction_max_degree
    tire_angle=-tire_single_direction_max_degree;
end
if tire_angle>tire_single_direction_max_degree
    tire_angle=tire_single_direction_max_degree;
end
 
% %滤波
% if t>0
%     if tire_angle-previous_tire_angle>filter_tire_angle
%         tire_angle=previous_tire_angle+filter_tire_angle;
%     end
%     if tire_angle-previous_tire_angle<-filter_tire_angle
%         tire_angle=previous_tire_angle-filter_tire_angle;
%     end
% end
%%
%输出
sys(1)=linear_v;
sys(2)=tire_angle;%tire_angle;

sys(3)=current_path_x;
sys(4)=current_path_y;
sys(5)=current_path_heading;

% sys(6)=target_path_x;
% sys(7)=target_path_y;
sys(6)=t;
sys(7)=0;

sys(8)=lateral_error;
sys(9)=lateral_error_rate;
sys(10)=heading_error;
sys(11)=heading_error_rate;

sys(12)=matrix_K(1,1);
sys(13)=matrix_K(1,2);
sys(14)=matrix_K(1,3);
sys(15)=matrix_K(1,4);

sys(16)=t;
sys(17)=index_min;
sys(18)=tire_angle_feedback;
sys(19)=tire_angle_feedforward;

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%


function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
