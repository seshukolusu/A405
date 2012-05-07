function answer_entrain
    filename='littlerock.nc';
    fprintf('reading file: %s\n',filename);
    file_struct=nc_info(filename);
    c=constants;
    %
    % grab the March 2 12Z sounding
    %
    sound_var = file_struct.Dataset(4).Name;
    fprintf('found sounding: %s\n',sound_var);
    press=nc_varget(filename,sound_var,[0,0],[Inf,1]);
    height=nc_varget(filename,sound_var,[0,1],[Inf,1]);
    temp=nc_varget(filename,sound_var,[0,2],[Inf,1]);
    dewpoint=nc_varget(filename,sound_var,[0,3],[Inf,1]);
    envHeight=nudgeHeight(height);
    interpTenv=@(zVals) interp1(envHeight,temp,zVals);
    interpTdEnv=@(zVals) interp1(envHeight,dewpoint,zVals);
    interpPress=@(zVals) interp1(envHeight,press,zVals);
    p900_level=find(abs(900 - press) < 2.);
    p800_level=find(abs(800 - press) < 7.);
    thetaeVal=thetaep(dewpoint(p900_level) + c.Tc,temp(p900_level) + c.Tc,press(p900_level)*100.);
    height_800=height(p800_level);
    wTcloud=wsat(dewpoint(p900_level) + c.Tc,press(p900_level)*100.);
    yinit=[0.5,height_800,thetaeVal,wTcloud];
    tspan=0:10:2500;
    lambda=2.e-4;
    derivs=@(t,y) F(t,y,lambda,interpTenv,interpTdEnv,interpPress);
    options=odeset('Events',@stopIt);
    [t,y]=ode45(derivs,tspan,yinit,options);
    wvel=y(:,1);
    cloud_height=y(:,2); 
    thetae_cloud=y(:,3);
    wT_cloud=y(:,4);

    %limit the plot to the updraft
    figure(1);
    clf;
    plot(wvel,cloud_height)
    xlabel('vertical velocity');
    ylabel('height above surface');
    title('wvel vs. height')
    grid on;
    
    for i=1:numel(cloud_height)
      the_press=interpPress(cloud_height(i))*100.;
      [Tcloud(i),wvCloud(i),wlCloud(i)]= ...
         tinvert_thetae(thetae_cloud(i),wT_cloud(i),the_press);
    end
    for i=1:numel(cloud_height)
      the_press=interpPress(cloud_height(i))*100.;
      [Tadia(i),wvAdia(i),wlAdia(i)]= ...
         tinvert_thetae(thetae_cloud(1),wT_cloud(1),the_press);
    end

    figure(2);
    clf;
    plot(Tcloud - c.Tc,cloud_height,'r-')
    hold on;
    plot(temp,envHeight,'g-');
    plot(Tadia - c.Tc,cloud_height,'b-');
    xlabel('vertical velocity');
    ylabel('height above surface');
    title('wvel vs. height')
    grid on;
    hold off;
    
    function [value,isterminal,direction] = stopIt(t,y)
        % Locate the time when height passes through stopHeight in a decreasing direction
        % and stop integration.  Here we use a nested function to avoid
        % passing the additional parameter stopHeight as an input argument.
        value = y(1);     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = -1;   % only when approaching stopHeight from above
   end

end

function yp=F(t,y,lambda,interpTenv,interpTdEnv,interpPress)
  c=constants;
  yp=zeros(4,1); % since output must be a column vector
  %first variable is velocity in m/s
  height=y(2);
  thetae_cloud=y(3);
  wT_cloud=y(4);
  yp(1)=calcBuoy(height,thetae_cloud,wT_cloud,interpTenv,interpTdEnv,interpPress);
  press=interpPress(height)*100.;%hPa
  Tdenv=interpTdEnv(height) + c.Tc;
  Tenv=interpTenv(height) + c.Tc;
  wTenv=wsat(Tdenv,press); %kg/kg
  thetaeEnv=thetaep(Tdenv,Tenv,press);
  %second variable is height in m
  yp(2)= y(1);
  yp(3)=lambda*(thetaeEnv - y(3));
  yp(4)= lambda*(wTenv - y(4));
end

function Bout=calcBuoy(height,thetae0,wT0,interpTenv,interpTdEnv,interpPress)
    %calcBuoy(height,thetae0,wT0,interpTenv,interpTdenv,interpPress)
    %input: height (m), thetae0 (K), plus function handles for
    %T,Td, press soundings
    %output: Bout = buoyant acceleration in m/s^2
    %neglect liquid water loading in the virtual temperature
    c=constants;
    press=interpPress(height)*100.;%hPa
    [Tcloud,wvCloud,wlCloud]=tinvert_thetae(thetae0,wT0,press); %K
    %wvcloud=wsat(Tcloud,press); %kg/kg
    Tvcloud=Tcloud*(1. + c.eps*wvCloud - wlCloud);
    Tenv=interpTenv(height) + c.Tc;
    Tdenv=interpTdEnv(height) + c.Tc;
    wvenv=wsat(Tdenv,press); %kg/kg
    Tvenv=Tenv*(1. + c.eps*wvenv);
    TvDiff=Tvcloud - Tvenv;
    fprintf('%10.3f %10.3f %10.3f\n',press*0.01,height,TvDiff);
    Bout=c.g0*(TvDiff/Tvenv);
end


function envHeight=nudgeHeight(zVec)
    %if two balloon pressure levels are idential
    %add a factor of 0.1% to the second one
    %so interpolation will work
    envHeight=zVec;
    hit=find(diff(envHeight) < 0);
    envHeight(hit+1)=zVec(hit) + 1.e-3*zVec(hit);
end
