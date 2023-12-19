clear all
for beta=0.2:0.01:1
    % Set up an initial magnetic configuration
    N=30;
    S = ones(N,N);
    % Choose the number of MCMC steps
    itermax=10^8;
    iter=0;
    % Magnetization m of the initial configuration
    m=mag(S);
    mu=m;
    var=0;
    for iter=1: itermax
        % Randomly pick a site (i,j) and propose to flip the spin at it
        i=randi(N);
        j=randi(N);

        % Calculate the energy difference Delta H 
        Delta_H=2*S(i,j)*(S(mod(i,N)+1,j)+S(mod(i-2,N)+1,j)+S(i,mod(j,N)+1)+S(i,mod(j-2,N)+1));
        if Delta_H <0
            Accept=true;
        else
            % Generate random variable u of uniformly
            u=rand;
            if u<exp(-beta*Delta_H)
                Accept=true;
            else
                Accept=false;
            end
        end
        if Accept==true
            % flip the spin that was proposed to flip
            S(i,j)=-S(i,j);
            % Calculate the magnetization m of the new state
            m=mag(S);
        end
        % Update Variance and the mean magnetization
        var=(iter-1)/iter*var+1/(iter+1)*(m-mu)^2;
        mu=(iter*mu+m)/(iter+1);
    end
    plot(beta,mu,'b.');
    hold on;
    plot(beta,mu+sqrt(var),'k*');
    hold on;
    plot(beta,mu-sqrt(var),'k*');
    hold on;
    xlabel('Beta');
    ylabel('Mean Megatization');

end

fplot(@(x) f(x),[0.2 1],'r');


%% function of calculating the magnetization
function m = mag(S)
    N=size(S,1);
    m=sum(S,"all")/N^2;
end
%% function of the analytic expression of mean magnetization
function mu=f(beta)
    if beta>0.4408
        mu=(1-(sinh(2*beta))^(-4))^(1/8);
    end
    if beta< 0.4408
        mu=0;
    end
end



