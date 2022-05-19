clear all;


%%%%%OUR APPROACH (Safety QP with Stability Penalty)
%%%%%%%%%%%%%%%%%

epsilons=[0.05];

thetas=linspace(0,2*pi,10);
r=8;

% x0s=[5,2,0.1];
% y0s=[1,6,8];
x0s=r*cos(thetas);
y0s=r*sin(thetas);

x0s=[x0s,0,0];
y0s=[y0s,9,1.5];

npoints=length(x0s);

hold on;

for l=1:npoints
	l
	for epsilon=epsilons

		x0=x0s(l);
		y0=y0s(l);

		[t,xvalues]=ode45(@(t,x)closedloop(t,x,epsilon),[0 10],[x0;y0]);

		if(l~=1)
			s='off';
		else
			s='on';
		end
		plot(xvalues(:,1),xvalues(:,2),'Linewidth',2,'Color','blue');

	end
end


% % %%%%%%QP AMES (CLF-CBF QP)
% % %%%%%%%%%%%%%

ps=[1];

thetas=linspace(0,2*pi,10);
r=8;

% x0s=[5,2,0.1];
% y0s=[1,6,8];

x0s=r*cos(thetas);
y0s=r*sin(thetas);

x0s=[x0s,0,0];
y0s=[y0s,9,1.5];

npoints=length(x0s);

hold on;

for l=1:npoints
	l
	for p=ps

		x0=x0s(l);
		y0=y0s(l);

		[t,xvalues]=ode45(@(t,x)qpcl(t,x,p),[0 50],[x0;y0]);

		if(l~=1)
			s='off';
		else
			s='on';
		end

		plot(xvalues(:,1),xvalues(:,2),'Linewidth',2,'Color','red');
	end
end

%%%%Modified QP Ames - Tan & Dimarogonas (M-CLF-CBF QP)
%%%%
%%%%Note: to use this one you need to change the function f

ps=[1];

thetas=linspace(0,2*pi,10);
r=8;

% x0s=[5,2,0.1];
% y0s=[1,6,8];

x0s=r*cos(thetas);
y0s=r*sin(thetas);

x0s=[x0s,0,0];
y0s=[y0s,9,1.5];

npoints=length(x0s);

hold on;

for l=1:npoints
	l
	for p=ps

		x0=x0s(l);
		y0=y0s(l);

		[t,xvalues]=ode45(@(t,x)qpcl(t,x,p),[0 50],[x0;y0]);

		if(l~=1)
			s='off';
		else
			s='on';
		end

		plot(xvalues(:,1),xvalues(:,2),'Linewidth',2,'Color','cyan','HandleVisibility',s);
	end
end


%%%%%%%%%Univeral formula for smooth safe stabilization
%%% (To do)


%%%%%%PLOTS

% xs=-2:0.01:2;
% plot(xs,4+sqrt(4-xs.^2),'Color','black')
% plot(xs,4-sqrt(4-xs.^2),'Color','black')
circle(0,4,2,'green');
hold on;
plot(0,0,'x','Markersize',14,'Color','black','Linewidth',2);
scatter(x0s,y0s,50,'black','filled');
xlabel('$x_1$','Interpreter','Latex','Fontsize',16);
ylabel('$x_2$','Interpreter','Latex','Fontsize',16);
xlim([-8,17]);
ylim([-9,17]);
pbaspect([1 1 1]);
yticks([-8,-6,-4,-2,0,2,4,6,8,10,12,14])
xticks([-8,-6,-4,-2,0,2,4,6,8,10,12,14]);
%region of attraction
xroa=-2:0.01:2;
yroa=sqrt(4-xroa.^2);
plot(xroa,yroa,':','color',[0.8500 0.3250 0.0980],'Linewidth',2);
plot(xroa,-yroa,':','color',[0.8500 0.3250 0.0980],'Linewidth',2);
%
ax=gca;
ax.FontSize=14
%legend(strcat('$\epsilon$=',num2str(epsilons(1))),strcat('$\epsilon$=',num2str(epsilons(2))),strcat('$\epsilon$=',num2str(epsilons(3))),'Location','northeastoutside','Fontsize',14,'Interpreter','latex');
%%Plots to have the right legend
h1=plot(-1,-1,'color','blue','Linewidth',2);
h2=plot(-1,-1,'color','red','Linewidth',2);
h3=plot(-1,-1,'color','cyan','Linewidth',2);
h4=plot(-1,-1,'color','green','Linewidth',2);
h5=plot(-1,-1,':','color',[0.8500 0.3250 0.0980],'Linewidth',2);
%%
legend([h1,h2,h3,h4,h5],'$\textrm{Safety QP with Stability Penalty}$','$\textrm{CLF-CBF QP}$','\textrm{M-CLF-CBF QP}','\textrm{$R^2 \backslash \mathcal{C}$}','\textrm{Region of attraction}','Interpreter','Latex','Fontsize',14,'Location','northeast')
hold off;


%%%FUNCTIONS%%%%%%%%

function uqp=uqp(x,p)

	set1=((Fv(x)<0) && (Fh(x)>0));
	set2=((Fv(x)<0) && Fh(x)==0 && all(Lgh(x)==0));
	set3=((Fh(x)<=0) && (Fv(x)*Lgh(x)*Lgh(x)'-Fh(x)*LgV(x)*Lgh(x)'<0));
	set4=((Fv(x)>=0) && (Fv(x)*Lgh(x)*LgV(x)'-Fh(x)*(1/p+LgV(x)*LgV(x)')<0));
	set5=((Fv(x)>=0) && Fh(x)==0 && all(Lgh(x)==0));
	set6=((Fv(x)*Lgh(x)*Lgh(x)'-Fh(x)*LgV(x)*Lgh(x)'>=0) && (Fv(x)*LgV(x)*Lgh(x)'-Fh(x)*(1/p+LgV(x)*LgV(x)')>=0) && all(Lgh(x)~=0));

	% x
	% LgV(x)*Lgh(x)'
	% [1/p+LgV(x)*LgV(x)' -LgV(x)*Lgh(x)';-LgV(x)*Lgh(x)' Lgh(x)*Lgh(x)']

	if(set1 || set2)
		uqp=[0;0];
	elseif(set3)
		uqp=-Fh(x)/(Lgh(x)*Lgh(x)')*Lgh(x)';
	elseif(set4 || set5)
		uqp=-Fv(x)/(1/p+LgV(x)*LgV(x)')*LgV(x)';
	%elseif(set6)
	else
		%v=inv([1/p+LgV(x)*LgV(x)' -LgV(x)*Lgh(x)';-LgV(x)*Lgh(x)' Lgh(x)*Lgh(x)'])*[Fv(x);-Fh(x)];
		%v=linsolve([1/p+LgV(x)*LgV(x)' -LgV(x)*Lgh(x)';-LgV(x)*Lgh(x)' Lgh(x)*Lgh(x)'],[Fv(x);-Fh(x)]);
		v=[1/p+LgV(x)*LgV(x)' -LgV(x)*Lgh(x)';-LgV(x)*Lgh(x)' Lgh(x)*Lgh(x)']\[Fv(x);-Fh(x)];
		%determ=norm(Lgh(x))^2/p+norm(LgV(x))^2*norm(Lgh(x))^2-(LgV(x)*Lgh(x)')^2;
		%v=1/determ*[norm(Lgh(x))^2 LgV(x)*Lgh(x)';LgV(x)*Lgh(x)' 1/p+norm(LgV(x))^2]*[Fv(x);-Fh(x)];
		v1=v(1);
		v2=v(2);
		uqp=-v1*LgV(x)'+v2*Lgh(x)';
	end

end

%%%%%%%%%%%%%%%%%%%%%

function dxdt=qpcl(t,x,p)
	dxdt=f(x)+g(x)*uqp(x,p);
end

%%%%%Penalty

function dxdt=closedloop(t,x,epsilon)
	dxdt=f(x)+g(x)*usafe(x,epsilon);
end

function usafe=usafe(x,epsilon)
	A=dot(gradh(x),f(x))-1/epsilon*gradh(x)'*g(x)*g(x)'*gradV(x)+alph(h(x));
	if(A>=0)
		usafe=-1/epsilon*g(x)'*gradV(x);
	else
		usafe=-1/epsilon*g(x)'*gradV(x)-A/norm(g(x)'*gradh(x))^2 *g(x)'*gradh(x);
	end
end

%%%%%%

function Fh=Fh(x)
	Fh=Lfh(x)+alph(h(x));
end

function Fv=Fv(x)
	Fv=LfV(x)+gam(V(x));
end

function LfV=LfV(x)
	LfV=gradV(x)'*f(x);
end

function LgV=LgV(x)
	LgV=gradV(x)'*g(x);
end

function Lfh=Lfh(x)
	Lfh=gradh(x)'*f(x);
end

function Lgh=Lgh(x)
	Lgh=gradh(x)'*g(x);
end

function gam=gam(x)
	gam=x;
end

function alph=alph(x)
	alph=x;
end

function gradh=gradh(x)
	%gradh=[-2*x(1);-8*x(2)];
	gradh=[2*x(1);2*(x(2)-4)];
end

function gradV=gradV(x)
	gradV=x;
end

function V=V(x)
	V=1/2*(x(1)^2+x(2)^2);
end

function h=h(x)
	h=x(1)^2+(x(2)-4)^2-4;
end


function f=f(x)
	%%Use this one for our approach and Ames QP
	%f=[x(1);x(2)];

	%%Use this one for M-CLF-CBF QP
	f=[x(1);x(2)]+g(x)*[-2*x(1);-2*x(2)];
end

function g=g(x)
	%g=[0;1];
	g=[1 0;0 1];
end


%%%%VIZ

function circles = circle(x,y,r,c)
	hold on
	th = 0:pi/50:2*pi;
	x_circle = r * cos(th) + x;
	y_circle = r * sin(th) + y;
	circles = plot(x_circle, y_circle);
	fill(x_circle, y_circle, c)
	hold off
end




