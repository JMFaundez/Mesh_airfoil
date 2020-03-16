clear all
close all

nelx = 165;
nely = 40;
input = 'fringe_m30.f00008';
m3 = base_case(input,nelx,nely);
[var, ind] = min(m3.xx(1,:));

nelx = 200;
nely = 34;
input = 'fringe_m40.f00008';
m4 = base_case(input,nelx,nely);
[var, ind2] = min(m4.xx(1,:));


figure(2000)
hold on
plot(m3.Re_s(ind:end),m3.p(1,ind:end))
plot(m4.Re_s(ind2:end),m4.p(1,ind2:end))
xlabel('$Re_s$','Interpreter','latex')
ylabel('P')
legend('Mesh 3','Mesh 4')
saveas(gca,['./figures/pressure.png'])

figure(2001)
contourf(m3.xx,m3.yy,m3.ut)
xlabel('x')
ylabel('y')
colorbar()
saveas(gca,['./figures/tan_velocity_m3.png'])

figure(2001)
contourf(m4.xx,m4.yy,m4.ut)
xlabel('x')
ylabel('y')
colorbar()
saveas(gca,['./figures/tan_velocity_m4.png'])

figure(2002)
hold on
plot(m3.Re_s(ind:end), m3.Uv(ind:end))
plot(m4.Re_s(ind2:end), m4.Uv(ind2:end))
xlabel('$Re_s$','Interpreter','latex')
ylabel("$U_{99}$",'Interpreter','latex')
legend('Mesh 3','Mesh 4')
saveas(gca,['./figures/U_inf.png'])

figure(2003)
hold on
plot(m3.Re_s(ind:end),m3.dth(ind:end))
plot(m4.Re_s(ind2:end),m4.dth(ind2:end))
plot(m3.Re_s(ind:end)*0+30300, m3.dth(ind:end))
ylabel('$\delta^*$','Interpreter','latex')
xlabel('$Re_s$','Interpreter','latex')
legend('Mesh 3','Mesh 4')
saveas(gca,['./figures/displacement_thick.png'])

figure(2004)
hold on
plot(m3.Re_s(ind:end),m3.y99(ind:end))
plot(m4.Re_s(ind2:end),m4.y99(ind2:end))
xlabel('$Re_s$','Interpreter','latex')
ylabel('$\delta$','Interpreter','latex')
legend('Mesh 3','Mesh 4')
saveas(gca,['./figures/y99.png'])

figure(2005)
hold on
plot(m3.Re_s(ind:end),m3.dUTdn(ind:end))
plot(m4.Re_s(ind2:end),m4.dUTdn(ind2:end))
xlabel('$Re_s$','Interpreter','latex')
ylabel('$\partial U_t \partial n$','Interpreter','latex')
legend('Mesh 3','Mesh 4')
saveas(gca,['./figures/dUtdn.png'])

figure(2006)
hold on
plot(m3.xx(1,ind:end),m3.Re_s(ind:end))
plot(m4.xx(1,ind2:end),m4.Re_s(ind2:end))
plot(m3.xx(1,ind:end),30300+0*m3.xx(1,ind:end))
xlabel('chord','Interpreter','latex')
ylabel('$Re_s$','Interpreter','latex')
legend('Mesh 3','Mesh 4')
saveas(gca,['./figures/Res_vs_chord.png'])
