%Akul Kesarwani (akesarw3)
%MATLAB R2021b
close all;clear;clc;
rng('shuffle')

loaded_data = importdata('Coords.dat');
coords = loaded_data.data;
beads = loaded_data.textdata(2:end, 1);

[t1, ke1, pe1, temp1, configs ] = ld( coords, beads, 1, 0.003, 500001);
plot( t1, ke1 )
hold on
plot( t1, pe1 )
plot( t1, ke1+pe1 )
legend( 'Kinetic Energy', 'Potential Energy', 'Total internal Energy', Location='best')


function [time_var, KE, PE, temp_hist, totconfig] = ld( coords, beads, init_temp, dt, simulation_iter)

    %temp is in reduced units
    %coords should be in reduced units
    n = length(coords);
    totconfig = zeros( simulation_iter*n, 3);

    totconfig(1:n, :) = coords;
    
    time_var = [0:dt:(simulation_iter-1)*dt]';
    KE = zeros(simulation_iter, 1);
    PE = zeros(simulation_iter, 1);
    temp_hist = zeros(simulation_iter, 1);
    vel = init_vel( init_temp, n );
    
    [zeta, KE(1), temp_hist(1)] = calc_z( vel, n );
    [fc, PE(1)] = calc_fc( coords, beads, n );
    rand_force = randn(n, 3)* sqrt(2*temp_hist(1)*0.05/dt);
    total_f = fc + rand_force + zeta;

    for i = [2:simulation_iter]
        coords = update_pos( coords, vel, total_f, dt );
        cstart = n*i - n + 1;
        cend = n*i;
        totconfig(cstart:cend,:) = coords;

        [fc, PE(i)] = calc_fc( coords, beads, n );

        

        rand_force = randn(n, 3)* sqrt(2*temp_hist(i-1)*0.05/dt);
        vel = calc_vel( vel, dt, fc, rand_force, total_f);
        [zeta, KE(i), temp_hist(i)] = calc_z( vel, n );
        total_f = fc + rand_force + zeta;
    end
end

function vel = init_vel( temp, n )
    vel = randn( n, 3 ) * sqrt(temp);
    vel = vel - (sum(vel)/n);
end


function [zeta, cur_ke, cur_temp] = calc_z( vel, n )
    zeta = vel*-0.05;
    vel_sq = sum(vel .* vel, 2);
    cur_ke = (sum(vel_sq))/2;
    cur_temp = (sum(vel_sq))/((3*n-3));
end

function coords = update_pos(in_coords, vel, forces, dt)
    coords = in_coords + dt * vel + dt^2*forces*0.5;
end

function vel = calc_vel( in_vel, dt, fc, rand_force, total_f)
    const_m = 1/(1+dt*0.05/2);
    vel = in_vel + 0.5*dt*( total_f + fc + rand_force );
    vel = const_m*vel;
end


function [fc, pe] = calc_fc( c, beads, n )
    fc = zeros( n, 3 );
    pe = 0;
    for i = [1:n-1]
        %calculate bonded with i and i+1
        dx = c(i, 1)-c(i+1, 1);
        dy = c(i, 2)-c(i+1, 2);
        dz = c(i, 3)-c(i+1, 3);
        r = sqrt( dx.^2 + dy.^2 + dz.^2 );
        pe = pe + 10*((r-1).^2);
        f_f = -20*(r-1)/r * [dx, dy, dz];
        fc(i,:) = fc(i,:) + f_f;
        fc(i+1,:) = fc(i+1,:) - f_f;        
        for j = [i+2:n]
            %calculate nonbonded with i and j
            eps = calc_eps( i, j, beads );
            dx = c(i, 1)-c(j, 1);
            dy = c(i, 2)-c(j, 2);
            dz = c(i, 3)-c(j, 3);
            r_sqr = dx*dx + dy*dy + dz*dz;
            pe = pe + 4*eps*((1/r_sqr)^6 - (1/r_sqr)^3);

            f_dir = 4*eps*(12/(r_sqr^7) - 6/(r_sqr^4)) * [dx, dy, dz];
            fc(i, :) = fc(i,:) + f_dir;
            fc(j, :) = fc(j,:) - f_dir;
        end
    end
    cos_naut = -0.897934315;
    k_con = 40;
    fb = zeros( n, 3 );
    
    cos_arr = zeros(n, 1);

    for i = [2:n-1]
        h=i-1;
        j=i+1;
        r1 = c(h,:) - c(i,:);
        r2 = c(j,:) - c(i,:);
        l1 = sqrt(sum(r1 .* r1));
        l2 = sqrt(sum(r2.*r2));
        cos_theta = dot(r1, r2) ./ (l1*l2);
        cos_arr(i) = cos_theta;
        pe = pe + k_con/2 * ((cos_theta-cos_naut).^2);
        
        dcos = k_con * (cos_theta-cos_naut);

        d1 = r2 ./ (l1*l2) - cos_theta * r1 ./ (l1*l1);
        d2 = r1 ./ (l1*l2) - cos_theta * r2 ./ (l2*l2);

        fb(i,:) = fb(i,:) + dcos*d1 + dcos*d2;
        fb(h,:) = fb(h,:) - dcos*d1;
        fb(j,:) = fb(j,:) - dcos*d2;
    end
    fc = fc + fb;
    %disp(cos_arr);
end

function eps = calc_eps(i, j, beads )
    t1 = beads{i};
    t2 = beads{j};
    if length(t1) == 2
        t1 = t1(2);
    end
    if length(t2) == 2
        t2 = t2(2);
    end
    inp = [t1, t2];
    inp = sort(inp);
    eps = -1;
    if inp == 'HH'
        eps = 1;
    elseif inp ==  'HP' | inp == 'PP'
        eps = 2/3;
    elseif inp == '-H' | inp == '+H'
        eps = 1/2;
    elseif inp == '-P' | inp == '+P'
        eps = 1;
    elseif inp == '+-'
        eps = 1.5;
    elseif inp == '--' | inp=='++'
        eps = 0.25;
    else
        disp( 'error: no epsilon value')
    end
end