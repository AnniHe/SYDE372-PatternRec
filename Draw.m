classdef Draw
    methods(Static)
        function draw_circle(mu, colour)
            x = mu(1);
            y = mu(2);

            ang =0:0.01:2*pi;
            xp = 1*cos(ang);
            yp = 1*sin(ang);
            plot(x + xp, y + yp, colour);
            axis([-50 50 -50 50]);
            axis('square');
            hold on;
        end

        function draw_ellipse(x,y,theta,a,b, colour)

            if nargin<5, error('Too few arguments to Plot_Ellipse.'); end;

            np = 100;
            ang = [0:np]*2*pi/np;
            pts = [x;y]*ones(size(ang)) + [cos(theta) -sin(theta); sin(theta) cos(theta)]*[cos(ang)*a; sin(ang)*b];
            plot( pts(1,:), pts(2,:), colour);

        end
    end
    
end

