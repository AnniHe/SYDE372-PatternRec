classdef Parzen
    methods (Static)
        function p = plot(sample, sigma)

            if size(sample, 2) == 2
                n = 2;
                sigma = cov(sample);
                range = [min(sample) - 1 : .1 ; max(sample) +1];
                
            else
                n = 1;
            end

            N = length(sample); %number of samples
            parzen = [];
            K = 100 %optimize this parameter
            h = K/sqrt(N); %width of window, formula in notes

            % height = 1/(N) * 1/width * 1;
            range = [0:.1:10];

            for x = range
                sum = 0;
                %sum the sample constributions at each discrete step
                for i = 1:N
                    mean = sample(i,:);  
                    sum = sum + 1/h* evalGuassian(x, mean', sigma, n);
                end

                %keep track of the y (1D) z (2D) vals at each discrete step
                parzen = [parzen 1/((N))*sum];
            end
            figure();
            plot(range, parzen);
            hold on;
            
            
        end
        
        function p2d = plot2D(al,bl,cl,sigma)
            %Create a matrix window with a Gaussian shape
            %http://dali.feld.cvut.cz/ucebna/matlab/toolbox/images/fspecial.html
            window = fspecial('gaussian', [100 100], sigma);

            %Resolution (step size)
            %Determines the spatial step between PDF estimates
            lowx = min([min(al(:,1)), min(bl(:,1)), min(cl(:,1))]);
            lowy = min([min(al(:,2)), min(bl(:,2)), min(cl(:,2))]);
            highx = max([max(al(:,1)), max(bl(:,1)), max(cl(:,1))]);
            highy = max([max(al(:,2)), max(bl(:,2)), max(cl(:,2))]);
            step = 1;
            resolution = [step (lowx) (lowy) (highx) (highy)];

            [pdfA, xA, yA] = Parzen.parzen2D(al, resolution, window);
            [pdfB, xB, yB] = Parzen.parzen2D(bl, resolution, window);
            [pdfC, xC, yC] = Parzen.parzen2D(cl, resolution, window);

            %Create a grid
            grid = zeros(length(xA),length(yA));

            %Classify each square on grid using ML
            for i = 1:length(xA)-1
                for j = 1:(length(yA)-1)
                    grid(i,j) = ML.classify(i,j, pdfA, pdfB, pdfC);
                end
            end
            
            %Plot points and ML boundary
            figure();
            hold on;
            contour(xA,yA,grid');
            scatter(al(:,1),al(:,2));
            scatter(bl(:,1),bl(:,2));
            scatter(cl(:,1),cl(:,2)); 
            hold off;
        end
        
        function [p,x,y] = parzen2D(data, res, win)
            %Source: Paul Fieguth (Almighty God)
            if (size(data,2)>size(data,1)), data = data'; end;
            if (size(data,2)==2), data = [data ones(size(data))]; end;
            numpts = sum(data(:,3));

            dl = min(data(:,1:2));
            dh = max(data(:,1:2));
            if length(res)>1, dl = [res(2) res(3)]; dh = [res(4) res(5)]; res = res(1); end;

            if (nargin == 2), win = 10; end;
            if (max(dh-dl)/res>1000), 
              error('Excessive data range relative to resolution.');
            end;

            if length(win)==1, win = ones(1,win); end;
            if min(size(win))==1, win = win(:) * win(:)'; end;
            win = win / (res*res*sum(sum(win)));

            p = zeros(2+(dh(2)-dl(2))/res,2+(dh(1)-dl(1))/res);
            fdl1 = find(data(:,1)>dl(1));
            fdh1 = find(data(fdl1,1)<dh(1));
            fdl2 = find(data(fdl1(fdh1),2)>dl(2));
            fdh2 = find(data(fdl1(fdh1(fdl2)),2)<dh(2));

            for i=fdl1(fdh1(fdl2(fdh2)))',
              j1 = round(1+(data(i,1)-dl(1))/res);
              j2 = round(1+(data(i,2)-dl(2))/res);
              p(j2,j1) = p(j2,j1) + data(i,3);
            end;

            p = conv2(p,win,'same')/numpts;
            x = [dl(1):res:(dh(1)+res)]; x = x(1:size(p,2));
            y = [dl(2):res:(dh(2)+res)]; y = y(1:size(p,1));

        end
    end     
end
        
function result = evalGuassian(x, mean, sigma, n)
    if n == 2
        sigmaMagnitude = det(sigma)^(.5);
        result = (1/(sigmaMagnitude*sqrt(2*pi)))*exp(-1/2*(inv(sigma)*(x-mean).^2));
    else
        result = (1/(sigma*sqrt(2*pi)))*exp(-1/(2*(sigma^2))*((x-mean)^2));
    end
    
end


        

