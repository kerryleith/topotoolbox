function [x,y] = refmat2XY(R,siz, varargin)

    % Convert referencing matrix (refmat) to coordinate vectors
    %
    % Syntax
    %
    %     [x,y] = refmat2XY(R,siz)
    %
    % Input arguments
    %
    %     R     referencing matrix
    %     siz   size of the grid
    %
    % Output arguments
    %
    %     x     x-coordinates
    %     y     y-coordinates
    %
    %
    
    p = inputParser;
    
    addParameter(p,'mode','point');   
    
    parse(p,varargin{:});    
    params = p.Results;
    
    mode = params.mode;
    
    nrrows = siz(1);
    nrcols = siz(2);
    
    x = [ones(nrcols,1) (1:nrcols)' ones(nrcols,1)]*R;
    x = x(:,1)';
    if ~strcmp(mode, 'point')
        x = x-R(2,1); %kerry added -R(2,1)
    else
        1;
    end
    
    y = [(1:nrrows)' ones(nrrows,2)]*R;
    y = y(:,2); %kerry added -R(1,2)
    if ~strcmp(mode, 'point')
        y = y-R(1,2); %kerry added -R(2,1)
    end
    % fprintf('%0.1f%s%0.1f\n\n', R(3,1), ' ', R(3,2));
    % fprintf('%0.1f%s%0.1f\n', x(1), ' ', y(1));
    % fprintf('%0.1f%s%0.1f\n\n', x(end), ' ', y(end));
    
    
    