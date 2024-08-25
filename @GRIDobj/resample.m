function DEMr = resample(DEM,target,method,swapzone)

%RESAMPLE change spatial resolution of a GRIDobj
%
% Syntax
%
%     DEMr = resample(DEM,cellsize)
%     DEMr = resample(DEM,GRID)
%     DEMr = resample(...,method)
%     DEMr = resample(DEM,GRID,method,swapzone)
%
% Description
%
%     resample changes the cellsize of a grid. The function uses the MATLAB
%     function imtransform. If an instance of GRIDobj is supplied as
%     second argument, resample interpolates values in DEM to match the
%     spatial reference of GRID.
%
% Input arguments
%
%     DEM       grid object (GRIDobj)
%     cellsize  cellsize of resampled grid
%     GRID      other grid object
%     method    'bicubic', 'bilinear', or 'nearest' 
%     swapzone  true or false. If true and if DEM and GRID have different
%               projected coordinate systems, the function will attempt to
%               reproject and resample the DEM in one step. Note that this
%               requires the mapping toolbox. In case, the DEM is in a
%               geographic coordinate system, please use the function
%               reproject2utm(DEM,GRID).
%
% Output arguments
%
%     DEMr    grid object (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMr = resample(DEM,100);
%     imagesc(DEMr)
%
%
% See also: griddedInterpolant, imtransform
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. August, 2015 
% 22/08/24 kerry changed the refmat section

% check input arguments
narginchk(2,4)
validateattributes(target,{'double' 'GRIDobj'},{'scalar'})
if nargin == 2
    method = 'bilinear';
    swapzone = false;
elseif nargin == 3
    method = validatestring(method,{'bicubic', 'bilinear', 'nearest' });
    swapzone = false;
else
    method = validatestring(method,{'bicubic', 'bilinear', 'nearest' });
end

% check underlying class
if islogical(DEM.Z)
    method = 'nearest';
end

if swapzone && isa(target,'GRIDobj')
    if ~isequal(DEM.georef.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey,...
            target.georef.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
        warning('UTM zones differ. Will attempt to match.');
        swapzone = true;
    else
        swapzone = false;
    end
end

% Fillvalues
if isinteger(DEM.Z)
    fillval = 0;
elseif islogical(DEM.Z)
    fillval = 0;
else
    fillval = nan;
end


if isa(target,'GRIDobj')
    % the target is another GRIDobj
    % Create an affine2d object
    T = affine2d([1 0 0; 0 1 0; 0 0 1]);

    % check whether both grids have the same projection
    if swapzone
        mstructsource = DEM.georef.mstruct;
        mstructtarget = target.georef.mstruct;
        T = maketform('custom', 2, 2, ...
                @FWDTRANS, ...
                @INVTRANS, ...
                []);
    end
    rat = DEM.refmat ./target.refmat;
    rat = [rat(2,1) rat(1,2) 1];
    DEMr    = target;

    x_celldim = target.refmat(2,1);
    y_celldim = target.refmat(1,2);
    
    if ~isempty(DEM.Z)
        %T.T = T.T .* rat;
        T.T = T.T .* rat;
        outputSize = round(DEMr.size .* rat(1:2));
        DEMr.Z = imwarp(DEM.Z, T, 'OutputView', imref2d(outputSize), 'interp', method, ...
            'FillValues', fillval);
    end

    DEMr.name = [DEM.name ' (resampled)'];
        
else   
    DEMr    = GRIDobj([]);
    rat = DEM.cellsize/target;
    x_celldim =target *sign(DEM.refmat(2,1));
    y_celldim =target *sign(DEM.refmat(1,2));

    if ~isempty(DEM.Z)
        T = affine2d([1 0 0; 0 1 0; 0 0 1]);
        T.T = T.T .* [rat rat 1]'; 
        outputSize = round(DEM.size .* rat); 
        DEMr.Z = imwarp(DEM.Z, T, 'OutputView', imref2d(outputSize), 'interp', method, ...
            'FillValues', fillval); 
    end   
    
end



if ~isempty(DEM.georef) 

    DEMr.georef = DEM.georef;

%     [XWorldLimits, YWorldLimits,... 
%      x_offset_cells, y_offset_cells,...
%      width_cells, height_cells] = GRIDobj_alignDimensions(  DEM.georef.SpatialRef.XWorldLimits, DEM.georef.SpatialRef.YWorldLimits,...
%                                                             DEM.georef.SpatialRef.XWorldLimits, DEM.georef.SpatialRef.YWorldLimits,...
%                                                             x_celldim,      y_celldim,...
%                                                             DEM.georef.SpatialRef.RowsStartFrom, DEM.georef.SpatialRef.ColumnsStartFrom);

     [XWorldLimits, YWorldLimits,... 
     x_offset_cells, y_offset_cells,...
     width_cells, height_cells] = GRIDobj_alignDimensions( DEM.georef.SpatialRef.XWorldLimits, DEM.georef.SpatialRef.YWorldLimits,...
                                                            DEM.georef.SpatialRef.XWorldLimits, DEM.georef.SpatialRef.YWorldLimits,...
                                                            x_celldim,      y_celldim,...
                                                            DEM.refmat(2,1),      DEM.refmat(1,2),...
                                                            DEM.georef.SpatialRef.RowsStartFrom, DEM.georef.SpatialRef.ColumnsStartFrom,...
                                                            0);


    DEMr.georef.SpatialRef = maprefcells(XWorldLimits,YWorldLimits,size(DEMr.Z), ...
        'ColumnsStartFrom', DEM.georef.SpatialRef.ColumnsStartFrom, 'RowsStartFrom', DEM.georef.SpatialRef.RowsStartFrom);
    
    if strcmp(DEMr.georef.SpatialRef.ColumnsStartFrom, 'north')
        yn = 2;
    else
        yn = 1;
    end

    if strcmp(DEMr.georef.SpatialRef.RowsStartFrom, 'west')
        xn = 1;
    else
        xn = 2;
    end

    yr = YWorldLimits(yn);
    xr = XWorldLimits(xn);
    
    DEMr.refmat = zeros(3,2);
    DEMr.refmat(1,2) = y_celldim;
    DEMr.refmat(2,1) = x_celldim;
    DEMr.refmat(3,:) = [xr yr] - [DEMr.refmat(2,1) DEMr.refmat(1,2)]/2; 
    %DEMr.refmat(3,:) = [xr yr]; %22/08/24 kerry changed

    %DEMr.georef.RefMatrix = DEMr.refmat;
    % size of the resampled grid           
    DEMr.size    = size(DEMr.Z);
    DEMr.cellsize = abs(DEMr.georef.SpatialRef.CellExtentInWorldX);
    %DEMr.georef.Height = DEMr.size(1);
    %DEMr.georef.Width  = DEMr.size(2);
    
end
DEMr.name    = [DEM.name ' (resampled)'];

fprintf('%s%s%s%s%s\n', string(datetime(now,'ConvertFrom','datenum')), ' ', mfilename, ': ', 'Resampling DEM from:')
fprintf('%0.1f%s%0.1f\n', DEM.refmat(2,1), ' ', DEM.refmat(1,2));
fprintf('%0.1f%s%0.1f\n', DEM.size(1), ' ', DEM.size(2));
fprintf('%0.1f%s%0.1f\n\n', DEM.refmat(3,1), ' ', DEM.refmat(3,2));
fprintf('%0.1f%s%0.1f\n', DEM.georef.SpatialRef.XWorldLimits(1), ' ', DEM.georef.SpatialRef.YWorldLimits(1));
fprintf('%0.1f%s%0.1f\n\n\n', DEM.georef.SpatialRef.XWorldLimits(2), ' ', DEM.georef.SpatialRef.YWorldLimits(2));

fprintf('%s%s%s%s%s\n', string(datetime(now,'ConvertFrom','datenum')), ' ', mfilename, ': ', 'Resampled DEM to:')
fprintf('%0.1f%s%0.1f\n', DEMr.refmat(2,1), ' ', DEMr.refmat(1,2));
fprintf('%0.1f%s%0.1f\n', DEMr.size(1), ' ', DEMr.size(2));
fprintf('%0.1f%s%0.1f\n\n', DEMr.refmat(3,1), ' ', DEMr.refmat(3,2));
fprintf('%0.1f%s%0.1f\n', DEMr.georef.SpatialRef.XWorldLimits(1), ' ', DEMr.georef.SpatialRef.YWorldLimits(1));
fprintf('%0.1f%s%0.1f\n\n\n', DEMr.georef.SpatialRef.XWorldLimits(2), ' ', DEMr.georef.SpatialRef.YWorldLimits(2));




%%

    function x = FWDTRANS(u,~)
        % invtrans first
        [lati,long] = minvtran(mstructsource,u(:,1),u(:,2));
        [x,y] = mfwdtran(mstructtarget,lati,long);
        x = [x y];        
    end

    function u = INVTRANS(x,~)
        [lati,long] = minvtran(mstructtarget,x(:,1),x(:,2));
        [x,y] = mfwdtran(mstructsource,lati,long);
        u = [x y];   
        
    end
end






