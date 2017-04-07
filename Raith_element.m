classdef Raith_element < handle
%    
% obj=Raith_element('polygon',layer,uv,DF)
% obj=Raith_element('path',layer,uv,w,DF)
% obj=Raith_element('dot',layer,uv,DF)
% obj=Raith_element('arc',layer,uv_c,r,theta,angle,w,N,DF)
% obj=Raith_element('circle',layer,uv_c,r,w,N,DF)
% obj=Raith_element('ellipse',layer,uv_c,r,w,angle,N,DF)
% obj=Raith_element('text',layer,uv_0,h,angle,uv_align,textlabel,DF)
% obj=Raith_element('fbmspath',layer,uv,cvtr,w,DF)
% obj=Raith_element('fbmscircle',layer,uv_c,r,w,DF)
% obj=Raith_element('sref',name,uv_0,[mag,angle,reflect])
% obj=Raith_element('aref',name,uv_0,n_colrow,a_colrow,[mag,angle,reflect])
%
% Raith_element objects define low-level elements used to generate 
% GDSII hierarchies for Raith beamwriting tools.
%
%
% Arguments:
%
% type - type of element: must be either 'polygon', 'path', 'dot', 'arc',
%   'circle', 'ellipse', 'text', 'fbmspath', 'fbmscircle', 'sref', or
%   'aref'
%
% The number and identity of the remaining arguments depends on type:
%
% 'polygon' (closed, filled polygon)
%       layer - GDSII layer for element (0-63)
%       uv - polygon vertices; 2 x n matrix [u;v] (um)
%       DF - dose factor for polygon
%
%   'path' (path of line segments)
%       layer - GDSII layer for element (0-63)
%       uv - path vertices; 2 x n matrix [u;v] (um)
%       w - width of path (um); value of zero yields single-pixel line; a
%           negative value denotes an absolute width (not affected by
%           magnification of any parent structure)
%       DF - dose factor for path
%
%   'dot' (single-pixel dot(s))
%       layer - GDSII layer for element (0-63)
%       uv - dot position(s); 2 x n matrix [u;v] (um)
%       DF - dose factor(s) for dot(s); if scalar, DF is applied to all 
%           dots specified by uv; if vector, must be of same length as uv
%
%   'arc' (segment of circular or elliptical path; Raith curved element)
%       layer - GDSII layer for element (0-63)
%       uv_c - arc centre; 1 x 2 vector [u_c v_c] (um)
%       r - radius of arc; may be scalar for a circular arc, or a 1 x 2
%           vector, [semi-major semi-minor] axes of an elliptical arc (um)
%       theta - starting and ending angles of arc w.r.t. axis defined by
%           angle argument; 1 x 2 vector [theta_1 theta_2] (degrees)
%       angle - angle of rotation between positive u-axis and theta = 0 
%           axis (degrees)
%       w - width of arc (um); if empty, arc is a filled elliptical disk
%           segment; if zero, arc is a single-pixel line; if non-zero, arc
%           has a width; a negative value denotes an absolute width (not 
%           affected by magnification of any parent structure)
%       N - number of vertices
%       DF - dose factor for arc
%
%   'circle' (circular Raith curved element)
%       layer - GDSII layer for element (0-63)
%       uv_c - circle centre; 1 x 2 vector [u_c v_c] (um)
%       r - radius of circle (um)
%       w - width of circle (um); if empty, circle is filled (disk); if 
%           zero, circle is a single-pixel line; if greater than zero,
%           circle has a width; a negative value denotes an absolute width 
%           (not affected by magnification of any parent structure)
%       N - number of vertices
%       DF - dose factor for circle
%
%   'ellipse' (elliptical Raith curved element)
%       layer - GDSII layer for element (0-63)
%       uv_c - ellipse centre; 1 x 2 vector [u_c v_c] (um)
%       r - 1 x 2 vector, [semi-major semi-minor] axes of ellipse (um)
%       w - width of ellipse (um); if empty, ellipse is filled (elliptical
%           disk); if zero, ellipse is a single-pixel line; if greater than 
%           zero, ellipse has a width; a negative value denotes an absolute
%           width (not affected by magnification of any parent structure)
%       angle - angle between semi-major axis and u axis (degrees)
%       N - number of vertices
%       DF - dose factor for ellipse
%
%   'text' (text rendered as simply-connected polygons)
%       layer - GDSII layer for element (0-63)
%       uv_0 - text anchor point [u_0 v_0] (um)
%       h - height of capital letters (um)
%       angle - angle of rotation of text w.r.t. positive u-axis (degrees)
%       uv_align - alignment w.r.t. anchor point; 1 x 2 vector 
%           [u_align v_align]; allowed values are 0 (left/top), 1 (centre),
%           2 (right/bottom) 
%       textlabel - the text to be written (string); allowed characters are
%           `1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~
%           !@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:"ZXCVBNM<>?µ and [space]
%       DF - dose factor for text
%
%   'fbmspath' (fixed beam moving stage path of line or arc segments)
%       layer - GDSII layer for element (0-63)
%       uv - path vertices; 2 x n matrix [u;v] (um)
%       cvtr - curvature of path segments (um); if scalar and zero, the path 
%           comprises line segments (no curvature); if a 1 x n vector, 
%           cvtr(k) yields a circular arc with chord endpoints of uv(:,k-1) 
%           and uv(:,k) such that the radial distance between the arc and 
%           the chord centre is cvtr(k); a positive (negative) value of 
%           cvtr(k) corresponds to an arc to the left (right) of the chord;
%           the value of cvtr(1) is ignored if cvtr is 1 x n
%       w - width of path (um); if zero, path is a single-pixel line (no 
%           beam raster); if greater than zero, circle has a width (beam 
%           rastered during stage motion)
%       DF - dose factor for path
%
%   'fbmscircle' (fixed beam moving stage circle)
%       layer - GDSII layer for element (0-63)
%       uv_c - circle centre; 1 x 2 vector [u_c v_c] (um)
%       r - radius of circle (um)
%       w - width of circle (um); if zero, circle is a single-pixel line 
%           (no beam raster); if greater than zero, circle has a width
%           (beam rastered during stage motion)
%       DF - dose factor for circle
%
%   'sref' (structure reference)
%   N.B.!  Transformations are applied in the following order: 1. scaling, 
%   mirroring; 2. rotation; 3. insertion.
%       name - name of structure being referenced (string)
%       uv_0 - structure origin; 1 x 2 vector [u_0 v_0] (um)
%       mag - magnification factor [optional]; default is no magnification
%           (mag = 1)
%       angle - angle of rotation, counter-clockwise positive (degrees) 
%           [optional]; default is no rotation (angle = 0)
%       reflect - Boolean flag for reflecting about u axis before other 
%           transformations [optional]; default is no reflection 
%           (reflect = 0)
%
%   'aref' (array reference)
%   N.B.!  Raith interprets aref objects differently than the GDSII
%   specification (e.g., as displayed using KLayout).  Given the number and
%   spacing of rows and columns, a lattice of instance origins is
%   generated, then rotation is applied to this lattice (if specified).  At
%   each of these lattice points, a structure is placed, after first being
%   scaled and/or rotated.
%       name - name of structure being referenced (string)
%       uv_0 - structure origin; 1 x 2 vector [u_0 v_0] (um)
%       n_colrow - 1 x 2 vector indicating number of columns and rows,
%           respectively
%       a_colrow - 1 x 2 vector indicating lattice spacing in rows and
%           columns, respectively (um)
%       mag - magnification factor [optional]; default is no magnification
%           (mag = 1)
%       angle - angle of rotation, counter-clockwise positive (degrees) 
%           [optional]; default is no rotation (angle = 0)
%       reflect - Boolean flag for reflecting about u axis before other 
%           transformations [optional]; default is no reflection 
%           (reflect = 0)
%
%
% Properties:
%
%   type - type of element:  'polygon', 'path', 'dot', 'arc', 'circle', 
%       'ellipse', 'text', 'fbmspath', 'fbmscircle', 'sref', or 'aref'
%   data - remaining record data for element (depends on element type)
%
%
% Methods:
%
%   plot([M,scDF]) - plot element with Raith dose factor colouring (filled
%       polygons where applicable)
%           M - augmented transformation matrix for plot [optional]
%           scDF - overall multiplicative scaling factor for DF specified
%               in obj.data.DF (e.g., passed from a positionlist entry)
%               [optional]
%
%   plotedges([M,scDF]) - plot element with Raith dose factor colouring 
%       (edges of polygons where applicable)
%           M - augmented transformation matrix for plot [optional]
%           scDF - overall multiplicative scaling factor for DF specified
%               in obj.data.DF (e.g., passed from a positionlist entry)
%               [optional]
%
%
% Aaron Hryciw
% 2013-03-07
%
% Version 1.2
% 2014-10-07
%
%
% The Raith_GDSII MATLAB toolbox was developed at the National Institute 
% for Nanotechnology (NINT), a joint initiative between the Government of 
% Canada, the Government of Alberta, the National Research Council (NRC), 
% and the University of Alberta.  If is currently maintained by the
% University of Alberta nanoFAB facility.
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0.  If a copy of the MPL was not distributed with this 
% file, you can obtain one at http://mozilla.org/MPL/2.0/.
%
    
    
	properties
        type
        data
    end
   
    
    methods
        
        function obj=Raith_element(Type,varargin)
            if nargin>0
                switch lower(Type)
                    
                    case {'polygon','dot'} % Arguments:  layer, uv, DF
                        if nargin~=4  % Check for correct number of arguments
                            error(['Construct a Raith_element ' Type ' as:  Raith_element(''' Type ''',layer,uv,DF).'])
                        else 
                            Data.layer=varargin{1};
                            Data.uv=varargin{2};  % Round uv to nearest nm (1 nm data grid for Raith)
                            Data.DF=varargin{3}; % Round DF to nearest thousandth (precision expected by Raith GDSII)
                        end % polygon
                        
                    case 'path' % Arguments:  layer, uv, w, DF
                        if nargin~=5
                            error('Construct a Raith_element path as:  Raith_element(''path'',layer,uv,w,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv=varargin{2};  
                            Data.w=varargin{3};  
                            Data.DF=varargin{4}; 
                        end % path
                        
                    case 'arc' % Arguments:  layer, uv_c, r, theta, angle, w, N, DF
                        if nargin~=9
                            error('Construct a Raith_element arc as:  Raith_element(''arc'',layer,uv_c,r,theta,angle,w,N,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv_c=varargin{2};
                            Data.r=varargin{3};
                            Data.theta=varargin{4};
                            Data.angle=varargin{5};
                            Data.w=varargin{6};
                            Data.N=varargin{7}; 
                            Data.DF=varargin{8};
                        end % arc
                        
                    case 'circle' % Arguments:  layer, uv_c, r, w, N, DF
                        if nargin~=7
                            error('Construct a Raith_element circle as:  Raith_element(''circle'',layer,uv_c,r,w,N,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv_c=varargin{2};
                            Data.r=varargin{3};
                            Data.w=varargin{4};
                            Data.N=varargin{5}; 
                            Data.DF=varargin{6};
                        end % circle
                        
                    case 'ellipse' % Arguments:  layer, uv_c, r, w, angle, N, DF
                        if nargin~=8
                            error('Construct a Raith_element ellipse as:  Raith_element(''ellipse'',layer,uv_c,r,w,angle,N,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv_c=varargin{2};
                            Data.r=varargin{3};
                            Data.w=varargin{4};
                            Data.angle=varargin{5};
                            Data.N=varargin{6}; 
                            Data.DF=varargin{7};
                        end % ellipse     
                        
                    case 'text' % Arguments:  layer, uv_0, h, angle, uv_align, textlabel, DF
                        if nargin~=8
                            error('Construct a Raith_element text as:  Raith_element(''text'',layer,uv_0,h,angle,uv_align,textlabel,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv_0=varargin{2};
                            Data.h=varargin{3};
                            Data.angle=varargin{4};
                            Data.uv_align=varargin{5};
                            Data.textlabel=varargin{6}; 
                            Data.DF=varargin{7};
                        end % text  
                        
                    case 'fbmspath' % Arguments:  layer, uv, cvtr, w, DF
                        if nargin~=6
                            error('Construct a Raith_element fbmspath as:  Raith_element(''fbmspath'',layer,uv,cvtr,w,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv=varargin{2};  
                            Data.cvtr=varargin{3};  
                            Data.w=varargin{4};  
                            Data.DF=varargin{5}; 
                        end % fbmspath
                        
                    case 'fbmscircle' % Arguments:  layer, uv_c, r, w, DF
                        if nargin~=6
                            error('Construct a Raith_element fbmscircle as:  Raith_element(''fbmscircle'',layer,uv_c,r,w,DF).')
                        else 
                            Data.layer=varargin{1};
                            Data.uv_c=varargin{2};  
                            Data.r=varargin{3};  
                            Data.w=varargin{4};  
                            Data.DF=varargin{5}; 
                        end % fbmscircle                        
                        
                    case 'sref' % Arguments:  name, uv_0, [mag, angle, reflect]
                        if nargin<3 || nargin>6
                            error('Construct a Raith_element structure reference as:  Raith_element(''sref'',name,uv_0,[mag,angle,reflect]).')
                        else 
                            Data.name=varargin{1};
                            Data.uv_0=varargin{2};
                            % Defaults for transformations
                            Data.mag=1;
                            Data.angle=0;
                            Data.reflect=0;
                            if nargin>3 % Mag specified
                                if ~isempty(varargin{3})
                                    Data.mag=varargin{3};
                                end
                            end
                            if nargin>4 % Mag and angle specified
                                if ~isempty(varargin{4})
                                    Data.angle=varargin{4};
                                end
                            end
                            if nargin>5 % Mag, angle, and reflect specified
                                if ~isempty(varargin{5})
                                    Data.reflect=varargin{5};
                                end
                            end
                        end % sref
                        
                    case 'aref' % Arguments:  name, uv_0, n_colrow, a_colrow, [mag, angle, reflect]
                        if nargin<5 && nargin>8
                            error('Construct a Raith_element array reference as:  Raith_element(''aref'',name,uv_0,n_colrow,a_colrow,[mag,angle,reflect]).')
                        else 
                            Data.name=varargin{1};
                            Data.uv_0=varargin{2};
                            Data.n_colrow=varargin{3};
                            Data.a_colrow=varargin{4};
                            % Defaults for transformations
                            Data.mag=1;
                            Data.angle=0;
                            Data.reflect=0;
                            if nargin>5 % Mag specified
                                if ~isempty(varargin{5})
                                    Data.mag=varargin{5};
                                end
                            end
                            if nargin>6 % Mag and angle specified
                                if ~isempty(varargin{6})
                                    Data.angle=varargin{6};
                                end
                            end
                            if nargin>7 % Mag, angle, and reflect specified
                                if ~isempty(varargin{7})
                                    Data.reflect=varargin{7};
                                end
                            end
                        end % aref
                        
                        
                    otherwise
                        error('Raith_element type must be ''polygon'', ''path'', ''dot'', ''circle'', ''ellipse'', ''text'', ''fbmspath'', ''fbmscircle'', ''sref'', or ''aref''.');
                        
                end
                
                obj.type=Type;
                obj.data=Data;
                        
            end
            
        end  % Constructor
        
        
        function set.type(obj,Type)

            Type=lower(Type);
            
            global checkdata;  % Global variable for argument checking

            if checkdata==false
                
                obj.type=Type;
            
            else

                if ~any(strcmp(Type,{'polygon','path','dot','arc','circle','ellipse','text','fbmspath','fbmscircle','sref','aref'}))
                    error('Raith_element type must be ''polygon'', ''path'', ''dot'', ''circle'', ''ellipse'', ''text'', ''fbmspath'', ''fbmscircle'', ''sref'', or ''aref''.')
                end

                obj.type=Type;

                % Fill data with default fields, depending on type
                switch obj.type
                    case {'polygon','dot'}
                        Data.layer=[];
                        Data.uv=[];
                        Data.DF=[];

                    case 'path'
                        Data.layer=[];
                        Data.uv=[];
                        Data.w=[];
                        Data.DF=[];                    

                    case 'arc'
                        Data.layer=[];
                        Data.uv_c=[];
                        Data.r=[];
                        Data.theta=[];
                        Data.angle=[];
                        Data.w=[];
                        Data.N=[];
                        Data.DF=[];                    

                    case 'circle'
                        Data.layer=[];
                        Data.uv_c=[];
                        Data.r=[];
                        Data.w=[];
                        Data.N=[];
                        Data.DF=[];                    

                    case 'ellipse'
                        Data.layer=[];
                        Data.uv_c=[];
                        Data.r=[];
                        Data.w=[];
                        Data.angle=[];
                        Data.N=[];
                        Data.DF=[];   

                    case 'text'
                        Data.layer=[];
                        Data.uv_0=[];
                        Data.h=[];
                        Data.angle=[];
                        Data.uv_align=[];
                        Data.textlabel=[];
                        Data.DF=[];   

                    case 'fbmspath'
                        Data.layer=[];
                        Data.uv=[];
                        Data.cvtr=[];
                        Data.w=[];
                        Data.DF=[];       
                        
                    case 'fbmscircle'
                        Data.layer=[];
                        Data.uv_c=[];
                        Data.r=[];
                        Data.w=[];
                        Data.DF=[];                          

                    case 'sref'
                        Data.name=[];
                        Data.uv_0=[];
                        Data.mag=[];
                        Data.angle=[];
                        Data.reflect=[];

                    case 'aref'
                        Data.name=[];
                        Data.uv_0=[];
                        Data.n_colrow=[];
                        Data.a_colrow=[];
                        Data.mag=[];
                        Data.angle=[];
                        Data.reflect=[];

                end

                obj.data=Data;
                
            end
                    
        end % set.type
        
        
        function set.data(obj,Data)
            
            % Check that type is already assigned
            if isempty(obj.type)
                error('Raith_element:  type must be set before data.');
            end
            
            if isempty(Data) % Clear all fields from obj.data
                obj.data=struct;
                return;
            end
            
            datafields=fieldnames(Data);  % Cell array of all fields
            
            global checkdata;  % Global variable for argument checking
            
            if checkdata==false  % Round data as necessary, but do not otherwise check

                for k=1:length(datafields)
                    
                    switch datafields{k}

                        case 'uv'
                            Data.uv=round(Data.uv*1000)/1000;  % Round uv to nearest nm (1 nm data grid for Raith)

                        case 'DF'
                            Data.DF=round(Data.DF*1000)/1000; % Round DF to nearest thousandth (precision expected by Raith GDSII)

                        case 'uv_c'
                            Data.uv_c=round(Data.uv_c*1000)/1000;  % Round uv_c to nearest nm (1 nm data grid for Raith)

                        case 'r'
                            Data.r=round(Data.r*1000)/1000;  % Round r to nearest nm (1 nm data grid for Raith)
                            
                        case 'cvtr'
                            Data.cvtr=round(Data.cvtr*1000)/1000;  % Round cvtr to nearest nm (1 nm data grid for Raith)

                        case 'w'
                            Data.w=round(Data.w*1000)/1000;  % Round w to nearest nm (1 nm data grid for Raith)

                        case 'h'
                            Data.h=round(Data.h*1000)/1000;  % Round h to nearest nm (1 nm data grid for Raith)

                        case 'uv_0'
                            Data.uv_0=round(Data.uv_0*1000)/1000;  % Round uv_0 to nearest nm (1 nm data grid for Raith)

                        case 'a_colrow'
                            Data.a_colrow=round(Data.a_colrow*1000)/1000;  % Round a_colrow to nearest nm (1 nm data grid for Raith)

                    end

                end
                
            else  % Check all data
            
                switch lower(obj.type)  % Check for illegal field(s)

                        case 'polygon' % data fields:  layer, uv, DF
                            if ~isempty(setdiff(datafields,{'layer','uv','DF'})) % Trying to set illegal field(s) for polygon
                                error('Raith_element:  allowed data fields for ''polygon'' element type are ''layer'', ''uv'', and ''DF''.');
                            end

                        case 'path' % data fields:  layer, uv, w, DF
                            if ~isempty(setdiff(datafields,{'layer','uv','w','DF'})) % Trying to set illegal field(s) for path
                                error('Raith_element:  allowed data fields for ''path'' element type are ''layer'', ''uv'', ''w'', and ''DF''.');
                            end

                        case 'dot' % data fields:  layer, uv, DF
                            if ~isempty(setdiff(datafields,{'layer','uv','DF'})) % Trying to set illegal field(s) for dot
                                error('Raith_element:  allowed data fields for ''dot'' element type are ''layer'', ''uv'', and ''DF''.');
                            end

                        case 'arc' % data fields:  layer, uv_c, r, theta, angle, w, N, DF
                            if ~isempty(setdiff(datafields,{'layer','uv_c','r','theta','angle','w','N','DF'})) % Trying to set illegal field(s) for arc
                                error('Raith_element:  allowed data fields for ''arc'' element type are ''layer'', ''uv_c'', ''r'', ''theta'', ''angle'', ''w'', ''N'', and ''DF''.');
                            end

                        case 'circle' % data fields:  layer, uv_c, r, w, N, DF
                            if ~isempty(setdiff(datafields,{'layer','uv_c','r','w','N','DF'})) % Trying to set illegal field(s) for circle
                                error('Raith_element:  allowed data fields for ''circle'' element type are ''layer'', ''uv_c'', ''r'', ''w'', ''N'', and ''DF''.');
                            end

                        case 'ellipse' % data fields:  layer, uv_c, r, w, angle, N, DF
                            if ~isempty(setdiff(datafields,{'layer','uv_c','r','w','angle','N','DF'})) % Trying to set illegal field(s) for ellipse
                                error('Raith_element:  allowed data fields for ''ellipse'' element type are ''layer'', ''uv_c'', ''r'', ''w'', ''angle'', ''N'', and ''DF''.');
                            end

                        case 'text' % data fields:  layer, uv_0, h, angle, uv_align, textlabel, DF
                            if ~isempty(setdiff(datafields,{'layer','uv_0','h','angle','uv_align','textlabel','DF'})) % Trying to set illegal field(s) for text
                                error('Raith_element:  allowed data fields for ''text'' element type are ''layer'', ''uv_0'', ''h'', ''angle'', ''uv_align'', ''textlabel'', and ''DF''.');
                            end
                            
                        case 'fbmspath' % data fields:  layer, uv, cvtr, w, DF
                            if ~isempty(setdiff(datafields,{'layer','uv','cvtr','w','DF'})) % Trying to set illegal field(s) for fbmspath
                                error('Raith_element:  allowed data fields for ''fbmspath'' element type are ''layer'', ''uv'', ''cvtr'', ''w'', and ''DF''.');
                            end
                            
                        case 'fbmscircle' % data fields:  layer, uv_c, r, w, DF
                            if ~isempty(setdiff(datafields,{'layer','uv_c','r','w','DF'})) % Trying to set illegal field(s) for fbmscircle
                                error('Raith_element:  allowed data fields for ''fbmscircle'' element type are ''layer'', ''uv_c'', ''r'', ''w'', and ''DF''.');
                            end  
                            
                        case 'sref' % data fields:  name, uv_0, mag, angle, reflect
                            if ~isempty(setdiff(datafields,{'name','uv_0','mag','angle','reflect','DF'})) % Trying to set illegal field(s) for sref
                                error('Raith_element:  allowed data fields for ''sref'' element type are ''layer'', ''uv_0'', ''mag'', ''angle'', and ''reflect''.');
                            end

                        case 'aref' % data fields:  name, uv_0, n_colrow, a_colrow, mag, angle, reflect
                            if ~isempty(setdiff(datafields,{'name','uv_0','n_colrow','a_colrow','mag','angle','reflect','DF'})) % Trying to set illegal field(s) for aref
                                error('Raith_element:  allowed data fields for ''aref'' element type are ''layer'', ''uv_0'', ''n_colrow'', ''a_colrow'', ''mag'', ''angle'', and ''reflect''.');
                            end

                end


                % Argument checking for all data fields
                for k=1:length(datafields)
                    f=getfield(Data,datafields{k});
                    if ~isempty(f)

                        % Check for infinities and NaNs
                        if any(isnan(f(:))) || any(isinf(f(:)))
                           error(['Raith_element ' obj.type ':  ' datafields{k} ' must be finite.'])
                        end

                        switch datafields{k}

                            case 'layer'
                                % Check that layer is an integer between 0 and 63 (may be stored as a float, but must have no fraction part)
                                if ~isnumeric(Data.layer) || floor(Data.layer)~=Data.layer || Data.layer<0 || Data.layer>63
                                    error(['Raith_element ' obj.type ':  layer must be an integer between 0 and 63 (inclusive).']);
                                end

                            case 'uv'
                                % Check size of uv
                                if ~isnumeric(Data.uv) || size(Data.uv,1)~=2
                                    error(['Raith_element ' obj.type ':  uv must a 2 x n matrix.']);
                                end

                                Data.uv=round(Data.uv*1000)/1000;  % Round uv to nearest nm (1 nm data grid for Raith)

                                switch obj.type
                                    case 'polygon'
                                    % Check for closedness and correct if not (with warning)
                                        if any(Data.uv(:,1)~=Data.uv(:,end))
                                            warning('Raith_element:openPolygon','Raith_element polygon is open; closing polygon.')
                                            Data.uv(:,end+1)=Data.uv(:,1);
                                        end
                                    case {'path','fbmspath'}
                                        % Check that uv contains at least two vertices
                                        if size(Data.uv,2)<2
                                            error(['Raith_element ' obj.type ':  uv must contain at least two vertices.']);
                                        end
                                end

                            case 'DF'
                                if strcmp(obj.type,'dot')
                                    % Check that DF is either scalar or a vector of the same length as uv
                                    if ~isnumeric(Data.DF) || ~isvector(Data.DF) || (length(Data.DF)~=size(Data.uv,2) && length(Data.DF)~=1)
                                        error('Raith_element dot:  DF must either be scalar or a vector of length size(data.uv,2).')
                                    end
                                else 
                                    % Check that DF is scalar
                                    if ~isnumeric(Data.DF) || ~isscalar(Data.DF)
                                        error(['Raith_element ' obj.type ':  DF must be a scalar.'])
                                    end
                                end
                                if any(Data.DF<0)
                                    error(['Raith_element ' obj.type ':  DF must be non-negative.']);
                                end
                                Data.DF=round(Data.DF*1000)/1000; % Round DF to nearest thousandth (precision expected by Raith GDSII)

                            case 'uv_c'
                                % Check that uv_c is a vector of length 2
                                if ~isnumeric(Data.uv_c) || ~isvector(Data.uv_c) || length(Data.uv_c)~=2
                                    error(['Raith_element ' obj.type ':  uv_c must be a vector of length 2.'])
                                end
                                Data.uv_c=round(Data.uv_c*1000)/1000;  % Round uv_c to nearest nm (1 nm data grid for Raith)

                            case 'r'
                                if strcmp(obj.type,'circle') || strcmp(obj.type,'fbmscircle')
                                    % Check that r is scalar
                                    if ~isnumeric(Data.r) || ~isscalar(Data.r) || Data.r<=0
                                        error(['Raith_element ' obj.type ':  r must be a positive scalar.'])
                                    end                      
                                elseif strcmp(obj.type,'ellipse')
                                    % Check that r is a vector of length 2
                                    if ~isnumeric(Data.r) || ~isvector(Data.r) || length(Data.r)~=2 || any(Data.r<=0)
                                        error('Raith_element ellipse:  r must be a vector of length 2, with positive elements.')
                                    end
                                else % Arc
                                    % Check that r is either scalar or a vector of length 2
                                    if ~isnumeric(Data.r) || ~isvector(Data.r) || ~any(length(Data.r)==[1 2]) || any(Data.r<=0)
                                        error('Raith_element arc:  r must be a vector of length 1 or 2, with positive elements.')
                                    end
                                end
                                Data.r=round(Data.r*1000)/1000;  % Round r to nearest nm (1 nm data grid for Raith)
                                
                            case 'cvtr'                                
                                % Check that cvtr is either 0 or a vector of the same length as uv
                                if ~isnumeric(Data.cvtr) || ~isvector(Data.cvtr) || (length(Data.cvtr)~=size(Data.uv,2) && length(Data.cvtr)~=1) || (isscalar(Data.cvtr) && Data.cvtr~=0)
                                    error('Raith_element fbmspath:  cvtr must either be 0 or a vector of length size(data.uv,2).')
                                end
                                Data.cvtr=round(Data.cvtr*1000)/1000;  % Round cvtr to nearest nm (1 nm data grid for Raith)

                            case 'w'
                                % Check that w is scalar (may be negative)
                                if ~isnumeric(Data.w) || ~isscalar(Data.w)
                                    error(['Raith_element ' obj.type ':  w must be a numeric scalar.'])
                                end
                                Data.w=round(Data.w*1000)/1000;  % Round w to nearest nm (1 nm data grid for Raith)

                            case 'h'
                                % Check that h is scalar and positive
                                if ~isnumeric(Data.h) || ~isscalar(Data.h) || Data.h<=0
                                    error('Raith_element text:  h must be a positive scalar.')
                                end
                                Data.h=round(Data.h*1000)/1000;  % Round h to nearest nm (1 nm data grid for Raith)

                            case 'N'
                                % Check that N is an integer greater than 2
                                if ~isnumeric(Data.N) || ~isscalar(Data.N) || round(Data.N)~=Data.N || Data.N<3
                                    error(['Raith_element ' obj.type ':  N must be an integer greater than 2.'])
                                end

                            case 'angle'
                                % Check that angle is scalar
                                if ~isnumeric(f) || ~isscalar(f)
                                    error(['Raith_element ' obj.type ':  angle must be a numeric scalar.'])
                                end

                            case 'name'
                                % Check that name is a string
                                if ~ischar(Data.name)
                                    error(['Raith_element ' obj.type ':  name must be a string.'])
                                end

                            case 'textlabel'
                                % Check that textlabel is a string
                                if ~ischar(Data.textlabel)
                                    error('Raith_element text:  textlabel must be a string.')
                                end
                                % Check for illegal characters
                                illegals=setdiff(Data.textlabel,obj.chars);
                                if ~isempty(illegals) % Trying to use illegal characters in textlabel
                                    if length(illegals)==1  % For unnecessary elegance, format error message according to number of illegal characters.
                                        ch='character ';
                                        be=' is';
                                    else
                                        ch='characters ';
                                        be=' are';
                                        illegals=sprintf('%c,',illegals(1:end-1));
                                        illegals=[illegals ' and ' illegals(end)];
                                    end
                                    error(['Raith_element text:  ' ch illegals be ' not allowed in textlabel.'])
                                end

                            case 'uv_0'
                                % Check that uv_0 is a vector of length 2
                                if ~isnumeric(Data.uv_0) || ~isvector(Data.uv_0) || length(Data.uv_0)~=2
                                    error(['Raith_element ' obj.type ':  uv_0 must be a vector of length 2.'])
                                end
                                Data.uv_0=round(Data.uv_0*1000)/1000;  % Round uv_0 to nearest nm (1 nm data grid for Raith)

                            case 'theta'
                                % Check that theta is a vector of length 2
                                if ~isnumeric(Data.theta) || ~isvector(Data.theta) || length(Data.theta)~=2
                                    error('Raith_element arc:  theta must be a vector of length 2.')
                                end

                            case 'mag'
                                % Check that mag is a positive scalar
                                if ~isnumeric(Data.mag) || ~isscalar(Data.mag) || Data.mag<=0
                                    error(['Raith_element ' obj.type ':  mag must be a positive scalar.'])
                                end

                            case 'reflect'
                                % Check that reflect is either 0 or 1
                                if Data.reflect~=1 && Data.reflect~=0
                                    error(['Raith_element ' obj.type ':  reflect must be either 0 or 1.'])
                                end

                            case 'uv_align'
                                % Check that uv_align is a vector of length 2
                                if ~isnumeric(Data.uv_align) || ~isvector(Data.uv_align) || length(Data.uv_align)~=2
                                    error('Raith_element text:  uv_align must be a vector of length 2.')
                                end
                                % Check that uv_align is either 0, 1, or 2
                                if ~isempty(setdiff(Data.uv_align,[0 1 2]))
                                    error('Raith_element text:  allowed values for uv_align elements are 0, 1, and 2.')
                                end

                            case 'n_colrow'
                                % Check that n_colrow is a vector of length 2
                                if ~isnumeric(Data.n_colrow) || ~isvector(Data.n_colrow) || length(Data.n_colrow)~=2
                                    error('Raith_element aref:  n_colrow must be of the form [n_columns n_rows].')
                                end
                                % Check that n_colrow elements are integers (though they may be stored as floats)
                                if ~all(round(Data.n_colrow)==Data.n_colrow)
                                    error('Raith_element aref:  n_colrow elements must be integers.')
                                end

                            case 'a_colrow'
                                % Check that a_colrow is a vector of length 2
                                if ~isnumeric(Data.a_colrow) || ~isvector(Data.a_colrow) || length(Data.a_colrow)~=2
                                    error('Raith_element aref:  a_colrow must be of the form [column_spacing row_spacing].')
                                end
                                Data.a_colrow=round(Data.a_colrow*1000)/1000;  % Round a_colrow to nearest nm (1 nm data grid for Raith)

                        end

                    end

                end
                
            end
                
            obj.data=Data;
            
        end % set.data
        
        
        
         function plot(obj,varargin)
        %
        % Raith_element.plot([M,scDF]) 
        %
        % Plot element with Raith dose factor colouring (filled polygons where applicable)
        %
        % Argument:
        %
        %   M - augmented transformation matrix for plot [optional]
        %   scDF - overall multiplicative scaling factor for DF specified
        %       in obj.data.DF (e.g., passed from a positionlist entry)
        %       [optional]
        %
            
            if nargin==1 
                M=[];  
                scDF=1;  % Scaling for DF specified in obj.data.DF
            elseif nargin==2
                M=varargin{1};
                scDF=1;  
            elseif nargin==3
                M=varargin{1};
                scDF=varargin{2};  
            else
                error('Raith_element.plot:  Too many input arguments.');
            end
            
            if isempty(M)
                M=eye(3);  % Identity matrix (no transformations)
            end
            
            if isempty(scDF)
                scDF=1;  % Empty scDF defaults to unity
            end
            
            if ~all(size(M)==[3 3])
                    error('Raith_element.plot:  augmented transformation matrix must be 3 x 3.');
            end
            
            if scDF<0
                error('Raith_element.plot:  scDF must be non-negative.')
            end
            
            hold on
            
            obj.renderplot(M,scDF,1);
                       
        end % plot
        

        function plotedges(obj,varargin)
        %
        % Raith_element.plotedges([M,scDF]) 
        %
        % Plot element with Raith dose factor colouring (edges of polygons where applicable)
        %
        % Argument:
        %
        %   M - augmented transformation matrix for plot [optional]
        %   scDF - overall multiplicative scaling factor for DF specified
        %       in obj.data.DF (e.g., passed from a positionlist entry)
        %       [optional]
        %
            
            if nargin==1 
                M=[];  
                scDF=1;  % Scaling for DF specified in obj.data.DF
            elseif nargin==2
                M=varargin{1};
                scDF=1;  
            elseif nargin==3
                M=varargin{1};
                scDF=varargin{2};  
            else
                error('Raith_element.plotedges:  Too many input arguments.');
            end
            
            if isempty(M)
                M=eye(3);  % Identity matrix (no transformations)
            end
            
            if isempty(scDF)
                scDF=1;  % Empty scDF defaults to unity
            end
            
            if ~all(size(M)==[3 3])
                    error('Raith_element.plotedges:  augmented transformation matrix must be 3 x 3.');
            end
            
            if scDF<0
                error('Raith_element.plotedges:  scDF must be non-negative.')
            end
            
            obj.renderplot(M,scDF,0);
                       
        end % plotedges
                 
        
    end % methods
    
    
    
    methods(Hidden)
        
        function UV=renderplot(obj,M,scDF,plflag)
        %
        % UV=Raith_element.renderplot(obj,M,scDF,plflag) 
        %
        % Plot element with Raith dose factor colouring (called from
        % Raith_element.plot or Raith_element.plotedges)
        % 
        % If plflag==1, plot as filled polygons where applicable
        % (Raith_element.plot). If plflag==0, plot as polygon outlines where 
        %  applicable (Raith_element.plotedges).  If plflag==2, do not
        %  plot, but only return UV.
        %
        % Arguments:
        %
        %   M - augmented transformation matrix for plot
        %   scDF - overall multiplicative scaling factor for DF specified
        %       in obj.data.DF (e.g., passed from a positionlist entry)
        %   plflag - flag for type of plot (1 for .plot, 0 for .plotedges)
        %
        %
        % Return value:
        %
        %   UV - 2 x n matrix [u;v] of points being plotted (um); for paths 
        %       with nonzero width, returns points of underlying path; 
        %       required by Raith_library.writegds for plain GDSII export
            
            hold on
            
           
            switch obj.type 
                
                case 'polygon'
                    UV=M*[obj.data.uv;ones(1,size(obj.data.uv,2))];
                    if plflag==1
                        fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');
                    elseif plflag==0
                        plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                    else % plflag==2
                        % Do not plot
                    end
                        
                        
                case 'path'
                    mag=sqrt(abs(det(M)));  % Total magnification
                    UV=M*[obj.data.uv;ones(1,size(obj.data.uv,2))];
                    w=obj.data.w; 
                    if w>0
                        w=abs(w)*mag;  % Positive w, so scale with everything else
                    else  % Negative w:  default to zero (Raith software interpretation)
                        w=0;
                    end
                        
                    if plflag~=2
                        if w==0 % Single-pixel line
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));    
                        else  % Finite-thickness line
                            [outx,outy]=obj.plotpathwidth(UV(1,:),UV(2,:),w);
                            if plflag==1
                                fill(outx,outy,obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                            elseif plflag==0
                                plot(outx,outy,'color',obj.RaithDF(obj.data.DF*scDF));
                            end
                        end                
                    end
                    
                        
                case 'dot'
                    if isscalar(obj.data.DF)
                        DF=obj.data.DF*scDF*ones(1,size(obj.data.uv,2));
                    else
                        DF=obj.data.DF*scDF;
                    end
                    UV=M*[obj.data.uv;ones(1,size(obj.data.uv,2))];
                    if plflag~=2
                        for k=1:size(obj.data.uv,2)
                           plot(UV(1,k),UV(2,k),'.','color',obj.RaithDF(DF(k)));
                        end
                    end
                    
                    
                case 'arc'
                    
                    mag=sqrt(abs(det(M)));  % Total magnification
                    w=obj.data.w; 
                    if w>=0
                        w=w*mag;  % Positive w, so scale with everything else
                    else  % Negative w:  default to filled element (Raith behaviour)
                        w=[];
                    end
                    
                    if length(obj.data.r)==1 % Circular arc
                        r=obj.data.r*[1 1];
                    else % Elliptical arc
                        r=obj.data.r;
                    end

                    phi=obj.data.angle/180*pi;
                    vertices=obj.data.N;
                    t=linspace(obj.data.theta(1),obj.data.theta(2),vertices)/180*pi;    
                    
                    u_wr=obj.data.uv_c(1)+r(1)*cos(t)*cos(phi)-r(2)*sin(t)*sin(phi);
                    v_wr=obj.data.uv_c(2)+r(1)*cos(t)*sin(phi)+r(2)*sin(t)*cos(phi);
                        
                    if isempty(w)  % Filled arc segment
                        u_wr=[obj.data.uv_c(1) u_wr obj.data.uv_c(1)];
                        v_wr=[obj.data.uv_c(2) v_wr obj.data.uv_c(2)];
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag==1
                            fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');
                        elseif plflag==0
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    elseif w==0 % Single pixel line
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag~=2
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    else % Line with some width
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        [outx,outy]=obj.plotpathwidth(UV(1,:),UV(2,:),w); 
                        if plflag==1
                            fill(outx,outy,obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                        elseif plflag==0
                            plot(outx,outy,'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    end
                     
                        
                case 'circle' 
                    
                    mag=sqrt(abs(det(M)));  % Total magnification
                    w=obj.data.w; 
                    if w>=0
                        w=w*mag;  % Positive w, so scale with everything else
                    else  % Negative w:  default to filled element
                        w=[];
                    end
                    
                    r=obj.data.r;
                    
                    vertices=obj.data.N;
                    th=linspace(0,2*pi,vertices+1);                        
                        
                    if ~isempty(w) % Ring
                        if w==0  % Single-pixel line
                            u_wr=r*cos(th)+obj.data.uv_c(1);
                            v_wr=r*sin(th)+obj.data.uv_c(2);
                            UV=M*[u_wr;v_wr;ones(size(u_wr))];
                            if plflag~=2
                                plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                            end
                        else % Finite-width circle
                            u1=(r+w/2)*cos(th)+obj.data.uv_c(1);
                            v1=(r+w/2)*sin(th)+obj.data.uv_c(2);
                            u2=(r-w/2)*cos(th)+obj.data.uv_c(1);
                            v2=(r-w/2)*sin(th)+obj.data.uv_c(2);
                            u_wr=[u1 u2(end:-1:1) u1(1)];
                            v_wr=[v1 v2(end:-1:1) v1(1)];
                            UV=M*[u_wr;v_wr;ones(size(u_wr))];
                            if plflag==1
                                fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                            elseif plflag==0
                                plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                            end
                        end
                    else % Disk
                        u_wr=r*cos(th)+obj.data.uv_c(1);
                        v_wr=r*sin(th)+obj.data.uv_c(2);
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag==1
                            fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                        elseif plflag==0
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    end
                        
                        
                case 'ellipse' 
                    
                    mag=sqrt(abs(det(M)));  % Total magnification
                    w=obj.data.w; 
                    if w>=0
                        w=w*mag;  % Positive w, so scale with everything else
                    else  % Negative w:  default to filled element
                        w=[];
                    end
                    
                    r=obj.data.r;
                    
                    phi=obj.data.angle/180*pi;
                    vertices=obj.data.N;
                    t=linspace(0,2*pi,vertices+1);    
                        
                    if isempty(w)  % Filled ellipse
                        u_wr=obj.data.uv_c(1)+r(1)*cos(t)*cos(phi)-r(2)*sin(t)*sin(phi);
                        v_wr=obj.data.uv_c(2)+r(1)*cos(t)*sin(phi)+r(2)*sin(t)*cos(phi);
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag==1
                            fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                        elseif plflag==0
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    elseif w==0 % Single pixel line
                        u_wr=obj.data.uv_c(1)+r(1)*cos(t)*cos(phi)-r(2)*sin(t)*sin(phi);
                        v_wr=obj.data.uv_c(2)+r(1)*cos(t)*sin(phi)+r(2)*sin(t)*cos(phi);
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag~=2
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    else % Line with some width:  approximate plot as a simply-connected polygon
                        r1=r-w/2;
                        r2=r+w/2;
                        u_in=obj.data.uv_c(1)+r1(1)*cos(t)*cos(phi)-r1(2)*sin(t)*sin(phi);
                        v_in=obj.data.uv_c(2)+r1(1)*cos(t)*sin(phi)+r1(2)*sin(t)*cos(phi);
                        u_out=obj.data.uv_c(1)+r2(1)*cos(t)*cos(phi)-r2(2)*sin(t)*sin(phi);
                        v_out=obj.data.uv_c(2)+r2(1)*cos(t)*sin(phi)+r2(2)*sin(t)*cos(phi);
                        u_wr=[u_out u_in(end:-1:1) u_out(1)];
                        v_wr=[v_out v_in(end:-1:1) v_out(1)];
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag==1
                            fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                        elseif plflag==0
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    end
                
                    
                case 'text'  % Render text as a structure of polygons, then plot
                    UV=[obj.data.uv_0(1);obj.data.uv_0(2);1];
                    E=Raith_element;
                    T=E.rendertext(obj.data.layer,obj.data.uv_0,obj.data.h,obj.data.angle,obj.data.uv_align,obj.data.textlabel,obj.data.DF*scDF);
                    if plflag==1
                        T.plot(M);
                    elseif plflag==0
                        T.plotedges(M);
                    end
                    
                    
                case 'fbmspath'
                    mag=sqrt(abs(det(M)));  % Total magnification
                    UV=M*[obj.data.uv;ones(1,size(obj.data.uv,2))];
                    w=abs(obj.data.w)*mag; % Defaults to non-negative (always scales) 
                    
                    if isscalar(obj.data.cvtr) % cvtr=0, so make vector
                        cvtr=zeros(1,size(obj.data.uv,2)); 
                    else
                        cvtr=obj.data.cvtr*mag; 
                    end
                    
                    % Construct arcs/line segments
                    uwr=UV(1,1);
                    vwr=UV(2,1);
                    for k=2:length(cvtr)
                        if cvtr(k)==0 % Line segment
                            uwr=[uwr UV(1,k)];
                            vwr=[vwr UV(2,k)];
                        else % Arc
                            % Construct three points for circle fitting
                            P1=UV(1:2,k-1);
                            P2=UV(1:2,k);
                            Pmean=(P1+P2)/2; % Centre of chord                           
                            theta=atan2(P2(2)-P1(2),P2(1)-P1(1))+sign(cvtr(k))*pi/2;
                            P3=Pmean+abs(cvtr(k))*[cos(theta);sin(theta)];
                            C=[norm(P1)^2 P1(1) P1(2) 1;norm(P2)^2 P2(1) P2(2) 1;norm(P3)^2 P3(1) P3(2) 1];
                            UV_c=[det(C(:,[1 3 4]))/det(C(:,[2 3 4]))/2;-det(C(:,[1 2 4]))/det(C(:,[2 3 4]))/2]; % Centre of circle
                            R=sqrt(norm(UV_c)^2+det(C(:,[1 2 3]))/det(C(:,[2 3 4])));  % Radius of circle
                            th1=atan2(P1(2)-UV_c(2),P1(1)-UV_c(1)); % Angle from circle centre to P1
                            th2=atan2(P2(2)-UV_c(2),P2(1)-UV_c(1)); % Angle from circle centre to P2
                            th3=atan2(P3(2)-UV_c(2),P3(1)-UV_c(1)); % Angle from circle centre to P3
                            if ~issorted(-sign(cvtr(k))*[th1  th3 th2]) % Adjust th1 if arc crosses negative y-axis
                                th1=th1+sign(cvtr(k))*2*pi;
                            end
                            th=linspace(th1,th2,33);  % Chosen to match Raith software plotting
                            th(1)=[];
                            uwr=[uwr UV_c(1)+R*cos(th)];  
                            vwr=[vwr UV_c(2)+R*sin(th)];
                        end
                    end
                        
                    if plflag~=2
                        if w==0 % Single-pixel line
                            plot(uwr,vwr,'color',obj.RaithDF(obj.data.DF*scDF));    
                        else  % Finite-thickness line
                            [outx,outy]=obj.plotpathwidth(uwr,vwr,w);
                            if plflag==1
                                fill(outx,outy,obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                            elseif plflag==0
                                plot(outx,outy,'color',obj.RaithDF(obj.data.DF*scDF));
                            end
                        end                
                    end
                    
                    
                case 'fbmscircle' 
                    
                    mag=sqrt(abs(det(M)));  % Total magnification
                    w=abs(obj.data.w)*mag; % Defaults to positive (always scales) 
                    
                    r=obj.data.r;
                    
                    vertices=64;  % Value used in Raith software when converting FBMS circle to path
                    th=linspace(0,2*pi,vertices);                        
                        
                    if w==0  % Single-pixel line
                        u_wr=r*cos(th)+obj.data.uv_c(1);
                        v_wr=r*sin(th)+obj.data.uv_c(2);
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag~=2
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    else % Finite-width circle
                        u1=(r+w/2)*cos(th)+obj.data.uv_c(1);
                        v1=(r+w/2)*sin(th)+obj.data.uv_c(2);
                        u2=(r-w/2)*cos(th)+obj.data.uv_c(1);
                        v2=(r-w/2)*sin(th)+obj.data.uv_c(2);
                        u_wr=[u1 u2(end:-1:1) u1(1)];
                        v_wr=[v1 v2(end:-1:1) v1(1)];
                        UV=M*[u_wr;v_wr;ones(size(u_wr))];
                        if plflag==1
                            fill(UV(1,:),UV(2,:),obj.RaithDF(obj.data.DF*scDF),'EdgeColor','none');    
                        elseif plflag==0
                            plot(UV(1,:),UV(2,:),'color',obj.RaithDF(obj.data.DF*scDF));
                        end
                    end
                    
                    
                case 'sref'  % Cannot plot structure, since sref element does not contain structure data; instead, mark origin as +, with structure name
                    UV=M*[obj.data.uv_0(1);obj.data.uv_0(2);1];
                    if plflag~=2
                        plot(UV(1),UV(2),'.','MarkerEdgeColor','none');  % Plot point for proper axis display
                        text(UV(1),UV(2),'+','color','r','HorizontalAlignment','center');
                        text(UV(1),UV(2),sprintf(['\n\n[' obj.data.name ']']),'HorizontalAlignment','center','Interpreter','none');
                    end
                     
                    
                case 'aref'  % Cannot plot structure, since aref element does not contain structure data; instead, mark origins as +, with structure name
                    
                    % Construct lattice of origins for structures
                    [U,V]=meshgrid((0:obj.data.n_colrow(1)-1)*obj.data.a_colrow(1),(0:obj.data.n_colrow(2)-1)*obj.data.a_colrow(2));
                    U=reshape(U,1,prod(obj.data.n_colrow));
                    V=reshape(V,1,prod(obj.data.n_colrow));
                    
                    % Construct 4 points enclosing lattice for good axis limits
                    Ua=[-1 -1 obj.data.n_colrow(1) obj.data.n_colrow(1)]*obj.data.a_colrow(1);
                    Va=[-1 obj.data.n_colrow(2) -1 obj.data.n_colrow(2)]*obj.data.a_colrow(2);
                    
                    if ~isempty(obj.data.angle)  % Rotate lattice if specified
                        th=obj.data.angle;
                        UVrot=[cosd(th) -sind(th);sind(th) cosd(th)]*[U;V];
                        UaVarot=[cosd(th) -sind(th);sind(th) cosd(th)]*[Ua;Va];
                        U=UVrot(1,:);
                        V=UVrot(2,:);
                        Ua=UaVarot(1,:);
                        Va=UaVarot(2,:);
                    end
                    
                    UaVa=M*[Ua;Va;ones(size(Ua))];
                    UV=M*[U;V;ones(size(U))];
                    
                    if plflag~=2

                        plot(UaVa(1,:)+obj.data.uv_0(1),UaVa(2,:)+obj.data.uv_0(2),'.','MarkerEdgeColor','none');  % Plot points for proper axis display

                        for k=1:length(U)
                            text(UV(1,k)+obj.data.uv_0(1),UV(2,k)+obj.data.uv_0(2),'+','color','r','HorizontalAlignment','center');
                            text(UV(1,k)+obj.data.uv_0(1),UV(2,k)+obj.data.uv_0(2),sprintf(['\n\n[' obj.data.name ']']),'HorizontalAlignment','center','Interpreter','none');
                        end
                        
                    end

            end
            
            % Return points being plotted (for paths with nonzero width, returns points of underlying path); required by Raith_library.writegds for plain GDSII export
            UV(3,:)=[];
                       
        end % renderplot
        
        
        
        function S=rendertext(obj,layer,uv_0,h,theta,uv_align,textlabel,DF)
        %
        % S=rendertext(obj,layer,uv_0,h,theta,uv_align,textlabel,DF)
        % 
        % Return a Raith_structure object for stencil text (simply-connected 
        % polygons, to prevent regions from being released during an undercut); 
        % based on Geogrotesque Stencil.
        %
        % Arguments:
        %
        % u_0 - abscissa (u) for text anchor point (um)
        % v_0 - ordinate (v) for text anchor point (um)
        % h - height of letters (um)
        % theta - angle of rotation of text w.r.t. positive x-axis (degrees)
        % u_align - alignment w.r.t. u_c:  0 (left), 1 (centre), 2 (right) 
        % v_align - alignment w.r.t. v_c:  0 (top), 1 (centre), 2 (bottom)
        % textlabel - the text to be written (string)
        % layer - layer of GDSII file
        % DF - e-beam dose factor
        
        % Construct cell array of text of proper size, but with first character
        % centred at x=half of first character width.
        cursorx=0;
        C={};
        for k=1:length(textlabel)
            in=strfind(obj.chars,textlabel(k));
            if in==length(obj.chars)  % Space
                cursorx=cursorx+obj.textw(in)/2*h;
            else
                chr=obj.font{in};
                cursorx=cursorx+obj.textw(in)/2*h;
                for kk=1:length(chr)
                    ch=chr{kk};
                    x=ch(1,:)*h+cursorx;
                    y=ch(2,:)*h;
                    C{end+1}=[x;y];
                end
                
            end

            if k<length(textlabel)
                in2=strfind(obj.chars,textlabel(k+1));
                cursorx=cursorx+(obj.textw(in)/2+obj.kern(in,in2))*h;
            end

        end

        minx=0;
        maxx=cursorx+obj.textw(in)/2*h;

        switch uv_align(1)
            case 0
                xshift=minx;
            case 1
                xshift=mean([minx maxx]);
            case 2
                xshift=maxx;
        end

        switch uv_align(2)
            case 0
                yshift=h;
            case 1
                yshift=h/2;
            case 2
                yshift=0;
        end
        
        S=Raith_structure;
        S.name='renderedtext';
        S.elements(length(C))=Raith_element;

        for k=1:length(C)
            C{k}(1,:)=C{k}(1,:)-xshift;  % Correct centering
            C{k}(2,:)=C{k}(2,:)-yshift;
            C{k}=[cosd(theta) -sind(theta);sind(theta) cosd(theta)]*C{k};  % Correct rotation
            C{k}(1,:)=C{k}(1,:)+uv_0(1);  % Correct origin
            C{k}(2,:)=C{k}(2,:)+uv_0(2);
            S.elements(k)=Raith_element('polygon',layer,C{k},DF);
        end
            
        end % rendertext
        
    end % Hidden methods
    
    
    methods(Static,Hidden)
        
        function [outx,outy]=plotpathwidth(X,Y,w)
        %
        % [outx,outy]=plotpathwidth(x,y,w)
        %
        % Return polygon vertices for an outlined path (no mitering)
        %

        x=X;
        y=Y;
        th=0;
        in=0;  % Counter for number of rotations required to eliminate infinite slopes

        % Check for vertical segments (infinite slope)
        mch=(y(2:end)-y(1:(end-1)))./(x(2:end)-x(1:(end-1)));
        while any(abs(mch)>1e10)  % Effectively infinite
            in=in+1;
            th=pi/2^in;  % Rotate everything
            xy=[cos(th) -sin(th);sin(th) cos(th)]*[X;Y];
            x=xy(1,:);
            y=xy(2,:);
            mch=(y(2:end)-y(1:(end-1)))./(x(2:end)-x(1:(end-1)));
        end
        
        % Remove repeated vertices
        x(isnan(mch))=[];
        y(isnan(mch))=[];
        
        % Remove colinear vertices
        alph=atan2(diff(y),diff(x));  % Angles between adjacent line segments
        ins=[true abs(diff(alph))>1e-14 true];
        x=x(ins);
        y=y(ins);

        outx=[];
        outy=[];

        
        for k=1:2 % Do procedure for path forwards and backwards to obtain outline

            alph=repmat(atan2(diff(y),diff(x)),2,1);  % Angles between adjacent line segments (in columns)

            % Vertices of line segment pairs (in columns)
            qx=reshape([x(1) reshape(repmat(x(2:(end-1)),2,1),1,2*length(x(2:(end-1)))) x(end)],2,length(x)-1);  
            qy=reshape([y(1) reshape(repmat(y(2:(end-1)),2,1),1,2*length(y(2:(end-1)))) y(end)],2,length(y)-1);

            % Vertices of line segment pairs shifted by w/2, 90° clockwise from (qx,qy)
            px=qx+w/2*cos(alph-pi/2);
            py=qy+w/2*sin(alph-pi/2);

            m=(py(2:2:end)-py(1:2:end))./(px(2:2:end)-px(1:2:end)); % Slopes of line segments

            if size(px,2)==1  % Single line segment
            
                % Vertices of line segments for clockwise-shifted outline
                outx=[outx px'];
                outy=[outy py'];
                
            else
                
                % Vertices of clockwise corners
                cx=(py(3:2:end)-py(1:2:(end-2))+m(1:(end-1)).*px(1:2:(end-2))-m(2:end).*px(3:2:end))./(m(1:(end-1))-m(2:end));  
                cy=py(1:2:(end-2))+m(1:(end-1)).*(cx-px(1:2:(end-2)));

                % Vertices of line segments for clockwise-shifted outline
                outx=[outx px(1) cx px(end)];
                outy=[outy py(1) cy py(end)];
                
            end

            % Reverse path, then repeat
            x=x(end:-1:1);
            y=y(end:-1:1);  

        end

        % Close outline
        outx(end+1)=outx(1);
        outy(end+1)=outy(1);

        % Rotate everything back if necessary
        if th~=0
            outxy=[cos(-th) -sin(-th);sin(-th) cos(-th)]*[outx;outy];
            outx=outxy(1,:);
            outy=outxy(2,:);
        end

        end % plotpathwidth
        
        
        function rgb=RaithDF(DF)
        %
        % rgb=RaithDF(DF)
        %
        % Returns the RGB triple for the Raith dose factor colourmap.
        %
        % Argument:
        %
        % DF - dose factor for Raith EBL
        %
        % Return value:
        %
        % rgb - RGB triple for colour correponding to this dose factor.
        %
        % Aaron Hryciw, 2011.10.04
        %
        %

        x=255/780*[-4 -3 -2 -1 0]+1.5; % Turning points

        if DF<=x(1)
            rgb=[0 0 1];
        elseif DF<=x(2)
            rgb=[0 780*(DF-x(1)) 255]/255;
        elseif DF<=x(3)
            rgb=[0 255 255-780*(DF-x(2))]/255;
        elseif DF<=x(4)
            rgb=[780*(DF-x(3)) 255 0]/255;
        elseif DF<=x(5)
            rgb=[255 255-780*(DF-x(4)) 0]/255;
        else 
            rgb=[1 0 0];
        end

        end % RaithDF
       
        
    end % static, hidden methods
    
    
    properties (Constant,Hidden)
        
        % Allowed characters for text rendering
        chars='`1234567890-=qwertyuiop[]\asdfghjkl;''zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:"ZXCVBNM<>?µ ';
        
        % Stencil font polygons (based on Geogrotesque Stencil)
        font={{[0.1051 0.0352 -0.1051 -0.0206 0.1051;0.8214 0.8214 1.0255 1.0255 0.8214]},{[0.0336 0.0336 -0.3062 -0.3062 -0.0406 -0.0406 -0.3062 -0.3062 0.3062 0.3062 0.0336;0.0565 1.0027 0.9150 0.8450 0.9132 0.0565 0.0565 -0.0110 -0.0110 0.0565 0.0565]},{[-0.1990 -0.1990 0.2923 0.2923 -0.2800 -0.2800 -0.2646 -0.2227 -0.1654 -0.1026 -0.0445 0.0061 0.0579 0.1156 0.1678 0.2073 0.2073 0.2014 0.1694 0.1124 0.0619 -0.0034 -0.0682 -0.1188 -0.1666 -0.2041 -0.2156 -0.2156 -0.2923 -0.2923 -0.2887 -0.2670 -0.2279 -0.1812 -0.1160 -0.0607 0.0061 0.0729 0.1271 0.1903 0.2342 0.2694 0.2883 0.2883 0.2792 0.2480 0.2038 0.1575 0.0994 0.0279 -0.0255 -0.0860 -0.1421 -0.1860 -0.1990;0.2807 0.0565 0.0565 -0.0110 -0.0110 0.2982 0.3567 0.3997 0.4303 0.4521 0.4667 0.4794 0.4929 0.5114 0.5379 0.5755 0.7992 0.8554 0.8995 0.9262 0.9359 0.9394 0.9357 0.9257 0.9045 0.8615 0.8065 0.7007 0.7007 0.7992 0.8547 0.9134 0.9560 0.9822 1.0010 1.0085 1.0105 1.0080 0.9994 0.9777 0.9484 0.9026 0.8447 0.6043 0.5503 0.5044 0.4723 0.4503 0.4299 0.4106 0.3988 0.3807 0.3555 0.3195 0.2807]},{[0.2814 0.2814 0.2759 0.2557 0.2186 0.1692 0.2206 0.2621 0.2842 0.2905 0.2905 0.2818 0.2530 0.2063 0.1553 0.0874 0.0328 -0.0324 -0.0874 -0.1553 -0.2063 -0.2530 -0.2818 -0.2905 -0.2905 -0.2123 -0.2123 -0.1992 -0.1609 -0.1032 -0.0383 0.0237 0.0775 0.1328 0.1790 0.2126 0.2126 0.2012 0.1704 0.1265 0.0656 -0.1482 -0.1482 0.0510 0.1020 0.1514 0.1885 0.2032 0.2032 0.1968 0.1680 0.1221 0.0696 0.0194 -0.0372 -0.0976 -0.1530 -0.1913 -0.2051 -0.2051 -0.2834 -0.2834 -0.2759 -0.2498 -0.2055 -0.1553 -0.0877 -0.0328 0.0332 0.0877 0.1553 0.2047 0.2490 0.2743 0.2814;0.8147 0.6947 0.6400 0.5871 0.5492 0.5292 0.4997 0.4623 0.4161 0.3617 0.1627 0.1070 0.0545 0.0165 -0.0063 -0.0217 -0.0274 -0.0274 -0.0219 -0.0069 0.0155 0.0533 0.1068 0.1581 0.2267 0.2267 0.1669 0.1146 0.0753 0.0519 0.0416 0.0406 0.0469 0.0645 0.0955 0.1420 0.3617 0.4122 0.4530 0.4781 0.4906 0.4906 0.5586 0.5586 0.5664 0.5907 0.6330 0.6827 0.7904 0.8467 0.8916 0.9198 0.9342 0.9390 0.9380 0.9281 0.9055 0.8671 0.8154 0.7555 0.7555 0.8154 0.8723 0.9256 0.9641 0.9872 1.0027 1.0084 1.0084 1.0032 0.9886 0.9664 0.9284 0.8722 0.8147]},{[0.2927 0.2188 0.2188 0.2927 0.2927;0.9946 0.9946 -0.0107 -0.0107 0.9946],[-0.0638 0.0128 -0.2081 0.1650 0.1650 -0.2927 -0.2927 -0.0638;0.9946 0.9946 0.3885 0.3885 0.3143 0.3143 0.3885 0.9946]},{[-0.2573 -0.2573 0.2498 0.2498 -0.1830 -0.1830 0.0375 0.0957 0.1518 0.2000 0.2379 0.2688 0.2767 0.2767 0.2696 0.2431 0.2028 0.1403 0.0779 0.0000 -0.0545 -0.1221 -0.1747 -0.2249 -0.2621 -0.2759 -0.2767 -0.2024 -0.2024 -0.1814 -0.1364 -0.0794 -0.0206 0.0439 0.0980 0.1498 0.1893 0.2024 0.2024 0.1822 0.1344 0.0858 0.0328 -0.2573;0.5266 0.9946 0.9946 0.9266 0.9266 0.5927 0.5927 0.5893 0.5754 0.5509 0.5172 0.4648 0.4154 0.1576 0.1043 0.0522 0.0159 -0.0114 -0.0243 -0.0289 -0.0273 -0.0171 0.0002 0.0320 0.0797 0.1302 0.2303 0.2303 0.1494 0.0968 0.0641 0.0457 0.0391 0.0414 0.0510 0.0728 0.1104 0.1576 0.4154 0.4645 0.4992 0.5175 0.5266 0.5266]},{[0.2875 0.2875 0.2812 0.2543 0.2128 0.1674 0.1073 0.0306 -0.0306 -0.0824 -0.1480 -0.1982 -0.2456 -0.2729 -0.2868 -0.2875 -0.2824 -0.2638 -0.2290 -0.1773 -0.1243 -0.0587 -0.0061 0.0464 0.1121 0.1650 0.2168 0.2441 0.2670 0.2753 0.2753 0.2010 0.2010 0.1745 0.1243 0.0682 0.0168 -0.0374 -0.0895 -0.1460 -0.1891 -0.2132 -0.2132 -0.1966 -0.1528 -0.0947 -0.0350 0.0239 0.0789 0.1358 0.1836 0.2132 0.2132 0.1966 0.1536 0.0974 0.0393 -0.0128 -0.0698 -0.1255 -0.1619 -0.1619 -0.1136 -0.0575 0.0065 0.0670 0.1306 0.1899 0.2409 0.2757 0.2875;0.4156 0.1638 0.1090 0.0552 0.0176 -0.0050 -0.0210 -0.0289 -0.0289 -0.0235 -0.0086 0.0138 0.0514 0.0939 0.1457 0.8132 0.8644 0.9110 0.9506 0.9820 0.9984 1.0077 1.0091 1.0076 0.9981 0.9819 0.9519 0.9252 0.8861 0.8457 0.7638 0.7638 0.8359 0.8847 0.9154 0.9328 0.9390 0.9368 0.9272 0.9047 0.8727 0.8227 0.1638 0.1106 0.0740 0.0515 0.0400 0.0382 0.0463 0.0652 0.0960 0.1410 0.3946 0.4487 0.4856 0.5064 0.5152 0.5161 0.5114 0.4972 0.4728 0.5478 0.5678 0.5794 0.5833 0.5812 0.5721 0.5519 0.5187 0.4675 0.4156]},{[-0.3164 -0.3164 0.2283 -0.1702 -0.0836 0.3164 0.3164 -0.3164;0.9946 0.9242 0.9242 -0.0134 -0.0134 0.9242 0.9946 0.9946]},{[0.1541 0.2087 0.2470 0.2704 0.2810 0.2810 0.2755 0.2545 0.2174 0.1652 0.1047 0.0597 0.0253 0.0253 0.0854 0.1399 0.1814 0.2012 0.2016 0.1917 0.1573 0.1095 0.0581 -0.0253 -0.0755 -0.1265 -0.1707 -0.2020 -0.2020 -0.1901 -0.1573 -0.1095 -0.0573 -0.0253 -0.0253 -0.0905 -0.1431 -0.1941 -0.2364 -0.2652 -0.2794 -0.2794 -0.2739 -0.2541 -0.2170 -0.1541 -0.2190 -0.2605 -0.2822 -0.2885 -0.2885 -0.2814 -0.2589 -0.2237 -0.1680 -0.1059 -0.0684 -0.0253 -0.0253 -0.0802 -0.1281 -0.1767 -0.2111 -0.2111 -0.2000 -0.1692 -0.1257 -0.0652 0.0249 0.0787 0.1356 0.1814 0.2103 0.2103 0.1870 0.1466 0.0933 0.0514 0.0253 0.0253 0.0870 0.1383 0.1925 0.2395 0.2723 0.2885 0.2885 0.2818 0.2601 0.2190 0.1541;0.5288 0.5449 0.5784 0.6252 0.6750 0.8164 0.8686 0.9173 0.9582 0.9860 1.0018 1.0074 1.0091 0.9395 0.9320 0.9114 0.8749 0.8216 0.6957 0.6427 0.5973 0.5701 0.5596 0.5596 0.5616 0.5775 0.6106 0.6624 0.8051 0.8606 0.9007 0.9249 0.9369 0.9395 1.0091 1.0040 0.9932 0.9732 0.9416 0.8982 0.8485 0.6957 0.6411 0.5882 0.5502 0.5288 0.5007 0.4634 0.4173 0.3628 0.1639 0.1128 0.0654 0.0286 -0.0014 -0.0181 -0.0234 -0.0259 0.0424 0.0492 0.0640 0.0957 0.1415 0.3628 0.4133 0.4542 0.4793 0.4917 0.4917 0.4906 0.4754 0.4414 0.3913 0.1566 0.1074 0.0738 0.0520 0.0442 0.0424 -0.0259 -0.0208 -0.0105 0.0095 0.0428 0.0890 0.1380 0.3628 0.4173 0.4634 0.5007 0.5288]},{[-0.2875 -0.2875 -0.2808 -0.2539 -0.2128 -0.1670 -0.1069 -0.0306 0.0306 0.0824 0.1480 0.1986 0.2456 0.2729 0.2872 0.2875 0.2824 0.2638 0.2294 0.1773 0.1247 0.0587 0.0061 -0.0460 -0.1121 -0.1650 -0.2168 -0.2516 -0.2721 -0.2753 -0.2010 -0.2010 -0.1741 -0.1243 -0.0682 -0.0168 0.0377 0.0895 0.1460 0.1891 0.2136 0.2136 0.1966 0.1528 0.0951 0.0354 -0.0235 -0.0785 -0.1354 -0.1836 -0.2132 -0.2132 -0.1962 -0.1536 -0.0970 -0.0389 0.0128 0.0698 0.1188 0.1619 0.1623 0.1140 0.0575 -0.0065 -0.0666 -0.1302 -0.1899 -0.2405 -0.2753 -0.2875;0.5632 0.8152 0.8699 0.9237 0.9613 0.9839 0.9999 1.0079 1.0079 1.0024 0.9875 0.9651 0.9275 0.8851 0.8333 0.1657 0.1146 0.0679 0.0283 -0.0031 -0.0195 -0.0288 -0.0302 -0.0287 -0.0192 -0.0030 0.0270 0.0631 0.1131 0.2150 0.2150 0.1430 0.0942 0.0635 0.0461 0.0399 0.0421 0.0517 0.0741 0.1061 0.1562 0.8152 0.8684 0.9049 0.9274 0.9389 0.9407 0.9326 0.9137 0.8830 0.8379 0.5842 0.5302 0.4933 0.4725 0.4637 0.4628 0.4675 0.4790 0.5054 0.4311 0.4111 0.3995 0.3956 0.3977 0.4068 0.4269 0.4602 0.5113 0.5632]},{[-0.0233 -0.0783 -0.1316 -0.1767 -0.2146 -0.2146 -0.2012 -0.1628 -0.1146 -0.0652 -0.0233 -0.0233 -0.0897 -0.1455 -0.2051 -0.2474 -0.2743 -0.2889 -0.2889 -0.2806 -0.2514 -0.2103 -0.1506 -0.0929 -0.0233 -0.0233;0.0393 0.0432 0.0567 0.0816 0.1485 0.8142 0.8634 0.9012 0.9226 0.9317 0.9339 1.0039 0.9997 0.9877 0.9620 0.9266 0.8844 0.8262 0.1500 0.1077 0.0557 0.0191 -0.0085 -0.0216 -0.0263 0.0393],[0.0233 0.0933 0.1506 0.1968 0.2427 0.2763 0.2889 0.2889 0.2794 0.2478 0.2051 0.1455 0.0897 0.0233 0.0233 0.0787 0.1356 0.1790 0.2146 0.2146 0.2024 0.1605 0.1111 0.0561 0.0233 0.0233;-0.0263 -0.0216 -0.0085 0.0115 0.0459 0.0971 0.1500 0.8142 0.8732 0.9266 0.9620 0.9877 0.9997 1.0039 0.9339 0.9303 0.9151 0.8902 0.8233 0.1600 0.1103 0.0700 0.0500 0.0407 0.0393 -0.0263]},{[0.2180 -0.2180 -0.2180 0.2180 0.2180;0.3769 0.3769 0.4418 0.4418 0.3769]},{[0.3122 -0.3122 -0.3122 0.3122 0.3122;0.3134 0.3134 0.3782 0.3782 0.3134],[0.3122 -0.3122 -0.3122 0.3122 0.3122;0.5669 0.5669 0.6319 0.6319 0.5669]},{[0.1917 0.2659 0.2659 0.1917 0.1917;0.7192 0.7192 -0.3196 -0.3196 0.7192],[-0.2659 -0.2659 -0.2569 -0.2298 -0.1860 -0.1280 -0.0715 -0.0136 0.0377 0.0900 0.1419 0.1419 0.0893 0.0363 -0.0200 -0.0723 -0.1250 -0.1719 -0.1918 -0.1918 -0.1784 -0.1391 -0.0871 -0.0321 0.0195 0.0706 0.1169 0.1419 0.1419 0.0969 0.0482 -0.0020 -0.0620 -0.1220 -0.1723 -0.2209 -0.2522 -0.2659;0.5665 0.1351 0.0854 0.0401 0.0038 -0.0199 -0.0259 -0.0219 -0.0115 0.0055 0.0298 0.1011 0.0747 0.0575 0.0472 0.0467 0.0606 0.0900 0.1275 0.5510 0.6004 0.6390 0.6585 0.6623 0.6562 0.6409 0.6179 0.6021 0.6732 0.6993 0.7175 0.7287 0.7329 0.7283 0.7135 0.6808 0.6316 0.5665]},{[0.3768 0.2099 0.0479 -0.0329 -0.1996 -0.3706 -0.4493 -0.2488 -0.1660 0.0039 0.1697 0.2526 0.4493 0.3768;0.7181 0.1016 0.7181 0.7181 0.0882 0.7181 0.7181 -0.0112 -0.0112 0.5917 -0.0112 -0.0112 0.7181 0.7181]},{[0.2605 0.2605 0.2548 0.2284 0.1797 0.1213 0.0660 0.0000 -0.0660 -0.1213 -0.1797 -0.2210 -0.2514 -0.2605 -0.2605 -0.2527 -0.2248 -0.1851 -0.1265 -0.0696 0.0000 0.0525 0.1167 0.1665 0.2132 0.2473 0.2605 0.2605 0.1821 0.1821 0.1604 0.1142 0.0644 0.0107 -0.0434 -0.0964 -0.1423 -0.1745 -0.1862 -0.1862 -0.1736 -0.1369 -0.0896 -0.0279 0.0293 0.0795 0.1305 0.1713 0.1882 0.1882 -0.1406 -0.1406 0.2605;0.3457 0.5516 0.6041 0.6529 0.6937 0.7167 0.7275 0.7313 0.7280 0.7182 0.6959 0.6632 0.6091 0.5562 0.1589 0.0986 0.0459 0.0118 -0.0125 -0.0236 -0.0275 -0.0259 -0.0159 0.0009 0.0313 0.0762 0.1312 0.2116 0.2116 0.1354 0.0835 0.0559 0.0432 0.0382 0.0406 0.0522 0.0743 0.1092 0.1474 0.5415 0.5906 0.6284 0.6498 0.6603 0.6608 0.6563 0.6406 0.5980 0.5360 0.4065 0.4065 0.3457 0.3457]},{[0.1716 0.1176 0.0674 0.0096 -0.0423 -0.0973 -0.0973 -0.1716 -0.1716 -0.0973 -0.0973 -0.0484 -0.0038 0.0433 0.0913 0.1267 0.1716 0.1716;0.7267 0.7240 0.7156 0.6976 0.6718 0.6317 0.7192 0.7192 -0.0110 -0.0110 0.5515 0.5892 0.6143 0.6347 0.6495 0.6565 0.6597 0.7267]},{[-0.0210 -0.0210 0.1759 0.1759 -0.0210 -0.0210 -0.0952 -0.0952 -0.2066 -0.2066 -0.0952 -0.0952 -0.0879 -0.0602 -0.0177 0.0307 0.0880 0.2066 0.2066 0.1426 0.0914 0.0401 -0.0043 -0.0210;0.1605 0.6508 0.6508 0.7168 0.7168 0.8949 0.8949 0.7168 0.7168 0.6508 0.6508 0.1605 0.1087 0.0568 0.0226 0.0024 -0.0093 -0.0112 0.0538 0.0561 0.0628 0.0773 0.1087 0.1605]},{[-0.2900 -0.2093 0.0027 0.2134 0.2900 -0.0670 -0.1350 -0.0361 -0.2900;0.7233 0.7233 0.0963 0.7233 0.7233 -0.3162 -0.3162 -0.0155 0.7233]}, ...
              {[-0.0762 -0.0071 0.0446 0.0953 0.1458 0.1867 0.1867 0.2610 0.2610 0.1867 0.1867 0.1365 0.0849 0.0350 -0.0190 -0.0771 -0.1263 -0.1699 -0.1868 -0.1868 -0.2610 -0.2610 -0.2552 -0.2308 -0.1899 -0.1327 -0.0762;-0.0259 -0.0232 -0.0142 0.0027 0.0282 0.0563 -0.0165 -0.0165 0.7181 0.7181 0.1278 0.0950 0.0708 0.0547 0.0449 0.0459 0.0597 0.0952 0.1457 0.7181 0.7181 0.1566 0.0990 0.0435 0.0042 -0.0187 -0.0259]},{[0.0371 -0.0371 -0.0371 0.0371 0.0371;0.7181 0.7181 -0.0110 -0.0110 0.7181],[0.0371 -0.0371 -0.0371 0.0371 0.0371;1.0253 1.0253 0.8480 0.8480 1.0253]},{[-0.0234 -0.0758 -0.1320 -0.1734 -0.1925 -0.1925 -0.1826 -0.1451 -0.0887 -0.0517 -0.0234 -0.0234 -0.0887 -0.1420 -0.1967 -0.2336 -0.2595 -0.2667 -0.2667 -0.2607 -0.2374 -0.2020 -0.1472 -0.0922 -0.0234 -0.0234;0.0398 0.0437 0.0611 0.0953 0.1491 0.5432 0.5923 0.6353 0.6570 0.6620 0.6628 0.7329 0.7287 0.7167 0.6910 0.6556 0.6022 0.5432 0.1606 0.1082 0.0562 0.0196 -0.0080 -0.0211 -0.0259 0.0398],[0.0234 0.0922 0.1472 0.2020 0.2374 0.2607 0.2667 0.2667 0.2595 0.2336 0.1967 0.1420 0.0887 0.0234 0.0234 0.0770 0.1285 0.1699 0.1916 0.1924 0.1826 0.1477 0.0957 0.0655 0.0234 0.0234;-0.0259 -0.0211 -0.0080 0.0196 0.0562 0.1082 0.1606 0.5432 0.6022 0.6556 0.6910 0.7167 0.7287 0.7329 0.6628 0.6592 0.6440 0.6128 0.5609 0.1606 0.1109 0.0706 0.0479 0.0422 0.0398 -0.0259]},{[-0.1918 -0.2660 -0.2660 -0.1918 -0.1918;0.7192 0.7192 -0.3196 -0.3196 0.7192],[0.2660 0.2660 0.2570 0.2298 0.1860 0.1281 0.0715 0.0136 -0.0377 -0.0900 -0.1419 -0.1419 -0.0893 -0.0363 0.0201 0.0722 0.1249 0.1718 0.1917 0.1917 0.1783 0.1390 0.0871 0.0320 -0.0195 -0.0706 -0.1419 -0.1419 -0.0969 -0.0483 0.0020 0.0621 0.1220 0.1724 0.2209 0.2522 0.2648 0.2660;0.5510 0.1351 0.0854 0.0401 0.0038 -0.0199 -0.0259 -0.0219 -0.0115 0.0055 0.0298 0.1011 0.0747 0.0575 0.0472 0.0467 0.0606 0.0900 0.1351 0.5510 0.6004 0.6390 0.6585 0.6623 0.6562 0.6409 0.6021 0.6732 0.6993 0.7175 0.7287 0.7329 0.7283 0.7135 0.6808 0.6316 0.5809 0.5510]},{[-0.0709 0.1484 0.1484 -0.1484 -0.1484 0.1484 0.1484 -0.0709 -0.0709;0.9483 0.9483 1.0132 1.0132 -0.3352 -0.3352 -0.2671 -0.2671 0.9483]},{[0.0709 -0.1484 -0.1484 0.1484 0.1484 -0.1484 -0.1484 0.0709 0.0709;-0.2702 -0.2702 -0.3352 -0.3352 1.0132 1.0132 0.9451 0.9451 -0.2702]},{[-0.1306 -0.2093 0.1306 0.2093 -0.1306;1.0124 1.0124 -0.2435 -0.2435 1.0124]},{[0.2634 0.2634 0.1892 0.1892 0.1762 0.1323 0.0791 0.0217 -0.0422 -0.0928 -0.1442 -0.1746 -0.1768 -0.2510 -0.2510 -0.2420 -0.2119 -0.1711 -0.1134 -0.0589 0.0062 0.0713 0.1257 0.1835 0.2242 0.2544 0.2634;0.5546 -0.0081 -0.0081 0.5546 0.6089 0.6488 0.6662 0.6722 0.6697 0.6597 0.6335 0.5871 0.5107 0.5107 0.5690 0.6190 0.6669 0.6998 0.7244 0.7361 0.7402 0.7369 0.7272 0.7052 0.6729 0.6201 0.5546],[-0.0345 0.1407 0.1407 -0.0778 -0.1309 -0.1828 -0.2264 -0.2558 -0.2634 -0.2634 -0.2494 -0.2121 -0.1606 -0.1041 -0.0513 -0.0013 0.0552 0.1040 0.1407 0.1407 0.0962 0.0399 -0.0200 -0.0730 -0.1249 -0.1656 -0.1912 -0.1912 -0.1779 -0.1376 -0.0887 -0.0345;0.3402 0.3402 0.3980 0.3980 0.3922 0.3732 0.3387 0.2879 0.2331 0.1211 0.0653 0.0267 -0.0023 -0.0203 -0.0247 -0.0186 -0.0035 0.0153 0.0330 0.1011 0.0758 0.0550 0.0442 0.0445 0.0604 0.0931 0.1396 0.2331 0.2819 0.3188 0.3358 0.3402]},{[-0.0206 -0.0749 -0.1287 -0.1774 -0.2178 -0.2450 -0.2543 -0.2522 -0.2318 -0.1965 -0.1488 -0.0934 -0.0345 0.0246 0.0841 0.1405 0.1893 0.2329 0.2488 0.2488 0.1802 0.1802 0.1491 0.0992 0.0386 -0.0264 -0.0763 -0.1239 -0.1632 -0.1769 -0.1763 -0.1526 -0.0995 -0.0399 0.0202 0.0755 0.1245 0.1795 0.2296 0.2489 0.2612 0.2612 0.2498 0.2147 0.1703 0.1110 0.0578 -0.0032 -0.0608 -0.1119 -0.1696 -0.2138 -0.2494 -0.2612 -0.2612 -0.1881 -0.1874 -0.1703 -0.1222 -0.0665 -0.0067 0.0542 0.1034 0.1490 0.1731 0.1890 0.1890 0.1744 0.1238 0.0669 -0.0206;0.3311 0.3477 0.3636 0.3856 0.4185 0.4692 0.5226 0.5828 0.6414 0.6832 0.7095 0.7237 0.7298 0.7296 0.7226 0.7086 0.6857 0.6434 0.5840 0.5109 0.5109 0.5840 0.6306 0.6561 0.6690 0.6715 0.6660 0.6504 0.6146 0.5648 0.5107 0.4647 0.4365 0.4156 0.3972 0.3803 0.3641 0.3418 0.3092 0.2854 0.2466 0.1051 0.0689 0.0293 0.0037 -0.0146 -0.0229 -0.0259 -0.0230 -0.0147 0.0038 0.0306 0.0737 0.1260 0.1932 0.1932 0.1349 0.0844 0.0522 0.0375 0.0337 0.0371 0.0474 0.0680 0.0897 0.1285 0.2189 0.2521 0.2835 0.3026 0.3311]},{[0.2675 0.1933 0.1933 0.2675 0.2675;1.0264 1.0264 -0.0120 -0.0120 1.0264],[-0.2675 -0.2675 -0.2629 -0.2441 -0.2123 -0.1592 -0.1091 -0.0551 -0.0022 0.0557 0.1090 0.1426 0.1426 0.0942 0.0420 -0.0134 -0.0681 -0.1248 -0.1658 -0.1907 -0.1933 -0.1805 -0.1417 -0.0885 -0.0304 0.0249 0.0789 0.1071 0.1426 0.1426 0.0890 0.0373 -0.0117 -0.0647 -0.1173 -0.1763 -0.2243 -0.2510 -0.2675;0.1502 0.5369 0.5905 0.6466 0.6860 0.7163 0.7281 0.7313 0.7272 0.7128 0.6914 0.6729 0.6065 0.6295 0.6488 0.6591 0.6598 0.6442 0.6132 0.5622 0.1689 0.1135 0.0709 0.0494 0.0452 0.0519 0.0673 0.0789 0.0963 0.0297 0.0047 -0.0127 -0.0235 -0.0275 -0.0236 -0.0049 0.0308 0.0749 0.1502]},{[0.1310 0.2001 0.2001 0.1359 0.0782 0.0266 -0.0277 -0.0632 -0.0810 -0.0876 -0.0876 -0.2001 -0.2001 -0.0876 -0.0876 -0.0133 -0.0133 0.1779 0.1779 -0.0133 -0.0133 -0.0032 0.0310 0.0830 0.1310;0.9638 0.9638 1.0299 1.0299 1.0275 1.0137 0.9829 0.9404 0.8969 0.8621 0.7132 0.7132 0.6513 0.6513 -0.0120 -0.0120 0.6513 0.6513 0.7132 0.7132 0.8416 0.8925 0.9358 0.9591 0.9638]},{[0.2666 0.2666 0.1923 0.1923 0.1757 0.1342 0.0803 0.0264 -0.2133 -0.2133 0.0577 0.1118 0.1676 0.2176 0.2534 0.2666;-0.1618 0.7168 0.7168 -0.1466 -0.1970 -0.2310 -0.2502 -0.2563 -0.2563 -0.3182 -0.3182 -0.3142 -0.2961 -0.2634 -0.2140 -0.1618],[-0.2666 -0.2666 -0.2620 -0.2432 -0.2114 -0.1583 -0.1082 -0.0542 -0.0013 0.0566 0.1099 0.1435 0.1435 0.0951 0.0430 -0.0125 -0.0672 -0.1240 -0.1649 -0.1898 -0.1924 -0.1796 -0.1408 -0.0876 -0.0295 0.0258 0.0798 0.1080 0.1435 0.1435 0.0899 0.0382 -0.0108 -0.0638 -0.1164 -0.1754 -0.2234 -0.2501 -0.2666;0.1502 0.5369 0.5905 0.6466 0.6860 0.7163 0.7281 0.7313 0.7272 0.7128 0.6914 0.6729 0.6065 0.6295 0.6488 0.6591 0.6598 0.6442 0.6132 0.5622 0.1689 0.1135 0.0709 0.0494 0.0452 0.0519 0.0673 0.0789 0.0963 0.0297 0.0047 -0.0127 -0.0235 -0.0275 -0.0236 -0.0049 0.0308 0.0749 0.1502]},{[0.0485 -0.0120 -0.0671 -0.1167 -0.1918 -0.1918 -0.2660 -0.2660 -0.1918 -0.1918 -0.1437 -0.0974 -0.0429 0.0190 0.0794 0.1340 0.1754 0.1917 0.1917 0.2660 0.2660 0.2612 0.2391 0.1984 0.1496 0.0867 0.0485;0.7291 0.7253 0.7142 0.6971 0.6551 1.0253 1.0253 -0.0110 -0.0110 0.5708 0.6034 0.6266 0.6455 0.6562 0.6548 0.6380 0.6011 0.5501 -0.0110 -0.0110 0.5419 0.5979 0.6541 0.6956 0.7168 0.7281 0.7291]},{[0.1442 0.1442 0.0700 0.0700 0.0524 0.0091 -0.0436 -0.0977 -0.1442 -0.1442 -0.0717 -0.0172 0.0366 0.0836 0.1215 0.1411 0.1442;-0.1498 0.7181 0.7181 -0.1498 -0.1986 -0.2306 -0.2475 -0.2549 -0.2565 -0.3215 -0.3215 -0.3164 -0.3016 -0.2764 -0.2361 -0.1859 -0.1498],[0.1442 0.0700 0.0700 0.1442 0.1442;1.0253 1.0253 0.8480 0.8480 1.0253]},{[-0.0760 0.2250 0.2250 0.1570 -0.1399 -0.1896 -0.1896 -0.2639 -0.2639 -0.1896 -0.1896 -0.1399 0.1858 0.2639 0.2639 -0.0760;0.4006 0.7017 0.7181 0.7181 0.4274 0.4274 1.0253 1.0253 -0.0110 -0.0110 0.3624 0.3624 -0.0110 -0.0110 0.0032 0.4006]},{[0.0371 -0.0371 -0.0371 0.0371 0.0371;1.0253 1.0253 -0.0110 -0.0110 1.0253]},{[0.0838 0.0095 0.0095 0.0838 0.0838;0.7074 0.7074 0.5574 0.5574 0.7074],[0.0830 0.0087 0.0087 0.0040 -0.0099 -0.0526 -0.0838 -0.0522 -0.0063 0.0375 0.0700 0.0830 0.0830;0.1251 0.1251 -0.0053 -0.0262 -0.0477 -0.0798 -0.0982 -0.1565 -0.1326 -0.1011 -0.0585 -0.0122 0.1251]},{[0.0372 -0.0372 -0.0281 0.0281 0.0372;0.9887 0.9887 0.6934 0.6934 0.9887]},{[0.1572 -0.2356 -0.2356 0.2453 0.2453 -0.1712 0.2631 0.2631 -0.2631 -0.2631 0.1572;0.6568 0.6568 0.7248 0.7248 0.6568 0.0570 0.0570 -0.0110 -0.0110 0.0570 0.6568]}, ...
              {[0.2975 0.0404 0.2931 0.2144 0.0013 -0.2100 -0.2931 -0.0405 -0.2975 -0.2189 -0.0012 0.2144 0.2975;-0.0110 0.3609 0.7233 0.7233 0.4176 0.7233 0.7233 0.3577 -0.0110 -0.0110 0.3010 -0.0110 -0.0110]},{[-0.2605 -0.2605 -0.2548 -0.2358 -0.2019 -0.1530 -0.1046 -0.0459 0.0239 0.0862 0.1380 0.1916 0.2358 0.2548 0.2605 0.2605 0.1883 0.1883 0.1797 0.1537 0.1117 0.0557 0.0007 -0.0543 -0.1101 -0.1581 -0.1863 -0.1863 -0.1746 -0.1341 -0.0863 -0.0325 0.0215 0.0749 0.1230 0.1656 0.1821 0.1821 0.2605 0.2605 0.2544 0.2302 0.1927 0.1340 0.0747 -0.0000 -0.0696 -0.1265 -0.1851 -0.2248 -0.2526 -0.2605;0.1589 0.5415 0.5967 0.6433 0.6805 0.7085 0.7225 0.7302 0.7313 0.7250 0.7121 0.6864 0.6438 0.6041 0.5727 0.5017 0.5017 0.5202 0.5765 0.6231 0.6484 0.6593 0.6611 0.6575 0.6423 0.6111 0.5674 0.1589 0.1092 0.0689 0.0489 0.0396 0.0386 0.0450 0.0598 0.0898 0.1264 0.2116 0.2116 0.1478 0.0950 0.0485 0.0151 -0.0107 -0.0230 -0.0275 -0.0236 -0.0125 0.0118 0.0459 0.0986 0.1589]},{[-0.2867 -0.2081 0.0029 0.2142 0.2867 0.0455 -0.0373 -0.2867;0.7181 0.7181 0.0964 0.7181 0.7181 -0.0112 -0.0112 0.7181]},{[-0.1938 -0.2680 -0.2680 -0.1938 -0.1938;1.0299 1.0299 -0.0081 -0.0081 1.0299],[0.2680 0.2680 0.2623 0.2378 0.1964 0.1504 0.0946 0.0302 -0.0257 -0.0814 -0.1131 -0.1433 -0.1433 -0.0929 -0.0446 0.0084 0.0593 0.1089 0.1601 0.1847 0.1938 0.1938 0.1816 0.1440 0.0927 0.0407 -0.0094 -0.0672 -0.1189 -0.1433 -0.1433 -0.0943 -0.0441 0.0056 0.0624 0.1237 0.1813 0.2286 0.2554 0.2680;0.1478 0.5546 0.6082 0.6615 0.7006 0.7204 0.7311 0.7311 0.7233 0.7076 0.6939 0.6764 0.6123 0.6359 0.6538 0.6651 0.6652 0.6563 0.6306 0.5998 0.5679 0.1663 0.1174 0.0765 0.0548 0.0516 0.0546 0.0687 0.0906 0.1031 0.0371 0.0121 -0.0055 -0.0162 -0.0203 -0.0148 0.0048 0.0425 0.0872 0.1478]},{[0.0494 0.0070 -0.0446 -0.0954 -0.1459 -0.1867 -0.1867 -0.2609 -0.2609 -0.1867 -0.1867 -0.1366 -0.0850 -0.0351 0.0190 0.0770 0.1263 0.1699 0.1867 0.1867 0.2609 0.2609 0.2552 0.2308 0.1899 0.1326 0.0761 0.0494;0.7284 0.7258 0.7168 0.6999 0.6744 0.6463 0.7192 0.7192 -0.0110 -0.0110 0.5748 0.6076 0.6318 0.6478 0.6576 0.6566 0.6429 0.6075 0.5568 -0.0110 -0.0110 0.5460 0.6036 0.6591 0.6984 0.7213 0.7284 0.7284]},{[0.3018 0.2288 0.1733 0.1188 0.0643 0.0201 -0.0075 -0.0482 -0.1010 -0.1539 -0.2083 -0.2601 -0.3110 -0.3615 -0.4125 -0.4125 -0.4867 -0.4867 -0.4125 -0.4125 -0.3623 -0.3108 -0.2608 -0.2068 -0.1488 -0.0995 -0.0559 -0.0390 -0.0390 0.0352 0.0352 0.0341 0.0643 0.1166 0.1688 0.2176 0.2679 0.3198 0.3731 0.4037 0.4125 0.4125 0.4867 0.4867 0.4810 0.4565 0.4156 0.3584 0.3018;0.7283 0.7253 0.7147 0.6948 0.6647 0.6312 0.6745 0.7047 0.7225 0.7283 0.7267 0.7189 0.7041 0.6803 0.6463 0.7192 0.7192 -0.0110 -0.0110 0.5747 0.6075 0.6317 0.6478 0.6576 0.6566 0.6429 0.6075 0.5568 -0.0110 -0.0110 0.5460 0.5710 0.5926 0.6215 0.6417 0.6538 0.6588 0.6534 0.6306 0.5928 0.5568 -0.0110 -0.0110 0.5460 0.6036 0.6590 0.6983 0.7212 0.7283]},{[0.0832 0.0089 0.0089 0.0061 -0.0097 -0.0524 -0.0832 -0.0520 -0.0061 0.0377 0.0706 0.0800 0.0832 0.0832;0.1224 0.1224 -0.0019 -0.0241 -0.0504 -0.0824 -0.1009 -0.1591 -0.1352 -0.1038 -0.0612 -0.0346 -0.0149 0.1224]},{[0.0370 -0.0370 -0.0370 0.0370 0.0370;0.1224 0.1224 -0.0275 -0.0275 0.1224]},{[0.1304 0.2095 -0.1308 -0.2095 0.1304;1.0124 1.0124 -0.2435 -0.2435 1.0124]},{[-0.2314 -0.2298 -0.2188 -0.1939 -0.1709 -0.1425 -0.1196 -0.0879 -0.0397 0.0109 0.0551 0.1041 0.1350 0.1547 0.1678 0.1729 0.1765 0.2314 0.2302 0.2192 0.2081 0.1895 0.1721 0.1460 0.1192 0.0868 0.0354 -0.0093 -0.0528 -0.1030 -0.1342 -0.1540 -0.1670 -0.1729 -0.1757 -0.2314;0.3022 0.3555 0.4081 0.4518 0.4692 0.4766 0.4755 0.4670 0.4445 0.4156 0.3896 0.3671 0.3627 0.3737 0.3990 0.4224 0.4730 0.4730 0.4152 0.3616 0.3378 0.3151 0.3047 0.2987 0.3004 0.3098 0.3335 0.3585 0.3835 0.4051 0.4092 0.4000 0.3772 0.3482 0.3022 0.3022]},{[0.0393 -0.0393 -0.0393 0.0393 0.0393;-0.0165 -0.0165 0.1339 0.1339 -0.0165],[0.0302 -0.0302 -0.0393 -0.0393 0.0393 0.0393 0.0302;0.2164 0.2164 0.4803 0.9936 0.9936 0.4803 0.2164]},{[0.4709 0.4709 0.4646 0.4401 0.4049 0.3520 0.3049 0.2460 0.1785 0.0974 0.0014 -0.0496 -0.1358 -0.2124 -0.2757 -0.3267 -0.3871 -0.4287 -0.4607 -0.4709 -0.4709 -0.4654 -0.4433 -0.4113 -0.3630 -0.2958 -0.2389 -0.1698 -0.0923 0.2939 0.2939 -0.0460 -0.1109 -0.1666 -0.2370 -0.2923 -0.3449 -0.3765 -0.3927 -0.3927 -0.3828 -0.3464 -0.2848 -0.2152 -0.1571 -0.0852 0.0014 0.0868 0.1575 0.2152 0.2650 0.3196 0.3642 0.3923 0.3923 0.3824 0.3551 0.3121 0.2599 0.2239 0.2089 0.2528 0.1828 0.1342 0.1346 0.1567 0.1947 0.2528 0.3148 0.3698 0.4196 0.4508 0.4674 0.4709;0.1992 0.6556 0.7152 0.7720 0.8106 0.8440 0.8627 0.8776 0.8896 0.8969 0.8988 0.8988 0.8934 0.8848 0.8716 0.8545 0.8236 0.7869 0.7314 0.6711 0.0587 0.0040 -0.0509 -0.0895 -0.1238 -0.1518 -0.1660 -0.1755 -0.1815 -0.1815 -0.1166 -0.1166 -0.1153 -0.1102 -0.0971 -0.0782 -0.0455 -0.0065 0.0475 0.6556 0.7125 0.7607 0.7939 0.8133 0.8232 0.8292 0.8309 0.8292 0.8232 0.8134 0.8013 0.7783 0.7418 0.6928 0.2082 0.1573 0.1152 0.0895 0.1023 0.1434 0.1955 0.6495 0.6495 0.2079 0.1564 0.1018 0.0633 0.0395 0.0350 0.0466 0.0779 0.1216 0.1703 0.1992],[0.1318 0.1192 0.0717 0.0180 -0.0366 -0.0899 -0.1354 -0.1599 -0.1958 -0.1832 -0.1401 -0.0919 -0.0389 0.0128 0.0591 0.0852 0.1022 0.0607 0.0144 -0.0476 -0.1085 -0.1634 -0.2156 -0.2496 -0.2662 -0.2658 -0.2302 -0.2184 -0.1844 -0.1377 -0.0891 -0.0374 0.0176 0.0721 0.0955 0.1318;0.6050 0.5484 0.5727 0.5909 0.5938 0.5816 0.5492 0.4968 0.1867 0.1504 0.1178 0.1034 0.1065 0.1211 0.1407 0.1571 0.1016 0.0724 0.0521 0.0378 0.0350 0.0420 0.0640 0.1013 0.1533 0.2066 0.5205 0.5721 0.6184 0.6466 0.6610 0.6659 0.6595 0.6427 0.6312 0.6050]},{[0.1028 0.1202 -0.0198 -0.0368 0.1028;0.5594 0.6213 0.6213 0.5594 0.5594],[-0.1324 -0.0842 -0.0672 -0.0213 -0.0913 -0.1368 -0.2901 -0.2901 -0.1538 -0.2020 -0.3474 -0.3474 -0.2194 -0.2652 -0.1953 -0.1494 -0.1324;0.3848 0.5594 0.6213 0.7865 0.7865 0.6213 0.6213 0.5594 0.5594 0.3848 0.3848 0.3228 0.3228 0.1574 0.1574 0.3228 0.3848],[0.0375 0.0545 -0.0850 -0.1020 0.0375;0.3228 0.3848 0.3848 0.3228 0.3228],[0.3474 0.3474 0.2375 0.2830 0.2134 0.1676 0.1506 0.1024 0.0850 0.0395 0.1095 0.1549 0.2901 0.2901 0.1723 0.2206 0.3474;0.5594 0.6213 0.6213 0.7865 0.7865 0.6213 0.5594 0.3848 0.3228 0.1574 0.1574 0.3228 0.3228 0.3848 0.3848 0.5594 0.5594]},{[0.1816 0.1816 0.2484 0.2484 0.2389 0.2113 0.1690 0.1160 0.0674 0.0674 0.0239 0.0239 0.0761 0.1263 0.1508 0.1761 0.1816;0.6983 0.6407 0.6407 0.7114 0.7562 0.7993 0.8293 0.8487 0.8572 0.9956 0.9956 0.8034 0.7980 0.7841 0.7699 0.7370 0.6983],[-0.1891 -0.1899 -0.2611 -0.2611 -0.2484 -0.2113 -0.1630 -0.1124 -0.0686 -0.0686 -0.0239 -0.0239 -0.0844 -0.1393 -0.1781 -0.1891;0.2593 0.3301 0.3301 0.2570 0.2026 0.1603 0.1343 0.1205 0.1139 -0.0173 -0.0173 0.1687 0.1754 0.1935 0.2253 0.2593],[-0.0018 -0.0634 -0.1219 -0.1682 -0.1828 -0.1828 -0.1650 -0.1180 -0.0623 -0.0623 -0.0239 -0.0239 -0.0686 -0.0686 -0.1184 -0.1725 -0.2156 -0.2437 -0.2532 -0.2516 -0.2346 -0.2041 -0.1623 -0.1113 -0.0551 -0.0069 0.0417 0.0951 0.1441 0.1741 0.1864 0.1864 0.1698 0.1263 0.0769 0.0239 0.0239 0.0674 0.0674 0.1251 0.1729 0.2192 0.2441 0.2611 0.2611 0.2539 0.2219 0.1686 0.1097 0.0575 -0.0018;0.5327 0.5540 0.5774 0.6080 0.6565 0.7071 0.7590 0.7893 0.8015 0.8015 0.8037 0.9956 0.9956 0.8582 0.8510 0.8321 0.8003 0.7514 0.6981 0.6455 0.5884 0.5461 0.5161 0.4941 0.4758 0.4613 0.4464 0.4287 0.4081 0.3865 0.3612 0.2486 0.2158 0.1882 0.1744 0.1690 -0.0173 -0.0173 0.1144 0.1246 0.1409 0.1679 0.1939 0.2380 0.3743 0.4056 0.4450 0.4751 0.4966 0.5136 0.5327]},{[-0.2648 -0.3146 -0.3447 -0.3601 -0.3601 -0.3545 -0.3379 -0.3063 -0.2648 -0.2648 -0.3209 -0.3715 -0.4091 -0.4269 -0.4269 -0.4111 -0.3672 -0.3103 -0.2648 -0.2648;0.6241 0.6321 0.6510 0.6832 0.8852 0.9078 0.9275 0.9410 0.9448 1.0009 0.9958 0.9765 0.9430 0.8954 0.6791 0.6299 0.5907 0.5733 0.5700 0.6241],[-0.2071 -0.1482 -0.0968 -0.0613 -0.0451 -0.0451 -0.0632 -0.1004 -0.1510 -0.2071 -0.2071 -0.1534 -0.1217 -0.1123 -0.1123 -0.1190 -0.1387 -0.1870 -0.2071 -0.2071;0.5700 0.5755 0.5957 0.6299 0.6749 0.8904 0.9430 0.9765 0.9958 1.0009 0.9448 0.9374 0.9148 0.8890 0.6832 0.6619 0.6417 0.6253 0.6241 0.5700],[0.2652 0.1913 -0.2696 -0.1953 0.2652;0.9956 0.9956 -0.0173 -0.0173 0.9956],[0.2071 0.1573 0.1241 0.1119 0.1119 0.1174 0.1340 0.1885 0.2071 0.2071 0.1510 0.1004 0.0632 0.0451 0.0451 0.0609 0.1051 0.1617 0.2071 0.2071;0.0243 0.0322 0.0544 0.0833 0.2852 0.3078 0.3275 0.3443 0.3448 0.4009 0.3958 0.3766 0.3430 0.2954 0.0792 0.0300 -0.0093 -0.0267 -0.0300 0.0243],[0.2648 0.3237 0.3751 0.4150 0.4269 0.4269 0.4091 0.3715 0.3209 0.2648 0.2648 0.2913 0.3186 0.3502 0.3601 0.3601 0.3336 0.2850 0.2648 0.2648;-0.0300 -0.0245 -0.0043 0.0360 0.0750 0.2904 0.3430 0.3766 0.3958 0.4009 0.3448 0.3434 0.3374 0.3148 0.2852 0.0881 0.0418 0.0254 0.0243 -0.0300]},{[0.2911 0.0354 -0.0358 -0.2911 -0.2168 -0.0002 0.2168 0.2911;0.4195 0.9886 0.9886 0.4195 0.4195 0.9027 0.4195 0.4195]},{[0.1852 0.1792 0.1583 0.1207 0.0678 0.0065 -0.0559 -0.0737 -0.0737 -0.0128 0.0425 0.0840 0.1041 0.1045 0.1840 0.1852;0.8069 0.8591 0.9078 0.9488 0.9765 0.9923 0.9991 0.9997 0.9300 0.9225 0.9019 0.8654 0.8121 0.7490 0.7490 0.8069],[-0.1500 -0.2014 -0.2508 -0.2883 -0.3030 -0.3030 -0.3030 -0.2911 -0.2579 -0.2101 -0.1571 -0.1247 -0.1247 -0.1749 -0.2318 -0.2860 -0.3318 -0.3626 -0.3804 -0.3832 -0.3824 -0.3816 -0.3725 -0.3488 -0.3105 -0.2607 -0.3117 -0.3571 -0.3820 -0.3907 -0.3907 -0.3836 -0.3571 -0.3140 -0.2630 -0.2144 -0.1567 -0.0891 -0.0891 -0.1512 -0.2061 -0.2563 -0.2974 -0.3121 -0.3124 -0.3124 -0.3010 -0.2698 -0.2263 -0.1650 0.3907 0.3907 -0.1500;0.5503 0.5581 0.5823 0.6246 0.6744 0.7410 0.7957 0.8512 0.8912 0.9154 0.9275 0.9300 0.9997 0.9963 0.9869 0.9681 0.9378 0.8972 0.8469 0.7919 0.7316 0.6657 0.6158 0.5689 0.5355 0.5197 0.4968 0.4610 0.4163 0.3631 0.1544 0.1035 0.0552 0.0172 -0.0079 -0.0224 -0.0318 -0.0352 0.0329 0.0358 0.0498 0.0732 0.1114 0.1596 0.2109 0.3534 0.4039 0.4447 0.4698 0.4823 0.4823 0.5503 0.5503],[0.1812 0.1808 0.1646 0.1243 0.0717 0.0168 -0.0421 -0.0421 0.0219 0.0777 0.1259 0.1785 0.2243 0.2524 0.2595 0.2595 0.1812 0.1812;0.2184 0.1467 0.1141 0.0769 0.0502 0.0369 0.0329 -0.0352 -0.0301 -0.0199 -0.0054 0.0192 0.0560 0.1035 0.1544 0.4345 0.4345 0.2184]},{[0.2327 0.2064 0.0193 0.0263 -0.0263 -0.0193 -0.2064 -0.2327 -0.0385 -0.2327 -0.2064 -0.0193 -0.0263 0.0263 0.0193 0.2064 0.2327 0.0385 0.2327;0.8822 0.9277 0.8115 1.0317 1.0317 0.8115 0.9277 0.8822 0.7782 0.6743 0.6287 0.7450 0.5248 0.5248 0.7450 0.6287 0.6743 0.7782 0.8822]},{[-0.0656 -0.0494 -0.0047 0.0490 0.1012 0.1439 0.1439 0.0767 0.0182 -0.0312 -0.0838 -0.1198 -0.1320 -0.1439 -0.1439 -0.1364 -0.1123 -0.0719 -0.0162 0.0364 0.0976 0.1439 0.1439 0.0874 0.0372 -0.0142 -0.0545 -0.0656 -0.0656;0.8503 0.8986 0.9348 0.9549 0.9634 0.9654 1.0264 1.0238 1.0158 1.0024 0.9762 0.9396 0.9175 0.8653 -0.1825 -0.2226 -0.2670 -0.3007 -0.3251 -0.3367 -0.3428 -0.3437 -0.2825 -0.2792 -0.2689 -0.2471 -0.2087 -0.1764 0.8503]}, ...
              {[0.0654 0.0492 0.0045 -0.0492 -0.1014 -0.1441 -0.1441 -0.0769 -0.0184 0.0310 0.0836 0.1196 0.1318 0.1441 0.1441 0.1362 0.1121 0.0717 0.0160 -0.0362 -0.0978 -0.1441 -0.1441 -0.0875 -0.0374 0.0140 0.0543 0.0654 0.0654;0.8592 0.8986 0.9348 0.9549 0.9634 0.9654 1.0264 1.0238 1.0158 1.0024 0.9762 0.9396 0.9175 0.8653 -0.1825 -0.2226 -0.2670 -0.3007 -0.3251 -0.3367 -0.3428 -0.3437 -0.2825 -0.2792 -0.2689 -0.2471 -0.2087 -0.1764 0.8592]},{[0.2832 -0.2832 -0.2832 0.2832 0.2832;-0.2216 -0.2216 -0.1566 -0.1566 -0.2216]},{[0.3067 0.0352 0.0352 -0.0352 -0.0352 -0.3067 -0.3067 -0.0352 -0.0352 0.0352 0.0352 0.3067 0.3067;0.5145 0.5145 0.7860 0.7860 0.5145 0.5145 0.4441 0.4441 0.1727 0.1727 0.4441 0.4441 0.5145]},{[0.1107 0.1661 0.2201 0.2666 0.2952 0.3113 0.3113 0.3061 0.2850 0.2525 0.2027 0.1474 0.0929 0.0391 0.0236 0.0236 0.0827 0.1445 0.1881 0.2222 0.2371 0.2371 0.2222 0.1881 0.1445 0.0827 0.0236 0.0236 0.2637 0.3113 0.1107;-0.0163 -0.0032 0.0201 0.0566 0.0990 0.1674 0.8105 0.8610 0.9109 0.9490 0.9809 1.0003 1.0105 1.0146 1.0148 0.9460 0.9401 0.9216 0.8956 0.8583 0.8105 0.1812 0.1333 0.0961 0.0701 0.0516 0.0457 -0.0231 -0.2679 -0.2262 -0.0163],[-0.0236 -0.0827 -0.1445 -0.1881 -0.2223 -0.2371 -0.2371 -0.2223 -0.1881 -0.1445 -0.0827 -0.0236 -0.0236 -0.0761 -0.1289 -0.1846 -0.2369 -0.2795 -0.3030 -0.3113 -0.3113 -0.3060 -0.2850 -0.2525 -0.2027 -0.1474 -0.0930 -0.0391 -0.0236 -0.0236;0.0457 0.0516 0.0701 0.0961 0.1333 0.1812 0.8105 0.8583 0.8956 0.9216 0.9401 0.9460 1.0148 1.0124 1.0047 0.9886 0.9612 0.9195 0.8719 0.8105 0.1812 0.1307 0.0808 0.0426 0.0108 -0.0086 -0.0189 -0.0229 -0.0231 0.0457]},{[0.4804 0.2646 0.0414 -0.0354 -0.2512 -0.4744 -0.5507 -0.2907 -0.2124 0.0009 0.2251 0.3035 0.5507 0.4804;1.0066 0.1143 1.0066 1.0066 0.1143 1.0066 1.0066 -0.0023 -0.0023 0.8678 -0.0023 -0.0023 1.0066 1.0066]},{[-0.2060 -0.2060 0.2307 0.2307 -0.2060 -0.2060 0.2576 0.2576 -0.2802 -0.2802 0.2802 0.2802 -0.2060;0.0609 0.4955 0.4955 0.5564 0.5564 0.9346 0.9346 1.0043 1.0043 -0.0044 -0.0044 0.0609 0.0609]},{[-0.2443 -0.3186 -0.3186 -0.2443 -0.2443;1.0035 1.0035 -0.0054 -0.0054 1.0035],[0.1071 0.1605 0.2180 0.2649 0.2882 0.2969 0.2969 0.2890 0.2677 0.2241 0.1690 0.1111 0.0455 -0.1918 -0.1918 0.0350 0.0919 0.1428 0.1869 0.2111 0.2197 0.2197 0.1998 0.1607 0.1111 0.0429 -0.1918 -0.1918 0.0305 0.2347 0.3186 0.1071;0.4158 0.4325 0.4615 0.5056 0.5493 0.5836 0.8139 0.8493 0.8951 0.9423 0.9745 0.9936 1.0035 1.0035 0.9341 0.9341 0.9257 0.9058 0.8732 0.8377 0.8056 0.5929 0.5485 0.5142 0.4911 0.4770 0.4770 0.4076 0.4076 -0.0054 -0.0054 0.4158]},{[0.3264 -0.3264 -0.3264 -0.0372 -0.0372 0.0371 0.0371 0.3264 0.3264;1.0066 1.0066 0.9378 0.9378 -0.0023 -0.0023 0.9378 0.9378 1.0066]},{[-0.3501 -0.0378 -0.0378 0.0395 0.0395 0.3501 0.2776 -0.0022 -0.2728 -0.3501;1.0000 0.3939 -0.0020 -0.0020 0.4027 1.0000 1.0000 0.4822 1.0000 1.0000]},{[0.2371 0.2371 0.2223 0.1882 0.1445 0.0826 0.0232 -0.0443 -0.0999 -0.1566 -0.2047 -0.2300 -0.2371 -0.2371 -0.3113 -0.3113 -0.3060 -0.2849 -0.2524 -0.2027 -0.1475 -0.0930 -0.0390 0.0117 0.0680 0.1197 0.1754 0.2287 0.2734 0.2991 0.3113 0.3113 0.2371;1.0066 0.1875 0.1396 0.1024 0.0764 0.0580 0.0514 0.0529 0.0616 0.0825 0.1170 0.1545 0.1818 1.0066 1.0066 0.1875 0.1371 0.0872 0.0489 0.0171 -0.0024 -0.0127 -0.0166 -0.0168 -0.0152 -0.0085 0.0060 0.0313 0.0703 0.1157 0.1738 1.0066 1.0066]},{[0.0371 -0.0371 -0.0371 0.0371 0.0371;1.0030 1.0030 -0.0058 -0.0058 1.0030]},{[-0.0633 -0.1151 -0.1684 -0.2106 -0.2296 -0.2367 -0.2367 -0.2169 -0.1788 -0.1306 -0.0633 -0.0231 -0.0231 -0.0759 -0.1287 -0.1843 -0.2365 -0.2790 -0.3025 -0.3108 -0.3108 -0.3055 -0.2845 -0.2521 -0.2024 -0.1472 -0.0928 -0.0390 -0.0231 -0.0231 -0.0633;0.0675 0.0784 0.1007 0.1365 0.1665 0.1937 0.8333 0.8830 0.9197 0.9441 0.9595 0.9630 1.0317 1.0293 1.0215 1.0054 0.9781 0.9366 0.8891 0.8414 0.1857 0.1490 0.0992 0.0611 0.0293 0.0098 -0.0003 -0.0044 -0.0046 0.0641 0.0675],[0.3076 0.2900 0.2594 0.2112 0.1566 0.1016 0.0456 0.0232 0.0232 0.0825 0.1443 0.1878 0.2219 0.2367 0.2367 0.2169 0.1789 0.1306 0.0634 0.0232 0.0232 0.0760 0.1287 0.1844 0.2366 0.2791 0.3025 0.3108 0.3108 0.3076;0.8661 0.9191 0.9593 0.9934 1.0147 1.0262 1.0312 1.0317 0.9630 0.9570 0.9386 0.9126 0.8754 0.8333 0.1937 0.1440 0.1074 0.0830 0.0675 0.0641 -0.0046 -0.0022 0.0055 0.0216 0.0489 0.0904 0.1380 0.1857 0.8414 0.8661]},{[-0.2335 -0.3077 -0.3077 -0.2335 -0.2335;1.0035 1.0035 -0.0054 -0.0054 1.0035],[0.3077 0.3077 0.2999 0.2786 0.2350 0.1799 0.1220 0.0563 -0.1809 -0.1809 0.0459 0.1028 0.1537 0.1978 0.2220 0.2306 0.2306 0.2106 0.1716 0.1220 0.0538 -0.1809 -0.1809 0.0538 0.0998 0.1566 0.2150 0.2659 0.2935 0.3077;0.5657 0.8139 0.8493 0.8951 0.9423 0.9745 0.9936 1.0035 1.0035 0.9341 0.9341 0.9257 0.9058 0.8732 0.8377 0.8056 0.5752 0.5307 0.4964 0.4733 0.4591 0.4591 0.3897 0.3897 0.3950 0.4097 0.4352 0.4749 0.5152 0.5657]},{[0.0308 0.0308 0.0308 0.0308 0.0308 0.0455 0.0846 0.1379 0.2411 0.2411 0.1628 0.1032 0.0522 0.0103 -0.0277 -0.0466 -0.0478 -0.0478 -0.0478 -0.0478 -0.0577 -0.0866 -0.1320 -0.1818 -0.2411 -0.2411 -0.1818 -0.1320 -0.0866 -0.0577 -0.0478 -0.0478 -0.0478 -0.0478 -0.0478 -0.0415 -0.0146 0.0296 0.0767 0.1320 0.2411 0.2411 0.1628 0.1063 0.0601 0.0336 0.0308 0.0308 0.0308 0.0308 0.0308 0.0225 -0.0036 -0.0462 -0.1000 -0.0462 -0.0036 0.0225 0.0308;0.4931 0.5930 0.6781 0.7631 0.8482 0.8983 0.9365 0.9572 0.9594 1.0276 1.0276 1.0219 1.0045 0.9769 0.9310 0.8738 0.7631 0.6781 0.5930 0.5079 0.4585 0.4172 0.3888 0.3782 0.3782 0.3045 0.3045 0.2940 0.2656 0.2242 0.1748 0.0898 0.0047 -0.0803 -0.1654 -0.2150 -0.2679 -0.3092 -0.3313 -0.3435 -0.3447 -0.2767 -0.2767 -0.2652 -0.2337 -0.1879 -0.1228 -0.0378 0.0472 0.1323 0.1896 0.2420 0.2854 0.3184 0.3413 0.3644 0.3974 0.4407 0.4931]},{[-0.0306 -0.0310 -0.0310 -0.0310 -0.0310 -0.0457 -0.0844 -0.1377 -0.2413 -0.2413 -0.1626 -0.1034 -0.0520 -0.0101 0.0279 0.0468 0.0480 0.0480 0.0480 0.0480 0.0480 0.0611 0.0927 0.1409 0.1935 0.2413 0.2413 0.1820 0.1322 0.0864 0.0579 0.0480 0.0480 0.0480 0.0480 0.0480 0.0417 0.0144 -0.0298 -0.0769 -0.1322 -0.2413 -0.2413 -0.1626 -0.1061 -0.0599 -0.0338 -0.0310 -0.0310 -0.0310 -0.0310 -0.0306 -0.0227 0.0034 0.0460 0.1002 0.0460 0.0034 -0.0227 -0.0306;0.4931 0.5930 0.6781 0.7631 0.8482 0.8983 0.9365 0.9572 0.9594 1.0276 1.0276 1.0219 1.0045 0.9769 0.9310 0.8738 0.8057 0.7206 0.6355 0.5505 0.4993 0.4508 0.4112 0.3854 0.3782 0.3782 0.3045 0.3045 0.2940 0.2656 0.2242 0.1748 0.0898 0.0047 -0.0803 -0.1654 -0.2150 -0.2679 -0.3092 -0.3313 -0.3435 -0.3447 -0.2767 -0.2767 -0.2652 -0.2337 -0.1879 -0.1228 -0.0378 0.0472 0.1323 0.1896 0.2420 0.2854 0.3184 0.3413 0.3644 0.3974 0.4407 0.4931]},{[0.0372 -0.0372 -0.0372 0.0372 0.0372;1.0124 1.0124 -0.3332 -0.3332 1.0124]},{[-0.3525 -0.2783 -0.0041 0.0024 0.1587 -0.1113 -0.1319 0.1794 0.2763 0.3525 0.0413 -0.0391 -0.3525;0.0000 0.0000 0.8979 0.8979 0.3646 0.3646 0.3052 0.3052 0.0000 0.0000 1.0000 1.0000 0.0000]},{[-0.0247 -0.0775 -0.1277 -0.1738 -0.2270 -0.2673 -0.2872 -0.2969 -0.2969 -0.2863 -0.2644 -0.2202 -0.1617 -0.1110 -0.0558 0.0023 0.0603 0.1155 0.1661 0.2247 0.2689 0.2908 0.3014 0.3014 0.2271 0.2271 0.2228 0.2040 0.1806 0.1294 0.0741 0.0010 -0.0507 -0.1114 -0.1670 -0.2066 -0.2226 -0.2237 -0.2196 -0.1939 -0.1422 -0.0921 -0.0247 0.0280 0.0952 0.1533 0.2023 0.2532 0.2868 0.3013 0.3062 0.3062 0.3013 0.2806 0.2425 0.1873 0.1355 0.0746 0.0047 -0.0657 -0.1276 -0.1809 -0.2384 -0.2787 -0.3009 -0.3062 -0.3062 -0.2319 -0.2309 -0.2238 -0.1912 -0.1396 -0.0796 -0.0223 0.0426 0.0964 0.1518 0.1993 0.2247 0.2320 0.2330 0.2172 0.1802 0.1334 0.0729 0.0174 -0.0247;0.4975 0.5154 0.5327 0.5521 0.5858 0.6350 0.6859 0.7520 0.8059 0.8722 0.9231 0.9702 0.9979 1.0093 1.0147 1.0161 1.0141 1.0076 0.9952 0.9679 0.9246 0.8801 0.8240 0.7213 0.7213 0.8017 0.8429 0.8837 0.9070 0.9321 0.9445 0.9490 0.9475 0.9376 0.9137 0.8727 0.8152 0.7651 0.7147 0.6685 0.6331 0.6093 0.5815 0.5611 0.5347 0.5103 0.4862 0.4519 0.4115 0.3755 0.3482 0.1579 0.1333 0.0870 0.0467 0.0149 -0.0028 -0.0140 -0.0179 -0.0145 -0.0044 0.0122 0.0437 0.0865 0.1398 0.1698 0.2882 0.2882 0.1804 0.1457 0.0995 0.0702 0.0552 0.0501 0.0513 0.0587 0.0766 0.1076 0.1441 0.1740 0.3316 0.3697 0.4051 0.4326 0.4605 0.4823 0.4975]},{[-0.3103 -0.2299 -0.2299 -0.3103 -0.3103;0.0000 0.0000 1.0000 1.0000 0.0000],[0.3103 0.3103 0.3016 0.2783 0.2314 0.1737 0.1155 0.0649 -0.1804 -0.1804 0.0464 0.1034 0.1595 0.2000 0.2242 0.2319 0.2319 0.2169 0.1780 0.1277 0.0777 -0.1804 -0.1804 0.0464 0.1007 0.1580 0.2168 0.2682 0.2960 0.3103;0.8227 0.1804 0.1408 0.0922 0.0462 0.0188 0.0052 0.0002 0.0000 0.0619 0.0619 0.0682 0.0928 0.1278 0.1672 0.1975 0.8131 0.8569 0.8970 0.9220 0.9354 0.9381 1.0000 1.0000 0.9966 0.9845 0.9611 0.9214 0.8785 0.8227]}, ...
              {[-0.2147 0.2621 0.2621 -0.2147 -0.2147 0.2889 0.2889 -0.2889 -0.2889 -0.2147 -0.2147;0.4958 0.4958 0.5567 0.5567 0.9350 0.9350 1.0048 1.0048 -0.0041 -0.0041 0.4958]},{[0.2286 0.2988 0.2988 0.2917 0.2652 0.2197 0.1684 0.0987 0.0422 -0.0259 -0.0832 -0.1339 -0.1933 -0.2368 -0.2737 -0.2931 -0.2991 -0.2991 -0.2817 -0.2491 -0.2066 -0.1452 -0.0907 -0.0279 0.0484 0.1102 0.1633 0.2217 0.2689 0.2923 0.2991 0.2988 -0.0001 -0.0001 0.2349 0.2349 0.2178 0.1730 0.1265 0.0654 0.0102 -0.0453 -0.1070 -0.1544 -0.1950 -0.2154 -0.2186 -0.2186 -0.2065 -0.1685 -0.1098 -0.0426 0.0215 0.0772 0.1356 0.1864 0.2196 0.2286 0.2286;0.7134 0.7134 0.8280 0.8765 0.9305 0.9703 0.9943 1.0105 1.0164 1.0164 1.0107 1.0016 0.9803 0.9530 0.9116 0.8695 0.8322 0.1608 0.1047 0.0611 0.0292 0.0030 -0.0089 -0.0164 -0.0164 -0.0096 0.0012 0.0259 0.0678 0.1128 0.1698 0.4679 0.4679 0.4040 0.4040 0.1814 0.1284 0.0879 0.0672 0.0547 0.0515 0.0539 0.0640 0.0831 0.1191 0.1693 0.1930 0.8186 0.8616 0.9047 0.9323 0.9451 0.9463 0.9423 0.9305 0.9072 0.8703 0.8354 0.7134]},{[0.2320 0.2320 -0.2319 -0.2319 -0.3061 -0.3061 -0.2319 -0.2319 0.2320 0.2320 0.3061 0.3061 0.2320;1.0048 0.5690 0.5690 1.0048 1.0048 -0.0041 -0.0041 0.4999 0.4999 -0.0041 -0.0041 1.0048 1.0048]},{[0.2708 0.1965 0.1965 0.1829 0.1454 0.0928 0.0360 -0.0166 -0.0773 -0.1333 -0.1746 -0.1924 -0.1924 -0.2708 -0.2708 -0.2658 -0.2431 -0.2058 -0.1630 -0.1053 -0.0303 0.0303 0.0813 0.1447 0.1923 0.2357 0.2626 0.2708 0.2708;1.0030 1.0030 0.1674 0.1223 0.0824 0.0586 0.0481 0.0471 0.0546 0.0755 0.1132 0.1668 0.2656 0.2656 0.1821 0.1300 0.0751 0.0322 0.0051 -0.0149 -0.0251 -0.0251 -0.0191 -0.0033 0.0200 0.0572 0.1053 0.1480 1.0030]},{[-0.1198 -0.1198 0.2773 0.1853 -0.1940 -0.2392 -0.2392 -0.3135 -0.3135 -0.2392 -0.2392 -0.1940 0.2284 0.3135 -0.1198;0.4986 0.5120 1.0030 1.0030 0.5253 0.5253 1.0030 1.0030 -0.0058 -0.0058 0.4661 0.4661 -0.0058 -0.0058 0.4986]},{[-0.1931 -0.1931 -0.2673 -0.2673 0.2673 0.2673 -0.1931;0.0574 1.0030 1.0030 -0.0058 -0.0058 0.0574 0.0574]},{[0.0372 -0.0372 -0.0372 0.0372 0.0372;0.7074 0.7074 0.5574 0.5574 0.7074],[0.0372 -0.0372 -0.0372 0.0372 0.0372;0.1251 0.1251 -0.0247 -0.0247 0.1251]},{[-0.0577 -0.1320 -0.1229 -0.0664 -0.0577;0.9887 0.9887 0.6934 0.6934 0.9887],[0.1320 0.0577 0.0668 0.1233 0.1320;0.9887 0.9887 0.6934 0.6934 0.9887]},{[0.1705 -0.2714 -0.2714 0.2756 0.2756 -0.1937 0.2988 0.2988 -0.2988 -0.2988 0.1705;0.9319 0.9319 1.0000 1.0000 0.9319 0.0606 0.0606 -0.0074 -0.0074 0.0606 0.9319]},{[0.3386 0.2552 0.0000 -0.2551 -0.3386 -0.0417 -0.3386 -0.2551 0.0000 0.2552 0.3386 0.0418 0.3386;1.0000 1.0000 0.5695 1.0000 1.0000 0.4990 -0.0020 -0.0020 0.4286 -0.0020 -0.0020 0.4990 1.0000]},{[0.2285 0.2988 0.2988 0.2917 0.2652 0.2197 0.1684 0.0988 0.0421 -0.0259 -0.0832 -0.1339 -0.1933 -0.2369 -0.2736 -0.2931 -0.2991 -0.2991 -0.2937 -0.2817 -0.2491 -0.2066 -0.1452 -0.0907 -0.0280 0.0483 0.1102 0.1633 0.2217 0.2689 0.2924 0.2991 0.2991 0.2349 0.2349 0.2178 0.1730 0.1265 0.0654 0.0102 -0.0452 -0.1070 -0.1544 -0.1950 -0.2153 -0.2187 -0.2187 -0.2065 -0.1685 -0.1243 -0.0609 -0.0001 0.0592 0.1226 0.1777 0.2156 0.2285 0.2285;0.7134 0.7134 0.8280 0.8765 0.9306 0.9703 0.9943 1.0105 1.0164 1.0164 1.0107 1.0016 0.9804 0.9530 0.9117 0.8695 0.8322 0.1608 0.1331 0.1047 0.0611 0.0293 0.0030 -0.0089 -0.0164 -0.0164 -0.0096 0.0012 0.0259 0.0678 0.1128 0.1698 0.3052 0.3052 0.1740 0.1284 0.0879 0.0672 0.0547 0.0515 0.0539 0.0640 0.0831 0.1191 0.1693 0.1930 0.8186 0.8616 0.9047 0.9280 0.9426 0.9464 0.9439 0.9341 0.9124 0.8769 0.8351 0.7134]},{[-0.3484 -0.2721 0.0044 0.2782 0.3484 0.0433 -0.0350 -0.3484;1.0066 1.0066 0.1143 1.0066 1.0066 -0.0023 -0.0023 1.0066]},{[-0.3082 -0.2278 -0.2278 -0.3082 -0.3082;0.0000 0.0000 1.0000 1.0000 0.0000],[0.3082 0.3082 0.3041 0.2875 0.2529 0.2098 0.1299 0.1917 0.2347 0.2612 0.2711 0.2711 0.2624 0.2389 0.1925 0.1368 0.0833 0.0330 -0.1804 -0.1804 0.0278 0.0816 0.1307 0.1737 0.1969 0.1969 0.1840 0.1504 0.1047 0.0428 -0.1804 -0.1804 0.0353 0.0852 0.1423 0.1846 0.2175 0.2340 0.2340 0.2251 0.1921 0.1438 0.0930 0.0278 -0.1804 -0.1804 0.0278 0.0807 0.1410 0.1889 0.2342 0.2723 0.2985 0.3082;0.1968 0.3294 0.3725 0.4231 0.4676 0.4975 0.5248 0.5521 0.5909 0.6381 0.7030 0.7958 0.8603 0.9095 0.9557 0.9825 0.9953 1.0000 1.0000 0.9361 0.9361 0.9265 0.9058 0.8683 0.8108 0.7030 0.6467 0.5986 0.5680 0.5526 0.5526 0.4887 0.4887 0.4834 0.4663 0.4390 0.3976 0.3374 0.2258 0.1688 0.1153 0.0842 0.0692 0.0639 0.0639 0.0000 0.0000 0.0011 0.0104 0.0259 0.0524 0.0926 0.1495 0.1968]},{[0.2335 0.3077 0.3077 0.2335 0.2335 -0.2335 -0.3077 -0.3077 -0.2335 -0.2335 0.2335;-0.0058 -0.0058 1.0030 1.0030 0.1410 1.0030 1.0030 -0.0058 -0.0058 0.8497 -0.0058]},{[0.0000 0.3133 0.3874 0.3874 0.3133 0.3133 0.0323 -0.0323 -0.3132 -0.3132 -0.3874 -0.3874 -0.3132 0.0000;0.2295 1.0030 1.0030 -0.0058 -0.0058 0.8270 0.1491 0.1491 0.8270 -0.0058 -0.0058 1.0030 1.0030 0.2295]},{[0.2777 -0.2777 -0.2777 0.2777 0.2777 -0.1939 0.2777 0.2777;0.7741 0.5188 0.4474 0.1920 0.2663 0.4831 0.7000 0.7741]},{[-0.2775 0.2775 0.2775 -0.2775 -0.2775 0.1937 -0.2775 -0.2775;0.1920 0.4474 0.5188 0.7741 0.7000 0.4831 0.2663 0.1920]},{[0.0249 -0.0494 -0.0494 0.0249 0.0249;0.1224 0.1224 -0.0275 -0.0275 0.1224],[0.2984 0.2984 0.2877 0.2585 0.2170 0.1688 0.1198 0.0755 0.0332 0.0249 0.0249 -0.0494 -0.0494 -0.0375 -0.0067 0.0372 0.0874 0.1379 0.1818 0.2126 0.2245 0.2245 0.2067 0.1684 0.1126 0.0585 -0.0012 -0.0609 -0.1142 -0.1692 -0.2067 -0.2241 -0.2241 -0.2984 -0.2984 -0.2929 -0.2680 -0.2269 -0.1798 -0.1158 -0.0632 0.0000 0.0632 0.1162 0.1798 0.2273 0.2684 0.2814 0.2933 0.2984;0.8111 0.6309 0.5935 0.5463 0.5084 0.4771 0.4494 0.4225 0.3829 0.3591 0.2091 0.2091 0.3790 0.4116 0.4522 0.4846 0.5125 0.5396 0.5696 0.6062 0.6353 0.8129 0.8559 0.8899 0.9135 0.9246 0.9282 0.9251 0.9149 0.8923 0.8584 0.8136 0.6855 0.6855 0.8120 0.8533 0.9073 0.9469 0.9713 0.9890 0.9959 0.9978 0.9958 0.9884 0.9699 0.9446 0.9039 0.8824 0.8502 0.8111]},{[-0.07620,-0.007100,0.04460,0.09530,0.1458,0.1867,0.1867,0.2610,0.2610,0.1867,0.1867,0.1365,0.08490,0.03500,-0.01900,-0.07710,-0.1263,-0.1699,-0.1868,-0.1868,-0.2610,-0.2610,-0.1868,-0.1868,-0.1327,-0.07620;-0.02590,-0.02320,-0.01420,0.002700,0.02820,0.05630,-0.01650,-0.01650,0.7181,0.7181,0.1278,0.09500,0.07080,0.05470,0.04490,0.04590,0.05970,0.09520,0.1457,0.7181,0.7181,-0.3183,-0.3183,0.004200,-0.01870,-0.02590]}};
              
        % Width of characters for 1 um height capitals
        textw=[0.2103 0.6123 0.5846 0.5810 0.5854 0.5534 0.5751 0.6328 0.5771 0.5751 0.5779 0.4360 0.6245 0.5319 0.8985 0.5209 0.3432 0.4131 0.5800 0.5219 0.0743 0.5334 0.5320 0.2968 0.2968 0.4186 0.5268 0.5224 0.5350 0.4002 0.5332 0.5319 0.2885 0.5277 0.0742 0.1676 0.0743 0.5262 0.5951 0.5209 0.5735 0.5360 0.5219 0.9734 0.1664 0.0739 0.4190 0.4636 0.0787 0.9419 0.6949 0.5221 0.8537 0.5822 0.7814 0.4654 0.2877 0.2881 0.5664 0.6134 0.6226 1.1014 0.5604 0.6371 0.6529 0.7003 0.6227 0.0742 0.6216 0.6154 0.4822 0.4826 0.0743 0.7051 0.6123 0.6206 0.5779 0.5981 0.6123 0.5415 0.6269 0.5347 0.0743 0.2640 0.5977 0.6772 0.5982 0.6969 0.6164 0.6154 0.7749 0.5553 0.5549 0.5968 0.5220 0.5100];
        
        % Kerning lookup table for text
        kern=[0.275 0.157 0.235 0.216 0.196 0.231 0.235 0.216 0.235 0.235 0.235 0.255 0.255 0.235 0.157 0.212 0.250 0.173 0.154 0.231 0.250 0.231 0.250 0.250 0.212 0.192 0.216 0.216 0.235 0.176 0.235 0.255 0.039 0.255 0.255 0.235 0.255 0.196 0.176 0.235 0.157 0.255 0.255 0.255 0.157 0.235 0.196 0.255 0.275 0.275 0.255 0.255 0.255 0.255 0.255 0.269 0.235 0.196 0.137 0.235 0.235 0.157 0.255 0.255 0.137 0.157 0.255 0.275 0.255 0.275 0.216 0.216 0.275 0.176 0.235 0.275 0.255 0.235 0.255 0.157 0.255 0.255 0.235 0.235 0.196 0.196 0.250 0.157 0.255 0.255 0.255 0.235 0.235 0.235 0.231 0.000 ; 0.216 0.173 0.176 0.154 0.118 0.173 0.173 0.154 0.173 0.173 0.173 0.192 0.192 0.173 0.096 0.173 0.192 0.115 0.098 0.154 0.173 0.154 0.173 0.173 0.154 0.115 0.154 0.135 0.154 0.096 0.154 0.173 0.192 0.077 0.192 0.192 0.192 0.135 0.096 0.154 0.096 0.173 0.173 0.173 0.096 0.154 0.118 0.176 0.196 0.196 0.176 0.176 0.176 0.176 0.173 0.192 0.173 0.154 0.077 0.173 0.173 0.096 0.212 0.212 0.096 0.096 0.192 0.212 0.192 0.216 0.135 0.135 0.192 0.096 0.154 0.192 0.192 0.173 0.216 0.098 0.196 0.196 0.157 0.176 0.137 0.118 0.192 0.098 0.216 0.196 0.196 0.176 0.176 0.176 0.154 0.000 ; 0.235 0.173 0.196 0.154 0.135 0.176 0.157 0.137 0.157 0.157 0.157 0.176 0.192 0.192 0.115 0.192 0.212 0.135 0.115 0.192 0.212 0.196 0.196 0.196 0.176 0.137 0.176 0.157 0.176 0.118 0.192 0.192 0.192 0.096 0.192 0.192 0.192 0.154 0.135 0.176 0.098 0.196 0.196 0.196 0.098 0.176 0.137 0.192 0.212 0.212 0.192 0.192 0.192 0.192 0.192 0.231 0.192 0.154 0.096 0.192 0.192 0.115 0.212 0.212 0.077 0.096 0.173 0.192 0.173 0.192 0.135 0.135 0.235 0.118 0.196 0.235 0.235 0.216 0.235 0.137 0.235 0.216 0.176 0.196 0.157 0.137 0.196 0.118 0.216 0.216 0.192 0.192 0.192 0.192 0.192 0.000 ; 0.255 0.212 0.192 0.173 0.154 0.192 0.212 0.173 0.192 0.192 0.192 0.212 0.212 0.192 0.115 0.212 0.231 0.154 0.135 0.212 0.231 0.212 0.231 0.212 0.173 0.154 0.173 0.173 0.192 0.135 0.192 0.212 0.020 0.060 0.149 0.134 0.149 0.119 0.090 0.104 0.075 0.212 0.212 0.212 0.135 0.212 0.154 0.212 0.231 0.250 0.231 0.231 0.231 0.231 0.231 0.250 0.231 0.196 0.137 0.235 0.235 0.157 0.255 0.255 0.137 0.115 0.212 0.231 0.212 0.231 0.173 0.192 0.231 0.135 0.216 0.216 0.216 0.216 0.216 0.137 0.216 0.216 0.196 0.212 0.173 0.154 0.212 0.135 0.231 0.231 0.231 0.231 0.212 0.212 0.212 0.000 ; 0.250 0.212 0.196 0.192 0.173 0.212 0.192 0.173 0.192 0.212 0.212 0.231 0.231 0.212 0.135 0.212 0.231 0.135 0.115 0.192 0.212 0.192 0.212 0.212 0.192 0.176 0.192 0.192 0.212 0.154 0.212 0.231 0.019 0.231 0.212 0.192 0.212 0.154 0.135 0.192 0.115 0.212 0.212 0.212 0.135 0.212 0.173 0.212 0.231 0.231 0.212 0.212 0.235 0.212 0.212 0.231 0.212 0.173 0.135 0.212 0.212 0.135 0.231 0.231 0.115 0.135 0.212 0.231 0.212 0.255 0.196 0.196 0.255 0.157 0.216 0.255 0.255 0.235 0.231 0.135 0.231 0.231 0.212 0.212 0.173 0.154 0.231 0.135 0.250 0.250 0.250 0.231 0.231 0.231 0.192 0.000 ; 0.255 0.212 0.196 0.192 0.173 0.192 0.192 0.173 0.192 0.192 0.192 0.212 0.212 0.192 0.135 0.192 0.212 0.135 0.135 0.192 0.212 0.216 0.212 0.231 0.192 0.154 0.192 0.173 0.192 0.135 0.212 0.212 0.019 0.212 0.212 0.192 0.212 0.173 0.154 0.212 0.118 0.216 0.216 0.216 0.118 0.196 0.157 0.216 0.212 0.231 0.212 0.212 0.212 0.212 0.212 0.231 0.212 0.192 0.135 0.231 0.231 0.154 0.250 0.250 0.135 0.135 0.212 0.231 0.212 0.231 0.173 0.192 0.231 0.137 0.212 0.250 0.250 0.231 0.250 0.154 0.250 0.250 0.212 0.212 0.173 0.173 0.231 0.135 0.231 0.231 0.231 0.231 0.231 0.231 0.192 0.000 ; 0.235 0.192 0.176 0.173 0.154 0.173 0.192 0.173 0.192 0.192 0.192 0.196 0.196 0.176 0.118 0.176 0.196 0.118 0.118 0.196 0.192 0.173 0.192 0.192 0.173 0.135 0.154 0.154 0.192 0.137 0.196 0.216 0.216 0.098 0.216 0.196 0.216 0.154 0.135 0.192 0.115 0.212 0.212 0.212 0.115 0.196 0.154 0.212 0.231 0.231 0.212 0.212 0.212 0.212 0.212 0.231 0.212 0.173 0.115 0.212 0.212 0.135 0.212 0.212 0.096 0.096 0.192 0.212 0.192 0.212 0.154 0.154 0.212 0.115 0.173 0.212 0.212 0.192 0.235 0.115 0.212 0.212 0.173 0.192 0.154 0.135 0.212 0.115 0.235 0.231 0.231 0.212 0.212 0.212 0.196 0.000 ; 0.216 0.173 0.157 0.154 0.135 0.154 0.173 0.154 0.154 0.154 0.154 0.173 0.173 0.154 0.077 0.154 0.173 0.096 0.096 0.154 0.173 0.154 0.173 0.173 0.135 0.115 0.135 0.154 0.176 0.118 0.176 0.196 0.196 0.078 0.196 0.176 0.173 0.115 0.096 0.154 0.077 0.173 0.173 0.173 0.196 0.019 0.115 0.173 0.192 0.192 0.173 0.173 0.173 0.192 0.192 0.212 0.192 0.154 0.096 0.192 0.192 0.115 0.216 0.216 0.098 0.098 0.196 0.216 0.196 0.216 0.135 0.135 0.173 0.077 0.154 0.173 0.173 0.173 0.216 0.096 0.192 0.192 0.173 0.173 0.135 0.135 0.173 0.115 0.192 0.192 0.192 0.173 0.173 0.173 0.154 0.000 ; 0.235 0.192 0.176 0.173 0.154 0.173 0.192 0.173 0.173 0.192 0.173 0.192 0.192 0.173 0.096 0.173 0.192 0.115 0.115 0.192 0.192 0.173 0.192 0.212 0.173 0.154 0.173 0.173 0.173 0.115 0.173 0.192 0.192 0.096 0.192 0.192 0.212 0.154 0.135 0.192 0.115 0.212 0.212 0.212 0.115 0.173 0.154 0.192 0.212 0.212 0.192 0.192 0.192 0.212 0.196 0.216 0.196 0.157 0.098 0.196 0.196 0.118 0.231 0.231 0.115 0.115 0.212 0.231 0.212 0.231 0.154 0.173 0.212 0.115 0.173 0.212 0.212 0.192 0.235 0.135 0.231 0.231 0.192 0.212 0.173 0.154 0.212 0.135 0.235 0.212 0.212 0.212 0.212 0.212 0.192 0.000 ; 0.235 0.192 0.176 0.192 0.154 0.173 0.192 0.154 0.173 0.173 0.173 0.192 0.192 0.173 0.115 0.173 0.192 0.115 0.115 0.192 0.212 0.192 0.212 0.212 0.173 0.154 0.173 0.173 0.173 0.115 0.173 0.192 0.192 0.077 0.192 0.173 0.192 0.154 0.135 0.173 0.115 0.192 0.192 0.192 0.115 0.192 0.154 0.212 0.231 0.231 0.212 0.212 0.212 0.216 0.216 0.235 0.216 0.176 0.118 0.216 0.216 0.137 0.231 0.231 0.115 0.115 0.212 0.231 0.212 0.231 0.173 0.176 0.231 0.135 0.192 0.231 0.231 0.212 0.231 0.135 0.235 0.231 0.192 0.212 0.173 0.154 0.212 0.115 0.231 0.231 0.235 0.192 0.192 0.173 0.192 0.000 ; 0.235 0.192 0.176 0.173 0.154 0.173 0.176 0.154 0.173 0.173 0.173 0.192 0.192 0.173 0.096 0.176 0.196 0.118 0.118 0.196 0.196 0.176 0.196 0.192 0.176 0.157 0.176 0.176 0.196 0.137 0.196 0.216 -0.020 0.061 0.152 0.136 0.152 0.121 0.091 0.106 0.076 0.212 0.212 0.212 0.115 0.192 0.154 0.212 0.231 0.231 0.192 0.192 0.192 0.192 0.212 0.212 0.212 0.176 0.096 0.192 0.192 0.115 0.192 0.192 0.096 0.096 0.212 0.231 0.212 0.231 0.173 0.173 0.231 0.115 0.196 0.212 0.212 0.192 0.212 0.115 0.212 0.212 0.173 0.212 0.173 0.154 0.192 0.115 0.231 0.231 0.231 0.212 0.196 0.192 0.196 0.000 ; 0.235 0.212 0.192 0.192 0.154 0.192 0.196 0.154 0.192 0.173 0.192 0.212 0.212 0.192 0.077 0.192 0.212 0.135 0.077 0.196 0.216 0.196 0.216 0.216 0.196 0.176 0.176 0.173 0.192 0.135 0.192 0.212 0.019 0.212 0.212 0.192 0.231 0.077 0.058 0.192 0.058 0.212 0.212 0.212 0.115 0.196 0.154 0.212 0.231 0.231 0.212 0.212 0.212 0.212 0.212 0.212 0.192 0.173 0.115 0.192 0.192 0.058 0.212 0.235 0.019 0.038 0.212 0.231 0.192 0.173 0.173 0.173 0.231 0.096 0.192 0.231 0.231 0.212 0.231 0.096 0.235 0.255 0.216 0.235 0.098 0.059 0.235 0.059 0.255 0.235 0.212 0.212 0.212 0.192 0.196 0.000 ; 0.235 0.192 0.192 0.173 0.154 0.192 0.176 0.154 0.173 0.173 0.173 0.192 0.212 0.192 0.115 0.173 0.192 0.115 0.096 0.173 0.192 0.173 0.192 0.192 0.173 0.212 0.212 0.173 0.212 0.135 0.212 0.212 0.212 0.096 0.212 0.192 0.212 0.173 0.135 0.154 0.115 0.212 0.212 0.212 0.115 0.192 0.154 0.212 0.192 0.235 0.216 0.216 0.216 0.216 0.216 0.235 0.192 0.176 0.115 0.192 0.192 0.115 0.212 0.212 0.096 0.118 0.216 0.235 0.216 0.235 0.176 0.176 0.235 0.118 0.196 0.235 0.235 0.216 0.235 0.137 0.235 0.231 0.192 0.212 0.173 0.154 0.192 0.135 0.231 0.231 0.231 0.212 0.212 0.212 0.173 0.000 ; 0.235 0.192 0.192 0.173 0.154 0.192 0.176 0.154 0.173 0.173 0.173 0.192 0.212 0.192 0.115 0.192 0.212 0.135 0.115 0.192 0.216 0.196 0.216 0.216 0.176 0.157 0.176 0.176 0.192 0.137 0.196 0.216 0.020 0.216 0.216 0.196 0.216 0.154 0.118 0.176 0.098 0.196 0.196 0.196 0.098 0.176 0.157 0.212 0.212 0.212 0.212 0.212 0.212 0.212 0.212 0.231 0.212 0.173 0.115 0.212 0.212 0.135 0.212 0.235 0.098 0.118 0.196 0.216 0.196 0.216 0.157 0.157 0.235 0.115 0.173 0.212 0.212 0.192 0.212 0.115 0.212 0.231 0.192 0.212 0.173 0.154 0.212 0.135 0.231 0.231 0.235 0.212 0.212 0.212 0.192 0.000 ; 0.135 0.115 0.115 0.096 0.077 0.115 0.098 0.077 0.096 0.096 0.098 0.096 0.115 0.096 0.038 0.096 0.135 0.058 0.038 0.115 0.135 0.096 0.135 0.135 0.096 0.077 0.077 0.096 0.096 0.058 0.096 0.135 0.135 0.019 0.135 0.118 0.137 0.078 0.059 0.098 0.039 0.137 0.137 0.135 0.078 0.059 0.078 0.137 0.157 0.137 0.137 0.137 0.135 0.135 0.135 0.154 0.135 0.096 0.038 0.135 0.115 0.038 0.135 0.135 0.019 0.019 0.135 0.157 0.137 0.157 0.098 0.098 0.157 0.039 0.118 0.135 0.135 0.115 0.135 0.038 0.135 0.154 0.115 0.135 0.096 0.077 0.115 0.038 0.154 0.154 0.157 0.115 0.115 0.115 0.115 0.000 ; 0.212 0.173 0.173 0.154 0.135 0.173 0.157 0.135 0.154 0.154 0.157 0.173 0.173 0.173 0.096 0.154 0.192 0.096 0.096 0.173 0.192 0.173 0.192 0.192 0.154 0.135 0.154 0.154 0.173 0.115 0.173 0.192 0.000 0.216 0.212 0.192 0.212 0.135 0.115 0.173 0.096 0.192 0.192 0.192 0.096 0.176 0.135 0.192 0.212 0.212 0.192 0.192 0.192 0.192 0.192 0.196 0.176 0.137 0.098 0.176 0.176 0.098 0.196 0.216 0.096 0.096 0.192 0.212 0.192 0.212 0.154 0.154 0.216 0.098 0.157 0.196 0.196 0.176 0.196 0.098 0.196 0.212 0.173 0.192 0.154 0.135 0.173 0.115 0.212 0.212 0.212 0.192 0.192 0.192 0.173 0.000 ; ...
              0.154 0.115 0.115 0.096 0.058 0.115 0.098 0.077 0.096 0.096 0.098 0.058 0.115 0.115 0.058 0.096 0.135 0.059 0.058 0.115 0.135 0.115 0.135 0.137 0.098 0.098 0.098 0.118 0.078 0.118 0.115 0.135 0.231 0.115 0.231 0.212 0.231 0.192 0.154 0.115 0.058 0.135 0.135 0.135 0.212 0.212 0.173 0.137 0.157 0.157 0.137 0.137 0.137 0.137 0.137 0.173 0.135 0.096 0.038 0.135 0.135 0.058 0.154 0.154 0.039 0.038 0.135 0.154 0.135 0.154 0.096 0.096 0.154 0.039 0.115 0.154 0.154 0.135 0.154 0.058 0.154 0.154 0.118 0.137 0.098 0.078 0.115 0.039 0.157 0.157 0.157 0.137 0.137 0.137 0.115 0.000 ; 0.173 0.137 0.135 0.115 0.077 0.135 0.118 0.096 0.115 0.115 0.118 0.154 0.135 0.135 0.058 0.115 0.154 0.078 0.059 0.137 0.157 0.118 0.157 0.157 0.115 0.096 0.115 0.096 0.135 0.077 0.115 0.154 0.154 0.030 0.119 0.104 0.119 0.090 0.060 0.075 0.045 0.090 0.104 0.104 0.058 0.135 0.096 0.154 0.173 0.173 0.154 0.154 0.154 0.154 0.154 0.173 0.154 0.115 0.058 0.154 0.154 0.077 0.173 0.173 0.058 0.058 0.154 0.173 0.154 0.173 0.115 0.115 0.173 0.058 0.137 0.173 0.173 0.154 0.173 0.077 0.173 0.173 0.135 0.154 0.115 0.098 0.157 0.059 0.176 0.176 0.176 0.157 0.157 0.157 0.137 0.000 ; 0.154 0.096 0.115 0.096 0.058 0.115 0.118 0.077 0.115 0.115 0.098 0.077 0.115 0.096 0.038 0.096 0.135 0.059 0.038 0.115 0.135 0.096 0.135 0.135 0.096 0.078 0.098 0.098 0.098 0.059 0.118 0.137 0.137 0.020 0.137 0.118 0.137 0.078 0.059 0.098 0.039 0.137 0.137 0.137 0.059 0.078 0.059 0.115 0.135 0.135 0.115 0.115 0.115 0.115 0.135 0.135 0.135 0.096 0.038 0.135 0.135 0.058 0.154 0.154 0.038 0.039 0.135 0.154 0.135 0.154 0.096 0.096 0.154 0.038 0.115 0.154 0.154 0.135 0.154 0.058 0.154 0.154 0.115 0.135 0.096 0.078 0.115 0.039 0.157 0.157 0.157 0.137 0.137 0.137 0.115 0.000 ; 0.231 0.173 0.192 0.173 0.135 0.192 0.173 0.173 0.173 0.192 0.176 0.196 0.192 0.196 0.115 0.192 0.212 0.137 0.115 0.192 0.212 0.192 0.212 0.212 0.173 0.154 0.173 0.176 0.196 0.137 0.196 0.216 0.216 0.098 0.216 0.196 0.212 0.137 0.118 0.176 0.098 0.196 0.196 0.196 0.098 0.176 0.157 0.196 0.216 0.216 0.196 0.196 0.196 0.196 0.196 0.231 0.212 0.173 0.115 0.212 0.212 0.135 0.231 0.235 0.096 0.115 0.192 0.212 0.192 0.212 0.154 0.154 0.212 0.115 0.173 0.212 0.212 0.192 0.212 0.115 0.212 0.231 0.192 0.212 0.173 0.154 0.192 0.135 0.231 0.231 0.231 0.212 0.212 0.212 0.192 0.000 ; 0.250 0.192 0.216 0.192 0.154 0.212 0.192 0.173 0.192 0.212 0.196 0.216 0.212 0.216 0.135 0.212 0.231 0.157 0.135 0.212 0.231 0.212 0.231 0.231 0.192 0.173 0.192 0.192 0.212 0.135 0.173 0.212 0.019 0.212 0.212 0.192 0.212 0.154 0.154 0.192 0.115 0.212 0.212 0.212 0.137 0.216 0.176 0.235 0.255 0.255 0.235 0.235 0.235 0.235 0.235 0.231 0.231 0.192 0.135 0.212 0.212 0.135 0.231 0.231 0.115 0.135 0.212 0.231 0.231 0.255 0.196 0.196 0.255 0.137 0.216 0.255 0.255 0.235 0.255 0.157 0.255 0.255 0.216 0.235 0.196 0.176 0.235 0.157 0.255 0.255 0.255 0.235 0.235 0.235 0.212 0.000 ; 0.231 0.173 0.176 0.173 0.135 0.196 0.192 0.154 0.173 0.192 0.176 0.196 0.192 0.196 0.115 0.192 0.212 0.137 0.115 0.192 0.212 0.196 0.216 0.216 0.176 0.157 0.176 0.173 0.196 0.137 0.196 0.216 0.216 0.098 0.216 0.196 0.216 0.154 0.115 0.173 0.077 0.192 0.192 0.192 0.115 0.173 0.135 0.216 0.235 0.235 0.216 0.216 0.216 0.216 0.216 0.235 0.216 0.176 0.118 0.216 0.216 0.137 0.235 0.235 0.098 0.118 0.196 0.216 0.196 0.216 0.157 0.157 0.216 0.137 0.196 0.235 0.235 0.216 0.235 0.137 0.235 0.212 0.192 0.192 0.154 0.154 0.192 0.115 0.212 0.212 0.212 0.212 0.212 0.212 0.192 0.000 ; 0.231 0.154 0.176 0.135 0.135 0.173 0.173 0.154 0.154 0.173 0.176 0.196 0.192 0.196 0.096 0.192 0.216 0.137 0.096 0.192 0.212 0.196 0.212 0.212 0.173 0.154 0.173 0.173 0.176 0.118 0.176 0.196 0.216 0.196 0.196 0.176 0.196 0.137 0.096 0.173 0.096 0.192 0.192 0.192 0.096 0.173 0.135 0.212 0.231 0.231 0.212 0.212 0.212 0.212 0.212 0.235 0.216 0.176 0.118 0.216 0.216 0.137 0.235 0.235 0.096 0.096 0.192 0.212 0.192 0.212 0.154 0.154 0.212 0.118 0.196 0.235 0.235 0.216 0.235 0.137 0.235 0.212 0.176 0.196 0.157 0.137 0.173 0.098 0.196 0.196 0.196 0.176 0.173 0.154 0.192 0.000 ; 0.231 0.173 0.176 0.154 0.135 0.192 0.173 0.154 0.173 0.192 0.173 0.196 0.192 0.196 0.115 0.192 0.216 0.135 0.115 0.192 0.212 0.196 0.212 0.192 0.154 0.135 0.154 0.154 0.192 0.135 0.192 0.212 0.038 0.212 0.212 0.192 0.212 0.154 0.135 0.192 0.118 0.216 0.216 0.216 0.118 0.196 0.157 0.216 0.235 0.235 0.212 0.212 0.212 0.212 0.212 0.231 0.212 0.173 0.115 0.212 0.192 0.137 0.235 0.235 0.118 0.118 0.216 0.235 0.216 0.235 0.176 0.157 0.216 0.118 0.176 0.216 0.216 0.196 0.216 0.118 0.216 0.212 0.192 0.192 0.154 0.154 0.216 0.115 0.212 0.212 0.212 0.192 0.212 0.212 0.192 0.000 ; 0.250 0.212 0.216 0.192 0.173 0.231 0.212 0.192 0.212 0.231 0.235 0.235 0.250 0.235 0.154 0.231 0.255 0.173 0.157 0.231 0.250 0.235 0.250 0.231 0.212 0.192 0.212 0.212 0.231 0.173 0.231 0.255 0.173 0.058 0.173 0.154 0.173 0.135 0.096 0.115 0.077 0.135 0.154 0.231 0.135 0.212 0.173 0.231 0.231 0.250 0.231 0.231 0.231 0.255 0.255 0.275 0.255 0.216 0.157 0.255 0.255 0.176 0.269 0.269 0.154 0.154 0.250 0.269 0.250 0.269 0.212 0.212 0.269 0.157 0.212 0.231 0.231 0.231 0.231 0.154 0.231 0.231 0.212 0.231 0.196 0.173 0.255 0.154 0.231 0.231 0.231 0.231 0.231 0.212 0.231 0.000 ; 0.216 0.154 0.157 0.135 0.115 0.173 0.154 0.135 0.154 0.154 0.157 0.176 0.173 0.176 0.096 0.173 0.216 0.115 0.098 0.154 0.192 0.157 0.192 0.173 0.154 0.157 0.154 0.157 0.098 0.157 0.176 0.039 0.176 0.176 0.157 0.176 0.192 0.118 0.098 0.157 0.078 0.176 0.176 0.176 0.154 0.115 0.154 0.192 0.212 0.212 0.192 0.192 0.192 0.192 0.192 0.212 0.173 0.154 0.096 0.192 0.192 0.096 0.212 0.212 0.096 0.096 0.173 0.196 0.176 0.196 0.137 0.137 0.196 0.078 0.157 0.196 0.196 0.173 0.192 0.096 0.192 0.192 0.154 0.173 0.135 0.115 0.176 0.096 0.212 0.212 0.212 0.192 0.192 0.192 0.154 0.000 ; 0.235 0.173 0.176 0.173 0.154 0.192 0.192 0.173 0.173 0.192 0.196 0.196 0.173 0.196 0.115 0.192 0.176 0.135 0.118 0.196 0.212 0.196 0.212 0.192 0.173 0.118 0.173 0.173 0.192 0.115 0.173 0.192 0.192 0.077 0.192 0.173 0.192 0.154 0.118 0.192 0.115 0.212 0.212 0.212 0.115 0.192 0.154 0.212 0.231 0.231 0.212 0.212 0.212 0.212 0.212 0.231 0.212 0.173 0.115 0.212 0.212 0.135 0.231 0.231 0.020 0.118 0.216 0.235 0.216 0.235 0.176 0.176 0.235 0.118 0.173 0.212 0.212 0.192 0.212 0.115 0.212 0.212 0.196 0.212 0.173 0.154 0.216 0.135 0.231 0.231 0.231 0.212 0.216 0.212 0.196 0.000 ; 0.196 0.135 0.137 0.115 0.096 0.154 0.154 0.137 0.154 0.154 0.157 0.173 0.192 0.157 0.077 0.154 0.216 0.096 0.078 0.157 0.173 0.154 0.157 0.173 0.135 0.157 0.135 0.135 0.154 0.096 0.154 0.173 0.000 0.173 0.173 0.154 0.173 0.115 0.096 0.154 0.058 0.176 0.176 0.176 0.078 0.157 0.118 0.176 0.196 0.196 0.173 0.176 0.176 0.176 0.176 0.196 0.176 0.137 0.078 0.176 0.173 0.096 0.192 0.192 0.077 0.077 0.173 0.192 0.173 0.176 0.118 0.137 0.176 0.078 0.137 0.176 0.176 0.157 0.154 0.077 0.154 0.154 0.135 0.154 0.115 0.096 0.176 0.077 0.196 0.192 0.192 0.173 0.173 0.173 0.157 0.000 ; 0.235 0.192 0.176 0.173 0.154 0.192 0.196 0.176 0.192 0.173 0.196 0.192 0.135 0.192 0.115 0.192 0.137 0.135 0.118 0.196 0.192 0.196 0.196 0.212 0.173 0.098 0.173 0.173 0.173 0.115 0.173 0.192 0.192 0.077 0.192 0.173 0.192 0.135 0.137 0.192 0.115 0.212 0.212 0.212 0.115 0.192 0.154 0.212 0.235 0.235 0.216 0.216 0.216 0.216 0.216 0.235 0.216 0.176 0.118 0.216 0.216 0.137 0.235 0.235 0.137 0.115 0.212 0.231 0.212 0.231 0.173 0.173 0.231 0.135 0.192 0.231 0.231 0.212 0.231 0.135 0.231 0.231 0.196 0.216 0.176 0.157 0.216 0.118 0.235 0.235 0.235 0.196 0.192 0.192 0.196 0.000 ; 0.157 0.096 0.115 0.077 0.058 0.115 0.118 0.098 0.096 0.096 0.118 0.058 0.192 0.098 0.038 0.096 0.216 0.058 0.039 0.118 0.070 0.098 0.118 0.135 0.096 0.157 0.077 0.096 0.077 0.000 0.096 0.115 -0.039 0.115 0.000 0.096 0.135 0.077 0.059 0.077 0.038 0.115 0.115 0.115 0.039 0.058 0.019 0.115 0.173 0.135 0.137 0.137 0.137 0.137 0.137 0.196 0.137 0.098 0.039 0.137 0.115 0.038 0.135 0.135 0.019 0.038 0.135 0.154 0.135 0.154 0.096 0.096 0.154 0.038 0.115 0.154 0.154 0.135 0.154 0.058 0.154 0.154 0.115 0.154 0.096 0.077 0.137 0.058 0.154 0.154 0.154 0.135 0.135 0.154 0.118 0.000 ; 0.235 0.192 0.192 0.173 0.154 0.192 0.196 0.176 0.192 0.173 0.196 0.192 0.212 0.196 0.115 0.192 0.212 0.135 0.118 0.196 0.192 0.196 0.196 0.212 0.196 0.176 0.173 0.173 0.173 0.135 0.173 0.192 0.000 0.192 0.192 0.173 0.192 0.135 0.115 0.173 0.096 0.216 0.216 0.216 0.118 0.196 0.157 0.216 0.235 0.235 0.212 0.212 0.212 0.212 0.212 0.231 0.212 0.173 0.118 0.216 0.216 0.137 0.235 0.235 0.118 0.118 0.212 0.212 0.192 0.212 0.154 0.173 0.212 0.115 0.173 0.212 0.231 0.212 0.231 0.135 0.231 0.231 0.192 0.212 0.173 0.154 0.216 0.135 0.231 0.231 0.231 0.212 0.212 0.212 0.196 0.000 ; 0.255 0.192 0.192 0.196 0.154 0.212 0.216 0.196 0.192 0.173 0.216 0.192 0.192 0.216 0.135 0.212 0.231 0.154 0.137 0.216 0.212 0.216 0.216 0.231 0.154 0.176 0.192 0.192 0.173 0.135 0.192 0.235 0.020 0.235 0.235 0.216 0.235 0.176 0.157 0.212 0.115 0.212 0.212 0.212 0.115 0.192 0.154 0.212 0.231 0.255 0.235 0.235 0.235 0.235 0.235 0.255 0.235 0.176 0.118 0.196 0.196 0.118 0.216 0.216 0.098 0.135 0.212 0.231 0.212 0.231 0.173 0.173 0.231 0.115 0.216 0.231 0.231 0.212 0.231 0.135 0.231 0.231 0.192 0.192 0.173 0.154 0.235 0.115 0.231 0.231 0.231 0.212 0.212 0.235 0.216 0.000 ; ...
              0.255 0.173 0.192 0.179 0.154 0.212 0.216 0.196 0.192 0.192 0.216 0.212 0.192 0.216 0.135 0.192 0.135 0.119 0.137 0.216 0.212 0.216 0.000 0.231 0.212 0.059 0.192 0.192 0.192 0.135 0.192 0.235 0.039 0.231 0.235 0.216 0.235 0.176 0.157 0.216 0.137 0.235 0.235 0.235 0.137 0.216 0.176 0.235 0.255 0.255 0.235 0.235 0.235 0.235 0.231 0.250 0.231 0.192 0.135 0.231 0.231 0.154 0.250 0.250 0.135 0.135 0.231 0.250 0.231 0.250 0.192 0.192 0.250 0.154 0.212 0.250 0.250 0.231 0.182 0.154 0.250 0.250 0.212 0.231 0.192 0.173 0.182 0.154 0.250 0.250 0.250 0.231 0.231 0.231 0.216 0.000 ; 0.137 0.173 0.192 0.179 0.058 0.096 0.216 0.196 0.192 0.192 0.152 0.096 0.192 0.098 0.135 0.192 0.135 0.119 0.137 0.216 0.096 0.216 0.098 0.115 0.212 0.176 0.192 0.058 0.192 0.019 0.077 0.118 -0.078 0.115 0.115 0.096 0.115 0.058 0.038 0.096 0.019 0.115 0.118 0.118 0.020 0.098 0.059 0.118 0.137 0.137 0.118 0.118 0.115 0.115 0.115 0.135 0.115 0.077 0.019 0.115 0.115 0.039 0.137 0.137 0.020 0.020 0.118 0.137 0.118 0.137 0.058 0.058 0.115 0.019 0.077 0.115 0.115 0.096 0.115 0.020 0.118 0.118 0.078 0.098 0.059 0.039 0.098 0.020 0.115 0.115 0.115 0.096 0.096 0.096 0.216 0.000 ; 0.255 0.154 0.173 0.164 0.154 0.212 0.196 0.176 0.173 0.173 0.136 0.212 0.173 0.216 0.137 0.173 0.115 0.104 0.118 0.196 0.212 0.196 0.216 0.231 0.192 0.157 0.173 0.192 0.173 0.115 0.192 0.235 0.020 0.231 0.212 0.192 0.212 0.173 0.157 0.196 0.137 0.216 0.216 0.216 0.118 0.196 0.157 0.216 0.235 0.235 0.231 0.231 0.231 0.231 0.231 0.250 0.231 0.192 0.135 0.231 0.231 0.154 0.212 0.212 0.115 0.115 0.212 0.212 0.212 0.212 0.173 0.173 0.212 0.135 0.216 0.231 0.231 0.212 0.231 0.135 0.231 0.231 0.212 0.212 0.173 0.154 0.182 0.137 0.255 0.255 0.255 0.235 0.235 0.235 0.196 0.000 ; 0.235 0.173 0.192 0.179 0.135 0.192 0.216 0.154 0.192 0.173 0.152 0.212 0.192 0.196 0.118 0.192 0.135 0.119 0.118 0.216 0.192 0.216 0.196 0.212 0.212 0.176 0.192 0.173 0.192 0.115 0.173 0.216 0.020 0.212 0.192 0.196 0.216 0.157 0.137 0.192 0.115 0.212 0.212 0.212 0.115 0.192 0.154 0.212 0.231 0.231 0.212 0.196 0.196 0.196 0.196 0.216 0.196 0.157 0.098 0.196 0.196 0.115 0.216 0.216 0.098 0.098 0.196 0.216 0.196 0.216 0.157 0.157 0.216 0.115 0.192 0.231 0.231 0.212 0.231 0.135 0.231 0.231 0.192 0.212 0.173 0.157 0.167 0.135 0.231 0.231 0.231 0.212 0.212 0.212 0.216 0.000 ; 0.255 0.192 0.135 0.134 0.154 0.212 0.212 0.173 0.212 0.192 0.121 0.212 0.154 0.212 0.137 0.212 0.077 0.075 0.137 0.212 0.212 0.157 0.216 0.231 0.154 0.173 0.135 0.192 0.154 0.135 0.192 0.235 0.020 0.231 0.212 0.216 0.231 0.173 0.154 0.212 0.135 0.231 0.231 0.231 0.135 0.212 0.173 0.231 0.250 0.250 0.231 0.235 0.235 0.235 0.235 0.255 0.235 0.196 0.137 0.235 0.235 0.157 0.250 0.255 0.137 0.137 0.235 0.255 0.235 0.255 0.196 0.196 0.255 0.059 0.216 0.250 0.250 0.231 0.250 0.154 0.250 0.250 0.212 0.231 0.192 0.173 0.182 0.154 0.250 0.250 0.250 0.231 0.231 0.231 0.212 0.000 ; 0.216 0.154 0.115 0.119 0.115 0.173 0.173 0.135 0.173 0.154 0.106 0.096 0.135 0.137 0.098 0.154 0.077 0.060 0.098 0.157 0.173 0.154 0.192 0.192 0.135 0.137 0.137 0.135 0.135 0.115 0.135 0.196 0.000 0.192 0.176 0.176 0.192 0.115 0.115 0.154 0.096 0.192 0.192 0.192 0.096 0.173 0.135 0.192 0.212 0.216 0.196 0.196 0.196 0.196 0.196 0.216 0.196 0.157 0.096 0.192 0.192 0.115 0.212 0.212 0.096 0.096 0.192 0.192 0.173 0.192 0.135 0.135 0.192 0.096 0.154 0.192 0.216 0.196 0.216 0.118 0.216 0.216 0.176 0.196 0.157 0.137 0.152 0.098 0.196 0.196 0.196 0.176 0.176 0.176 0.157 0.000 ; 0.176 0.115 0.137 0.164 0.096 0.135 0.135 0.096 0.135 0.115 0.136 0.077 0.173 0.098 0.059 0.115 0.154 0.090 0.059 0.118 0.154 0.096 0.154 0.154 0.192 0.098 0.096 0.096 0.078 0.078 0.077 0.157 -0.039 0.154 0.137 0.135 0.154 0.096 0.078 0.098 0.059 0.157 0.157 0.157 0.059 0.137 0.098 0.157 0.176 0.176 0.157 0.157 0.157 0.157 0.157 0.173 0.154 0.115 0.058 0.154 0.154 0.077 0.173 0.173 0.058 0.058 0.154 0.173 0.154 0.173 0.115 0.115 0.173 0.058 0.137 0.154 0.154 0.135 0.154 0.058 0.154 0.154 0.115 0.135 0.115 0.096 0.121 0.058 0.173 0.173 0.173 0.154 0.154 0.154 0.118 0.000 ; 0.196 0.135 0.157 0.104 0.096 0.157 0.154 0.115 0.154 0.135 0.091 0.135 0.154 0.137 0.078 0.154 0.173 0.045 0.078 0.137 0.154 0.135 0.173 0.176 0.115 0.118 0.135 0.115 0.154 0.077 0.135 0.173 -0.020 0.173 0.157 0.154 0.173 0.115 0.098 0.154 0.077 0.173 0.173 0.173 0.077 0.157 0.096 0.154 0.154 0.154 0.154 0.154 0.154 0.154 0.173 0.192 0.173 0.135 0.077 0.173 0.173 0.096 0.192 0.196 0.059 0.078 0.157 0.176 0.157 0.176 0.118 0.118 0.176 0.078 0.154 0.192 0.192 0.173 0.192 0.096 0.192 0.192 0.154 0.176 0.137 0.118 0.136 0.078 0.196 0.196 0.196 0.176 0.176 0.154 0.137 0.000 ; 0.157 0.096 0.118 0.096 0.058 0.118 0.115 0.077 0.115 0.096 0.115 0.077 0.115 0.078 0.039 0.096 0.135 0.119 0.039 0.098 0.115 0.096 0.135 0.137 0.212 0.078 0.096 0.098 0.115 0.038 0.118 0.135 -0.059 0.135 0.118 0.115 0.135 0.077 0.059 0.115 0.019 0.135 0.115 0.115 0.000 0.058 0.058 0.115 0.154 0.154 0.135 0.135 0.135 0.135 0.135 0.154 0.135 0.098 0.039 0.118 0.118 0.039 0.137 0.137 0.020 0.039 0.135 0.154 0.135 0.154 0.096 0.096 0.154 0.058 0.115 0.154 0.157 0.137 0.157 0.059 0.157 0.157 0.118 0.137 0.098 0.078 0.106 0.058 0.154 0.154 0.154 0.135 0.135 0.135 0.098 0.000 ; 0.235 0.154 0.196 0.173 0.135 0.196 0.173 0.154 0.192 0.173 0.173 0.192 0.173 0.176 0.118 0.173 0.192 0.119 0.118 0.176 0.192 0.154 0.192 0.216 0.212 0.137 0.173 0.176 0.192 0.115 0.196 0.192 0.020 0.216 0.176 0.192 0.192 0.154 0.118 0.192 0.077 0.192 0.216 0.196 0.118 0.196 0.157 0.216 0.235 0.235 0.216 0.216 0.216 0.192 0.192 0.212 0.192 0.154 0.096 0.192 0.192 0.115 0.212 0.212 0.096 0.096 0.192 0.212 0.192 0.212 0.173 0.173 0.231 0.115 0.192 0.231 0.231 0.212 0.152 0.135 0.231 0.231 0.192 0.212 0.173 0.154 0.167 0.135 0.235 0.235 0.235 0.216 0.216 0.216 0.176 0.000 ; 0.255 0.173 0.216 0.192 0.154 0.216 0.192 0.173 0.212 0.192 0.192 0.212 0.192 0.196 0.137 0.192 0.212 0.119 0.137 0.196 0.212 0.173 0.212 0.235 0.173 0.157 0.192 0.196 0.212 0.135 0.216 0.212 0.020 0.235 0.196 0.212 0.212 0.173 0.157 0.212 0.096 0.212 0.235 0.216 0.137 0.216 0.157 0.216 0.235 0.235 0.216 0.216 0.216 0.216 0.231 0.255 0.235 0.196 0.137 0.235 0.235 0.157 0.255 0.231 0.115 0.115 0.212 0.231 0.212 0.231 0.173 0.196 0.250 0.154 0.212 0.250 0.250 0.231 0.250 0.154 0.255 0.231 0.192 0.212 0.173 0.154 0.182 0.135 0.231 0.231 0.255 0.212 0.212 0.212 0.196 0.000 ; 0.235 0.173 0.196 0.173 0.154 0.196 0.192 0.176 0.192 0.173 0.192 0.192 0.192 0.176 0.096 0.192 0.212 0.135 0.118 0.176 0.192 0.173 0.212 0.216 0.154 0.216 0.173 0.176 0.192 0.115 0.196 0.212 0.020 0.216 0.196 0.192 0.212 0.154 0.137 0.192 0.096 0.192 0.216 0.216 0.118 0.196 0.157 0.216 0.235 0.235 0.216 0.216 0.216 0.216 0.216 0.231 0.212 0.173 0.115 0.212 0.212 0.115 0.216 0.216 0.098 0.098 0.196 0.216 0.212 0.231 0.173 0.173 0.231 0.115 0.192 0.231 0.231 0.212 0.231 0.135 0.231 0.231 0.176 0.196 0.157 0.137 0.212 0.098 0.196 0.235 0.235 0.216 0.216 0.216 0.176 0.000 ; 0.235 0.173 0.196 0.173 0.154 0.196 0.192 0.154 0.118 0.192 0.192 0.192 0.192 0.176 0.078 0.192 0.058 0.135 0.059 0.176 0.216 0.173 0.212 0.216 0.154 0.154 0.173 0.176 0.192 0.115 0.196 0.212 0.020 0.216 0.196 0.192 0.212 0.154 0.137 0.192 0.038 0.216 0.216 0.216 0.115 0.192 0.154 0.212 0.231 0.231 0.212 0.212 0.212 0.212 0.212 0.154 0.192 0.154 0.096 0.192 0.192 0.019 0.077 0.077 0.096 0.135 0.216 0.231 0.212 0.231 0.096 0.096 0.154 0.038 0.115 0.115 0.212 0.231 0.167 0.115 0.212 0.212 0.192 0.192 0.154 0.154 0.212 0.019 0.115 0.154 0.154 0.135 0.216 0.216 0.176 0.000 ; 0.235 0.176 0.196 0.173 0.154 0.196 0.196 0.154 0.173 0.192 0.192 0.216 0.192 0.196 0.078 0.196 0.077 0.135 0.059 0.176 0.216 0.173 0.212 0.216 0.154 0.154 0.173 0.176 0.192 0.115 0.196 0.212 0.020 0.216 0.196 0.192 0.212 0.154 0.137 0.196 0.038 0.216 0.196 0.216 0.115 0.196 0.157 0.216 0.235 0.235 0.216 0.212 0.212 0.212 0.212 0.154 0.212 0.173 0.115 0.212 0.212 0.038 0.235 0.212 0.019 0.192 0.212 0.192 0.212 0.154 0.173 0.212 0.115 0.173 0.212 0.231 0.235 0.216 0.235 0.137 0.235 0.235 0.196 0.216 0.176 0.157 0.212 0.020 0.235 0.235 0.235 0.216 0.216 0.216 0.176 0.000 ; 0.196 0.137 0.154 0.154 0.115 0.176 0.154 0.154 0.135 0.154 0.154 0.173 0.173 0.154 0.137 0.173 0.176 0.096 0.078 0.157 0.176 0.157 0.173 0.176 0.135 0.115 0.135 0.137 0.173 0.096 0.157 0.173 -0.020 0.196 0.157 0.154 0.173 0.115 0.098 0.154 0.077 0.176 0.176 0.176 0.096 0.157 0.118 0.176 0.196 0.192 0.173 0.173 0.173 0.173 0.173 0.192 0.173 0.135 0.077 0.192 0.176 0.098 0.196 0.196 0.078 0.078 0.176 0.196 0.176 0.196 0.135 0.135 0.192 0.038 0.154 0.192 0.192 0.173 0.192 0.096 0.176 0.176 0.137 0.157 0.118 0.098 0.192 0.078 0.137 0.176 0.196 0.154 0.154 0.154 0.157 0.000 ; 0.255 0.196 0.212 0.192 0.173 0.192 0.212 0.173 0.192 0.212 0.212 0.231 0.231 0.212 0.157 0.212 0.235 0.154 0.135 0.216 0.235 0.216 0.231 0.235 0.173 0.173 0.192 0.196 0.216 0.135 0.216 0.231 0.020 0.235 0.216 0.212 0.231 0.173 0.157 0.212 0.135 0.235 0.235 0.235 0.135 0.216 0.176 0.212 0.231 0.231 0.212 0.235 0.231 0.231 0.231 0.250 0.231 0.192 0.135 0.231 0.212 0.157 0.235 0.235 0.137 0.137 0.235 0.235 0.235 0.255 0.192 0.192 0.250 0.135 0.212 0.250 0.250 0.231 0.182 0.154 0.250 0.250 0.212 0.231 0.192 0.173 0.231 0.154 0.250 0.212 0.212 0.212 0.212 0.192 0.216 0.000 ; ...
              0.275 0.216 0.231 0.212 0.192 0.231 0.231 0.192 0.212 0.231 0.231 0.231 0.231 0.212 0.157 0.231 0.255 0.173 0.154 0.216 0.255 0.235 0.250 0.255 0.192 0.192 0.212 0.216 0.235 0.154 0.235 0.255 0.039 0.255 0.235 0.231 0.250 0.196 0.176 0.212 0.154 0.255 0.235 0.255 0.154 0.235 0.196 0.231 0.231 0.231 0.250 0.250 0.250 0.250 0.250 0.269 0.250 0.212 0.154 0.250 0.231 0.154 0.250 0.250 0.135 0.154 0.231 0.250 0.231 0.250 0.192 0.212 0.275 0.173 0.231 0.269 0.269 0.250 0.269 0.173 0.269 0.269 0.231 0.250 0.212 0.192 0.250 0.173 0.269 0.269 0.269 0.250 0.250 0.250 0.216 0.000 ; 0.255 0.196 0.212 0.192 0.173 0.212 0.212 0.173 0.192 0.212 0.212 0.212 0.235 0.192 0.137 0.212 0.235 0.154 0.135 0.196 0.235 0.216 0.231 0.231 0.173 0.173 0.192 0.196 0.216 0.135 0.212 0.235 0.020 0.235 0.216 0.212 0.231 0.176 0.157 0.192 0.135 0.235 0.235 0.235 0.135 0.216 0.173 0.212 0.250 0.250 0.231 0.231 0.231 0.231 0.231 0.250 0.231 0.192 0.135 0.231 0.212 0.154 0.250 0.250 0.135 0.135 0.231 0.250 0.231 0.250 0.192 0.192 0.250 0.154 0.212 0.250 0.250 0.231 0.250 0.154 0.231 0.231 0.192 0.212 0.173 0.154 0.231 0.135 0.250 0.250 0.250 0.231 0.231 0.231 0.196 0.000 ; 0.255 0.196 0.212 0.192 0.173 0.212 0.212 0.173 0.192 0.212 0.192 0.231 0.235 0.212 0.137 0.212 0.235 0.154 0.135 0.216 0.235 0.216 0.231 0.231 0.192 0.173 0.192 0.192 0.216 0.157 0.212 0.235 0.020 0.235 0.231 0.216 0.231 0.176 0.157 0.212 0.135 0.235 0.235 0.235 0.135 0.212 0.173 0.212 0.250 0.250 0.231 0.231 0.231 0.231 0.235 0.255 0.235 0.196 0.137 0.235 0.235 0.154 0.255 0.255 0.137 0.137 0.235 0.255 0.235 0.250 0.192 0.192 0.250 0.154 0.212 0.250 0.255 0.231 0.250 0.154 0.250 0.250 0.212 0.231 0.176 0.176 0.231 0.137 0.255 0.255 0.255 0.235 0.235 0.235 0.216 0.000 ; 0.255 0.196 0.212 0.192 0.176 0.212 0.212 0.173 0.192 0.212 0.192 0.212 0.216 0.192 0.135 0.212 0.235 0.154 0.135 0.196 0.235 0.216 0.231 0.231 0.173 0.173 0.192 0.196 0.216 0.157 0.212 0.235 0.020 0.235 0.231 0.196 0.235 0.176 0.157 0.192 0.135 0.235 0.216 0.235 0.135 0.212 0.173 0.235 0.250 0.250 0.231 0.235 0.235 0.235 0.231 0.235 0.216 0.176 0.118 0.216 0.216 0.137 0.235 0.231 0.115 0.115 0.212 0.231 0.212 0.231 0.173 0.196 0.231 0.135 0.192 0.231 0.231 0.212 0.231 0.135 0.255 0.250 0.212 0.231 0.192 0.173 0.231 0.135 0.250 0.250 0.255 0.235 0.235 0.235 0.196 0.000 ; 0.255 0.196 0.212 0.192 0.173 0.212 0.212 0.192 0.192 0.212 0.192 0.212 0.216 0.192 0.135 0.212 0.235 0.154 0.135 0.196 0.235 0.216 0.231 0.231 0.196 0.173 0.192 0.196 0.216 0.157 0.212 0.235 0.020 0.231 0.231 0.196 0.235 0.176 0.157 0.192 0.135 0.235 0.216 0.235 0.135 0.212 0.173 0.231 0.250 0.250 0.231 0.235 0.235 0.235 0.231 0.250 0.231 0.192 0.135 0.231 0.231 0.157 0.255 0.255 0.137 0.137 0.235 0.231 0.231 0.250 0.192 0.192 0.250 0.154 0.212 0.250 0.250 0.231 0.250 0.154 0.250 0.231 0.212 0.231 0.192 0.173 0.235 0.154 0.250 0.250 0.250 0.231 0.231 0.231 0.196 0.000 ; 0.255 0.192 0.212 0.192 0.173 0.212 0.212 0.192 0.212 0.216 0.192 0.231 0.235 0.212 0.135 0.212 0.235 0.154 0.135 0.216 0.235 0.216 0.231 0.231 0.196 0.173 0.192 0.196 0.216 0.157 0.212 0.235 0.019 0.231 0.231 0.196 0.235 0.176 0.157 0.212 0.135 0.231 0.235 0.235 0.135 0.212 0.173 0.231 0.250 0.250 0.235 0.235 0.231 0.235 0.235 0.255 0.235 0.196 0.137 0.235 0.235 0.157 0.255 0.255 0.137 0.137 0.235 0.255 0.235 0.250 0.192 0.192 0.250 0.135 0.212 0.250 0.250 0.231 0.250 0.154 0.250 0.250 0.212 0.231 0.192 0.173 0.235 0.154 0.250 0.250 0.250 0.231 0.231 0.231 0.216 0.000 ; 0.192 0.135 0.154 0.135 0.115 0.154 0.154 0.135 0.137 0.157 0.135 0.154 0.157 0.135 0.077 0.154 0.176 0.096 0.077 0.154 0.176 0.157 0.176 0.173 0.137 0.115 0.135 0.137 0.157 0.098 0.154 0.176 -0.020 0.173 0.173 0.137 0.176 0.118 0.098 0.154 0.077 0.173 0.173 0.173 0.077 0.154 0.115 0.173 0.192 0.192 0.176 0.173 0.173 0.176 0.173 0.192 0.173 0.137 0.077 0.173 0.173 0.096 0.192 0.192 0.077 0.077 0.173 0.192 0.173 0.192 0.135 0.135 0.196 0.077 0.135 0.173 0.173 0.154 0.173 0.077 0.196 0.192 0.154 0.173 0.135 0.115 0.176 0.077 0.192 0.192 0.192 0.173 0.173 0.173 0.154 0.000 ; 0.255 0.212 0.212 0.192 0.192 0.231 0.231 0.212 0.216 0.235 0.212 0.231 0.231 0.192 0.154 0.196 0.212 0.154 0.154 0.231 0.231 0.216 0.255 0.250 0.216 0.192 0.192 0.137 0.216 0.176 0.212 0.235 0.038 0.250 0.250 0.216 0.255 0.196 0.173 0.212 0.154 0.250 0.235 0.212 0.135 0.231 0.173 0.231 0.269 0.269 0.255 0.235 0.250 0.255 0.250 0.250 0.231 0.192 0.154 0.231 0.231 0.154 0.255 0.255 0.137 0.137 0.235 0.255 0.235 0.255 0.216 0.212 0.250 0.154 0.173 0.250 0.250 0.231 0.250 0.154 0.255 0.255 0.216 0.235 0.196 0.176 0.255 0.137 0.255 0.255 0.231 0.231 0.231 0.212 0.231 0.000 ; 0.196 0.154 0.154 0.157 0.135 0.173 0.173 0.154 0.157 0.176 0.154 0.173 0.173 0.154 0.096 0.157 0.192 0.115 0.096 0.173 0.173 0.176 0.196 0.192 0.157 0.115 0.154 0.157 0.176 0.118 0.173 0.176 0.019 0.192 0.192 0.157 0.196 0.137 0.115 0.173 0.098 0.192 0.196 0.192 0.077 0.173 0.135 0.192 0.212 0.212 0.196 0.176 0.192 0.196 0.192 0.192 0.173 0.135 0.077 0.192 0.196 0.098 0.216 0.216 0.098 0.098 0.196 0.216 0.196 0.216 0.157 0.157 0.216 0.098 0.176 0.216 0.216 0.196 0.216 0.118 0.216 0.212 0.173 0.192 0.154 0.135 0.196 0.096 0.212 0.212 0.212 0.192 0.196 0.196 0.173 0.000 ; 0.235 0.192 0.192 0.196 0.173 0.212 0.212 0.192 0.196 0.216 0.216 0.212 0.212 0.192 0.135 0.196 0.231 0.154 0.135 0.212 0.212 0.216 0.235 0.231 0.196 0.173 0.192 0.196 0.216 0.157 0.216 0.216 0.019 0.231 0.231 0.196 0.235 0.176 0.154 0.212 0.137 0.231 0.235 0.231 0.115 0.212 0.173 0.231 0.250 0.250 0.235 0.216 0.231 0.235 0.235 0.231 0.212 0.196 0.137 0.235 0.235 0.135 0.255 0.255 0.137 0.137 0.235 0.255 0.235 0.255 0.196 0.196 0.231 0.135 0.192 0.231 0.231 0.212 0.231 0.135 0.231 0.231 0.216 0.216 0.176 0.157 0.235 0.137 0.235 0.235 0.235 0.216 0.216 0.216 0.212 0.000 ; 0.157 0.115 0.115 0.118 0.096 0.135 0.135 0.115 0.118 0.137 0.115 0.115 0.135 0.115 0.058 0.118 0.154 0.058 0.058 0.135 0.135 0.137 0.157 0.154 0.118 0.096 0.115 0.118 0.137 0.078 0.137 0.118 -0.039 0.154 0.154 0.118 0.157 0.096 0.077 0.135 0.059 0.154 0.157 0.154 0.038 0.135 0.096 0.154 0.173 0.173 0.157 0.137 0.154 0.157 0.154 0.154 0.135 0.118 0.059 0.137 0.137 0.059 0.157 0.157 0.039 0.059 0.137 0.173 0.154 0.173 0.115 0.115 0.173 0.058 0.135 0.173 0.173 0.157 0.176 0.078 0.176 0.176 0.137 0.157 0.118 0.098 0.154 0.077 0.173 0.173 0.173 0.154 0.154 0.154 0.135 0.000 ; 0.235 0.192 0.192 0.196 0.173 0.212 0.212 0.192 0.196 0.216 0.192 0.212 0.212 0.192 0.135 0.196 0.231 0.154 0.135 0.212 0.212 0.216 0.235 0.212 0.196 0.173 0.192 0.192 0.216 0.157 0.216 0.216 0.019 0.231 0.231 0.196 0.235 0.173 0.154 0.212 0.137 0.231 0.235 0.231 0.115 0.212 0.173 0.231 0.250 0.250 0.235 0.216 0.231 0.235 0.231 0.231 0.231 0.196 0.137 0.231 0.231 0.154 0.250 0.250 0.039 0.137 0.235 0.255 0.235 0.255 0.196 0.196 0.255 0.137 0.196 0.235 0.235 0.216 0.235 0.137 0.235 0.231 0.212 0.231 0.192 0.173 0.235 0.154 0.250 0.250 0.250 0.212 0.212 0.212 0.212 0.000 ; 0.216 0.173 0.192 0.176 0.154 0.192 0.192 0.173 0.176 0.196 0.173 0.192 0.192 0.173 0.115 0.196 0.212 0.135 0.115 0.192 0.192 0.196 0.216 0.216 0.176 0.154 0.173 0.173 0.196 0.115 0.196 0.196 -0.020 0.216 0.212 0.173 0.216 0.154 0.135 0.192 0.118 0.212 0.216 0.212 0.096 0.173 0.157 0.192 0.212 0.212 0.216 0.196 0.216 0.216 0.212 0.231 0.216 0.176 0.118 0.212 0.212 0.115 0.231 0.231 0.096 0.096 0.212 0.231 0.212 0.231 0.173 0.173 0.231 0.115 0.192 0.231 0.231 0.216 0.235 0.137 0.235 0.235 0.196 0.216 0.157 0.137 0.216 0.118 0.235 0.235 0.235 0.216 0.216 0.216 0.192 0.000 ; 0.157 0.115 0.115 0.098 0.077 0.115 0.115 0.098 0.115 0.115 0.115 0.077 0.115 0.077 0.058 0.078 0.135 0.058 0.058 0.115 0.115 0.078 0.137 0.137 0.098 0.096 0.077 0.115 0.078 0.058 0.098 0.118 -0.039 0.137 0.135 0.118 0.137 0.077 0.058 0.077 0.039 0.135 0.118 0.115 0.212 0.038 0.039 0.137 0.135 0.154 0.135 0.137 0.137 0.137 0.135 0.154 0.137 0.096 0.039 0.135 0.115 0.058 0.154 0.154 0.038 0.038 0.135 0.154 0.135 0.157 0.098 0.098 0.157 0.000 0.118 0.157 0.173 0.154 0.173 0.019 0.173 0.173 0.115 0.137 0.098 0.078 0.118 0.039 0.157 0.157 0.154 0.135 0.135 0.135 0.115 0.000 ; 0.196 0.154 0.154 0.157 0.115 0.173 0.154 0.157 0.173 0.154 0.154 0.176 0.173 0.154 0.077 0.157 0.173 0.096 0.096 0.173 0.173 0.176 0.196 0.196 0.154 0.135 0.154 0.135 0.176 0.096 0.157 0.176 0.000 0.196 0.173 0.157 0.192 0.135 0.096 0.173 0.078 0.192 0.173 0.176 0.212 0.173 0.137 0.176 0.192 0.212 0.196 0.173 0.196 0.196 0.192 0.196 0.196 0.157 0.098 0.192 0.192 0.096 0.216 0.216 0.098 0.098 0.196 0.216 0.196 0.216 0.157 0.154 0.212 0.096 0.173 0.212 0.212 0.192 0.212 0.115 0.176 0.176 0.157 0.176 0.137 0.118 0.196 0.078 0.176 0.216 0.216 0.196 0.196 0.196 0.173 0.000 ; 0.196 0.173 0.154 0.157 0.115 0.173 0.173 0.157 0.173 0.173 0.173 0.115 0.173 0.157 0.096 0.176 0.196 0.115 0.096 0.173 0.173 0.176 0.196 0.196 0.154 0.135 0.154 0.154 0.157 0.096 0.176 0.176 0.000 0.196 0.192 0.157 0.196 0.135 0.115 0.157 0.098 0.192 0.192 0.176 0.192 0.173 0.137 0.196 0.192 0.212 0.196 0.173 0.196 0.196 0.192 0.216 0.196 0.157 0.098 0.192 0.192 0.115 0.216 0.216 0.098 0.078 0.196 0.216 0.192 0.212 0.154 0.154 0.212 0.096 0.173 0.212 0.212 0.192 0.212 0.096 0.212 0.212 0.173 0.192 0.157 0.137 0.196 0.078 0.216 0.216 0.216 0.196 0.196 0.176 0.173 0.000 ; ...
              0.157 0.115 0.096 0.096 0.077 0.115 0.115 0.118 0.135 0.115 0.115 0.019 0.157 0.020 0.038 0.038 0.058 0.058 0.000 0.038 0.135 0.000 0.135 0.157 0.115 0.096 0.118 0.038 0.118 0.058 0.039 0.118 -0.039 0.157 0.135 0.118 0.157 0.019 0.077 0.020 0.000 0.154 0.077 0.059 0.212 0.019 0.020 0.137 0.154 0.173 0.157 0.115 0.157 0.157 0.135 0.176 0.157 0.098 0.059 0.137 0.135 0.058 0.176 0.176 0.039 0.039 0.137 0.154 0.096 0.135 0.077 0.096 0.135 -0.019 0.077 0.135 0.135 0.137 0.157 -0.039 0.157 0.157 0.118 0.137 0.098 0.078 0.157 0.039 0.157 0.157 0.157 0.137 0.137 0.137 0.038 0.000 ; 0.157 0.115 0.096 0.115 0.077 0.135 0.115 0.098 0.115 0.115 0.115 0.212 0.137 0.039 0.038 0.058 0.096 0.058 0.039 0.096 0.115 0.039 0.077 0.137 0.115 0.077 0.039 0.038 0.077 0.058 0.059 0.154 -0.039 0.137 0.135 0.118 0.137 0.077 0.077 0.039 0.020 0.135 0.096 0.098 0.192 0.192 0.020 0.137 0.135 0.154 0.137 0.135 0.137 0.137 0.135 0.157 0.137 0.098 0.039 0.157 0.115 0.058 0.157 0.157 0.137 0.058 0.154 0.173 0.135 0.173 0.115 0.115 0.157 -0.020 0.098 0.157 0.157 0.118 0.157 -0.059 0.157 0.173 0.135 0.154 0.115 0.096 0.118 0.077 0.173 0.173 0.157 0.154 0.154 0.154 0.096 0.000 ; 0.235 0.212 0.192 0.192 0.154 0.212 0.192 0.196 0.212 0.192 0.212 0.231 0.235 0.196 0.137 0.212 0.231 0.154 0.135 0.212 0.212 0.196 0.231 0.235 0.192 0.173 0.196 0.192 0.212 0.154 0.212 0.231 0.019 0.235 0.212 0.196 0.235 0.173 0.135 0.196 0.135 0.231 0.231 0.216 0.096 0.231 0.176 0.216 0.231 0.250 0.235 0.212 0.212 0.235 0.231 0.235 0.235 0.196 0.135 0.235 0.231 0.135 0.255 0.255 0.157 0.135 0.235 0.255 0.235 0.255 0.196 0.196 0.255 0.118 0.216 0.255 0.255 0.235 0.255 0.157 0.255 0.255 0.216 0.231 0.231 0.173 0.235 0.173 0.135 0.212 0.250 0.250 0.235 0.212 0.212 0.000 ; 0.255 0.231 0.212 0.212 0.173 0.231 0.212 0.216 0.231 0.212 0.231 0.212 0.255 0.216 0.157 0.231 0.250 0.173 0.154 0.231 0.250 0.216 0.250 0.255 0.212 0.196 0.216 0.212 0.231 0.173 0.212 0.250 0.038 0.255 0.231 0.216 0.255 0.173 0.173 0.216 0.154 0.250 0.250 0.235 0.154 0.154 0.196 0.235 0.250 0.269 0.255 0.231 0.250 0.255 0.250 0.255 0.255 0.216 0.154 0.255 0.250 0.154 0.275 0.269 0.135 0.154 0.255 0.250 0.231 0.250 0.192 0.192 0.250 0.154 0.212 0.275 0.275 0.255 0.275 0.176 0.275 0.275 0.235 0.255 0.216 0.196 0.255 0.157 0.269 0.269 0.269 0.250 0.250 0.250 0.231 0.000 ; 0.235 0.196 0.173 0.173 0.157 0.192 0.192 0.176 0.192 0.192 0.192 0.231 0.216 0.176 0.118 0.192 0.212 0.135 0.115 0.192 0.216 0.176 0.212 0.216 0.173 0.157 0.176 0.173 0.192 0.135 0.173 0.212 0.020 0.216 0.192 0.176 0.216 0.135 0.135 0.176 0.115 0.212 0.212 0.212 0.077 0.038 0.157 0.196 0.212 0.231 0.212 0.192 0.212 0.212 0.212 0.216 0.216 0.176 0.115 0.216 0.212 0.135 0.235 0.231 0.096 0.096 0.216 0.212 0.212 0.231 0.173 0.173 0.231 0.135 0.192 0.212 0.212 0.192 0.212 0.115 0.212 0.212 0.173 0.216 0.157 0.137 0.216 0.118 0.235 0.235 0.235 0.216 0.216 0.216 0.192 0.000 ; 0.196 0.154 0.154 0.154 0.118 0.173 0.154 0.135 0.154 0.154 0.173 0.173 0.196 0.157 0.098 0.173 0.192 0.115 0.096 0.173 0.196 0.157 0.192 0.196 0.154 0.137 0.157 0.137 0.154 0.115 0.154 0.192 0.000 0.196 0.173 0.157 0.196 0.115 0.115 0.157 0.096 0.192 0.173 0.192 0.173 0.173 0.059 0.196 0.192 0.212 0.192 0.173 0.192 0.192 0.192 0.196 0.196 0.157 0.096 0.196 0.192 0.118 0.216 0.212 0.077 0.096 0.196 0.192 0.192 0.216 0.157 0.157 0.216 0.039 0.176 0.216 0.216 0.196 0.216 0.250 0.154 0.192 0.231 0.250 0.173 0.192 0.173 0.154 0.216 0.216 0.216 0.196 0.196 0.196 0.173 0.000 ; 0.216 0.154 0.154 0.154 0.137 0.173 0.154 0.135 0.154 0.176 0.173 0.192 0.196 0.157 0.098 0.173 0.192 0.115 0.096 0.173 0.196 0.157 0.192 0.176 0.154 0.137 0.157 0.137 0.173 0.115 0.154 0.192 0.019 0.173 0.173 0.157 0.196 0.115 0.115 0.157 0.096 0.192 0.192 0.192 0.173 0.173 0.135 0.192 0.192 0.212 0.192 0.173 0.192 0.192 0.192 0.216 0.196 0.157 0.096 0.196 0.192 0.118 0.216 0.212 0.077 0.096 0.196 0.192 0.192 0.216 0.157 0.157 0.216 0.098 0.176 0.216 0.216 0.196 0.216 0.118 0.216 0.216 0.176 0.196 0.157 0.137 0.173 0.098 0.216 0.216 0.216 0.196 0.196 0.196 0.173 0.000 ; 0.196 0.154 0.176 0.154 0.118 0.173 0.154 0.135 0.154 0.173 0.173 0.173 0.196 0.157 0.098 0.173 0.192 0.115 0.096 0.173 0.196 0.157 0.192 0.176 0.154 0.137 0.157 0.137 0.173 0.115 0.154 0.192 0.000 0.173 0.173 0.157 0.196 0.115 0.115 0.157 0.096 0.192 0.196 0.192 0.231 0.231 0.135 0.192 0.192 0.212 0.192 0.196 0.192 0.192 0.192 0.192 0.196 0.135 0.096 0.196 0.192 0.118 0.212 0.212 0.077 0.098 0.196 0.192 0.192 0.216 0.157 0.157 0.216 0.098 0.176 0.216 0.216 0.196 0.216 0.118 0.196 0.196 0.157 0.176 0.137 0.118 0.173 0.098 0.196 0.196 0.216 0.196 0.196 0.196 0.173 0.000 ; 0.255 0.212 0.235 0.212 0.176 0.235 0.212 0.192 0.212 0.231 0.231 0.231 0.255 0.235 0.157 0.235 0.250 0.173 0.154 0.212 0.255 0.216 0.255 0.235 0.212 0.196 0.216 0.196 0.231 0.173 0.212 0.250 0.038 0.231 0.231 0.216 0.255 0.173 0.173 0.216 0.154 0.250 0.250 0.250 0.115 0.115 0.192 0.250 0.275 0.269 0.250 0.231 0.250 0.250 0.255 0.250 0.255 0.192 0.154 0.255 0.250 0.176 0.269 0.269 0.135 0.157 0.255 0.250 0.250 0.275 0.216 0.216 0.275 0.157 0.235 0.275 0.275 0.255 0.275 0.176 0.275 0.275 0.235 0.255 0.216 0.196 0.231 0.157 0.275 0.275 0.275 0.255 0.255 0.255 0.212 0.000 ; 0.157 0.096 0.118 0.098 0.078 0.115 0.096 0.077 0.096 0.115 0.115 0.096 0.137 0.096 0.118 0.098 0.137 0.058 0.000 0.096 0.137 0.098 0.137 0.118 0.098 0.078 0.077 0.078 0.115 0.058 0.096 0.135 -0.039 0.115 0.135 0.115 0.078 0.058 0.058 0.118 0.020 0.135 0.135 0.135 0.192 0.192 0.077 0.135 0.154 0.154 0.135 0.115 0.135 0.135 0.115 0.058 0.137 0.077 0.038 0.137 0.115 0.000 0.154 0.154 -0.019 -0.020 0.118 0.135 0.115 0.157 0.098 0.098 0.157 0.039 0.078 0.157 0.157 0.118 0.157 0.058 0.137 0.137 0.098 0.039 0.078 0.059 0.096 -0.020 0.137 0.137 0.157 0.137 0.137 0.059 0.096 0.000 ; 0.235 0.173 0.196 0.157 0.137 0.192 0.173 0.154 0.173 0.192 0.196 0.192 0.216 0.173 0.118 0.176 0.212 0.137 0.115 0.173 0.216 0.176 0.216 0.196 0.173 0.157 0.154 0.157 0.192 0.135 0.173 0.216 0.020 0.192 0.216 0.192 0.216 0.135 0.137 0.192 0.096 0.212 0.212 0.212 0.231 0.192 0.154 0.212 0.231 0.231 0.212 0.192 0.212 0.212 0.192 0.212 0.216 0.154 0.115 0.196 0.212 0.118 0.231 0.231 0.096 0.118 0.216 0.235 0.192 0.235 0.176 0.176 0.235 0.118 0.196 0.235 0.235 0.216 0.231 0.135 0.231 0.231 0.192 0.212 0.173 0.135 0.192 0.118 0.235 0.235 0.235 0.216 0.216 0.216 0.173 0.000 ; 0.235 0.192 0.216 0.176 0.157 0.212 0.192 0.173 0.192 0.212 0.192 0.212 0.235 0.192 0.135 0.196 0.231 0.154 0.137 0.192 0.235 0.196 0.235 0.216 0.173 0.176 0.173 0.176 0.212 0.154 0.212 0.231 0.019 0.212 0.212 0.212 0.231 0.176 0.154 0.212 0.135 0.231 0.231 0.231 0.231 0.192 0.173 0.231 0.250 0.250 0.231 0.212 0.231 0.231 0.212 0.231 0.235 0.173 0.135 0.216 0.231 0.154 0.250 0.250 0.115 0.118 0.235 0.255 0.212 0.255 0.196 0.196 0.255 0.137 0.216 0.250 0.250 0.231 0.250 0.154 0.250 0.250 0.212 0.231 0.173 0.154 0.212 0.135 0.250 0.255 0.255 0.235 0.235 0.235 0.192 0.000 ; 0.176 0.135 0.157 0.118 0.118 0.154 0.135 0.115 0.135 0.154 0.135 0.096 0.176 0.096 0.077 0.098 0.154 0.096 0.058 0.096 0.176 0.098 0.176 0.157 0.115 0.115 0.077 0.118 0.115 0.077 0.115 0.173 -0.020 0.154 0.154 0.154 0.173 0.118 0.096 0.115 0.078 0.173 0.154 0.154 0.115 0.020 0.115 0.173 0.192 0.192 0.176 0.154 0.173 0.173 0.154 0.173 0.176 0.115 0.078 0.157 0.157 0.096 0.192 0.192 0.058 0.078 0.176 0.196 0.135 0.196 0.137 0.137 0.196 0.020 0.137 0.192 0.192 0.154 0.192 0.250 0.154 0.192 0.231 0.250 0.196 0.135 0.154 0.096 0.192 0.192 0.192 0.173 0.173 0.173 0.096 0.000 ; 0.216 0.196 0.196 0.157 0.135 0.192 0.173 0.154 0.173 0.192 0.173 0.192 0.216 0.173 0.115 0.176 0.212 0.135 0.115 0.173 0.216 0.196 0.216 0.196 0.154 0.154 0.154 0.154 0.192 0.135 0.192 0.212 0.020 0.192 0.192 0.192 0.212 0.157 0.135 0.192 0.118 0.212 0.212 0.212 0.154 0.196 0.154 0.212 0.231 0.231 0.192 0.192 0.212 0.212 0.192 0.212 0.216 0.154 0.118 0.196 0.216 0.135 0.231 0.231 0.118 0.118 0.216 0.235 0.192 0.235 0.176 0.176 0.235 0.118 0.196 0.231 0.231 0.216 0.235 0.137 0.235 0.235 0.196 0.216 0.176 0.157 0.192 0.098 0.235 0.235 0.235 0.216 0.216 0.216 0.173 0.000 ; 0.255 0.216 0.235 0.196 0.173 0.231 0.235 0.216 0.235 0.231 0.212 0.231 0.255 0.212 0.154 0.216 0.250 0.173 0.154 0.212 0.255 0.216 0.255 0.235 0.192 0.192 0.192 0.192 0.231 0.173 0.231 0.250 0.038 0.235 0.231 0.231 0.250 0.196 0.173 0.231 0.157 0.197 0.250 0.250 0.121 0.235 0.192 0.197 0.269 0.269 0.231 0.231 0.250 0.250 0.231 0.250 0.255 0.192 0.157 0.235 0.255 0.173 0.269 0.269 0.157 0.157 0.255 0.275 0.231 0.269 0.216 0.216 0.275 0.157 0.231 0.269 0.269 0.255 0.197 0.136 0.197 0.197 0.182 0.231 0.212 0.192 0.250 0.154 0.269 0.269 0.269 0.250 0.255 0.250 0.212 0.000 ; 0.235 0.196 0.216 0.176 0.154 0.212 0.192 0.173 0.212 0.212 0.192 0.212 0.235 0.192 0.135 0.196 0.231 0.154 0.135 0.192 0.235 0.196 0.235 0.216 0.173 0.173 0.173 0.173 0.212 0.154 0.212 0.231 0.019 0.216 0.212 0.212 0.231 0.176 0.154 0.212 0.137 0.231 0.231 0.231 0.135 0.216 0.157 0.231 0.250 0.231 0.212 0.212 0.231 0.231 0.212 0.255 0.235 0.173 0.137 0.216 0.235 0.154 0.235 0.250 0.137 0.137 0.235 0.255 0.212 0.212 0.196 0.176 0.255 0.135 0.212 0.250 0.192 0.235 0.182 0.154 0.250 0.250 0.212 0.231 0.192 0.157 0.231 0.137 0.235 0.235 0.235 0.216 0.216 0.216 0.192 0.000 ; ...
              0.137 0.098 0.098 0.098 0.058 0.115 0.115 0.077 0.115 0.115 0.096 0.000 0.135 0.077 0.135 0.078 0.115 0.058 0.000 0.077 0.137 0.098 0.137 0.154 0.096 0.077 0.058 0.077 0.096 0.058 0.096 0.135 -0.059 0.118 0.135 0.115 0.135 0.078 0.058 0.096 0.000 0.135 0.118 0.115 0.038 0.118 0.059 0.135 0.154 0.135 0.135 0.137 0.135 0.135 0.137 0.157 0.137 0.077 0.039 0.118 0.118 0.058 0.137 0.154 0.039 0.059 0.118 0.157 0.096 0.212 0.098 0.078 0.157 0.059 0.096 0.154 0.192 0.118 0.106 0.038 0.157 0.157 0.115 0.135 0.096 0.096 0.115 0.077 0.154 0.154 0.154 0.135 0.135 0.135 0.077 0.000 ; 0.176 0.137 0.137 0.118 0.096 0.154 0.135 0.115 0.154 0.154 0.135 0.020 0.154 0.135 0.077 0.154 0.176 0.096 0.077 0.154 0.176 0.135 0.154 0.173 0.115 0.115 0.115 0.115 0.137 0.096 0.135 0.173 -0.020 0.157 0.154 0.154 0.173 0.118 0.096 0.135 0.078 0.173 0.173 0.157 0.077 0.157 0.098 0.173 0.192 0.173 0.154 0.173 0.154 0.173 0.173 0.039 0.173 0.115 0.078 0.154 0.157 0.019 0.176 0.192 -0.020 -0.039 0.176 0.196 0.135 0.173 0.137 0.118 0.196 0.078 0.135 0.192 0.154 0.157 0.136 0.077 0.196 0.173 0.135 0.154 0.115 0.096 0.154 -0.020 0.196 0.196 0.196 0.176 0.176 0.020 0.154 0.000 ; 0.216 0.176 0.176 0.173 0.135 0.192 0.173 0.154 0.192 0.192 0.192 0.216 0.192 0.192 0.115 0.192 0.216 0.135 0.115 0.192 0.216 0.173 0.196 0.212 0.154 0.154 0.176 0.154 0.196 0.135 0.192 0.192 0.020 0.196 0.192 0.192 0.212 0.157 0.135 0.196 0.118 0.212 0.212 0.196 0.115 0.196 0.137 0.212 0.231 0.212 0.192 0.212 0.212 0.212 0.212 0.235 0.212 0.176 0.118 0.212 0.216 0.137 0.216 0.231 0.118 0.115 0.212 0.235 0.216 0.192 0.176 0.157 0.235 0.118 0.192 0.231 0.173 0.216 0.167 0.135 0.231 0.212 0.192 0.212 0.173 0.154 0.212 0.115 0.231 0.231 0.231 0.212 0.212 0.216 0.192 0.000 ; 0.235 0.196 0.196 0.192 0.154 0.212 0.192 0.173 0.212 0.212 0.212 0.235 0.212 0.212 0.135 0.212 0.235 0.154 0.135 0.212 0.235 0.192 0.216 0.231 0.098 0.173 0.192 0.173 0.216 0.154 0.212 0.231 0.019 0.216 0.212 0.212 0.231 0.176 0.154 0.216 0.137 0.231 0.231 0.216 0.135 0.216 0.157 0.231 0.250 0.231 0.212 0.231 0.231 0.231 0.231 0.255 0.231 0.176 0.137 0.231 0.235 0.157 0.235 0.255 0.137 0.135 0.096 0.255 0.235 0.135 0.196 0.176 0.255 0.078 0.212 0.250 0.135 0.235 0.231 0.154 0.250 0.231 0.212 0.212 0.173 0.154 0.231 0.135 0.231 0.231 0.231 0.212 0.235 0.235 0.212 0.000 ; 0.196 0.137 0.137 0.135 0.096 0.154 0.135 0.115 0.154 0.135 0.135 0.078 0.154 0.154 0.038 0.154 0.176 0.077 0.058 0.154 0.176 0.135 0.157 0.173 0.115 0.096 0.135 0.115 0.157 0.096 0.135 0.154 -0.020 0.157 0.154 0.135 0.173 0.118 0.096 0.137 0.059 0.173 0.154 0.137 0.058 0.176 0.098 0.173 0.192 0.173 0.216 0.154 0.173 0.154 0.154 0.216 0.154 0.118 0.078 0.154 0.137 0.098 0.157 0.196 0.059 0.077 0.192 0.196 0.157 0.115 0.137 0.118 0.196 0.078 0.154 0.173 0.196 0.157 0.173 0.059 0.192 0.173 0.135 0.135 0.137 0.118 0.154 0.137 0.098 0.176 0.196 0.176 0.176 0.176 0.154 0.000 ; 0.192 0.137 0.137 0.135 0.096 0.154 0.135 0.115 0.154 0.154 0.154 0.059 0.154 0.135 0.038 0.135 0.157 0.098 0.039 0.135 0.176 0.115 0.157 0.173 0.115 0.115 0.115 0.115 0.137 0.096 0.154 0.173 -0.020 0.157 0.154 0.157 0.173 0.118 0.096 0.157 0.039 0.173 0.173 0.137 0.077 0.157 0.098 0.173 0.192 0.173 0.176 0.173 0.173 0.173 0.173 0.176 0.173 0.118 0.078 0.173 0.157 0.098 0.157 0.196 0.078 0.077 0.173 0.196 0.137 0.096 0.137 0.118 0.196 0.078 0.135 0.192 0.192 0.157 0.192 0.059 0.192 0.173 0.154 0.154 0.118 0.115 0.154 0.096 0.192 0.192 0.192 0.173 0.173 0.176 0.135 0.000 ; 0.173 0.154 0.137 0.135 0.115 0.154 0.154 0.135 0.154 0.154 0.135 0.173 0.154 0.154 0.058 0.135 0.176 0.098 0.077 0.154 0.176 0.135 0.157 0.176 0.137 0.118 0.137 0.137 0.157 0.098 0.157 0.176 -0.020 0.154 0.136 0.121 0.136 0.091 0.076 0.121 0.061 0.136 0.176 0.173 0.077 0.154 0.115 0.173 0.192 0.192 0.173 0.173 0.176 0.176 0.176 0.196 0.176 0.137 0.077 0.176 0.176 0.098 0.196 0.196 0.078 0.078 0.176 0.196 0.176 0.173 0.115 0.115 0.173 0.077 0.135 0.173 0.173 0.154 0.192 0.096 0.192 0.192 0.154 0.173 0.135 0.115 0.173 0.096 0.192 0.192 0.196 0.173 0.173 0.176 0.154 0.000 ; 0.137 0.118 0.098 0.096 0.077 0.115 0.115 0.077 0.115 0.115 0.115 0.059 0.135 0.077 0.038 0.077 0.118 0.059 0.039 0.096 0.137 0.058 0.118 0.135 0.077 0.077 0.058 0.058 0.078 0.058 0.077 0.135 -0.059 0.137 0.137 0.115 0.135 0.059 0.058 0.078 0.038 0.135 0.096 0.098 0.212 0.020 0.176 0.135 0.154 0.154 0.137 0.135 0.135 0.135 0.135 0.157 0.135 0.078 0.038 0.135 0.118 0.059 0.137 0.157 0.039 0.038 0.135 0.154 0.118 0.157 0.098 0.078 0.157 -0.020 0.098 0.154 0.154 0.118 0.154 -0.020 0.154 0.157 0.115 0.115 0.078 0.077 0.115 0.039 0.157 0.157 0.157 0.137 0.137 0.137 0.096 0.000 ; 0.216 0.176 0.196 0.173 0.154 0.192 0.196 0.154 0.196 0.192 0.192 0.216 0.192 0.192 0.096 0.192 0.216 0.137 0.118 0.192 0.216 0.173 0.196 0.212 0.154 0.154 0.173 0.176 0.196 0.135 0.192 0.212 0.020 0.192 0.216 0.192 0.212 0.137 0.135 0.196 0.096 0.216 0.212 0.196 0.212 0.196 0.176 0.192 0.231 0.231 0.216 0.212 0.212 0.212 0.212 0.235 0.212 0.157 0.115 0.212 0.216 0.118 0.216 0.235 0.098 0.096 0.250 0.231 0.216 0.235 0.176 0.157 0.235 0.118 0.196 0.235 0.231 0.216 0.231 0.118 0.231 0.235 0.192 0.192 0.196 0.135 0.212 0.098 0.235 0.235 0.235 0.196 0.216 0.216 0.192 0.000 ; 0.255 0.216 0.212 0.212 0.192 0.231 0.231 0.192 0.212 0.235 0.231 0.255 0.231 0.231 0.154 0.231 0.255 0.176 0.157 0.231 0.255 0.212 0.235 0.250 0.192 0.192 0.212 0.212 0.235 0.173 0.231 0.250 0.038 0.231 0.255 0.231 0.250 0.176 0.173 0.235 0.154 0.255 0.250 0.255 0.212 0.235 0.176 0.231 0.269 0.269 0.255 0.250 0.250 0.250 0.250 0.250 0.250 0.196 0.154 0.250 0.255 0.176 0.275 0.275 0.157 0.154 0.250 0.269 0.255 0.275 0.216 0.196 0.275 0.157 0.235 0.275 0.269 0.255 0.269 0.157 0.269 0.275 0.231 0.231 0.216 0.192 0.250 0.157 0.275 0.275 0.275 0.255 0.255 0.255 0.231 0.000 ; 0.255 0.216 0.212 0.212 0.192 0.231 0.231 0.192 0.212 0.212 0.231 0.231 0.255 0.235 0.157 0.231 0.255 0.176 0.157 0.231 0.255 0.212 0.235 0.250 0.192 0.192 0.212 0.212 0.235 0.173 0.231 0.250 0.038 0.231 0.255 0.231 0.250 0.176 0.173 0.235 0.154 0.255 0.255 0.255 0.192 0.235 0.196 0.231 0.269 0.269 0.255 0.255 0.250 0.250 0.250 0.250 0.250 0.196 0.154 0.231 0.255 0.173 0.275 0.275 0.157 0.157 0.250 0.269 0.255 0.275 0.216 0.216 0.275 0.157 0.235 0.275 0.269 0.255 0.269 0.157 0.269 0.275 0.231 0.231 0.216 0.192 0.250 0.157 0.275 0.275 0.275 0.235 0.250 0.250 0.231 0.000 ; 0.235 0.196 0.192 0.173 0.173 0.212 0.212 0.173 0.192 0.192 0.196 0.212 0.231 0.212 0.115 0.212 0.235 0.157 0.137 0.212 0.235 0.212 0.212 0.231 0.173 0.173 0.192 0.192 0.196 0.154 0.212 0.231 0.019 0.212 0.235 0.212 0.231 0.157 0.154 0.216 0.135 0.235 0.212 0.235 0.137 0.216 0.173 0.212 0.250 0.250 0.235 0.235 0.231 0.231 0.231 0.231 0.231 0.176 0.135 0.212 0.235 0.154 0.255 0.255 0.137 0.135 0.231 0.250 0.235 0.255 0.196 0.196 0.255 0.137 0.216 0.255 0.250 0.235 0.250 0.137 0.250 0.255 0.212 0.235 0.196 0.173 0.231 0.137 0.235 0.255 0.235 0.216 0.216 0.216 0.212 0.000 ; 0.235 0.196 0.192 0.173 0.173 0.212 0.212 0.173 0.192 0.192 0.192 0.212 0.231 0.212 0.115 0.212 0.216 0.157 0.137 0.212 0.235 0.212 0.212 0.231 0.173 0.173 0.196 0.192 0.192 0.154 0.212 0.235 0.019 0.212 0.235 0.212 0.231 0.157 0.154 0.192 0.135 0.235 0.212 0.235 0.137 0.216 0.173 0.212 0.250 0.250 0.235 0.235 0.231 0.231 0.231 0.231 0.235 0.176 0.135 0.212 0.235 0.154 0.255 0.255 0.137 0.135 0.235 0.250 0.235 0.255 0.196 0.196 0.255 0.137 0.216 0.255 0.250 0.235 0.255 0.137 0.250 0.255 0.212 0.235 0.196 0.173 0.231 0.137 0.255 0.255 0.250 0.216 0.235 0.235 0.212 0.000 ; 0.235 0.196 0.192 0.173 0.173 0.212 0.212 0.173 0.192 0.192 0.192 0.212 0.231 0.212 0.115 0.212 0.235 0.157 0.137 0.212 0.235 0.212 0.212 0.231 0.173 0.173 0.192 0.192 0.192 0.154 0.216 0.235 0.019 0.212 0.235 0.212 0.231 0.157 0.154 0.192 0.135 0.235 0.212 0.235 0.137 0.216 0.173 0.212 0.250 0.250 0.235 0.235 0.231 0.231 0.231 0.231 0.235 0.196 0.135 0.212 0.212 0.154 0.255 0.255 0.137 0.135 0.212 0.250 0.235 0.255 0.196 0.196 0.255 0.137 0.216 0.255 0.255 0.235 0.250 0.137 0.250 0.255 0.216 0.235 0.196 0.176 0.231 0.137 0.255 0.255 0.250 0.216 0.235 0.235 0.212 0.000 ; 0.231 0.173 0.192 0.173 0.135 0.192 0.173 0.173 0.173 0.192 0.176 0.196 0.192 0.196 0.115 0.192 0.212 0.137 0.115 0.192 0.212 0.192 0.212 0.212 0.173 0.154 0.173 0.176 0.196 0.137 0.196 0.216 0.216 0.098 0.216 0.196 0.212 0.137 0.118 0.176 0.098 0.196 0.196 0.196 0.098 0.176 0.157 0.196 0.216 0.216 0.196 0.196 0.196 0.196 0.196 0.231 0.212 0.173 0.115 0.212 0.212 0.135 0.231 0.235 0.096 0.115 0.192 0.212 0.192 0.212 0.154 0.154 0.212 0.115 0.173 0.212 0.212 0.192 0.212 0.115 0.212 0.231 0.192 0.212 0.173 0.154 0.192 0.135 0.231 0.231 0.231 0.212 0.212 0.212 0.192 0.000 ; 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
        
    end % constant properties
    
end