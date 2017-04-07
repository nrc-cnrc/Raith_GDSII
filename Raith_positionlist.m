classdef Raith_positionlist < handle
%
% obj=Raith_positionlist(library,csf_path,WF,chipUV)
%
% Raith_positionlist objects define .pls files (positionlists) to be used 
% with an associated Raith GDSII hierarchy created by Raith_library
% objects, for Raith beamwriting tools.
%
%
% Arguments:
%
%	library - Raith_library object containing all structures to be placed 
%       in positionposlist
%   csf_path - full path of GDSII hierarchy file, as found on Raith
%       computer
%   WF - writefield size [size_u size_v] (um)
%   chipUV - size of rectangular chip [size_u size_v] (mm)
%
%
% Properties:
%
%   library - Raith_library object containing all structures
%   csf_path - full path of GDSII hierarchy file, as found on Raith
%       computer
%   WF - writefield size [size_u size_v] (um)
%   chipUV - size of rectangular chip [size_u size_v] (mm)
%   poslist  - structure array of positionlist entries; has fields name, 
%       uv_c, DF, WA, and layers
%
%
% Methods:
%
%   append(name,uv_c,DF,WA,[layers]) - append Raith_structure object to positionlist
%       name - Raith_structure object name (as found in GDSII hierarchy 
%           file)
%       uv_c - centre of first writefield of structure:  [u_c v_c] (mm)
%       DF - overall dose factor scaling for entire structure
%       WA - working area of structure:  [u_min v_min u_max v_max] (um)
%       layers - vector of layers to expose [optional]; defaults to all
%           layers present in structure
%
%   plot - plot all structures in positionlist with Raith dose factor
%       colouring (filled polygons where applicable)
%
%   plotedges - plot all structures in positionlist with Raith dose factor
%       colouring (edges of polygons where applicable)
%
%   plotWA - plot working area of all structures in positionlist in dotted 
%       blue lines
%
%   plotWF - plot writefields of all structures in positionlist in dotted 
%       green lines; writefield centres are marked with a '+'
%
%   centre([mbyn]) - centre current positionlist entries on the chip, with 
%       the option of matrix-copying them.  If called with no argument, the 
%       current positionlist entries are shifted such that the overall 
%       pattern (as defined by the working areas) are centred both
%       vertically and horizontally on the chip.  If called with an
%       argument, the chip is divided into a matrix of equal-sized sub-chips,
%       and the positionlist entries are centred as described above in each
%       sub-chip.
%           mbyn - vector [m n] of the number of rows and columns of the
%               matrix of sub-chips [optional]
%
%   shift(uv_sh) - shift current positionlist entries on the chip.
%       uv_sh - vector [u_sh v_sh] specifying the shift of the overall 
%           positionlist with respect to their current positions (mm)
% 
%   writepls([filepath]) - write positionlist
%       filepath - full path of positionlist file, including .pls
%           extension [optional]; if none is specified, a 
%           [library name].pls file is written to the current directory.
%
%
% Dependencies:  Raith_element.m, Raith_structure.m, Raith_library.m
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
        library
        csf_path
        WF
        chipUV
        poslist 
    end 
    
	properties (Hidden)
        hChip=-1; % Handle for chip outlines being plotted
    end % set access private properties
   
    
    methods
    
        function obj=Raith_positionlist(library,csf_path,WF,chipUV)
            if nargin>0 
                obj.library=library;
                obj.csf_path=csf_path;
                obj.WF=WF;
                obj.chipUV=chipUV;
            end
        end % Constructor
        
        
        function set.library(obj,library)
            if ~isa(library,'Raith_library')
                error('Raith_positionlist:  library must be a Raith_library object.')
            end
            obj.library=library;
        end % set.library
        
        
        function set.csf_path(obj,csf_path)
            if ~ischar(csf_path)
                error('Raith_positionlist:  csf_path must be a string.')
            end
            if ~strcmp(csf_path((end-3):end),'.csf')
                error('Raith_positionlist:  csf_path must include .csf extension.')
            end
            obj.csf_path=csf_path;
        end % set.csf_path
        
        
        function set.WF(obj,WF)
            if ~isnumeric(WF) || ~isvector(WF) || length(WF)~=2
                error('Raith_positionlist:  WF must be of form [size_u size_v].')
            end
            obj.WF=WF;
        end % set.WF
        
        
        function set.chipUV(obj,chipUV)
            if ~isnumeric(chipUV) || ~isvector(chipUV) || length(chipUV)~=2
                error('Raith_positionlist:  chipUV must be of form [size_u size_v].')
            end
            obj.chipUV=chipUV;
        end % set.chipUV

        
        function append(obj,structname,uv_c,DF,WA,varargin)
        %
        %   Raith_positionlist.append(structname,uv_c,DF,WA,[layers]) - append 
        %       Raith_structure object to positionlist
        %
        %   Arguments:
        %
        %       structname - Raith_structure object name (as found in GDSII 
        %           hierarchy file)
        %       uv_c - centre of first writefield of structure, [u_c v_c] (mm)
        %       DF - overall dose factor scaling for entire structure
        %       WA - working area of structure: [u_min v_min u_max v_max] (um)            
        %       layers - vector of layers to expose [optional]; defaults to all
        %           layers present in structure
        %
        
        if nargin==6
            layers=varargin{1};
        elseif nargin==5 % Determine all layers present in structure
            layers=obj.alllayers(structname,[]);
        else
            error('Raith_positionlist:  too many arguments.')
        end
        
        % Append entry to positionlist
            obj.poslist(end+1).name=structname;
            obj.poslist(end).uv_c=uv_c;
            obj.poslist(end).DF=DF;
            obj.poslist(end).WA=WA;
            obj.poslist(end).layers=layers;
        end % append
        
        
        function plot(obj)
        %   Raith_positionlist.plot - plot all structures in positionlist with 
        %       Raith dose factor colouring (filled polygons where applicable)
        %            
            obj.plotchip;
            
            for k=1:length(obj.poslist)
                uv_0=obj.poslist(k).uv_c*1000-obj.WF/2-obj.poslist(k).WA(1:2);  % Origin of structure in positionlist (um); not to be confused with centre of first WF
                obj.library.plot(obj.poslist(k).name,obj.library.trans(uv_0),obj.poslist(k).DF);  
            end
            
        end % plot
        
        
        function plotedges(obj)
        %   Raith_positionlist.plotedges - plot all structures in positionlist with 
        %       Raith dose factor colouring (edges of polygons where applicable)
        %            
            obj.plotchip;
        
            for k=1:length(obj.poslist)
                uv_0=obj.poslist(k).uv_c*1000-obj.WF/2-obj.poslist(k).WA(1:2);  % Origin of structure in positionlist (um); not to be confused with centre of first WF
                obj.library.plotedges(obj.poslist(k).name,obj.library.trans(uv_0),obj.poslist(k).DF);  
            end
            
        end % plot

        
        function plotWA(obj)
        %   Raith_positionlist.plotWA - plot working area of all structures in 
        %       positionlist as blue dotted line
        %
            obj.plotchip;
            
            for k=1:length(obj.poslist)
                uv_0=obj.poslist(k).uv_c*1000-obj.WF/2-obj.poslist(k).WA(1:2);  % Origin of structure in positionlist (um); not to be confused with centre of first WF
                hold on
                plot(obj.poslist(k).WA([1 3 3 1 1])+uv_0(1),obj.poslist(k).WA([2 2 4 4 2])+uv_0(2),'LineStyle',':','Color',[39 170 255]/255);
            end
            
        end % plotWA
        
        
        function plotWF(obj)
        %   Raith_positionlist.plotWF - plot first writefield of all structures in 
        %       positionlist as green dotted line; writefield centre is marked with a '+'.
        %
            obj.plotchip;
            
            for k=1:length(obj.poslist)
                hold on;
                plot(obj.WF(1)/2*[-1 1 1 -1 -1]+obj.poslist(k).uv_c(1)*1000,obj.WF(2)/2*[-1 -1 1 1 -1]+obj.poslist(k).uv_c(2)*1000,'LineStyle',':','Color',[0 166 81]/255);
                plot(obj.poslist(k).uv_c(1)*1000,obj.poslist(k).uv_c(2)*1000,'Marker','+','Color',[0 166 81]/255);
            end
            
        end % plotWF
        
        
        function centre(obj,varargin)
        %   Raith_positionlist.centre([mbyn]) - centre current positionlist
        %   entries on the chip, with the option of matrix-copying them.
        %   If called with no argument, the current positionlist entries
        %   are shifted such that the overall pattern (as defined by the 
        %   working areas) are centred both vertically and horizontally on
        %   the chip.
        %
        %   Argument:
        %
        %   mbyn - vector [m n] of the number of rows and columns of the
        %       matrix of instances of the current positionlist to be
        %       copied [optional]
        %
        
        obj.checkposlist;  % Check that all poslist entries are of the correct format
        
        if nargin==1
            mbyn=[1 1];
        elseif nargin==2
            mbyn=varargin{1};
        else
            error('Raith_positionlist.centre:  too many arguments.');
        end
            
        pltemp=obj.poslist;  % Current poslist
        
        % Find boundaries of structure working areas (in mm)
        minU=Inf;
        minV=Inf;
        maxU=-Inf;
        maxV=-Inf;
        for k=1:length(pltemp)
            minU=min(minU,pltemp(k).WA(1)/1000+pltemp(k).uv_c(1));
            minV=min(minV,pltemp(k).WA(2)/1000+pltemp(k).uv_c(2));
            maxU=max(maxU,pltemp(k).WA(3)/1000+pltemp(k).uv_c(1));
            maxV=max(maxV,pltemp(k).WA(4)/1000+pltemp(k).uv_c(2));
        end
        wU=maxU-minU;  % Width of pattern in U (mm)
        wV=maxV-minV;  % Width of pattern in V (mm)
        
        % Centre pattern of pltemp entries about origin (preserving spacing)
        for k=1:length(pltemp)
            pltemp(k).uv_c=pltemp(k).uv_c-[minU+wU/2 minV+wV/2];
        end
        
        % Construct centres of instances
        uc=(1:2:(2*mbyn(2)))/(2*mbyn(2))*obj.chipUV(1);  % Centres of instances in U (mm)
        vc=(1:2:(2*mbyn(1)))/(2*mbyn(1))*obj.chipUV(2);  % Centres of instances in V (mm)
        [UC,VC]=meshgrid(uc,vc);
        UC=reshape(UC,prod(mbyn),1);
        VC=reshape(VC,prod(mbyn),1);
        
        obj.poslist=[];  % Clear current poslist
        
        % Reload positionlist
        for kC=1:length(UC)
            for k=1:length(pltemp)
               obj.append(pltemp(k).name,pltemp(k).uv_c+[UC(kC) VC(kC)],pltemp(k).DF,pltemp(k).WA,pltemp(k).layers);
            end
        end
        
        end % centre
        
        
        function shift(obj,uv_sh)
        %   Raith_positionlist.shift(uv_sh) - shift current positionlist
        %   entries on the chip.
        %
        %   Argument:
        %
        %   uv_sh - vector [u_sh v_sh] specifying the relative shift of the
        %       overall positionlist with respect to their current positions (mm)
        %
        
        obj.checkposlist;  % Check that all poslist entries are of the correct format
                   
        pltemp=obj.poslist;  % Current poslist
        
        obj.poslist=[];  % Clear current poslist
        
        % Reload positionlist
        for k=1:length(pltemp)
           obj.append(pltemp(k).name,pltemp(k).uv_c+uv_sh,pltemp(k).DF,pltemp(k).WA,pltemp(k).layers);
        end
        
        end % shift
        
        
        function writepls(obj,varargin)
        %   Raith_positionlist.writepls([filepath]) - write positionlist
        %
        %   Argument:
        %
        %       filepath - full path of positionlist file, including .pls
        %           extension [optional]; if none is specified, a 
        %           [library name].pls file is written to the current 
        %           directory
        %
        
        obj.checkposlist;  % Check that all poslist entries are of the correct format

            if isunix
                slsh='/';
            else
                slsh='\';
            end
        
            if nargin==1 % No filepath given
                filepath=[pwd slsh obj.library.name '.pls'];
            elseif nargin==2
                filepath=varargin{1};
            else
                error('Raith_positionlist.writepls:  too many arguments.');
            end
            
            if ~strcmp(filepath((end-3):end),'.pls')
                error('Raith_positionlist:  filepath must include .pls extension.')
            end
            
            disp(' ');
            disp(['Writing ' filepath '...']);
        
            fid=fopen(filepath,'w');
 
            % Write header
            headline=sprintf('\n[HEADER]\nFORMAT=IXYZRTUVWATC,Options,0,Type,0,Size-U,0,Size-V,0,Points-U,0,Points-V,0,Dir,0,Avg,0,Pos1,0,Pos2,0,Pos3,0,Link,0,File,0,Layer,0,Area,0,DoseFactor,0,Dwelltime,0,Stepsize,0,SplDwell,0,SplStep,0,CurveStep,0,CurveDwell,0,DotDwell,0,FBMSArea,0,FBMSLines,0,SplDot,0,Time,0,Timestamp,0,Method,0,Dot,0,StepsizeU,0,StepsizeV,0,CurveLine,0\nWAFERLAYOUT=DEFAULT.wlo\nLotID=\nWaferID=\nSlot=\nMinimizeWin=FALSE\n');
            colline=sprintf('[COLUMNS]\nNo.=W:25,!VISIBLE,!SHOWDIM\nID=W:25,VISIBLE,!SHOWDIM\nX=W:50,!VISIBLE,SHOWDIM\nY=W:50,!VISIBLE,SHOWDIM\nZ=W:50,!VISIBLE,SHOWDIM\nR=W:50,!VISIBLE,SHOWDIM\nT=W:50,!VISIBLE,SHOWDIM\nU=W:50,VISIBLE,SHOWDIM\nV=W:50,VISIBLE,SHOWDIM\nW=W:50,!VISIBLE,SHOWDIM\nAttribute=W:55,VISIBLE,DEFAULT:A,!SHOWDIM\nTemplate=W:55,VISIBLE,DEFAULT:UV,!SHOWDIM\nComment=W:100,VISIBLE,!SHOWDIM\nOptions=W:165,VISIBLE,!SHOWDIM\nType=W:85,VISIBLE,!SHOWDIM\nSize-U=W:55,!VISIBLE,DIM:um,SHOWDIM\nSize-V=W:55,!VISIBLE,DIM:um,SHOWDIM\nPoints-U=W:55,!VISIBLE,DIM:px,SHOWDIM\nPoints-V=W:55,!VISIBLE,DIM:px,SHOWDIM\nDir=W:15,!VISIBLE,!SHOWDIM\nAvg=W:25,!VISIBLE,!SHOWDIM\nPos1=W:85,VISIBLE,DIM:um,SHOWDIM\nPos2=W:85,VISIBLE,DIM:um,SHOWDIM\nPos3=W:85,VISIBLE,DIM:um,SHOWDIM\nLink=W:25,VISIBLE,!SHOWDIM\nFile=W:160,VISIBLE,!SHOWDIM\nLayer=W:80,VISIBLE,!SHOWDIM\nArea=W:245,!VISIBLE,!SHOWDIM\nDoseFactor=W:55,VISIBLE,!SHOWDIM\nDwelltime=W:55,!VISIBLE,DIM:ms,SHOWDIM\nStepsize=W:55,!VISIBLE,DIM:um,SHOWDIM\nSplDwell=W:55,!VISIBLE,DIM:ms,SHOWDIM\nSplStep=W:55,!VISIBLE,DIM:um,SHOWDIM\nCurveStep=W:69,!VISIBLE,DIM:um,SHOWDIM\nCurveDwell=W:55,!VISIBLE,DIM:ms,SHOWDIM\nDotDwell=W:55,!VISIBLE,DIM:ms,SHOWDIM\nFBMSArea=W:88,VISIBLE,DIM:mm/s,SHOWDIM\nFBMSLines=W:81,VISIBLE,DIM:mm/s,SHOWDIM\nSplDot=W:10,!VISIBLE,!SHOWDIM\nTime=W:85,VISIBLE,!SHOWDIM\nTimestamp=W:85,!VISIBLE,!SHOWDIM\nMethod=W:85,!VISIBLE,!SHOWDIM\nDot=W:20,!VISIBLE,!SHOWDIM\nStepsizeU=W:50,VISIBLE,!SHOWDIM\nStepsizeV=W:50,VISIBLE,!SHOWDIM\nCurveLine=W:50,VISIBLE,!SHOWDIM\n\n[DATA]');
            fprintf(fid,'%s\r\n',[headline colline]);
            disp('     Header information');
            
            % Write positionlist entries
                            
            for k=1:length(obj.poslist)
                % Construct string for layers to expose:  multiple values are separated by '/'
                layerstr=regexprep(num2str(obj.poslist(k).layers),'\s+','/');
                plsline=sprintf('%d,0.000000,0.000000,0.000000,0.000000,0.000000,%.6f,%.6f,0.000000,XN,UV,%s,,EXPOSURE,%.3f,%.3f,,,,,%.3f,%.3f,,,%s,%s,%.3f;%.3f;%.3f;%.3f,%.3f,,,,,,,,,,,,,,',k-1,obj.poslist(k).uv_c(1),obj.poslist(k).uv_c(2),obj.poslist(k).name,obj.poslist(k).WA(3)-obj.poslist(k).WA(1),obj.poslist(k).WA(4)-obj.poslist(k).WA(2),obj.WF(1)/2,obj.WF(2)/2,obj.csf_path,layerstr,obj.poslist(k).WA(1),obj.poslist(k).WA(2),obj.poslist(k).WA(3),obj.poslist(k).WA(4),obj.poslist(k).DF);
                fprintf(fid,'%s\r\n',plsline);
                if ~isscalar(obj.poslist(k).layers) % Formatting for plural
                    layerstr=['s ' layerstr];
                else
                    layerstr=[' ' layerstr];
                end
                fprintf(1,'     Positionlist entry %d/%d:  structure %s, layer%s\n',k,length(obj.poslist),obj.poslist(k).name,layerstr);
            end
            
            fclose(fid);
            
            slin=strfind(filepath,slsh)+1;
            
            disp(['Positionlist ' filepath(slin(end):end) ' successfully written.']);
            disp(' ');
            
        end % writepls
        
        
    end % methods
    
    
    
    methods(Hidden)
        
        function plotchip(obj)
            if ~ishandle(obj.hChip)  % Only replot if it isn't already in figure
                obj.hChip=fill([0 1 1 0 0]*obj.chipUV(1)*1000,[0 0 1 1 0]*obj.chipUV(2)*1000,'w','EdgeColor','k');  % Plot chip outlines
                axis equal
            end
        end
        
        
        function checkposlist(obj)
        % Check that all positionlist data is of the correct format (i.e,.
        % before writing .pls file)
        
        errtxt=[];  % Record of errors
        warntxt=[]; % Record of warnings
            
            for kp=1:length(obj.poslist)
                Poslist=obj.poslist(kp);
            
                % Check for illegal fields
                posfields=fieldnames(Poslist);  % Cell array of all fields
                if ~isempty(setdiff(posfields,{'name','uv_c','DF','WA','layers'})) % Trying to set illegal field(s) for positionlist entry
                    errtxt=[errtxt sprintf('     Entry %d:  illegal data fields (allowed fields are ''name'', ''uv_c'', ''DF'', ''WA'', and ''layers''.\n',kp)];
                end

                % Argument checking for all data fields
                for k=1:length(posfields)
                    f=getfield(Poslist,posfields{k});
                    if ~isempty(f)

                        % Check for infinities and NaNs
                        if any(isnan(f(:))) || any(isinf(f(:)))
                           errtxt=[errtxt sprintf('     Entry %d:  ''%s'' must be finite.\n',kp,posfields{k})];
                        end

                        switch posfields{k}

                            case 'name'
                                % Check that name is a string
                                if ~ischar(Poslist.name)
                                    errtxt=[errtxt sprintf('     Entry %d:  ''name'' must be a string.\n',kp)];
                                end
                                % Check whether name contains any commas; this would be caught in the next check unless data checking is turned off (since Raith_structure objects don't allow commas in their names)
                                illchars=regexp(Poslist.name,'[^a-zA-Z_0-9\.\$\?-]');
                                if ~isempty(illchars)
                                    errtxt=[errtxt sprintf('     Entry %d:  ''name'' cannot contain these characters:  %s.\n',kp,Poslist.name(illchars))];                                    
                                end
                                % Check whether structname is in library
                                if ~any(strcmp(obj.library.structlist,Poslist.name))
                                    errtxt=[errtxt sprintf('     Entry %d:  structure ''%s'' is not found in GDSII library ''%s''.\n',kp,Poslist.name,obj.library.name)];
                                end

                            case 'uv_c'
                                if ~isnumeric(Poslist.uv_c) || ~isvector(Poslist.uv_c) || length(Poslist.uv_c)~=2
                                    errtxt=[errtxt sprintf('     Entry %d:  ''uv_c'' must be of the form [u_c v_c].\n',kp)];
                                end
                                if ~inpolygon(Poslist.uv_c(1),Poslist.uv_c(2),[0 1 1 0 0]*obj.chipUV(1),[0 0 1 1 0]*obj.chipUV(2))
                                    warntxt=[warntxt sprintf('     Entry %d:  ''uv_c'' for structure ''%s'' is outside chip boundaries.\n',kp,Poslist.name)];
                                end

                            case 'DF'
                                if ~isnumeric(Poslist.DF) || ~isscalar(Poslist.DF) || Poslist.DF<0
                                    errtxt=[errtxt sprintf('     Entry %d:  ''DF'' must be a non-negative scalar.\n',kp)];
                                end

                            case 'WA'
                                if ~isnumeric(Poslist.WA) || ~isvector(Poslist.WA) || length(Poslist.WA)~=4
                                    errtxt=[errtxt sprintf('     Entry %d:  ''WA'' must be of the form [u_min v_min u_max v_max].\n',kp)];
                                end
                                if Poslist.WA(1)>=Poslist.WA(3)
                                    errtxt=[errtxt sprintf('     Entry %d:  in ''WA'', u_min must be smaller than u_max.\n',kp)];
                                end
                                if Poslist.WA(2)>=Poslist.WA(4)
                                    errtxt=[errtxt sprintf('     Entry %d:  in ''WA'', v_min must be smaller than v_max.\n',kp)];
                                end
                                
                            case 'layers'
                                if ~isnumeric(Poslist.layers) || any(floor(Poslist.layers)~=Poslist.layers) || any(Poslist.layers<0) || any(Poslist.layers>63)
                                    errtxt=[errtxt sprintf('     Entry %d:  in ''layers'', all values must be integers between 0 and 63 (inclusive).\n',kp)];
                                end
                                
                        end

                    end

                end
                
            end

            if ~isempty(warntxt)
                warntxt=[sprintf('\nRaith_positionlist:  warnings in poslist entries\n') warntxt sprintf('\n')];
                warning(warntxt);
            end

            if ~isempty(errtxt)
                errtxt=[sprintf('\nRaith_positionlist:  errors in poslist entries\n') errtxt sprintf('\n')];
                error(errtxt);
            end
            
        end % checkposlist
        
        
        function morelayers=alllayers(obj,structname,layers)
        % Recursively find all layers present in a given structure,
        % including 'sref' and 'aref' objects
        
            STR=obj.library.structures(strcmp(obj.library.structlist,structname)); % Structure to be appended
            C=cell(1,length(STR.elements));
            [C{:}]=deal(STR.elements.type);
            explayins=~strcmp(C,'sref')&~strcmp(C,'aref');  % Indices of all elements in structure which are neither 'sref' or 'aref' (i.e., those that explicitly have a layer property)
            refins=strcmp(C,'sref')|strcmp(C,'aref');  % Indices of all 'sref' or 'aref' elements
            [C{:}]=deal(STR.elements.data);
            explays=unique(cellfun(@(x)x.layer,C(explayins))); % All layers present in structure
            morelayers=unique([layers explays]);
            
            % Grab all layers from 'sref' and 'aref' structures
            refnames=unique(cellfun(@(x)x.name,C(refins),'UniformOutput',0));
            for k=1:length(refnames)
               morelayers=obj.alllayers(refnames{k},morelayers);
            end
            
        end % alllayers
        
    end % hidden methods
    
end % classdef