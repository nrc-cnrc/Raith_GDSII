classdef Raith_library < handle
%
% obj=Raith_library(name,structures)
%
% Raith_library objects define GDSII hierarchies for Raith beamwriting
% tools.
%
%
% Arguments:
%
% name - string specifying name of GDSII library, not including .csf
%   extension
% structures - array of Raith_structure objects in library
%
%
% Properties:
%
%   name - name of GDSII library
%   structures - array of Raith_structure objects in library
%   structlist - ordered cell array of all names of structures found in
%       library (private set access)
%
%
% Methods:
%
%   append(S) - append Raith_structure S (or array thereof) to library; 
%       structure names are checked for uniqueness
%
%   writegds([outdir],[dialect]) - write Raith GDSII hierarchy file of all
%       structures as [library.name].csf
%           outdir - string specifying directory in which to write .csf 
%               file; if called without arguments, file is written to 
%               working directory
%           dialect - string specifying dialect of GDSII to write
%               [optional]; may be 'Raith' (default) or 'plain'; if 'plain' 
%               is selected, Raith curved elements are converted to 
%               boundary (polygon) or path elements, as appropriate
%
%   plot(structname,[M,scDF]) - plot structure with Raith dose factor 
%       colouring (filled polygons where applicable)
%           structname - name of structure to be plotted (must be in
%              structlist); if structure contains sref or aref elements,
%              the entire hierarchy is plotted if the referenced structures
%              are in the library, and placeholder names if not
%           M - augmented transformation matrix for plot [optional]
%           scDF - overall multiplicative scaling factor for all elements 
%               in structure (e.g., passed from a positionlist entry)
%               [optional]
%
%   plotedges(structname,[M,scDF]) - plot structure with Raith dose factor 
%       colouring (edges of polygons where applicable)
%           structname - name of structure to be plotted (must be in
%              structlist); if structure contains sref or aref elements,
%              the entire hierarchy is plotted if the referenced structures
%              are in the library, and placeholder names if not
%           M - augmented transformation matrix for plot [optional]
%           scDF - overall multiplicative scaling factor for all elements 
%               in structure (e.g., passed from a positionlist entry)
%               [optional]
%
%
% Dependencies:  Raith_element.m, Raith_structure.m
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
        name
        structures=Raith_structure.empty;
    end % properties
    
    properties (SetAccess=private)
        structlist % Ordered cell array of all names of structures found in library 
    end % set access private properties
    
       
    methods
        
        function obj=Raith_library(name,structures)
            if nargin>0
                
                % Check that name is a string
                if ~ischar(name)
	                error('Raith_library:  name must be a string.');
                end
                
                obj.name=name;
                
                if nargin==2
                    
                    for k=1:numel(structures)
                        snames{k}=structures(k).name;
                    end
                    
                    % Check for name uniqueness
                    if length(snames)~=length(unique(snames))
                        error('Raith_library:  all structures in library must have unique names.');
                    end
                    
                    obj.structlist=snames;
                    obj.structures=structures;
                end
                
            end
            
        end % Constructor
        
        
        function set.name(obj,Name)
            
            if ~ischar(Name)
                error('Raith_library:  name must be a string.');
            end
            
            if length(Name)>3 && (strcmpi(Name((end-3):end),'.csf') || strcmpi(Name((end-3):end),'.gds'))
                xtn=Name((end-3):end);
                Name((end-3):end)=[];  % Remove extension
                warning(['Raith_library:  removing ' xtn ' from library name.'])
            end
            
            obj.name=Name;
            
        end % set.name
        
        
        function set.structures(obj,S)
            
            global checkdata
            
            if checkdata==false
                % Nothing to do
            else
                % Type checking
                if ~isa(S,'Raith_structure')
                  error('Raith_library:  only Raith_structure objects may be appended to library.');
                end
            end

            obj.structures=S;

            % Populate structure list
            for k=1:numel(obj.structures)
               snames{k}=obj.structures(k).name;
            end

            % Check for name uniqueness
            if length(snames)~=length(unique(snames))
                error('Raith_library:  all structures in library must have unique names.');
            end

            obj.structlist=snames;
           
        end % set.structures
        
        
        function append(obj,S)
        %
        %   Raith_library.append(S) - append Raith_structure S (or array thereof) to library
        %
        %   Argument:
        %
        %       S - Raith_structure object to append to library
        %


            for k=1:numel(S)
                obj.structlist{end+1}=S(k).name;  % Add name to structlist
                obj.structures(end+1)=S(k);  % Add structure to library
            end
            
        end % append
        
        
        
        function writegds(obj,varargin)
        %
        %   Raith_library.writegds([outdir],[dialect]) - write Raith GDSII hierarchy file of 
        %       all structures as [library.name].csf
        %
        %   Argument:
        %   
        %       outdir - string specifying directory in which to write .csf file [optional]; 
        %           if called without arguments, file is written to working directory
        %
        %       dialect - string specifying dialect of GDSII to write
        %           [optional]; may be 'Raith' (default) or 'plain'; if
        %           'plain' is selected, Raith curved elements are
        %           converted to boundary (polygon) or path elements, as
        %           appropriate
        %
            
            % See http://www.rulabinsky.com/cavd/text/chapc.html

            if nargin==1 % No path for output directory or dialect choice
                outdir=pwd;
                dialect='raith';
            elseif nargin==2  % Either output directory or dialect is given
                if any(strcmpi(varargin{1},{'Raith','plain'}))
                    dialect=lower(varargin{1});
                    outdir=pwd;
                elseif isdir(varargin{1})
                    outdir=varargin{1};
                    dialect='raith';
                else
                    error('Raith_library.writegds:  invalid output directory or GDSII dialect.');
                end
            elseif nargin==3 % Both output directory and dialect are given, in that order
                if isdir(varargin{1})
                    outdir=varargin{1};
                else
                    error('Raith_library.writegds:  invalid output directory.');
                end
                if any(strcmpi(varargin{2},{'Raith','plain'}))
                    dialect=lower(varargin{2});
                else
                   error('Raith_library.writegds:  invalid GDSII dialect.');
                end 
            else
                error('Raith_library.writegds:  too many arguments.');
            end
            
            if strcmpi(dialect,'plain')
                ext = '.gds';
            elseif strcmpi(dialect,'Raith')
                ext = '.csf';
            end
            
            % Kludgy way to ensure structure name uniqueness (i.e., if
            % names were changed after assigning obj.structures); rewrite
            % using listeners in a later version.
            obj.structures=obj.structures;  
            
            global checkdata;
            
            if checkdata==false
                
                fprintf(1,'\nSkipping all data checking.\n');
                
            else

                % Check whether all objects referenced in sref and aref elements are contained in the library
                fprintf(1,'\nChecking for missing structures...');
                C=cell(1,numel(obj.structures));
                [C{:}]=deal(obj.structures.reflist);
                ein=~cellfun('isempty',C);  % Indices of non-empty cells
                allrefs=unique([C{ein}]);  % List of all structures referenced in 'sref' and 'aref' elements in library
				matlabver = version;
				% setdiff behaviour was changed in R2013a (Matlab version 8)
				if matlabver(1) <= '8'
					missing=setdiff(allrefs,obj.structlist);  % Missing structures (those in allrefs that aren't in structlist)
				else
					missing=setdiff(allrefs,obj.structlist,'Legacy');  % Use legacy option for Matlab R2013a and later
				end

                if ~isempty(missing)
                    fprintf(2,'fail.\n');
                    fprintf(2,'\n The following referenced structures are not present in the library:\n\n');
                    for kk=1:length(missing)
                        fprintf(2,'     %s\n',missing{kk});
                    end
                    fprintf('\n');
                    error('Raith_library.writegds:  missing structures.');
                else
                    fprintf(1,'OK.\n');
                end
                
            end
            
            % Write all structures to a Raith-readable GDSII file
            if isunix
                slsh='/';
            else
                slsh='\';
            end
            
            FileID=fopen([outdir slsh obj.name ext],'w');            
            disp(['Writing ' outdir slsh obj.name ext '...']);
            
            disp('     Header information');
            % Write all header information
            obj.writehead(FileID,obj.name);
            
            % Write all structures in library
            for ks=1:numel(obj.structures)
                
                STR=obj.structures(ks);  % Current structure
                fprintf(1,'     Structure %d/%d:  %s\n',ks,numel(obj.structures),STR.name);
                
                obj.writebeginstruct(FileID,STR.name);  % Begin a structure
                
                % Write all elements in structure
                for ke=1:numel(STR.elements)
                    ELE=STR.elements(ke);
                    if strcmp(dialect,'plain') && any(strcmp(ELE.type,{'arc','circle','ellipse','fbmspath','fbmscircle'}))
                        
                        if isempty(ELE.data.w) ||  (ELE.data.w~=0 && any(strcmp(ELE.type,{'circle','ellipse','fbmscircle'}))) % Filled object:  convert to polygon
                            UV=ELE.renderplot(eye(3),1,2);  % Use polygon as rendered in Raith_element.plot
                            wELE=Raith_element('polygon',ELE.data.layer,UV,ELE.data.DF);
                        else  % Convert to path
                            w=ELE.data.w;  % Original path width
                            ELE.data.w=0;  % Set to single-pixel line to get vertices of underlying path
                            UV=ELE.renderplot(eye(3),1,2);  % Use path as rendered in Raith_element.plot
                            ELE.data.w=w;  % Restore original width
                            wELE=Raith_element('path',ELE.data.layer,UV,ELE.data.w,ELE.data.DF);
                        end
                        
                        disp(['          Converting ' ELE.type ' element to ' wELE.type]); 
                    else
                        wELE=ELE;
                    end
                    obj.writeelement(FileID,wELE);
                end
                
                obj.writeendstruct(FileID);  % End the structure
                         
            end
            
            obj.writeendlib(FileID);  % End the library
            
            fclose(FileID);
            fprintf(1,['GDSII library ' obj.name '.csf successfully written.\n\n']);
            
        end % writegds
        
        
        function plot(obj,structname,varargin)
        %
        %   Raith_library.plot(structname,[M,scDF]) - plot structure with Raith 
        %       dose factor colouring (filled polygons where applicable)
        %
        %   Arguments:
        %
        %       structname - name of structure to be plotted (must be in
        %       	structlist); if structure contains sref or aref elements,
        %           the entire hierarchy is plotted if the referenced structures
        %           are in the library, and placeholder names if not
        %
        %       M - augmented transformation matrix for plot [optional]
        %       scDF - overall multiplicative scaling factor for all elements 
        %           in structure (e.g., passed from a positionlist entry)
        %           [optional]
        %
        
            if nargin==1
                error('Raith_library.plot:  name of structure to be plotted is required as an argument.');
            elseif nargin==2 
                M=eye(3);  
                scDF=1;
            elseif nargin==3
                M=varargin{1};
                scDF=1;
            elseif nargin==4
                M=varargin{1};
                scDF=varargin{2};  
            else
                error('Raith_library.plot:  too many input arguments.');
            end
            
            if ~any(strcmp(obj.structlist,structname))  % Structure is not in library
                error(['Raith_library.plot:  structure ''' structname ''' is not in library.']);
            end
            
            STR=obj.structures(strcmp(obj.structlist,structname)); % Current structure
            
            for ke=1:numel(STR.elements)
                
                ELE=STR.elements(ke);  % Current element
                
                switch ELE.type

                    case {'polygon','path','dot','circle','ellipse','text','arc','fbmspath','fbmscircle'}
                        ELE.plot(M,scDF);
                        
                    case 'sref'
                        M2=M*obj.trans([ELE.data.uv_0(1) ELE.data.uv_0(2)])*obj.rot(ELE.data.angle)*obj.scale(ELE.data.mag)*obj.refl(ELE.data.reflect);
                        if ~any(strcmp(obj.structlist,ELE.data.name))
                            ELE.plot(M,scDF);
                        else
                            obj.plot(ELE.data.name,M2,scDF);
                        end
                        
                    case 'aref'
                        % Construct lattice of origins for structures
                        [U,V]=meshgrid((0:ELE.data.n_colrow(1)-1)*ELE.data.a_colrow(1),(0:ELE.data.n_colrow(2)-1)*ELE.data.a_colrow(2));
                        U=reshape(U,1,prod(ELE.data.n_colrow));
                        V=reshape(V,1,prod(ELE.data.n_colrow));
                        
                        UV=M*obj.trans([ELE.data.uv_0(1) ELE.data.uv_0(2)])*obj.rot(ELE.data.angle)*[U;V;ones(size(U))];  % Transformations applied to aref object (if coming from a parent sref or aref element)
                        M2=M;
                        M2(7:8)=0;  % Remove translation component
                        
                        for k=1:length(U)
                            Mk=obj.trans(UV(:,k))*M2*obj.rot(ELE.data.angle)*obj.scale(ELE.data.mag)*obj.refl(ELE.data.reflect);    
                            if ~any(strcmp(obj.structlist,ELE.data.name))
                                ELE.plot(M2,scDF);
                            else
                                obj.plot(ELE.data.name,Mk,scDF);
                            end
                        end
                        
                end
            
            end
           
        end % plot
        
        
         function plotedges(obj,structname,varargin)
        %
        %   Raith_library.plotedges(structname,[M,scDF]) - plot structure with Raith 
        %       dose factor colouring (edges of polygons where applicable).
        %
        %   Arguments:
        %
        %       structname - name of structure to be plotted (must be in
        %       	structlist); if structure contains sref or aref elements,
        %           the entire hierarchy is plotted if the referenced structures
        %           are in the library, and placeholder names if not
        %
        %       M - augmented transformation matrix for plot [optional]
        %       scDF - overall multiplicative scaling factor for all elements 
        %           in structure (e.g., passed from a positionlist entry)
        %           [optional]
        %
            
            if nargin==1
                error('Raith_library.plotedges:  name of structure to be plotted is required as an argument.');
            elseif nargin==2 
                M=eye(3);  
                scDF=1;
            elseif nargin==3
                M=varargin{1};
                scDF=1;
            elseif nargin==4
                M=varargin{1};
                scDF=varargin{2};  
            else
                error('Raith_library.plotedges:  too many input arguments.');
            end
            
            if ~any(strcmp(obj.structlist,structname))  % Structure is not in library
                error(['Raith_library.plotedges:  structure ''' structname ''' is not in library']);
            end
            
            STR=obj.structures(strcmp(obj.structlist,structname)); % Current structure
            
            for ke=1:numel(STR.elements)
                
                ELE=STR.elements(ke);  % Current element
                
                switch ELE.type

                    case {'polygon','path','dot','circle','ellipse','text','arc'}
                        ELE.plotedges(M,scDF);
                        
                    case 'sref'
                        M2=M*obj.trans([ELE.data.uv_0(1) ELE.data.uv_0(2)])*obj.rot(ELE.data.angle)*obj.scale(ELE.data.mag)*obj.refl(ELE.data.reflect);
                        if ~any(strcmp(obj.structlist,ELE.data.name))
                            ELE.plotedges(M,scDF);
                        else
                            obj.plotedges(ELE.data.name,M2,scDF);
                        end
                        
                    case 'aref'
                        
                        % Construct lattice of origins for structures
                        [U,V]=meshgrid((0:ELE.data.n_colrow(1)-1)*ELE.data.a_colrow(1),(0:ELE.data.n_colrow(2)-1)*ELE.data.a_colrow(2));
                        U=reshape(U,1,prod(ELE.data.n_colrow));
                        V=reshape(V,1,prod(ELE.data.n_colrow));
                        
                        UV=M*obj.trans([ELE.data.uv_0(1) ELE.data.uv_0(2)])*obj.rot(ELE.data.angle)*[U;V;ones(size(U))];  % Transformations applied to aref object (if coming from a parent sref or aref element)
                        M(7:8)=0;  % Remove translation component
                        
                        for k=1:length(U)
                            Mk=obj.trans(UV(:,k))*M*obj.rot(ELE.data.angle)*obj.scale(ELE.data.mag)*obj.refl(ELE.data.reflect);    
                            if ~any(strcmp(obj.structlist,ELE.data.name))
                                ELE.plotedges(M,scDF);
                            else
                                obj.plotedges(ELE.data.name,Mk,scDF);
                            end
                        end
                        
                end
            
            end
            
        end % plotedges
        
        
    end % methods
    
    
    
    methods(Static)
       
        
        function M=rot(theta)
        %
        % M=rot(theta)
        %
        % Return augemented 2D rotation matrix by an angle theta (in degrees).
        %

            M=[cosd(theta) -sind(theta) 0;sind(theta) cosd(theta) 0;0 0 1];

        end % rot
        
        
        function M=trans(p)
        %
        % M=trans(p)
        %
        % Return augemented translation matrix for a translation vector p.
        %

            M=[1 0 p(1);0 1 p(2);0 0 1];

        end % trans
        
        
        function M=refl(n)
        %
        % M=refl(n)
        %
        % Return augemented matrix for reflection about u-axis n times
        %

            M=[1 0 0;0 (-1)^n 0;0 0 1];

        end % refl
        
        
        function M=scale(mag)
        %
        % M=scale(mag)
        %
        % Return augemented matrix for scaling by a factor mag
        %

        M=[mag 0 0;0 mag 0;0 0 1];

        end % scale
        
        
        function writerec(FileID,rectype,datatype,parameters)  
                % Write GDSII record.  
                %   
                % Arguments:
                %   FileID - file identifier
                % 	rectype - 1-byte record type
                % 	datatype - 1-byte data type: 
                %       00 - no data present
                %       01 - bit array (2 bytes)
                %       02 - 2-byte signed integer
                %       03 - 4-byte signed integer
                %       05 - 8-byte float
                %       06 - ASCII string
                % 	parameters - record parameters, of type defined by datatype
                %
                % The length of the record is computed, padded to an even
                % number of bytes if necessary, and the appropriate 2-byte
                % length header prepended before writing.
                %
                
                switch datatype
                    case 0 % No data
                        lp=0; 
                        fmt='uint8'; % Ignored, since parameters should be empty
                    case 1 % 2-byte integer (bit array)
                        lp=2; 
                        fmt='uint16';
                    case 2 % 2-byte integer
                        lp=2; 
                        fmt='uint16';
                    case 3 % 4-byte signed integer
                        lp=4; 
                        fmt='int32';
                    case 5 % 8-byte float, but must convert to excess-64 format
                        lp=1;
                        fmt='uint8';
                        parms=parameters;
                        parameters=zeros(1,8*length(parms));  % Each float converted to 8 bytes
                        for kp=0:(length(parms)-1)
                            parameters((8*kp+1):(8*kp+8))=excess64(parms(kp+1));
                        end
                    case 6 % 1-byte ASCII character
                        lp=1;
                        fmt='uint8';
                    otherwise
                        error('Unknown datatype.');
                end
                
                l=length(parameters)*lp;  % Length of parameters in bytes
                
                fwrite(FileID,ceil(l/2)*2+4,'uint16','b'); % Record length
                fwrite(FileID,[rectype datatype],'uint8','b');  % Record and data types
                fwrite(FileID,parameters,fmt,'b');
                
                if mod(l,2)==1 % Odd-size parameter data
                    fwrite(FileID,0,'integer*1','b');
                end
                
            function d=excess64(N)
                % Convert number to GDSII 8-byte floating-point format
                % (excess-64), return as 8 bytes in decimal representation.
                % Based on gdsii_excess64enc.c by Ulf Griesmann, NIST (Jan. 2008)
                
                if N==0 % Special case of zero input
                    d=zeros(1,8);
                    return;
                end
                   
                % Extract sign
                if N<0
                    sgn=-1;
                    N=-N;
                else
                    sgn=1;
                end
                
                [mn,ex]=log2(N);  % N=mn*2^ex
                ex16=ex/4;  % Scale exponent to base 16
                E=ceil(ex16);  % Base 16 exponent
                M=mn*16^(ex16-E);  % Base 16 mantissa
                
                % Normalise representation such that 1/16 <= M < 1
                while M>=1
                    M=M/16;
                    E=E+1;
                end
                
                while M<(1/16)
                    M=M*16;
                    E=E-1;
                end
                
                E=E+64;  % Exponent for excess-64 format
                if sgn<0
                   E=bitor(E,128);  % Leading bit=1 for negative number
                end
                
                % Acquire each byte of mantissa
                b=zeros(1,7);
                for k=1:7
                    M=M*256;  % Shift 8 bits to the right
                    b(k)=floor(M);  % Grab integer part
                    M=M-b(k);  
                end
                
                d=[E b];  % Sign bit, 7-bit exponent, 7 bytes of mantissa
                
            end % excess64
             
        end % writerec
            
        
        function writehead(FileID,name)
            
            % Write GDSII library header records
            
            % Preliminary info
            Raith_library.writerec(FileID,0,2,3); % HEADER (0002); release info
            d=round(datevec(now));
            Raith_library.writerec(FileID,1,2,[d d]);  % BGNLIB (0102); last modified and accessed timestamps
            Raith_library.writerec(FileID,2,6,name);  % LIBNAME (0206)
            Raith_library.writerec(FileID,3,5,[1e-3 1e-9]);  % UNITS (0305)
            
        end % writehead
        
        
        function writeelement(FileID,ELE)
            
            % Write GDSII element records
            
            switch lower(ELE.type)
                case 'polygon'
                    Raith_library.writerec(FileID,8,0,[]);  % BOUNDARY (0800)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    XY=reshape([ELE.data.uv]*1000,1,2*size(ELE.data.uv,2));  % Reshape to 1D array of sequential XY pairs (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

                case 'path'
                    Raith_library.writerec(FileID,9,0,[]);  % PATH (0900)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    Raith_library.writerec(FileID,15,3,1000*ELE.data.w);  % WIDTH (0f03); in nm
                    XY=reshape([ELE.data.uv]*1000,1,2*size(ELE.data.uv,2));  % Reshape to 1D array of sequential XY pairs (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

                case 'dot'  % Single-pixel dots are defined as paths with two identical vertices, no width
                    if isscalar(ELE.data.DF)
                        DF=ELE.data.DF*ones(1,size(ELE.data.uv,2));
                    else
                        DF=ELE.data.DF;
                    end
                    for kd=1:size(ELE.data.uv,2)
                        Raith_library.writerec(FileID,9,0,[]);  % PATH (0900)
                        Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                        Raith_library.writerec(FileID,14,2,1000*DF(kd));  % DATATYPE (0e02); 1000*DF
                        XY=[ELE.data.uv(1,kd) ELE.data.uv(2,kd) ELE.data.uv(1,kd) ELE.data.uv(2,kd)]*1000;  % XY pairs of "path" defining dot (in nm)
                        Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                        Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)
                    end

                case 'arc'
                    Raith_library.writerec(FileID,86,0,[]);  % Raith curved element (5600)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    if isscalar(ELE.data.r)
                        fflag=4; % Flag for unfilled circular disk segment
                        r=[1 1]*ELE.data.r;
                    else
                        fflag=5; % Flag for unfilled elliptical disk segment
                        r=ELE.data.r;
                    end
                    if isempty(ELE.data.w)  % Filled disk segment
                        fflag=fflag+2;  % Add 2 for filled segment
                    end
                    if ELE.data.w~=0  % Width specified
                        Raith_library.writerec(FileID,15,3,1000*ELE.data.w);  % WIDTH (0f03); in nm
                    end
                    if ELE.data.angle~=0
                        Raith_library.writerec(FileID,28,5,ELE.data.angle);  % ANGLE (1c05); in degrees w.r.t. positive u axis
                    end
                    XY=[ELE.data.uv_c(1)*1000 ELE.data.uv_c(2)*1000 r*1000 round(ELE.data.theta*785398/45) ELE.data.N fflag];  % XY pairs defining arc properties (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

                case 'circle'
                    Raith_library.writerec(FileID,86,0,[]);  % Raith curved element (5600)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    if isempty(ELE.data.w) 
                        fflag=2; % Flag for filled circle
                    else
                        fflag=0;  % Flag for unfilled circle
                        if ELE.data.w~=0  % Width specified
                            Raith_library.writerec(FileID,15,3,1000*ELE.data.w);  % WIDTH (0f03); in nm
                        end
                    end
                    XY=[ELE.data.uv_c(1)*1000 ELE.data.uv_c(2)*1000 ELE.data.r*[1 1]*1000 0 0 ELE.data.N fflag];  % XY pairs defining circle properties (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

                case 'ellipse'
                    Raith_library.writerec(FileID,86,0,[]);  % Raith curved element (5600)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    if isempty(ELE.data.w) 
                        fflag=3; % Flag for filled ellipse
                    else
                        fflag=1;  % Flag for unfilled ellipse
                        if ELE.data.w~=0  % Width specified
                            Raith_library.writerec(FileID,15,3,1000*ELE.data.w);  % WIDTH (0f03); in nm
                        end
                    end
                    Raith_library.writerec(FileID,28,5,ELE.data.angle);  % ANGLE (1c05); in degrees w.r.t. positive u axis
                    XY=[ELE.data.uv_c(1)*1000 ELE.data.uv_c(2)*1000 ELE.data.r*1000 0 0 ELE.data.N fflag];  % XY pairs defining ellipse properties (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

                case 'text'  % Text - must first render as a series of polygon elements, then write as such
                    T=ELE.rendertext(ELE.data.layer,ELE.data.uv_0,ELE.data.h,ELE.data.angle,ELE.data.uv_align,ELE.data.textlabel,ELE.data.DF);
                    for kT=1:length(T.elements)
                        Raith_library.writerec(FileID,8,0,[]);  % BOUNDARY (0800)
                        Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                        Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                        XY=reshape([T.elements(kT).data.uv]*1000,1,2*size(T.elements(kT).data.uv,2));  % Reshape to 1D array of sequential XY pairs (in nm)
                        Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                        Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)
                    end
                    
                case 'fbmspath'  % FBMS path
                    Raith_library.writerec(FileID,88,0,[]);  % Raith curved element (5800)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    Raith_library.writerec(FileID,15,3,1000*ELE.data.w);  % WIDTH (0f03); in nm
                    if isscalar(ELE.data.cvtr) % Zero curvature for all segments
                        cvtr=zeros(1,size(ELE.data.uv,2));
                    else
                        cvtr=ELE.data.cvtr;
                    end
                    cvXY=[(cvtr~=0)+1;ELE.data.uv*1000;cvtr*1000]; % Each vertex specified by a quadruple:  [curvature_flag u_i v_i curvature_i], where curvature_flag is 1 (2) for line segment (arc)...
                    cvXY(1)=0;  % ...except curvature_flag for the first vertex, which is always zero
                    XY=[0 0 0 0 reshape(cvXY,1,numel(cvXY))];  % Reshape to 1D array of sequential XY pairs
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)
                    
                case 'fbmscircle'  % FBMS circle
                    Raith_library.writerec(FileID,88,0,[]);  % Raith curved element (5800)
                    Raith_library.writerec(FileID,13,2,ELE.data.layer);  % LAYER (0d02)
                    Raith_library.writerec(FileID,14,2,1000*ELE.data.DF);  % DATATYPE (0e02); 1000*DF
                    Raith_library.writerec(FileID,15,3,1000*ELE.data.w);  % WIDTH (0f03); in nm
                    XY=[0 0 0 0 0 ELE.data.uv_c*1000 ELE.data.r*1000];  % Reshape to 1D array of sequential XY pairs
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)
                    
                case 'sref'  % Structure reference
                    Raith_library.writerec(FileID,10,0,[]);  % SREF (0A00)
                    Raith_library.writerec(FileID,18,6,ELE.data.name);  % SNAME (1206)            
                    if ELE.data.reflect
                        Raith_library.writerec(FileID,26,1,32768);  % STRANS (1A01); reflect in u
                    end
                    if ELE.data.mag~=1 || ELE.data.angle~=0
                        if ~ELE.data.reflect
                            Raith_library.writerec(FileID,26,1,0);  % STRANS (1A01); no reflection in u
                        end
                        Raith_library.writerec(FileID,27,5,ELE.data.mag);  % MAG (1B05)
                        Raith_library.writerec(FileID,28,5,ELE.data.angle);  % ANGLE (1C01)
                    end
                    XY=[ELE.data.uv_0(1)*1000 ELE.data.uv_0(2)*1000];  % Single origin for structure reference (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

                case 'aref'  % Array reference  
                    Raith_library.writerec(FileID,11,0,[]);  % AREF (0B00)
                    Raith_library.writerec(FileID,18,6,ELE.data.name);  % SNAME (1206)            
                    if ELE.data.reflect
                        Raith_library.writerec(FileID,26,1,32768);  % STRANS (1A01); reflect in u
                    end
                    if ELE.data.mag~=1 || ELE.data.angle~=0
                        if ~ELE.data.reflect
                            Raith_library.writerec(FileID,26,1,0);  % STRANS (1A01); no reflection in u
                        end
                        Raith_library.writerec(FileID,27,5,ELE.data.mag);  % MAG (1B05)
                        Raith_library.writerec(FileID,28,5,ELE.data.angle);  % ANGLE (1C01)
                    end
                    Raith_library.writerec(FileID,19,2,ELE.data.n_colrow);  % COLROW (1302)
                    XY=[ELE.data.uv_0(1) ELE.data.uv_0(2) ELE.data.uv_0(1)+ELE.data.n_colrow(1)*ELE.data.a_colrow(1) ELE.data.uv_0(2) ELE.data.uv_0(1) ELE.data.uv_0(2)+ELE.data.n_colrow(2)*ELE.data.a_colrow(2)]*1000;  % Corner instance origin, n_columns*column_spacing+ref_x, n_rows*row_spacing+ref_x (in nm)
                    Raith_library.writerec(FileID,16,3,XY);  % XY (1003)
                    Raith_library.writerec(FileID,17,0,[]);  % ENDEL (1100)

            end
            
        end % writeelement
        
        
        function writebeginstruct(FileID,name)
            
            % Write records to begin a structure (BGNSTR and STRNAME)
            
            d=round(datevec(now));
            Raith_library.writerec(FileID,5,2,[d d]);  % BGNSTR (0502); last modified and accessed timestamps
            Raith_library.writerec(FileID,6,6,name);  % STRNAME (0606)
            
        end  % writebeginstruct
        
        
        function writeendstruct(FileID)
            
            % Write ENDSTR record
            
            Raith_library.writerec(FileID,7,0,[]);  % ENDSTR (0700)
            
        end  % writeendstruct
        
        
        function writeendlib(FileID)
            
            % Write ENDLIB record
            
            Raith_library.writerec(FileID,4,0,[]);  % ENDLIB (0400)
            
        end  % writeendlib

        
        
    end % Static, hidden methods
    
   
    
end
