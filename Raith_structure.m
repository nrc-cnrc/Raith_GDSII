classdef Raith_structure < handle
%
% obj=Raith_structure(name,elements)
%
% Raith_structure objects define named structures composed of low-level 
% elements, used to create GDSII hierarchies for Raith beamwriting tools.
%
%
% Arguments:
%
% name - string specifying name of structure; may be up to 127 characters
%   long; allowed characters are A-Z, a-z, 0-9, _, ., $ ,?, and -.
% elements - array of Raith_element objects in structure
%
%
% Properties:
%
%   name - name of structure
%   elements - Raith_element array of low-level elements in structure
%   reflist - cell array of structure names referenced by 'sref' or 'aref'
%       elements within the structure
%
%
% Methods:
%
%   plot([M,scDF]) - plot structure with Raith dose factor colouring (filled
%       polygons where applicable)
%           M - augmented transformation matrix for plot [optional]
%           scDF - overall multiplicative scaling factor for all elements 
%               in structure (e.g., passed from a positionlist entry)
%               [optional]
%
%   plotedges([M,scDF]) - plot structure with Raith dose factor colouring 
%       (edges of polygons where applicable)
%           M - augmented transformation matrix for plot [optional]
%           scDF - overall multiplicative scaling factor for all elements 
%               in structure (e.g., passed from a positionlist entry)
%               [optional]
%
%
% Dependencies:  Raith_element.m
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
        elements=Raith_element.empty;
    end % properties
    
    properties (SetAccess=private)
        reflist % cell array of structure names referenced by 'sref' or 'aref' elements within the structure
    end % set access private properties

    methods
        
        function obj=Raith_structure(name,elements)
            if nargin>0
                obj.name=name;
                if nargin==2
                    obj.elements=elements;
                end
            end
        end % Constructor
        
        
        function set.name(obj,Name)
            
            if ~ischar(Name)
                error('Raith_structure:  name must be a string.');
            end
            
            if length(Name)>127
                error('Raith_structure:  maximum name length is 127.');
            end
            
            if ~isempty(regexp(Name,'[^a-zA-Z_0-9\.\$\?-]','once'))
                Name=regexprep(Name,'[^a-zA-Z_0-9\.\$\?-]','_');  % Replace illegal characters with underscores
                warning('Raith_structure:illegalCharacters','Illegal characters in structure name:  changing name to %s',Name);
            end
                
            obj.name=Name;
            
        end % set.name
        
        
        function set.elements(obj,E)
            
            global checkdata;
            
            if checkdata==false

                obj.elements=E;
            
            else
                
                if ~isa(E,'Raith_element')
                    error('Raith_structure:  all elements must be Raith_element objects.');
                end
            
                obj.elements=E;
            
                % Populate reflist
                C=cell(1,numel(E));
                [C{:}]=deal(E.type);
                ins=[find(strcmp(C,'sref')) find(strcmp(C,'aref'))];  % Indices of 'sref' or 'aref' elements
                Reflist=cell(1,length(ins));
                for k=1:length(ins)
                    Reflist{k}=E(ins(k)).data.name;
                end
                obj.reflist=unique(Reflist);
            
            end

        end % set.name
        
        
        
        
        function plot(obj,varargin)
        %
        % Raith_structure.plot([M,scDF]) 
        %
        % Plot structure with Raith dose factor colouring (filled polygons where applicable)
        %
        % Argument:
        %
        %   M - augmented transformation matrix for plot [optional]
        %   scDF - overall multiplicative scaling factor for all elements 
        %   	in structure (e.g., passed from a positionlist entry)
        %   	optional]
        %            
            
            if nargin==1 
                M=[];  
                scDF=1;  % Scaling for DF
            elseif nargin==2
                M=varargin{1};
                scDF=1;  
            elseif nargin==3
                M=varargin{1};
                scDF=varargin{2};  
            else
                error('Raith_structure.plot:  Too many input arguments.');
            end
        
            for k=1:numel(obj.elements)
                obj.elements(k).plot(M,scDF);
            end
            
        end % plot
        
        
        function plotedges(obj,varargin)
        %
        % Raith_structure.plotedges([M,scDF]) 
        %
        % Plot structure with Raith dose factor colouring (edges of polygons where applicable)
        %
        % Argument:
        %
        %   M - augmented transformation matrix for plot [optional]
        %   scDF - overall multiplicative scaling factor for all elements 
        %   	in structure (e.g., passed from a positionlist entry)
        %       [optional]
        %
            
            if nargin==1 
                M=[];  
                scDF=1;  % Scaling for DF
            elseif nargin==2
                M=varargin{1};
                scDF=1;  
            elseif nargin==3
                M=varargin{1};
                scDF=varargin{2};  
            else
                error('Raith_structure.plotedges:  Too many input arguments.');
            end
        
            for k=1:numel(obj.elements)
                obj.elements(k).plotedges(M,scDF);
            end
            
        end % plotedges
                
    end
    
end