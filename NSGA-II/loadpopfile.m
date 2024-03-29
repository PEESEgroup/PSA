function result = loadpopfile(fileName)
% Function: result = loadpopfile(fileName)
% Description: Load population file which generated by last optimization.
% Syntax:
%       oldresult = loadpopfile('populations.txt');
% Return: 
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************


%*************************************************************************
% Open the populations file
%*************************************************************************
if( ~ischar(fileName))
    error('NSGA2:LoadPopError', 'The fileName parameter should be string!');
end

fid = fopen(fileName, 'r');
if(fid==-1)
    error('NSGA2:LoadPopError', 'The populations file "%s" could not opened!', fileName);
end
%fprintf('...Loading the population file "%s"......\n', fileName);



%*************************************************************************
% Read file head
%*************************************************************************
popsize = 1;
maxGen  = 1;
nVar  = 0;
nObj  = 0;
nCons = 0;
fieldNames = {''};


strLine = fgetl(fid);
if( ~ischar(strLine) || strcmp(strLine, '#NSGA2')==0 )
    error('NSGA2:PopFileError', ...
        'The population file "%s" is not a Nsga2 populations file! Line��\n%s\n', ...
        fileName, strLine);
end

strLine = fgetl(fid);
while( ischar(strLine) && strcmp(strLine, '#end')==0)
    token = textscan(strLine, '%s');
    keyword = strtrim(token{1}{1});
    switch keyword
        case 'popsize'
            popsize = str2double(token{1}{2});
        case 'maxGen'
            maxGen = str2double(token{1}{2});
        case 'numVar'
            nVar = str2double(token{1}{2});
        case 'numObj'
            nObj = str2double(token{1}{2});
        case 'numCons'
            nCons = str2double(token{1}{2});
        case 'stateFieldNames'
            nfield = length(token{1});
            fieldNames = cell(nfield - 1, 1);
            [fieldNames{:}] = token{1}{2:end};
        otherwise
            warning('NSGA2:PopFileError', 'No support state keyword: "%s"', token{1}{1});
    end
    
    strLine = fgetl(fid);
end


%*************************************************************************
% Initialize the result structure
%*************************************************************************
pop = repmat( struct(...
    'var', zeros(1,nVar), ...
    'obj', zeros(1,nObj), ...
    'cons', zeros(1,nCons)),...
    [1,popsize]);

% state: optimization state of one generation
state = struct();
for i = 1:length(fieldNames)
    state.(fieldNames{i}) = 0;
end

result.pops     = repmat(pop, [maxGen, 1]);     % each row is the population of one generation
result.states   = repmat(state, [maxGen, 1]);   % each row is the optimizaiton state of one generation
clear i fieldNames pop


%*************************************************************************
% Parse the populations file
%*************************************************************************
lastGen = 0;    % The last generation which has correct datas.
try
    strLine = fgetl(fid);
    while ischar(strLine)
        %*****************************************************************
        % 1. Skip empty lines
        strLine = strtrim(strLine);
        if( isempty(strLine)  )
            strLine = fgetl(fid);   % read new line
            continue;
        end
        
        %*****************************************************************
        % 2. The first line of one generation
        if( strcmp(strLine(1:11), '#Generation') == 1)
            % Only the 'ngen' is needed
            ngen = sscanf(strLine(12:end), ' %d');
        else
            error('NSGA2:PopFileError', 'The population file format Error��Line��\n%s\n', strLine);
        end


        %*****************************************************************
        % 3. Read optimization states
        strLine = fgetl(fid);
        while( ischar(strLine) && strcmp(strLine, '#end')==0 )
            token = textscan(strLine, '%s%f');
            result.states(ngen).(token{1}{1}) = token{1,2};
            
            strLine = fgetl(fid);
        end
            
        
        %*****************************************************************
        % 4. Read population
        strLine = fgetl(fid);
        val = fscanf(fid, '%f');
        ncols = nVar+nObj+nCons;
        nrows = popsize;
        if( length(val) ~= ncols*nrows )
            error('NSGA2:PopFileError', 'File error when read population datas! ');
        end

        val = reshape(val, ncols, nrows)';  % reshape the vector column-wise
        for i = 1:popsize
            result.pops(ngen, i).var = val(i, 1:nVar);
            result.pops(ngen, i).obj = val(i, (nVar+1):(nVar+nObj));
            result.pops(ngen, i).cons= val(i, (nVar+nObj+1):end);
        end
        
        %*****************************************************************
        % Read next line
        lastGen = ngen;
        strLine = fgetl(fid);

    end
catch exception
    id = exception.identifier;
    msg = [exception.message, 'File = ', fileName];
    warning(id, msg);       % If the file is wrong at the end, return the correct data.
end


%*************************************************************************
% Do some clean works
%*************************************************************************

% Delete unused datas
result.pops(lastGen+1:end, :) = [];
result.states(lastGen+1:end)  = [];
%fprintf('...Load population file success. Last generation is %d.\n', lastGen);


% Close file
fclose(fid);






