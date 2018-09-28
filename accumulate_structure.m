%% fill the structure array options

% in this file the structure options (some sort of global that contains the most
% important settings) is built, and common default values are set for
% the parameters.
% if different parameter values are specified in casename_parameters.m,
% then the ones in that file will be used.

%% accumulate options
% object=[];
% 
% voi={
%     
%     };
% 
% accumulate_object;

%% case information
object='case';

voi={
    'project',      [];
    };

accumulate_object;

%% physical properties
object='fluid';

voi={
        'Re',      [];
    };
    
accumulate_object;

%% grid parameters
object='grid';

voi={
    'Nx',      [];...
    'Ny',      [];...
    'x',       [];...
    'y',       [];...
    'x1',      [];...
    'x2',      [];...
    'y1',      [];...
    'y2',      [];...
    'deltax',  [];...
    'deltay',  [];...
    'sx',  [];...
    'sy',  [];        
    };

accumulate_object;

%% visualization settings
object = 'visualization';

voi={
    'plotgrid', [0];
    };

accumulate_object;