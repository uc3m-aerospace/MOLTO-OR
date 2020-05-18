function varargout = plotnsga_app(varargin)

%PLOTNSGA M-file for plotnsga.fig
%      PLOTNSGA, by itself, creates a new PLOTNSGA or raises the existing
%      singleton*.
%
%      H = PLOTNSGA returns the handle to a new PLOTNSGA or the handle to
%      the existing singleton*.
%
%      PLOTNSGA('Property','Value',...) creates a new PLOTNSGA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to plotnsga_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PLOTNSGA('CALLBACK') and PLOTNSGA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PLOTNSGA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%*************************************************************************
% Syntax:
%   plotnsga(result)
%     Plot the optimization result.
% 
%   plotnsga(result, curGen)
%     Plot the optimization result with specified generation.
% 
%   plotnsga('pops.txt')
%     Load from population file and plot the result.
%     Note: A global variable "oldresult" which contains the population loaded from file
%       would be created.
%*************************************************************************

% Edit the above text to modify the response to help plotnsga

% Last Modified by GUIDE v2.5 11-Jul-2011 15:53:30

% Begin initialization code - DO NOT EDIT 


app = varargin{1}.app;


if nargout
    [varargout{1:nargout}] = plotnsga_OpeningFcn(app,varargin);  %AQUI VOY A TENER QUE METER LA FUNCION QUE MANDE LOS PLOTS
else
    %gui_mainfcn(gui_State, varargin{:});
    plotnsga_OpeningFcn(app,varargin); 
end
% End initialization code - DO NOT EDIT


% --- Executes just before plotnsga is made visible.
function plotnsga_OpeningFcn(app,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)



% UIWAIT makes plotnsga wait for user response (see UIRESUME)
% uiwait(handles.figure1);

app.output = 0; 
%*************************************************************************
% 1. Save the result
%*************************************************************************
handles.bLoadFromFile = 0;          % Load from the population file.

        app.result  = varargin{1}{1};
        app.currentGen = varargin{1}{2};




%*************************************************************************
% 2. Initialize the population ID listbox
%*************************************************************************
popsize = size(app.result.pops, 1);
strList = repmat({''}, [1,popsize]);
for i = 1:popsize
    strList{i} = sprintf('%d', i);
end

 curSel = app.currentGen;   % the generation ID of population which would be ploted
 set(app.listPop, 'string', strList);
 set(app.listPop, 'value', curSel);

app.GenerationListBox.Value = num2str(curSel);
%*************************************************************************
% 3. Update data and plot population
%*************************************************************************
%dispState(handles, curSel);
plotPopulation( app, curSel );
%guidata(hObject,handles);








% --- Outputs from this function are returned to the command line.
function varargout = plotnsga_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;







function plotPopulation(app, gen)
% Function: plotPopulation(handles, gen)
% Description: Plot population with the first two objective values.
% Parameters: 
%   gen : the generation ID
%
%         LSSSSWC, NWPU
%    Revision: 1.3  Data: 2011-07-26
%*************************************************************************

cla(app.UIAxes_pareto)

%*************************************************************************
% Initialize data
%*************************************************************************
pop     = app.result.pops(gen, :);
obj     = vertcat(pop.obj);
numObj  = length(pop(1).obj);
strObj1 = 'objective 1';
strObj2 = 'objective 2';
strObj3 = 'objective 3';

opt = app.result.opt;
strObj1 = ['obj 1 : ', opt.nameObj{1}];
strObj2 = ['obj 2 : ', opt.nameObj{2}];
maxGen = size(app.result.pops, 1);



% Determin if reference points exist
refPoints = [];
refPlotStyle = {'kd', ...
        'LineWidth', 1,...
        'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', 'g',...
        'MarkerSize',6};
if( isfield(app.result, 'opt') && ~isempty(app.result.opt.refPoints) )
    refPoints = app.result.opt.refPoints;
end


%*************************************************************************
% Plot population with different methods for every "numObj" number
%*************************************************************************
if(numObj == 2)
      app.pareto = plot(app.UIAxes_pareto,obj(:,1), obj(:,2), 'ob');
%     app.UIAxes_pareto.YData = obj(:,2);
%     app.UIAxes_pareto.XData = obj(:,1);
    xlabel(app.UIAxes_pareto,strObj1, 'interpreter', 'latex');
    ylabel(app.UIAxes_pareto,strObj2, 'interpreter', 'latex'); 
    drawnow

    % plot reference points
    if(~isempty(refPoints))
        hold on
        plot(app.UIAxes_pareto,refPoints(:, 1), refPoints(:, 2), refPlotStyle{:});
    end
elseif(numObj == 3)
    [az, el] = view;    % Save the last view angles
    plot3(obj(:,1), obj(:,2), obj(:,3), 'ob');
    view(az,el);        % Restore the last view angles
    xlabel(strObj1, 'interpreter', 'none');
    ylabel(strObj2, 'interpreter', 'none');
    zlabel(strObj3, 'interpreter', 'none');

    % plot reference points
    if(~isempty(refPoints))
        hold on
        plot3(refPoints(:, 1), refPoints(:, 2), refPoints(:, 3), refPlotStyle{:});
    end
else
    plot(obj', 'b-');
    xlim([1,numObj]);
    set(gca, 'XGrid', 'on');
    xlabel('Objective number');
    ylabel('Objective value');

    % plot reference points
    if(~isempty(refPoints))
        hold on
        refPlotStyle{1} = 'gd-';
        plot(refPoints', refPlotStyle{:});
    end
end


%*************************************************************************
% Common operations
%*************************************************************************
% Title
strTitle = sprintf('Generation %d / %d', gen, maxGen);

title(app.UIAxes_pareto,strTitle, 'interpreter', 'latex');
drawnow






