function h = replay(Fem,varargin)

    FPS = 60; 
    [AX, cmap, cAX, BK] = deal([]);
    [Gif, Palidrome, topo] = deal(false);
    
    for ii = 1:2:length(varargin)
        if     strcmpi(varargin{ii},'axis'), AX = varargin{ii+1};
        elseif strcmpi(varargin{ii},'fps'),  FPS = varargin{ii+1};    
        elseif strcmpi(varargin{ii},'caxis'),  cAX = varargin{ii+1};    
        elseif strcmpi(varargin{ii},'palidrome'), Palidrome = true;    
        elseif strcmpi(varargin{ii},'gif'), Gif = true;    
        elseif strcmpi(varargin{ii},'background'), BK = varargin{ii+1};  
        elseif strcmpi(varargin{ii},'colormap'), cmap = varargin{ii+1};    
        elseif strcmpi(varargin{ii},'topo'), topo = true;
        end
    end
    
    t = Fem.solver.sol.tout;
    y = Fem.solver.sol.yout;

    framestep = fps(t(:),FPS);
    frameIteration = 1:framestep:numel(t);

    if Palidrome
        frameIteration = [frameIteration, ...
                          fliplr(frameIteration(2:end))]; 
    end
    
    % make clean new figure
    figure(101); clf;
    
    for ii = frameIteration

        % overwrite the current solutions
        Fem.solver.sol.x = y(ii,:).';

        Fem = Fem.compute();
        Fem.options.Display(Fem);
        
        background();
        if ~isempty(AX), axis(AX); end
        if ~isempty(cAX), caxis(cAX); end
        if ~isempty(cmap), colormap(cmap); end
        if ~isempty(BK), background(BK); end
        drawnow();
        
        if Gif
           if ii == 1, snapshotGif('Start',FPS);
           else, snapshotGif('',FPS);
           end
        end
    end   
end

% computes frame steps based on FPS and tout
function A = fps(time,FPS)
    dt = max(diff(time));
    A = max([round((1/FPS)/dt),1]);
end 

% produce gif 
function snapshotGif(Request,FPS)       
    switch(Request)
        case('Start')
            Name = 'fem';
            filename = string([Name,'_', char(datetime(now,...
                'ConvertFrom','datenum')),'.gif']);
            filename = erase(filename,[":"," "]);
            gif(char(filename),'frame',gcf,'nodither',...
            'Timestep',1/FPS);
        otherwise
            gif;
    end
end    