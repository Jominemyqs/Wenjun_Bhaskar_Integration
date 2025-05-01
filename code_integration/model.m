function pattern = model(modelpar)

path(path,'Models/Reaction_Diffusion/');
path(path,'Models/Spiral_Wave/');
path(path,'Models/');

% Generate solution
switch modelpar.model
    case 'Brusselator'
        [U,~,X,Y] = Brusselator(modelpar.a,modelpar.b);
    case 'SH'
        [U,X,Y]= SH_2D(modelpar.a,modelpar.b, modelpar.init);
     
    case 'SH_1D'
        [U,X]= SH_1D(modelpar.a,modelpar.b);
    case 'GS'
        [U,V,X,Y]= Gray_Scott(modelpar.a,modelpar.b);
    case 'Schnakenberg'
        [U,X,Y]= Schnakenberg(modelpar.a,modelpar.b);
    case 'Bullara'
        [U,X,Y] = bullara_code(modelpar.a,modelpar.b,100, 5e7, 0, 0);
        idx0 = U == 0;
        idx1 = U == 1;
        U(idx0) = 1;
        U(idx1) = 0;
    case 'Barkley'
        anew = num2str(modelpar.a,'%1.6f');
        bnew = num2str(modelpar.b,'%1.6f');
        
        
        
        fid = fopen('Models/Spiral_Wave/Barkley/task.txt');
        C=textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        C{1}{1} = anew; C{1}{2} = bnew;
        
        writecell(C{1},'Models/Spiral_Wave/Barkley/task.txt','QuoteStrings',0)
%                 if b > 0.05
%                     C{1}{11}(1:5) = num2str(100000);
%                 end
        cd Models/Spiral_Wave/Barkley
        system('rm tip.txt')
        system('touch tip.txt')
        system('./ezspiral')
        
        cd ../../..
        tip_pts = readtable('Models/Spiral_Wave/Barkley/tip.txt');
        tip_pts = table2array(tip_pts);
        U = readtable('Models/Spiral_Wave/Barkley/fc.txt');
        U = U.Var2;
        U = reshape(U(5:end),501,501);
        X = linspace(0,100,501); Y = X;
        [X,Y] = meshgrid(X,Y);
        if size(tip_pts,1) > 1000 || size(tip_pts,1)==1
            tip_pts = tip_pts(1000:end,:);

        end
        idx = find( (tip_pts(:,2)>98)+(tip_pts(:,2)<2) +...
            + (tip_pts(:,3)>98) + (tip_pts(:,3)<2)  );
        idx
        if length(idx) > 0
            idx = [1; idx;size(tip_pts,1) ];
            d = diff(idx);
            max_diff = [idx( find(d == max(d))):idx( find(d == max(d))+1)  ] ;
            tip_pts = tip_pts(max_diff,:);
        % if ~isempty(idx)
        %     tip_pts = tip_pts(1:idx(1),:);
        end
        %U
                        

     case 'Bar-Eiswirth'
        if modelpar.a < 0
            modelpar.a = 1E-14;
        end
        anew = num2str(1/(modelpar.a),'%1.6f');
        bnew = num2str(modelpar.b,'%1.6f');
        
        length(anew)
        length(bnew)
        
        fid = fopen('Models/Spiral_Wave/Bar-Eiswirth/task.txt');
        C=textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        C{1}{2} = bnew; C{1}{3} = anew;
        writecell(C{1},'Models/Spiral_Wave/Bar-Eiswirth/task.txt','QuoteStrings',0)
        
        cd Models/Spiral_Wave/Bar-Eiswirth
        system('rm tip.txt')
        system('touch tip.txt')
        system('./ezspiral')
        
        cd ../../../
        tip_pts = readtable('Models/Spiral_Wave/Bar-Eiswirth/tip.txt');
        tip_pts = table2array(tip_pts);

        U = readtable('Models/Spiral_Wave/Bar-Eiswirth/fc.txt');
        U = U.Var2;
        U = reshape(U(5:end),501,501);

        % if strcmp(featpar.feature,'turbulence') == 0
        if size(tip_pts,1) > 2000 || size(tip_pts,1)==1
            tip_pts = tip_pts(2000:end,:);

        end
        idx = find( (tip_pts(:,2)>49)+(tip_pts(:,2)<1) +...
            + (tip_pts(:,3)>49) + (tip_pts(:,3)<1)  );
        % else
        %     idx = [];
        % end
        if ~isempty(idx)
            tip_pts = tip_pts(1:idx(1),:);
        end
                    
    

    case 'Rossler'
                
                
        cnew = num2str(modelpar.a,'%1.6f');
        anew = num2str(modelpar.b,'%1.6f');
        
        fid = fopen('Models/Spiral_Wave/Rossler/task.txt');
        C=textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        C{1}{3} = cnew; C{1}{1} = anew;
        writecell(C{1},'Models/Spiral_Wave/Rossler/task.txt','QuoteStrings',0)

        cd Models/Spiral_Wave/Rossler

        system('./ezspiral')
        cd ../../..

        U = readtable('Models/Spiral_Wave/Rossler/fc.txt');

        
        U = U.Var2; 
        U = reshape(U(5:end),526,526);

        X = linspace(0,250,526); Y = X;
        [X,Y] = meshgrid(X,Y);
%                 cd ../
        tip_pts = [];

    case 'Bhaskar'
        % If needed, add path to Bhaskar code:
        % addpath('path_to_Bhaskar_Code');

        % If ParamSweep_1.mat is used:
        if isfile('ParamSweep_1.mat')
            load('ParamSweep_1.mat','n');  % typically has 'n'
        end
        
        % Overwrite if needed
        RR = modelpar.RR; % Orange-Orange
        RG = modelpar.RG; % Blue-Orange
        GG = modelpar.GG; % Blue-Blue
        pol_val = modelpar.pol_val; 
        cell_pop_prop = modelpar.cell_pop_prop;
        save('ParamSweep_1.mat','RR','RG','GG','pol_val','cell_pop_prop','-append');

        % Run the Bhaskar code (task_id=1)
        coculture_model(1);

        % Example: load time-step 500000
        pos_file = 'ParamSweep_1_Output/Pos_0500000.dat';
        type_file = 'ParamSweep_1_Output/Types_0500000.dat';
        
        if ~exist(pos_file,'file') || ~exist(type_file,'file')
            warning('Bhaskar output not found at time 500000. pattern=[]');
            pattern = [];
            return;
        end
        
        % Load positions & types
        posdata = load(pos_file);
        posdata = posdata(:);  % ensure column
        typedata = load(type_file);
        typedata = typedata(:);        
        
        if length(posdata) ~= length(typedata)
            warning('Mismatch: positions vs. types. pattern=[]');
            pattern = [];
            return;
        end

        % Create final pattern array [X, Y, type]
        X = real(posdata);
        Y = imag(posdata);
        pattern = [X, Y, typedata];
        
        % Return immediately, skipping 'select set' logic
        return
    % ============ END BHASKAR CASE ============

end % end switch modelpar.model

% Select set
switch modelpar.sets
    case 'pos'
        if max(U(:)) - min(U(:)) < 1E-5
            X = []; Y = [];
        else
            X = X(U > min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  quantile(U(:), 0.7) );
            if modelpar.xdim == 2 
                Y = Y(U > min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  quantile(U(:), 0.7));
            end
        end
    case 'neg'
        if max(U(:)) - min(U(:)) < 1E-5
            X = []; Y = [];
        else
            X = X(U < min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  < quantile(U(:), 0.3));
            if modelpar.xdim == 2 
                Y = Y(U < min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  < quantile(U(:), 0.3));
            end
        end
    % case 'tip_points'
    %     if ~isempty(tip_pts)
    %         X = tip_pts(:,2);
    %         Y = tip_pts(:,3);
    %     else
    %         X = [];
    %         Y = [];
    %     end
    case 'pos-filter'
        % This is for fish
        if max(U(:)) - min(U(:)) < 1E-5
            X = []; Y = [];
        else
            X = X(U > min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  quantile(U(:), 0.7) );
            if modelpar.xdim == 2 
                Y = Y(U > min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  quantile(U(:), 0.7));
            end
        end
        pattern = [X,Y];
        Dists = sqrt( (X-X').^2 + (Y - Y').^2 );
        Dists_nn = zeros(1, size(X,1));
        for i = 1:length(Dists_nn)
            Dists(i,:) = sort( Dists(i,:) );
            Dists_nn(i) = sum(Dists(i,2:11).^2/10);
        end
        idx = Dists_nn < quantile(Dists_nn,0.9);
        X = X(idx);
        Y = Y(idx);
    case 'neg-filter'
        if max(U(:)) - min(U(:)) < 1E-5
            X = []; Y = [];
        else
            X = X(U < min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  < quantile(U(:), 0.3));
            if modelpar.xdim == 2 
                Y = Y(U < min(U(:)) + modelpar.threshold*( max(U(:))-min(U(:))) );%  < quantile(U(:), 0.3));
            end
        end
        pattern = [X,Y];
        Dists = sqrt( (X-X').^2 + (Y - Y').^2 );
        Dists_nn = zeros(1, size(X,1));
        for i = 1:length(Dists_nn)
            Dists(i,:) = sort( Dists(i,:) );
            Dists_nn(i) = sum(Dists(i,2:11).^2/10);
        end
        idx = Dists_nn < quantile(Dists_nn,0.9);
        X = X(idx);
        Y = Y(idx);
end

if strcmp( modelpar.sets, 'tip_points')
    pattern = tip_pts;
elseif strcmp( modelpar.sets, 'homogeneous')
    pattern = max(U(:))-min(U(:)) <= 1E-10

else
    if modelpar.xdim == 2 
        pattern = [X,Y];
    else
        pattern = X;
    end
end

% if size(pattern,1) > 1
%     if strcmp( modelpar.sets, 'tip_points')
%         figure(29);plot(pattern(:,2),pattern(:,3))
%     else
%         figure(29);scatter(pattern(:,1),pattern(:,2))
%     end
% end







