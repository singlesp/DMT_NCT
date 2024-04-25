function Data_clean = regress_confound(Data, confounds, varargin)

    args = inputParser;
    args.addParameter('fitrowmask',[]);
    args.addParameter('fittype','lsq');
    args.addParameter('addconstant',true);

    args.parse(varargin{:});
    args = args.Results;
    
    if(iscell(confounds))
        confounds=cell_to_mat3d(size(Data),confounds);
    end
    
    if(ndims(confounds)==3)
        %Subj x ROI x confounds
        Data_clean=zeros(size(Data));
        %loop over ROIs
        for i = 1:size(confounds,2)
            Data_clean(:,i)=regress_confound(Data(:,i), squeeze(confounds(:,i,:)), varargin{:});
        end
        return;
    end


    fittype=args.fittype;

    %normpool=which rows to use for fitting
    fitrowmask=true(size(Data,1),1);
    if(~isempty(args.fitrowmask))
        fitrowmask=args.fitrowmask(:);
    end

    fitcolmask=~all(isnan(Data),1);

    %exclude rows with nans in any columns (except columns that are ALL nans)
    fitrowmask=fitrowmask & ~any(isnan(Data(:,fitcolmask)),2);
    fitrowmask=fitrowmask & ~any(isnan(confounds),2);

    %add a constant term if user requested it and there isn't one already 
    confounds_const=all(confounds==repmat(confounds(1,:),size(confounds,1),1));
    if(args.addconstant && ~any(confounds_const))
        confounds(:,end+1)=1;
        confounds_const(end+1)=true;
    end

    %make sure there's ONLY ONE constant term in the final confounds
    if(any(confounds_const))
        confounds=confounds(:,~confounds_const);
        confounds(:,end+1)=1;
    end
    
    confounds_const=all(confounds==repmat(confounds(1,:),size(confounds,1),1));
    
    B=nan(size(confounds,2),size(Data,2));
    if(strcmpi(fittype,'lsq'))
        B(:,fitcolmask)=confounds(fitrowmask,:)\Data(fitrowmask,fitcolmask);
    elseif(strcmpi(fittype,'robust'))

        for i = 1:size(Data,2)
            if(~fitcolmask(i))
                continue;
            end

            if(any(confounds_const))
                %if confound included >=1 constant column, remove it from the
                %input to robustfit, and then split the output since constant
                %goes first
                btmp=robustfit(confounds(fitrowmask,~confounds_const),Data(fitrowmask,i));
                B(confounds_const,i)=0;
                B(find(confounds_const,1),i)=btmp(1);
                B(~confounds_const,i)=btmp(2:end);
            else
                btmp=robustfit(confounds(fitrowmask,~confounds_const),Data(fitrowmask,i),[],[],'off');
                B(:,i)=btmp;
            end
        end
    end

    Data_clean=Data-confounds*B;
end

function output3d = cell_to_mat3d(datasize,inputcells)
    if(~iscell(inputcells) || isempty(inputcells))
        output3d=inputcells;
        return;
    end
    output3d=zeros(datasize(1),datasize(2),numel(inputcells));
    for i = 1:numel(inputcells)
        c=inputcells{i};
        if(all(size(c)==datasize))
            %donothing
        elseif(size(c,2)==1)
            c=repmat(c,1,datasize(2));
        else
            error('invalid confounds{} size. Cell items must be either single vector or full Data-size matrix');
        end
        output3d(:,:,i)=c;
    end
end