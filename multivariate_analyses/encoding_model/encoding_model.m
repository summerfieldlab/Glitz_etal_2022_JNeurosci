%% Encoding model script 
% script used to calculate the cross-entropy loss 
% between dependent and independent blocks when predicting probability bins
% in crossvalidation - either across sessions and contexts or just across 
% sessions.

% c Christopher Summerfield, 2020, University of Oxford




% settings used to generate results presented in Figure 3f:
% varx = 0;  % variance of population codes (0 = one hot)
% dolog = 1; % 1 = use log odds; 0 = use probabilities in unit range.
% doz = 1;   % z transform values
% doH = 1;   % use crossentropy losss
% np = i;    % number of input nodes


clear all;
close all;

%%  crossvalidation

%loop through number of input nodes (i) for one-hot case
for i = 4
 
    
load('BOLD.mat');
load('pValues.mat'); %load model-derived subjective transition probability values (p heist state|door presented)

submat = [1:20 22:27 29:31]; % subject matrix

pvals = p;  %rename for convenience

varx = 0;  % variance of population codes (0 = one hot)
dolog = 1; % 1 = use log odds; 0 = use probabilities in unit range.
doz = 1;   % z transform values
doH = 1;   % use crossentropy losss
np = i;    % number of input nodes

if dolog==0;vv = linspace(0,1,np); else; vv = linspace(-2,2,np);end
 clear log_pvals logp p_list_test Xt X;


for n = 1:length(submat);
    
    sub = submat(n);
    
    for nss = 1:4 %session
        
        ss = setdiff(1:4,nss); % sessions, excluding the one we are currently on
        
        for g = 1:4;
            
            p_list_train = [];
            Y = [];
            
            for sess = 1:length(ss); %loop through the three cval sessions and parse pvals and BOLD into a list
                
                p_list_train = [p_list_train; pvals(sub).sess(ss(sess)).gem(g).data];
                
                if doz %zscore
                Y = [Y; ztransf(BOLD(sub).sess(ss(sess)).gem(g).data)];
                else
                Y = [Y; (BOLD(sub).sess(ss(sess)).gem(g).data)];
                end
            end
            
            % same as above
            p_list_train(p_list_train==0) = 0.01;
            p_list_train(p_list_train==1) = 0.99;
            
            %log transform
            if dolog;p_list_train = logodds(p_list_train);end
            
            % create input matrix (inputs x trials)
            clear X;
            for p = 1:length(p_list_train)
                if varx > 0;
                    xp  = normpdf(vv,p_list_train(p),varx);
                    X(:,p) = [xp./sum(xp) 1];
                else % one hot
                    X(:,p) = hist(p_list_train(p),vv);
                end
            end
            
            W{g} = (pinv(X')*Y)';
            
        end
        
        for g = 1:4;
            
            % get test data
            p_list_test = [pvals(sub).sess(nss).gem(g).data];
            
            % same as above
            p_list_test(p_list_test==0) = 0.01;
            p_list_test(p_list_test==1) = 0.99;
            if dolog;p_list_test = logodds(p_list_test);end
            
            if doz 
            Ynew = ztransf([BOLD(sub).sess(nss).gem(g).data]);
            else
            Ynew = ([BOLD(sub).sess(nss).gem(g).data]);
            end                
            
            for wg = 1:4;
                
                % decoding model
                Xhat = pinv(W{wg})*Ynew';
                Xhat = scaler(Xhat,[0.01 0.99]);
                Xhat = Xhat./sum(Xhat);
                
                % max (note: not used in paper)
                for p = 1:length(p_list_test);
                    va = abs(vv-p_list_test(p));
                    ii(p) = first(find(va==min(va)));
                    logp(p) = log(Xhat(ii(p),p));
                end
                
                % this computes crossentropy
                for p = 1:length(p_list_test)
                    
                    if varx > 0;
                        xp  = normpdf(vv,p_list_test(p),varx);
                        Xt(:,p) = [xp./sum(xp) 1];
                    else % one hot
                        Xt(:,p) = hist(p_list_test(p),vv);
                    end
                    
                    % compute x-entropy
                    H(p) = crossentropy(Xt(1:end-1,p),Xhat(1:end-1,p));
                    
                end
                
            
            if doH
            cross_log_pvals(n,g,wg) = mean(H);  % switch in logp to use max
            else
            cross_log_pvals(n,g,wg) = mean(logp);  % switch in logp to use max
            end
            end
        end
        
    end
    
end

cross_log_pvals = mean(cross_log_pvals,4);


log12(:,1) = cross_log_pvals(:,1,2) + cross_log_pvals(:,2,1);
log12(:,2) = cross_log_pvals(:,3,4) + cross_log_pvals(:,4,3);
[tval, pval] = masst(log12(:,1)-log12(:,2)); %one-sample ttest

[num2str(i), ' one-hot bins, crossval across sessions and contexts t = ',num2str(tval),' p = ',num2str(pval), ', crossentropy loss difference is', num2str(nanmean(log12(:,1)-log12(:,2)))] %this gives us the cval across sessions and ontexts, plotted in Figure 3f


log12b(:,1) = cross_log_pvals(:,1,2) + cross_log_pvals(:,2,1) + cross_log_pvals(:,1,1) + cross_log_pvals(:,2,2);
log12b(:,2) = cross_log_pvals(:,3,4) + cross_log_pvals(:,4,3) + cross_log_pvals(:,3,3) + cross_log_pvals(:,4,4);
[tval,pval] = masst(log12b(:,1)-log12b(:,2)); %one-sample ttest

disp([num2str(i), ' one-hot bins, crossval across sessions t = ',num2str(tval),' p = ',num2str(pval)]); %this gives us the crossval across sessions, including same context to same context (from different sessions)

end


%% helper functions, c Christopher Summerfield

function o = logodds(pb,e);

if nargin < 2;
    e=2.718281828;
end

if pb < 1;
    o = slog(pb./(1-pb),e);
else
    error('probability must be less than one');
end
end

function output = scaler(input,minmaxval)
% output = Scale(input)
% Perform an affine scaling to put data in range [0-1].
r = ranger(input);
minval=r(1); maxval=r(2);



output = (input - minval) ./ (maxval-minval);

if nargin>1;
    output=minmaxval(1)+(output.*(minmaxval(2)-minmaxval(1)));
end
end

function r = ranger(x)

if ~iscell(x)
x=x(:);
r(1)=min(x);
r(2)=max(x);
else
% cell
for n=1:length(x);
    xx=x{n};
    if ~iscell(xx)    
    rr(n,1)=min(xx(:));
    rr(n,2)=max(xx(:));
    else
        % cell within celll
        for n1=1:length(xx)        
        xxx=xx{n1};
        rrr(n,n1,1)=min(xxx(:));
        rrr(n,n1,2)=max(xxx(:));
        end
        rr(n,1)=min(rrr(n,:,1));
        rr(n,2)=max(rrr(n,:,2));
    end
        
end

r(1)=min(rr(:,1));
r(2)=max(rr(:,2));
end
end 

function z=ztransf(x,trim);
%z transforms a matrix
a=size(x);
numdim=ndims(x);

stdlx=std(x(:));
meanlx=mean(x(:));
z=(x(:)-meanlx)./stdlx;

if nargin>1
z(find(z>trim))=trim;
z(find(z<trim*-1))=trim*-1;
end


if numdim==2;
z=reshape(z,a(1),a(2));    
elseif numdim==3;
z=reshape(z,a(1),a(2),a(3));
elseif numdim==4;
z=reshape(z,a(1),a(2),a(3),a(4));
elseif numdim==5;
z=reshape(z,a(1),a(2),a(3),a(4),a(5));
end
end

function y=slog(x,base);
% computes the logarithm of any number for a given base.
y=log(x)/log(base);
end

function out=first(in,k,default);

if nargin<2;
    k=1;
end

if nargin<3;
    default = Inf;
end

if isempty(in);
    out = default;
    return
end


try
    out=in(1:k);
catch
    out=in{1:k};
end
end 

function [tval pval]= masst(x,m,alpha,tail)


if nargin < 1, 
    error('Requires at least one input argument.'); 
end

if nargin < 3, 
    m=0;
end

samplesize  = size(x,1);
xmean = mean(x);
ser = std(x) ./ sqrt(samplesize);
tval = (xmean - m) ./ ser;

if (nargout > 1)
if isnan(tval)
   pval = NaN;
else
   pval = tcdf(tval,samplesize - 1);
end

pval=squeeze(pval);
end

tval=squeeze(tval);

return;
end

