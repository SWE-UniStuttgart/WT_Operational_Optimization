function Out = BinData(inData,V_bin,Price_bin_size)

% binning the data to wsp and price bins and evaluating their
% contributions to Damage and  Revenue based on the input time series.
% These will be fed (sorted) to the optimizer to assign values per bin

bin.V.edge = 4:V_bin:24;         %  only consider conditions where the turbine operates
bin.V.dbin = bin.V.edge(2)-bin.V.edge(1);
bin.V.center = bin.V.edge(1:end-1)+bin.V.dbin/2;

%remove negative prices
inData.Price(inData.Price<=0) = NaN;

if rem(ceil(max(inData.Price))-floor(min(inData.Price)),Price_bin_size)~=0
    remy = Price_bin_size-rem(ceil(max(inData.Price))-floor(min(inData.Price)),Price_bin_size);
    bin.Price.edge = floor(min(inData.Price))-remy/2:Price_bin_size:ceil(max(inData.Price))+remy/2;
else
    bin.Price.edge = floor(min(inData.Price))-Price_bin_size:Price_bin_size:ceil(max(inData.Price))+Price_bin_size;
end
clear remy
bin.Price.center = bin.Price.edge(1:end-1)+Price_bin_size/2;

% if mod(floor(min(inData.Price)),2) == 1 && mod(ceil(max(inData.Price)),2) ==1
%     bin.Price.edge = floor(min(inData.Price))-1:Price_bin_size:ceil(max(inData.Price))+1;
% elseif mod(floor(min(inData.Price)),2) == 1
%     bin.Price.edge = floor(min(inData.Price))-1:Price_bin_size:ceil(max(inData.Price));
% elseif mod(ceil(max(inData.Price)),2) ==1
%     bin.Price.edge = floor(min(inData.Price)):Price_bin_size:ceil(max(inData.Price))+1;
% elseif mod(floor(min(inData.Price)),2) ==0 && mod(ceil(max(inData.Price)),2) ==0
%     bin.Price.edge = floor(min(inData.Price)):Price_bin_size:ceil(max(inData.Price))+1;
% end
% bin.Price.center = bin.Price.edge(1:end-1)+Price_bin_size/2;


% discretize the data and keep indices and values
discr.V = discretize(inData.V,bin.V.edge); % discretize all velocities to bins
for i = 1:length(bin.V.center)
    [Vbin{i,1},~] = find(discr.V==i);  % indices corresponding to each of the bins
    Vbin{i,2} = inData.V(Vbin{i,1});   %#ok<*AGROW,*SAGROW>  Values corresponding to each of the bins
    Vbin{i,3} = length(Vbin{i,1})/nnz(~isnan(inData.Price)); % probability of the bin
    Pricediscr = discretize (inData.Price(Vbin{i,1}),bin.Price.edge); % for this velocity bin discretize prices in their bins
    for  ii = 1:length(bin.Price.center)
        if isempty(find(Pricediscr==ii)) %#ok<*EFIND> 
            Pbin{i,ii}(:,1) = 0; %#ok<*FNDSB> % indices of the price bin based on the total vector of pricecoming from the dicretized bin of V !!!
            Pbin{i,ii}(:,2) = 0;   % values
            Pbin{i,ii}(1,3) = 0; % rows are WSPs and colums are Prices. In each sell col 1 is global index, col2 value,col3 first row the global probability

        else
            Pbin{i,ii}(:,1) = Vbin{i,1}(find(Pricediscr==ii)); %#ok<*FNDSB> % indices of the price bin based on the total vector of pricecoming from the dicretized bin of V !!!
            Pbin{i,ii}(:,2) = inData.Price(Pbin{i,ii}(:,1));   % values
            Pbin{i,ii}(1,3) = length(Pbin{i,ii}(:,1))/nnz(~isnan(inData.Price)); % rows are WSPs and colums are Prices. In each sell col 1 is global index, col2 value,col3 first row the global probability
        end
        Prob(i,ii) =  Pbin{i,ii}(1,3);
    end
    clear Pricediscr
end
Prob = Prob/sum(sum(Prob));
[bin.mesh.V,bin.mesh.price] =  ndgrid(bin.V.center,bin.Price.center);
Out.Vdisc  = bin.mesh.V;
Out.Prdisc = bin.mesh.price;
Out.Prob   = Prob;
Out.Pbin   = Pbin;
Out.Vbin   = Vbin;
Out.bin    = bin;
