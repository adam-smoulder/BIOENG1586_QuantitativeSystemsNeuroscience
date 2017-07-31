function [ rasters ] = MultiNeuronRaster( trialLength, neuronCells, color )
%FOURNEURONRASTER Summary of this function goes here
%   Detailed explanation goes here
    rasters = zeros(2*length(neuronCells),trialLength);
    for i=1:length(neuronCells)
        for j=1:length(neuronCells(i).spikeTimes)
            
            timept = neuronCells(i).spikeTimes(j);
            if j  < trialLength
                rasters(2*i-1:2*i,j) = timept;
            end
            
            if color
            plot([timept timept], [i-1 i], 'k-',...
                'Color',ColorSelection2(i),'LineWidth',1)
            else
                plot([timept timept], [i-1 i], 'k-','LineWidth',1)
            end
        end
    end
    axis([0,trialLength,0,length(neuronCells)]);
    set(gca,'FontSize',14)
end

