number_of_cells = size(neuron.A , 2); 
spatial_footprints =  zeros(number_of_cells, neuron.options.d1, ...
    neuron.options.d2);

for n=1:number_of_cells
   spatial_footprints(n,:,:) = reshape(neuron.A(:,n), ...
       neuron.options.d1 , neuron.options.d2);
end