function filtered_data = convolve(data,filter)

% discrete Fourier transform of the data and of the filter.
ft_data = fft(data);
ft_filter = fft(filter);

% perform convolutions via fast Fourier transform.

if size(ft_filter,1)==size(ft_data,1)
    filtered_data = fftshift(ifft(ft_filter.*ft_data));
elseif size(ft_filter,1)==size(ft_data,2)
    filtered_data = fftshift(ifft(ft_filter'.*ft_data));
else
    disp('convolve error. Array sizes to not match.')
    filtered_data = NaN;
end

% try
%     filtered_data = fftshift(ifft(ft_filter'.*ft_data));
% catch
%     filtered_data = fftshift(ifft(ft_filter.*ft_data));
% end

end