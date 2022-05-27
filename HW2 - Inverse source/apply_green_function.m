function received_signal = apply_green_function(signal,c,h,ns,xr,zr,xs,zs)
recieved_signals_at_transmitter=zeros(1,6000);
t=linspace(1,10,6000);
% for j =-(ns/2)+d:ns/2
for j =0:ns-1    
if mod(j,2)==0 % j is even
    r=sqrt((xr-xs)^2+(zs-zr+j*h)^2);
else         % j is odd
    r=sqrt((xr-xs)^2+(zs+zr-(j+1)*h)^2);
end        
delay0=r/c;
delay=uint64((delay0/2.5)*length(signal));

recieved_signals_at_transmitter(1+delay:length(signal)+delay)=(signal/r)...
+recieved_signals_at_transmitter(1+delay:length(signal)+delay);

end
received_signal=recieved_signals_at_transmitter;
end


