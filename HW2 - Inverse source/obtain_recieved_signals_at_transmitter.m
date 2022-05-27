function recieved_signals_at_transmitter=obtain_recieved_signals_at_transmitter(signal,x,z)
% r=100;
c=1500;
r=sqrt(x^2+z^2);

delay0=r/c;
delay=round((delay0/2.5)*length(signal));
recieved_signals_at_transmitter=zeros(1,2*length(signal));
recieved_signals_at_transmitter(1+delay:length(signal)+delay)=signal/r;

end
