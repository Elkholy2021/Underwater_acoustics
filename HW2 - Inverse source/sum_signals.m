function final_recieved_signal_at_transmitter=sum_signals(recieved_signals_at_transmitter)
sum=zeros(1,size(recieved_signals_at_transmitter,2));
for i =1:size(recieved_signals_at_transmitter,1)
sum=sum+recieved_signals_at_transmitter(i,:);
end
final_recieved_signal_at_transmitter=sum;
end