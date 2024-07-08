figure
semilogy(avg_rel_state_errors_intrusive,'d-','DisplayName',"intrusive ROM")
hold on
semilogy(avg_rel_state_errors(:,1,1),'x-','DisplayName',"standard operator inference")
semilogy(avg_rel_state_errors(:,1,2),'+-','DisplayName',"standard operator inference- closure-clean data")
semilogy(avg_rel_state_errors(:,3,1),'s-','DisplayName',"block-skew-symm operator inference")

xlabel("number of ROM modes")
ylabel("average relative state error")