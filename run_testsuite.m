function [ output_args ] = run_testsuite( casefiles )
%RUN_TESTSUITE run test cases given by casefiles 
% casefiles can be a string of cases to run
% which should be present as folders in the testsuite folder

testsuite_folder = 'testsuite';

if (~iscell(casefiles))
    casefiles = {casefiles};
end

Ncases = length(casefiles);

errorV = zeros(Ncases,1);
errorp = zeros(Ncases,1);

for i=1:Ncases
    
    case_file = casefiles{i};

    disp(['running case: ' case_file]);
    [V,p] = main(case_file,testsuite_folder);

    ref_data_path = [testsuite_folder '/' case_file '/' case_file '_refdata.mat'];
    disp(['loading benchmark data: ' ref_data_path]);
    ref_data = load(ref_data_path);
    
    errorV(i) = max(abs(ref_data.V(:) - V(:)))
    errorp(i) = max(abs(ref_data.p(:) - p(:)))
   
    
end

%%
figure
semilogy(errorV,'x')
hold on
semilogy(errorp,'o');
grid
legend('velocity error','pressure error')
xlabel('test case number')

end

